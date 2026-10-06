#!/usr/bin/env python3
"""Deterministic end-to-end entrypoint for isolate screening."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Any

from analyze_similarity_distribution import run_diagnostics
from expand_cluster_members import discover_members
from refinement_expansion import assess_progress, member_coverage, schedule_focused_batch, validate_expansion_policy
from collect_provenance import collect
from classify_with_mlst import classify, route_for_organism
from common import CONFIG_DIR, WorkflowError, effective_database_root, load_json, sha256_file, write_json
from fetch_reference_genomes import fetch_genomes
from generate_report import generate_report
from inspect_input import inspect
from interpret_snp_resolution import interpret as interpret_snp_resolution
from parse_mashpit_results import interpret, load_candidates, locate_candidate_file, locate_tree_image
from render_snp_tree import render_from_interpretation
from run_assembly_workflow import run_workflow
from run_mashpit import run_mashpit
from run_ska import run_ska
from run_mummer import run_mummer
from run_minimap import run_minimap
from interpret_assembly_comparison import interpret as interpret_assembly_comparison
from select_snp_targets import select_targets
from validate_assembly import assess


def sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fasta", ".fna", ".fastq", ".fq", ".fa"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    for mate in ("_R1", "_1", ".R1", ".1", "-R1", "-1"):
        if name.endswith(mate):
            name = name[: -len(mate)]
            break
    return name


def summary_text(result: dict[str, Any]) -> str:
    lines = [f"Sample: {result['sample']}", f"Status: {result['status']}"]
    if result.get("input_type"):
        lines.append(f"Input: {result['input_type']}")
    if result.get("organism"):
        lines.append(f"Organism: {result['organism']}")
    if result.get("assembly_qc"):
        lines.append(f"Assembly QC: {result['assembly_qc']['status']}")
    if result.get("database"):
        lines.append(f"Mashpit database: {result['database']['name']} {result['database']['version']}")
    match = result.get("mashpit_result") or {}
    best = match.get("best_candidate")
    if best:
        lines.extend([
            f"Best candidate cluster: {best['cluster']} (score {best['score']:.3f})",
            f"Screening result: {match['screening_result']}",
            "Interpretation: This isolate is most similar to representatives associated with the candidate cluster.",
            "This screening result does not by itself establish outbreak relatedness.",
            "Recommended next step: Run a validated high-resolution SNP comparison against representative isolates.",
        ])
    snp = result.get("snp_resolution")
    if snp and snp.get("ranked") and snp.get("nearest_snp_distance") is not None:
        lines.append(
            f"{snp.get('backend', 'ska2').upper()} comparison: minimum observed distance {snp['nearest_snp_distance']:.2f} SNPs; "
            f"{len(snp.get('nearest_samples', []))} equally nearest reference(s). "
            "This is a comparison among qualifying references, not strain or outbreak confirmation."
        )
        lines.append(
            "Nearest SNP cluster(s) among examined references: "
            + (", ".join(snp.get("nearest_clusters", [])) or "unresolved")
            + ". Correct cluster concordance does not establish recovery of the closest genome."
        )
        if snp.get("cluster_status"):
            lines.append("Cluster resolution: " + snp["cluster_status"] + "; supported stored label(s): "
                         + (", ".join(snp.get("cluster_candidates", [])) or "none")
                         + ". Genome resolution: " + snp["genome_status"] + ".")
        uncertainty = (snp.get("expansion") or {}).get("remaining_uncertainty") or {}
        if uncertainty.get("budget_or_round_limit_prevented_required_exploration"):
            lines.append("Focused or alternative exploration remained incomplete at the resource or round limit.")
    if snp and snp.get("ranking_basis") == "common_finalist_regions":
        lines.append("Common-region finalist comparison: cluster " + snp.get("cluster_status", "INSUFFICIENT_DATA") +
                     "; genome " + snp.get("genome_status", "INSUFFICIENT_DATA") + ".")
        lines.append("Supported stored cluster: " + (snp.get("nearest_cluster") or "unresolved") + ".")
        lines.append(snp.get("confidence", {}).get("statement", ""))
        uncertainty = (snp.get("expansion") or {}).get("remaining_uncertainty") or {}
        if uncertainty.get("budget_or_round_limit_prevented_required_exploration"):
            lines.append("Focused or alternative exploration remained incomplete at the resource or round limit.")
    if result.get("stop_reason"):
        lines.append(f"Analysis stopped: {result['stop_reason']}")
    warnings = result.get("warnings", [])
    if warnings:
        lines.append("Warnings: " + "; ".join(warnings))
    return "\n".join(lines)


def run_snp_resolution(
    mashpit_output_dir: Path, assembly: Path, best_cluster: str, output_dir: Path,
    database_metadata: dict[str, Any] | None = None, expand: bool = False,
    query_reads: list[str] | None = None, selection_mode: str | None = None,
    query_accession: str | None = None,
    refinement_policy_path: Path | None = None,
    snp_backend: str | None = None, assembly_manifest: Path | None = None,
    comparison_cache: Path | None = None, assembly_policy_path: Path | None = None,
    comparison_workers: int | None = None,
) -> dict[str, Any]:
    policy = load_json(CONFIG_DIR / "snp-resolution-policy.json")
    workflow = load_json(CONFIG_DIR / "workflow.json")
    refinement_policy_path = refinement_policy_path or CONFIG_DIR / "refinement-policy.json"
    refinement_policy = load_json(refinement_policy_path)
    expansion_policy = refinement_policy["expansion"]
    validate_expansion_policy(expansion_policy)
    kmer_size = workflow["snp_resolution"]["kmer_size"]
    backend = snp_backend or ("ska" if query_reads else workflow["snp_resolution"].get("assembly_backend","minimap2"))
    if backend not in {"ska","mummer","minimap2"}: raise WorkflowError("Unknown SNP backend")
    if backend != "ska" and query_reads:
        raise WorkflowError("Assembly-comparison backends accept assemblies; paired-read refinement uses SKA")
    comparison_cache = comparison_cache or output_dir/"comparison_cache"
    comparison_runs = []
    commands: list[list[str]] = []
    targets_result = select_targets(
        mashpit_output_dir, policy, database_metadata,
        workflow["mashpit"] if database_metadata is not None else None,
        mode=selection_mode, query_accession=query_accession,
    )
    write_json(output_dir / "targets.json", targets_result)
    if targets_result["status"] != "SELECTED":
        return {"snp_resolution": targets_result, "commands": commands}
    if database_metadata is None:
        database_metadata = load_json(mashpit_output_dir / "mashpit_run.json")["database"]
    targets = {item["accession"]: item for item in targets_result["targets"]}
    initial_accessions = set(targets)
    batch = list(targets.values())
    genomes: dict[str, Any] = {"QUERY": query_reads if query_reads else str(assembly)}
    attempted: set[str] = set()
    unavailable: set[str] = set()
    rounds = []
    membership = None
    members: list[dict[str, Any]] = []
    previous = None
    stable_rounds = 0
    allocation_state: dict[str, Any] = {}
    stop_reason = "representatives_only"
    warnings = list(targets_result.get("warnings", []))
    # The total cap includes attempted downloads, even when some fail.
    cap = expansion_policy["max_total_genomes"]
    expansion_cap = expansion_policy["max_additional_reference_attempts"]
    if not isinstance(expansion_cap, int) or isinstance(expansion_cap, bool) or expansion_cap <= 0:
        raise WorkflowError("Expansion max_additional_reference_attempts must be positive.")
    if len(batch) > cap:
        raise WorkflowError("Representative selection exceeds the configured total reference-attempt budget.")
    interpretation = {"status": "ERROR", "ranked": [], "warnings": []}
    last_compared_count = 0
    for round_index in range(expansion_policy["max_rounds"] if expand else 1):
        round_dir = output_dir / f"round_{round_index:02d}"
        accessions = [item["accession"] for item in batch]
        attempted.update(accessions)
        targets.update({item["accession"]: item for item in batch})
        write_json(round_dir / "requested.json", {"targets": batch})
        try:
            fetch_args = {"local_manifest":assembly_manifest} if assembly_manifest else {}
            fetch_result = fetch_genomes(accessions, round_dir / "genomes", policy["download_attempts"], policy["download_retry_delay_seconds"], **fetch_args)
            commands.extend(fetch_result["commands"])
            genomes.update(fetch_result["verified"])
            unavailable.update(fetch_result["unavailable"])
            if len(genomes) < 2:
                interpretation = {"status": "INSUFFICIENT_DATA", "ranked": [], "warnings": ["No usable reference genomes could be downloaded."]}
            elif backend != "ska":
                comparison = (run_minimap if backend=="minimap2" else run_mummer)(genomes,round_dir/backend,comparison_cache,assembly_policy_path,comparison_workers)
                commands.extend(comparison["commands"])
                comparison_runs.append({key:comparison[key] for key in ("alignment_jobs","cache_hits","elapsed_seconds","status","tools","policy","policy_sha256","implementation_sha256","cache_directory")})
                interpretation = interpret_assembly_comparison(comparison["comparisons"],list(targets.values()),best_cluster,comparison["policy"])
            else:
                ska_result = run_ska(genomes, round_dir / "ska", kmer_size)
                commands.extend(ska_result["commands"])
                interpretation = interpret_snp_resolution(
                    ska_result["distances"], list(targets.values()), best_cluster,
                    refinement_policy["comparability"],
                )
        except (WorkflowError, OSError, ValueError) as error:
            stop_reason = "round_failed"
            if interpretation["status"] == "ERROR":
                interpretation["error"] = str(error)
            warnings.append(f"Refinement round {round_index} failed: {error}. Earlier successful comparisons are retained.")
            record = {"round": round_index, "status": "FAILED", "error": str(error),
                      "requested_accessions": accessions, "attempted_total": len(attempted),
                      "compared_total": last_compared_count, "feedback": {"action": "STOP", "reason": stop_reason}}
            rounds.append(record)
            write_json(round_dir / "round.json", record)
            break
        last_compared_count = len(genomes) - 1
        record = {
            "round": round_index, "requested_accessions": accessions,
            "downloaded_accessions": sorted(fetch_result["verified"]),
            "unavailable_accessions": fetch_result["unavailable"], "attempted_total": len(attempted),
            "compared_total": len(genomes) - 1,
            "nearest_samples": interpretation.get("nearest_samples", []),
            "nearest_clusters": interpretation.get("nearest_clusters", []),
            "nearest_snp_distance": interpretation.get("nearest_snp_distance"),
            "excluded_comparisons": interpretation.get("excluded_comparisons", []),
        }
        if not expand:
            record["feedback"] = {"action": "STOP", "reason": "representatives_only"}
            rounds.append(record)
            write_json(round_dir / "interpretation.json", interpretation)
            write_json(round_dir / "round.json", record)
            break
        if len(attempted) >= cap:
            stop_reason = "hard_resource_ceiling"
        elif len(attempted) - targets_result["selected_genomes"] >= expansion_cap:
            stop_reason = "additional_expansion_budget"
        elif round_index + 1 >= expansion_policy["max_rounds"]:
            stop_reason = "round_limit"
        if stop_reason in {"hard_resource_ceiling", "additional_expansion_budget", "round_limit"}:
            record["feedback"] = {"action": "STOP", "reason": stop_reason, "stable": False}
            record["allocation"] = {"strategy": expansion_policy["allocation_strategy"], "chosen": [], "not_requested_reason": stop_reason}
            rounds.append(record)
            write_json(round_dir / "interpretation.json", interpretation)
            write_json(round_dir / "round.json", record)
            break
        if membership is None:
            membership = discover_members(database_metadata, targets_result["selected_clusters"], output_dir / "membership", expansion_policy)
            if membership["status"] != "AVAILABLE":
                stop_reason = "membership_unavailable"
                warnings.append("Cluster expansion unavailable: " + membership.get("error", "unknown error"))
            else:
                # Include returned but unselected representatives in the same
                # finite pool as exact-release members; never reintroduce self.
                combined = {item["accession"]: item for item in membership["members"] if item["accession"] != query_accession}
                for item in targets_result["decisions"]:
                    if item["cluster"] not in targets_result["selected_clusters"] or item["accession"] == query_accession:
                        continue
                    if item["accession"] in combined and combined[item["accession"]]["cluster"] != item["cluster"]:
                        membership["status"] = "CONFLICT"
                        stop_reason = "membership_conflict"
                        warnings.append("Pinned membership and Mashpit disagree about a reference's cluster; expansion stopped and earlier comparisons are retained.")
                        write_json(output_dir / "membership" / "membership.json", membership)
                        break
                    combined[item["accession"]] = item
                if membership["status"] == "AVAILABLE":
                    members = list(combined.values())
        if stop_reason in {"membership_unavailable", "membership_conflict"}:
            record["feedback"] = {"action": "STOP", "reason": stop_reason, "stable": False}
            record["allocation"] = {"strategy": expansion_policy["allocation_strategy"], "chosen": [], "not_requested_reason": stop_reason}
            rounds.append(record)
            write_json(round_dir / "interpretation.json", interpretation)
            write_json(round_dir / "round.json", record)
            break
        coverage = member_coverage(members, initial_accessions, attempted,
                                   set(genomes) - {"QUERY"}, interpretation,
                                   targets_result["selected_clusters"])
        record["member_coverage"] = coverage
        if not any(row["unexamined_additional"] for row in coverage):
            stop_reason = "available_member_pool_exhausted" if interpretation.get("ranked") else "insufficient_evidence_pool_exhausted"
            record["feedback"] = {"action": "STOP", "reason": stop_reason, "stable": False}
            record["allocation"] = {"strategy": expansion_policy["allocation_strategy"], "chosen": [], "not_requested_reason": stop_reason}
            rounds.append(record)
            write_json(round_dir / "interpretation.json", interpretation)
            write_json(round_dir / "round.json", record)
            break
        remaining_budget = min(cap - len(attempted), expansion_cap - (len(attempted) - targets_result["selected_genomes"]))
        planned, next_state, allocation = schedule_focused_batch(
            members, attempted, targets_result["selected_clusters"], interpretation,
            expansion_policy, remaining_budget, expansion_policy["max_rounds"] - round_index - 1,
            allocation_state, query_accession,
        )
        decision = assess_progress(previous, interpretation, coverage, expansion_policy,
                                   allocation["finite_leader_completion_feasible"], stable_rounds)
        stable_rounds = decision["stable_rounds"]
        record["feedback"] = decision
        record["allocation"] = allocation if decision["action"] != "STOP" else {**allocation, "chosen": [], "not_requested_reason": "stable_sampled_neighborhood"}
        rounds.append(record)
        write_json(round_dir / "interpretation.json", interpretation)
        write_json(round_dir / "round.json", record)
        if decision["action"] == "STOP":
            stop_reason = decision["reason"]
            break
        if not planned:
            stop_reason = "insufficient_evidence" if not interpretation.get("ranked") else "available_member_pool_exhausted"
            break
        batch = planned
        allocation_state = next_state
        previous = interpretation
    pending_count = sum(item["accession"] not in attempted for item in members) if membership and membership["status"] == "AVAILABLE" else None
    final_member_coverage = (member_coverage(members, initial_accessions, attempted,
                             set(genomes) - {"QUERY"}, interpretation,
                             targets_result["selected_clusters"])
                             if membership and membership["status"] == "AVAILABLE" else None)
    final_leaders = set(interpretation.get("nearest_clusters", []))
    unresolved_clusters = ([row["cluster"] for row in final_member_coverage
                            if not any(item.get("cluster") == row["cluster"] for item in interpretation.get("ranked", []))]
                           if final_member_coverage is not None else None)
    focused_below_min = ([row["cluster"] for row in final_member_coverage
                          if row["cluster"] in final_leaders
                          and row["qualifying_additional"] < expansion_policy["min_focused_qualifying_per_leader"]
                          and row["unexamined_additional"]]
                         if final_member_coverage is not None else None)
    resource_stop = stop_reason in {"hard_resource_ceiling", "additional_expansion_budget", "round_limit"}
    expansion = {
        "enabled": expand, "policy": expansion_policy,
        "policy_version": refinement_policy["policy_version"],
        "policy_path": str(refinement_policy_path.resolve()),
        "policy_sha256": sha256_file(refinement_policy_path),
        "stop_reason": stop_reason,
        "rounds": rounds, "attempted_genomes": len(attempted), "pending_available_members": pending_count,
        "additional_reference_attempts": max(0, len(attempted) - targets_result["selected_genomes"]),
        "membership": membership, "all_attempted_targets": list(targets.values()),
        "interpretation": "Stability applies to sampled references only; it does not prove no closer unexamined genome exists.",
        "allocation_rule": "Supported cluster sets share a focused lane; reserved alternative slots rotate across selected clusters, with one configurable unresolved-cluster slot. Accession order is deterministic sampling, not similarity ranking.",
        "member_coverage_by_cluster": final_member_coverage,
        "remaining_uncertainty": {
            "unresolved_selected_clusters": unresolved_clusters,
            "leading_clusters_below_focused_minimum": focused_below_min,
            "unexamined_additional_pool": sum(row["unexamined_additional"] for row in final_member_coverage) if final_member_coverage is not None else None,
            "budget_or_round_limit_prevented_required_exploration": bool(resource_stop and (final_member_coverage is None or
                unresolved_clusters or focused_below_min or any(row["unexamined_additional"] for row in final_member_coverage))),
            "no_global_nearest_claim": True,
        },
    }
    if expand and (final_member_coverage is None or any(row["unexamined_additional"] for row in final_member_coverage)):
        warnings.append(f"Expansion stopped at {stop_reason}; the candidate neighborhood remains incompletely examined.")
    if membership and any(membership.get(field) for field in ("targets_without_metadata", "targets_without_assembly", "clusters_without_members")):
        warnings.append("The pinned membership snapshot contains isolates without usable assemblies or requested clusters without retrievable members.")
    if unavailable:
        warnings.append(f"{len(unavailable)} selected reference genomes could not be downloaded; the SNP comparison does not cover the entire selected set.")
    interpretation["warnings"] = list(dict.fromkeys(interpretation.get("warnings", []) + warnings))
    interpretation["selection"] = targets_result
    interpretation["unavailable"] = sorted(unavailable)
    interpretation["expansion"] = expansion
    qualifying = {row["sample"] for row in interpretation.get("ranked", [])}
    excluded = {row["sample"] for row in interpretation.get("excluded_comparisons", [])}
    returned_accessions = {row["accession"] for row in targets_result["decisions"]}
    member_accessions = {row["accession"] for row in membership["members"] if row["accession"] != query_accession} if membership and membership["status"] == "AVAILABLE" else set()
    cluster_coverage = []
    for cluster in dict.fromkeys(row["cluster"] for row in targets_result["decisions"]):
        returned = {row["accession"] for row in targets_result["decisions"] if row["cluster"] == cluster and row["reason"] != "query_self_excluded"}
        all_cluster_references = {accession for accession, row in targets.items() if row["cluster"] == cluster}
        cluster_members = {row["accession"] for row in membership["members"] if row["cluster"] == cluster and row["accession"] != query_accession} if membership and membership["status"] == "AVAILABLE" else set()
        cluster_coverage.append({
            "cluster": cluster, "returned": len(returned),
            "selected": len(returned & {row["accession"] for row in targets_result["targets"]}),
            "download_attempted": len(all_cluster_references & attempted),
            "download_verified": len(all_cluster_references & (set(genomes) - {"QUERY"})),
            "qualifying_query_comparisons": len(all_cluster_references & qualifying),
            "excluded_query_comparisons": len(all_cluster_references & excluded),
            "returned_unexamined": len(returned - attempted),
            "members_discovered": len(cluster_members) if membership and membership["status"] == "AVAILABLE" else None,
            "members_unexamined": len(cluster_members - attempted) if membership and membership["status"] == "AVAILABLE" else None,
        })
    interpretation["coverage"] = {
        "mashpit_returned_unique": targets_result["returned_genomes"],
        "returned_selected": targets_result["selected_genomes"],
        "returned_unexamined": sorted(returned_accessions - attempted - ({query_accession} if query_accession else set())),
        "reference_downloads_attempted": sorted(attempted),
        "reference_downloads_verified": sorted(set(genomes) - {"QUERY"}),
        "reference_downloads_unavailable": sorted(unavailable),
        "query_comparisons_qualifying": sorted(qualifying),
        "query_comparisons_excluded": sorted(excluded),
        "downloaded_without_recorded_query_comparison": sorted((set(genomes) - {"QUERY"}) - qualifying - excluded),
        "members_discovered": len(member_accessions) if membership and membership["status"] == "AVAILABLE" else None,
        "members_unexamined": sorted(member_accessions - attempted) if membership and membership["status"] == "AVAILABLE" else None,
        "member_coverage_unknown": not (membership and membership["status"] == "AVAILABLE"),
        "query_self_excluded_from_members": bool(query_accession and membership and membership["status"] == "AVAILABLE" and any(row["accession"] == query_accession for row in membership["members"])),
        "clusters": cluster_coverage,
        "omission_interpretation": "Unattempted or unavailable references are resource/availability gaps, not evidence of SNP dissimilarity.",
    }
    interpretation["search_scope"] = {
        "scope": "examined_reference_pool_only",
        "retrieval_limit_reached": targets_result.get("retrieval_limit_reached", False),
        "returned_unexamined_count": len(interpretation["coverage"]["returned_unexamined"]),
        "unavailable_count": len(unavailable),
        "global_nearest_established": False,
        "interpretation": "Expansion explores selected returned clusters. It cannot establish the absence of a closer cluster outside the retrieved pool.",
    }
    interpretation["query_input"] = {"type": "paired_reads" if query_reads else "assembly", "paths": query_reads if query_reads else [str(assembly)]}
    interpretation["file_checksums"] = [
        {"sample": sample, "path": path, "sha256": sha256_file(Path(path))}
        for sample, value in genomes.items()
        for path in ([value] if isinstance(value, str) else value)
        if Path(path).is_file()
    ]
    interpretation["backend"] = backend
    interpretation["comparison_runs"] = comparison_runs
    interpretation["tree_image"] = ({"status":"SKIPPED","reason":"Target-to-candidate comparisons do not define a phylogenetic tree."}
                                    if backend != "ska" else render_from_interpretation(interpretation, output_dir / "tree.png"))
    if backend != "ska":
        from plot_assembly_comparison import render
        interpretation["ranking_image"] = render(interpretation,output_dir/"candidate_distances.png")
    write_json(output_dir / "expansion.json", expansion)
    write_json(output_dir / "interpretation.json", interpretation)
    return {"snp_resolution": interpretation, "commands": commands}


def screen(
    inputs: list[Path], output_dir: Path, database_root: Path,
    organism: str | None = None, snp_resolve: bool = False, snp_expand: bool = False,
    snp_selection_mode: str | None = None, query_accession: str | None = None,
    refinement_policy_path: Path | None = None,
    snp_backend: str | None = None, assembly_manifest: Path | None = None,
    comparison_cache: Path | None = None, assembly_policy_path: Path | None = None,
    comparison_workers: int | None = None,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=False)
    commands: list[list[str]] = []
    assembly: Path | None = None
    query_reads: list[str] | None = None
    snp_resolve = snp_resolve or snp_expand or snp_selection_mode is not None
    database_metadata = None
    result: dict[str, Any] = {
        "schema_version": "1.0.0",
        "sample": sample_name(inputs[0]),
        "status": "FAILED",
        "warnings": [],
        "requested_organism": organism,
    }
    if query_accession is None and re.fullmatch(r"GC[AF]_\d+\.\d+", result["sample"]):
        query_accession = result["sample"]
    result["query_accession_excluded"] = query_accession
    try:
        detected = inspect(inputs)
        result["input_type"] = detected["input_type"]
        write_json(output_dir / "input.json", detected)
        if detected["input_type"] == "assembly":
            assembly = Path(detected["assembly_path"])
        else:
            assembly_run = run_workflow(
                Path(detected["r1_path"]), Path(detected["r2_path"]), output_dir / "assembly"
            )
            commands.extend(assembly_run["commands"])
            result["read_qc"] = assembly_run["read_qc"]
            result["warnings"].extend(
                ["Read QC is WARN; downstream interpretation may be affected."]
                if assembly_run["read_qc"]["status"] == "WARN" else []
            )
            assembly = Path(assembly_run["assembly_path"])
            query_reads = assembly_run["read_paths"]
        basic_qc = assess(assembly)
        write_json(output_dir / "assembly_qc.initial.json", basic_qc)
        if basic_qc["status"] == "FAIL":
            result["assembly_qc"] = basic_qc
            raise WorkflowError("Assembly QC failed: " + "; ".join(basic_qc["failures"]))
        if result.get("read_qc"):
            total_bases = result["read_qc"].get("total_bases")
            if total_bases:
                coverage = total_bases / basic_qc["metrics"]["total_length"]
                result["read_qc"]["estimated_coverage"] = coverage
                if coverage < 20.0:
                    result["read_qc"]["status"] = "FAIL"
                    raise WorkflowError("Estimated read coverage is below the fixed 20× minimum.")
        if organism:
            routing = route_for_organism(organism)
        else:
            routing = classify(assembly, output_dir / "mlst")
            if routing.get("command"):
                commands.append(routing["command"])
        write_json(output_dir / "routing.json", routing)
        result["routing"] = routing
        if routing["status"] == "UNSUPPORTED":
            raise WorkflowError(
                "The organism does not match one of the currently supported Mashpit databases."
            )
        if routing["status"] != "SUPPORTED":
            raise WorkflowError(
                "Local mlst could not confidently select a supported organism database; "
                "specify --organism only when the organism is known independently."
            )
        result["organism"] = routing["organism"]
        result["warnings"].extend(routing.get("warnings", []))
        routed_qc = assess(assembly, routing["organism_key"])
        write_json(output_dir / "assembly_qc.json", routed_qc)
        result["assembly_qc"] = routed_qc
        if routed_qc["status"] == "FAIL":
            raise WorkflowError("Assembly QC failed: " + "; ".join(routed_qc["failures"]))
        if routed_qc["status"] == "WARN":
            result["warnings"].extend(routed_qc["warnings"])
        mashpit_run = run_mashpit(
            assembly,
            database_root / routing["database_name"],
            routing["database_name"],
            output_dir / "mashpit",
        )
        commands.append(mashpit_run["command"])
        database_metadata = mashpit_run["database"]
        result["database"] = database_metadata
        candidate_file = locate_candidate_file(Path(mashpit_run["output_directory"]))
        mashpit_profile = load_json(CONFIG_DIR / "workflow.json")["mashpit"]
        parsed = interpret(load_candidates(candidate_file), threshold=mashpit_profile["threshold"])
        parsed["source_file"] = str(candidate_file)
        mashpit_tree_source = locate_tree_image(Path(mashpit_run["output_directory"]))
        if mashpit_tree_source:
            shutil.copy2(mashpit_tree_source, output_dir / "mashpit_tree.png")
            parsed["tree_image"] = {"status": "PASS", "path": str(output_dir / "mashpit_tree.png")}
        else:
            parsed["tree_image"] = {
                "status": "UNAVAILABLE",
                "reason": (
                    "Mashpit did not generate a tree for this query (top hit below --threshold, "
                    "or fewer than two candidates qualified)."
                ),
            }
        write_json(output_dir / "mashpit_result.json", parsed)
        result["mashpit_result"] = parsed
        result["warnings"].extend(parsed.get("warnings", []))
        diagnostics = run_diagnostics(
            Path(mashpit_run["output_directory"]), database_metadata,
            load_json(CONFIG_DIR / "workflow.json"),
            load_json(CONFIG_DIR / "snp-resolution-policy.json"),
            output_dir / "similarity_distribution",
            snp_selection_mode, query_accession,
        )
        result["similarity_distribution"] = diagnostics
        result["warnings"].extend(diagnostics.get("warnings", []))
        if snp_resolve and parsed.get("below_threshold"):
            result["warnings"].append(
                "SNP resolution skipped: Mashpit's top hit is below its own query threshold, "
                "so there is no real candidate to resolve."
            )
        elif snp_resolve and parsed.get("best_candidate"):
            try:
                snp_outcome = run_snp_resolution(
                    Path(mashpit_run["output_directory"]), assembly,
                    parsed["best_candidate"]["cluster"], output_dir / "snp_resolution",
                    database_metadata, snp_expand, query_reads, snp_selection_mode, query_accession,
                    refinement_policy_path, snp_backend, assembly_manifest, comparison_cache, assembly_policy_path, comparison_workers,
                )
                result["snp_resolution"] = snp_outcome["snp_resolution"]
                commands.extend(snp_outcome["commands"])
                result["warnings"].extend(snp_outcome["snp_resolution"].get("warnings", []))
                if snp_outcome["snp_resolution"].get("reason"):
                    result["warnings"].append(snp_outcome["snp_resolution"]["reason"])
            except (WorkflowError, OSError, ValueError) as error:
                result["snp_resolution"] = {"status": "ERROR", "error": str(error)}
                result["warnings"].append(f"SNP resolution failed: {error}")
        result["warnings"] = list(dict.fromkeys(result["warnings"]))
        result["status"] = "COMPLETED_WITH_WARNINGS" if result["warnings"] else "COMPLETED"
    except (WorkflowError, OSError, ValueError) as error:
        result["status"] = "STOPPED"
        result["stop_reason"] = str(error)
    finally:
        result["user_summary"] = summary_text(result)
        write_json(output_dir / "result.json", result)
        (output_dir / "report.md").write_text(generate_report(result), encoding="utf-8")
        provenance = collect(inputs, assembly, database_metadata, commands, result, refinement_policy_path)
        write_json(output_dir / "provenance.json", provenance)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--database-root")
    parser.add_argument(
        "--organism",
        choices=("salmonella", "ecoli_shigella", "listeria", "campylobacter", "cronobacter"),
        help="Known organism/database key. If omitted, local mlst auto-detects a PubMLST scheme.",
    )
    parser.add_argument(
        "--snp-resolve",
        action="store_true",
        help=(
            "After a Mashpit candidate is found, download representative genomes from NCBI "
            "and compare the target assembly with candidates using minimap2 (SKA for reads). This may reach "
            "out to NCBI, unlike the rest of the screen, which stays fully local."
        ),
    )
    parser.add_argument("--snp-expand", action="store_true", help="Enable SNP resolution plus bounded exploration of members from the pinned NCBI cluster release.")
    parser.add_argument("--snp-selection-mode", choices=("adaptive", "all_returned"), help="Representative selection mode; all_returned requires budget for every unique returned representative.")
    parser.add_argument("--query-accession", help="Known assembly accession of this query, excluded from SNP references.")
    parser.add_argument("--refinement-policy", help="Explicit versioned refinement policy JSON; defaults to config/refinement-policy.json.")
    parser.add_argument("--snp-backend",choices=("minimap2","mummer","ska"),help="Assembly default: minimap2 target-to-candidate; reads default: legacy SKA.")
    parser.add_argument("--assembly-manifest",help="Local accession/path/sha256 TSV; download only missing accessions.")
    parser.add_argument("--comparison-cache",help="Reusable content-addressed assembly alignment cache.")
    parser.add_argument("--assembly-comparison-policy",help="Explicit assembly comparison policy JSON.")
    parser.add_argument("--comparison-workers",type=int,help="Concurrent single-thread target-to-candidate alignments (default 4).")
    args = parser.parse_args()
    output_dir = Path(args.output).expanduser().resolve()
    if output_dir.exists():
        print(f"Refusing to overwrite existing output directory: {output_dir}")
        return 2
    if args.refinement_policy and not Path(args.refinement_policy).expanduser().is_file():
        print(f"Refinement policy file not found: {args.refinement_policy}")
        return 2
    database_root = effective_database_root(args.database_root)
    result = screen(
        [Path(value).expanduser().resolve() for value in args.inputs],
        output_dir,
        database_root,
        args.organism,
        args.snp_resolve,
        args.snp_expand,
        args.snp_selection_mode,
        args.query_accession,
        Path(args.refinement_policy).expanduser().resolve() if args.refinement_policy else None,
        args.snp_backend,
        Path(args.assembly_manifest).expanduser().resolve() if args.assembly_manifest else None,
        Path(args.comparison_cache).expanduser().resolve() if args.comparison_cache else None,
        Path(args.assembly_comparison_policy).expanduser().resolve() if args.assembly_comparison_policy else None,
        args.comparison_workers,
    )
    print(result["user_summary"])
    return 0 if result["status"].startswith("COMPLETED") else 2


if __name__ == "__main__":
    raise SystemExit(main())
