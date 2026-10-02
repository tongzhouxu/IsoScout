#!/usr/bin/env python3
"""Deterministic end-to-end entrypoint for isolate screening."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Any

from analyze_similarity_distribution import run_diagnostics
from expand_cluster_members import discover_members
from refinement_expansion import feedback, next_batch
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
    if snp and snp.get("ranked"):
        lines.append(
            f"SKA2 comparison: minimum observed distance {snp['nearest_snp_distance']:.2f} SNPs; "
            f"{len(snp.get('nearest_samples', []))} equally nearest reference(s). "
            "This is a comparison among qualifying references, not strain or outbreak confirmation."
        )
    if result.get("stop_reason"):
        lines.append(f"Analysis stopped: {result['stop_reason']}")
    warnings = result.get("warnings", [])
    if warnings:
        lines.append("Warnings: " + "; ".join(warnings))
    return "\n".join(lines)


def run_snp_resolution(
    mashpit_output_dir: Path, assembly: Path, best_cluster: str, output_dir: Path,
    database_metadata: dict[str, Any] | None = None, expand: bool = False,
    query_reads: list[str] | None = None,
) -> dict[str, Any]:
    policy = load_json(CONFIG_DIR / "snp-resolution-policy.json")
    workflow = load_json(CONFIG_DIR / "workflow.json")
    refinement_policy = load_json(CONFIG_DIR / "refinement-policy.json")
    expansion_policy = refinement_policy["expansion"]
    kmer_size = workflow["snp_resolution"]["kmer_size"]
    commands: list[list[str]] = []
    targets_result = select_targets(
        mashpit_output_dir, policy, database_metadata,
        workflow["mashpit"] if database_metadata is not None else None,
    )
    write_json(output_dir / "targets.json", targets_result)
    if targets_result["status"] != "SELECTED":
        return {"snp_resolution": targets_result, "commands": commands}
    if database_metadata is None:
        database_metadata = load_json(mashpit_output_dir / "mashpit_run.json")["database"]
    targets = {item["accession"]: item for item in targets_result["targets"]}
    batch = list(targets.values())
    genomes: dict[str, Any] = {"QUERY": query_reads if query_reads else str(assembly)}
    attempted: set[str] = set()
    unavailable: set[str] = set()
    rounds = []
    membership = None
    members: list[dict[str, Any]] = []
    previous = None
    stable_rounds = 0
    stop_reason = "representatives_only"
    warnings = list(targets_result.get("warnings", []))
    # The total cap includes attempted downloads, even when some fail.
    cap = expansion_policy["max_total_genomes"]
    if expand and len(batch) > cap:
        raise WorkflowError("Initial selection exceeds the configured expansion ceiling.")
    interpretation = {"status": "ERROR", "ranked": [], "warnings": []}
    last_compared_count = 0
    for round_index in range(expansion_policy["max_rounds"] if expand else 1):
        round_dir = output_dir / f"round_{round_index:02d}"
        accessions = [item["accession"] for item in batch]
        attempted.update(accessions)
        targets.update({item["accession"]: item for item in batch})
        write_json(round_dir / "requested.json", {"targets": batch})
        try:
            fetch_result = fetch_genomes(accessions, round_dir / "genomes", policy["download_attempts"], policy["download_retry_delay_seconds"])
            commands.extend(fetch_result["commands"])
            genomes.update(fetch_result["verified"])
            unavailable.update(fetch_result["unavailable"])
            if len(genomes) < 2:
                interpretation = {"status": "INSUFFICIENT_DATA", "ranked": [], "warnings": ["No usable reference genomes could be downloaded."]}
            else:
                ska_result = run_ska(genomes, round_dir / "ska", kmer_size)
                commands.extend(ska_result["commands"])
                interpretation = interpret_snp_resolution(ska_result["distances"], list(targets.values()), best_cluster)
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
        decision = feedback(previous, interpretation)
        if decision["action"] != "STOP" and previous is not None and len(genomes) - 1 <= last_compared_count:
            decision = {"action": "EXPAND", "reason": "no_new_comparisons", "stable": False}
        last_compared_count = len(genomes) - 1
        stable_rounds = stable_rounds + 1 if decision["stable"] else 0
        record = {
            "round": round_index, "requested_accessions": accessions,
            "downloaded_accessions": sorted(fetch_result["verified"]),
            "unavailable_accessions": fetch_result["unavailable"], "attempted_total": len(attempted),
            "compared_total": len(genomes) - 1, "feedback": decision,
            "nearest_samples": interpretation.get("nearest_samples", []),
            "nearest_snp_distance": interpretation.get("nearest_snp_distance"),
            "excluded_comparisons": interpretation.get("excluded_comparisons", []),
        }
        rounds.append(record)
        write_json(round_dir / "interpretation.json", interpretation)
        write_json(round_dir / "round.json", record)
        if not expand:
            break
        if decision["action"] == "STOP":
            stop_reason = decision["reason"]
            break
        if len(attempted) >= cap:
            stop_reason = "hard_resource_ceiling"
            break
        if membership is None:
            membership = discover_members(database_metadata, targets_result["relevant_clusters"], output_dir / "membership", expansion_policy)
            if membership["status"] != "AVAILABLE":
                stop_reason = "membership_unavailable"
                warnings.append("Cluster expansion unavailable: " + membership.get("error", "unknown error"))
                break
            # Include returned but unselected representatives as well as members.
            combined = {item["accession"]: item for item in membership["members"]}
            for item in targets_result["decisions"]:
                if item["cluster"] in targets_result["relevant_clusters"]:
                    if item["accession"] in combined and combined[item["accession"]]["cluster"] != item["cluster"]:
                        membership["status"] = "CONFLICT"
                        stop_reason = "membership_conflict"
                        warnings.append("Pinned membership and Mashpit disagree about a reference's cluster; expansion stopped and earlier comparisons are retained.")
                        write_json(output_dir / "membership" / "membership.json", membership)
                        break
                    combined[item["accession"]] = item
            if membership["status"] == "CONFLICT":
                break
            members = list(combined.values())
        pending = [item for item in members if item["accession"] not in attempted]
        if not pending:
            stop_reason = "available_member_pool_exhausted"
            break
        if stable_rounds >= expansion_policy["stable_rounds_to_stop"]:
            stop_reason = "stable_sampled_neighborhood"
            break
        if round_index + 1 >= expansion_policy["max_rounds"]:
            stop_reason = "round_limit"
            break
        priority = list(interpretation.get("nearest_clusters", []))
        priority.extend(cluster for cluster in targets_result["relevant_clusters"] if cluster not in priority)
        batch = next_batch(members, attempted, priority, expansion_policy["batch_size"], cap - len(attempted))
        previous = interpretation
    pending_count = sum(item["accession"] not in attempted for item in members) if membership and membership["status"] == "AVAILABLE" else None
    expansion = {
        "enabled": expand, "policy": expansion_policy, "stop_reason": stop_reason,
        "rounds": rounds, "attempted_genomes": len(attempted), "pending_available_members": pending_count,
        "membership": membership, "all_attempted_targets": list(targets.values()),
        "interpretation": "Stability applies to sampled references only; it does not prove no closer unexamined genome exists.",
    }
    if expand and (pending_count is None or pending_count > 0):
        warnings.append(f"Expansion stopped at {stop_reason}; the candidate neighborhood remains incompletely examined.")
    if membership and any(membership.get(field) for field in ("targets_without_metadata", "targets_without_assembly", "clusters_without_members")):
        warnings.append("The pinned membership snapshot contains isolates without usable assemblies or requested clusters without retrievable members.")
    if unavailable:
        warnings.append(f"{len(unavailable)} selected reference genomes could not be downloaded; the SNP comparison does not cover the entire selected set.")
    interpretation["warnings"] = list(dict.fromkeys(interpretation.get("warnings", []) + warnings))
    interpretation["selection"] = targets_result
    interpretation["unavailable"] = sorted(unavailable)
    interpretation["expansion"] = expansion
    interpretation["query_input"] = {"type": "paired_reads" if query_reads else "assembly", "paths": query_reads if query_reads else [str(assembly)]}
    interpretation["file_checksums"] = [
        {"sample": sample, "path": path, "sha256": sha256_file(Path(path))}
        for sample, value in genomes.items()
        for path in ([value] if isinstance(value, str) else value)
        if Path(path).is_file()
    ]
    interpretation["tree_image"] = render_from_interpretation(interpretation, output_dir / "tree.png")
    write_json(output_dir / "expansion.json", expansion)
    write_json(output_dir / "interpretation.json", interpretation)
    return {"snp_resolution": interpretation, "commands": commands}


def screen(
    inputs: list[Path], output_dir: Path, database_root: Path,
    organism: str | None = None, snp_resolve: bool = False, snp_expand: bool = False,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=False)
    commands: list[list[str]] = []
    assembly: Path | None = None
    query_reads: list[str] | None = None
    snp_resolve = snp_resolve or snp_expand
    database_metadata = None
    result: dict[str, Any] = {
        "schema_version": "1.0.0",
        "sample": sample_name(inputs[0]),
        "status": "FAILED",
        "warnings": [],
        "requested_organism": organism,
    }
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
                    database_metadata, snp_expand, query_reads,
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
        provenance = collect(inputs, assembly, database_metadata, commands, result)
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
            "and run ska2 to compute query-to-representative SNP distances. Opt-in: this reaches "
            "out to NCBI, unlike the rest of the screen, which stays fully local."
        ),
    )
    parser.add_argument("--snp-expand", action="store_true", help="Enable SNP resolution plus bounded exploration of members from the pinned NCBI cluster release.")
    args = parser.parse_args()
    output_dir = Path(args.output).expanduser().resolve()
    if output_dir.exists():
        print(f"Refusing to overwrite existing output directory: {output_dir}")
        return 2
    database_root = effective_database_root(args.database_root)
    result = screen(
        [Path(value).expanduser().resolve() for value in args.inputs],
        output_dir,
        database_root,
        args.organism,
        args.snp_resolve,
        args.snp_expand,
    )
    print(result["user_summary"])
    return 0 if result["status"].startswith("COMPLETED") else 2


if __name__ == "__main__":
    raise SystemExit(main())
