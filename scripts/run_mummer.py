#!/usr/bin/env python3
"""Incremental target-to-candidate MUMmer4 comparisons; never builds a matrix."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import fcntl
import hashlib
import json
import math
import os
import re
from pathlib import Path
import shutil
import subprocess
import time
from typing import Any
from common import CONFIG_DIR, WorkflowError, load_json, require_executable, sha256_file, write_json
from shared_target_regions import save_evidence


def key_for(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths = {}
    name = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields or fields[0] in lengths:
                    raise WorkflowError(f"Empty or duplicate FASTA identifier: {path}")
                name = fields[0]
                lengths[name] = 0
            elif name is None or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
                raise WorkflowError(f"Invalid assembly FASTA: {path}")
            else:
                lengths[name] += len(line)
    if not lengths or not all(lengths.values()):
        raise WorkflowError(f"Empty assembly or contig: {path}")
    return lengths


def union_length(intervals: list[tuple[int, int]]) -> int:
    total = 0
    end = 0
    for lo, hi in sorted(intervals):
        if not 1 <= lo <= hi:
            raise WorkflowError("Invalid alignment interval")
        total += max(0, hi - max(end, lo - 1))
        end = max(end, hi)
    return total


def parse_pair(directory: Path, target_lengths: dict, candidate_lengths: dict) -> dict:
    target_intervals = {name: [] for name in target_lengths}
    candidate_intervals = {name: [] for name in candidate_lengths}
    identities = []
    with (directory / "coords.tsv").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 13:
                raise WorkflowError("Unexpected show-coords -rclTH output")
            a,b,c,d = map(int, fields[:4])
            t,q = fields[11:13]
            if t not in target_lengths or q not in candidate_lengths:
                raise WorkflowError("Alignment contains unknown contig")
            if int(fields[7]) != target_lengths[t] or int(fields[8]) != candidate_lengths[q]:
                raise WorkflowError("Alignment contig length does not match input")
            if min(a,b,c,d) < 1 or max(a,b) > target_lengths[t] or max(c,d) > candidate_lengths[q]:
                raise WorkflowError("Alignment coordinate outside input")
            identity = float(fields[6])
            if not math.isfinite(identity) or not 0 <= identity <= 100:
                raise WorkflowError("Invalid alignment identity")
            target_intervals[t].append((min(a,b),max(a,b)))
            candidate_intervals[q].append((min(c,d),max(c,d)))
            identities.append((identity, abs(b-a)+1))
    snps, indels, ambiguous = set(), 0, 0
    with (directory / "snps.tsv").open() as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) != 12:
                raise WorkflowError("Unexpected show-snps -rlTHC output")
            pos, rbase, qbase, qpos = int(f[0]), f[1].upper(), f[2].upper(), int(f[3])
            t,q = f[10:12]
            if t not in target_lengths or q not in candidate_lengths:
                raise WorkflowError("SNP contains unknown contig")
            if rbase in "ACGT" and qbase in "ACGT" and len(rbase)==len(qbase)==1:
                if rbase == qbase or not 1 <= pos <= target_lengths[t] or not 1 <= qpos <= candidate_lengths[q]:
                    raise WorkflowError("Invalid SNP coordinate or substitution")
                site = (t,pos)
                if site in snps:
                    raise WorkflowError("Duplicate target SNP position in unique mapping")
                snps.add(site)
            elif rbase == "." or qbase == ".":
                indels += 1
            else:
                ambiguous += 1
    ta = sum(union_length(v) for v in target_intervals.values())
    ca = sum(union_length(v) for v in candidate_intervals.values())
    if snps and not ta:
        raise WorkflowError("SNPs without aligned bases")
    return {"snp_positions":sorted(snps), "target_unique_alignment_intervals":target_intervals, "snp_distance":len(snps), "indel_bases":indels, "ambiguous_difference_rows":ambiguous,
            "target_aligned_bases":ta, "candidate_aligned_bases":ca,
            "target_bases":sum(target_lengths.values()), "candidate_bases":sum(candidate_lengths.values()),
            "target_aligned_fraction":ta/sum(target_lengths.values()),
            "candidate_aligned_fraction":ca/sum(candidate_lengths.values()),
            "snps_per_target_aligned_mb":len(snps)*1e6/ta if ta else None,
            "alignment_identity_percent":sum(i*n for i,n in identities)/sum(n for _,n in identities) if identities else None,
            "rate_denominator":"union of target intervals in one-to-one alignments; not a callable-site count"}


def tool_identity(policy: dict) -> dict:
    paths = {name:require_executable(name) for name in ("nucmer","delta-filter","show-coords","show-snps")}
    version = subprocess.run([paths["nucmer"],"--version"], capture_output=True,text=True,check=True,timeout=30).stdout.strip()
    if version != policy["mummer_version"]:
        raise WorkflowError(f"Expected MUMmer {policy['mummer_version']}; found {version}")
    return {"version":version,"executables":{k:{"path":v,"sha256":sha256_file(Path(v))} for k,v in paths.items()}}


def command_run(command: list[str], directory: Path, output: str, timeout: int) -> None:
    env = {**os.environ,"OMP_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1"}
    with (directory/output).open("w") as stdout, (directory/(output+".stderr.log")).open("w") as stderr:
        completed = subprocess.run(command,cwd=directory,stdout=stdout,stderr=stderr,env=env,timeout=timeout)
    if completed.returncode:
        raise WorkflowError(f"{Path(command[0]).name} failed; inspect {directory/(output+'.stderr.log')}")


def compare_pair(target: Path, candidate: Path, cache: Path, policy: dict, tools: dict,
                 target_hash: str, target_lengths: dict) -> dict:
    candidate_hash = sha256_file(candidate)
    identity = {"schema":"mummer-one-to-one-v1","target_sha256":target_hash,"candidate_sha256":candidate_hash,
                "alignment":policy["alignment"],"tools":tools}
    cache_key = key_for(identity)
    directory = cache/cache_key
    directory.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    commands = []
    parser_hash = sha256_file(Path(__file__))
    with (directory/"pair.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        complete = directory/"alignment.json"
        hit = False
        if complete.exists():
            record = load_json(complete)
            hit = record.get("identity") == identity and all((directory/n).is_file() and sha256_file(directory/n)==h for n,h in record.get("artifacts",{}).items()) and set(record.get("artifacts",{})) == {"pair.delta","one.delta","coords.tsv","snps.tsv"}
            if hit:
                commands = record["commands"]
        if not hit:
            complete.unlink(missing_ok=True)
            # These are private, fixed-name inputs in a locked content-addressed directory.
            # Copy instead of linking so cached raw alignments remain interpretable after input relocation.
            inputs = cache/"inputs"
            inputs.mkdir(exist_ok=True)
            for name,source,digest in (("target.fna",target,target_hash),("candidate.fna",candidate,candidate_hash)):
                saved=inputs/(digest+".fna")
                with (inputs/(digest+".lock")).open("a") as input_lock:
                    fcntl.flock(input_lock,fcntl.LOCK_EX)
                    if not saved.exists():
                        temporary=saved.with_suffix(".tmp")
                        shutil.copyfile(source,temporary)
                        if sha256_file(temporary)!=digest: raise WorkflowError("Input changed while caching")
                        temporary.replace(saved)
                link=directory/name
                if link.exists() or link.is_symlink(): link.unlink()
                link.symlink_to(saved)
            a = policy["alignment"]
            exe = {n:v["path"] for n,v in tools["executables"].items()}
            commands = [
                [exe["nucmer"],"--maxmatch","--threads","1","--minmatch",str(a["min_match"]),"--mincluster",str(a["min_cluster"]),"--breaklen",str(a["break_length"]),"--maxgap",str(a["max_gap"]),"--prefix","pair","target.fna","candidate.fna"],
                [exe["delta-filter"],"-1","-i",str(a["min_alignment_identity"]),"-l",str(a["min_alignment_length"]),"pair.delta"],
                [exe["show-coords"],"-rclTH","one.delta"],
                [exe["show-snps"],"-rlTHC","one.delta"],
            ]
            for command,output in zip(commands,("nucmer.stdout.log","one.delta","coords.tsv","snps.tsv")):
                command_run(command,directory,output,policy["execution"]["command_timeout_seconds"])
            # Parse before accepting even an empty SNP table as a successful comparison.
            metrics = parse_pair(directory,target_lengths,fasta_lengths(candidate))
            write_json(complete,{"identity":identity,"commands":commands,"parser_sha256":parser_hash,"metrics":metrics,
                "artifacts":{name:sha256_file(directory/name) for name in ("pair.delta","one.delta","coords.tsv","snps.tsv")}})
        if hit:
            metrics = (record["metrics"] if record.get("parser_sha256")==parser_hash
                       else parse_pair(directory,target_lengths,fasta_lengths(candidate)))
        evidence=save_evidence(directory/"site_evidence.json",metrics,target_hash,candidate_hash,target_lengths)
    summary = {k: v for k, v in metrics.items() if k not in {"snp_positions","target_unique_alignment_intervals"}}
    return {**summary,"site_evidence":evidence,"status":"COMPARED","cache_key":cache_key,"cache_hit":hit,
            "target_sha256":target_hash,"candidate_sha256":candidate_hash,"artifact_directory":str(directory),
            "commands":commands,"elapsed_seconds":time.monotonic()-started}


def validate_policy(policy: dict) -> None:
    a=policy["alignment"]; e=policy["execution"]; c=policy["comparability"]
    if a["anchor_mode"] != "maxmatch":
        raise WorkflowError("Unsupported MUMmer anchor mode")
    for key in ("min_match","min_cluster","break_length","max_gap"):
        if type(a[key]) is not int or a[key] <= 0: raise WorkflowError(f"Invalid alignment {key}")
    if type(a["min_alignment_length"]) is not int or a["min_alignment_length"] < 0: raise WorkflowError("Invalid minimum alignment length")
    if not math.isfinite(a["min_alignment_identity"]) or not 0<=a["min_alignment_identity"]<=100: raise WorkflowError("Invalid minimum identity")
    for key in ("min_target_aligned_fraction","min_candidate_aligned_fraction"):
        if not math.isfinite(c[key]) or not 0<c[key]<=1: raise WorkflowError("Invalid coverage criterion")
    for key in ("workers","command_timeout_seconds"):
        if type(e[key]) is not int or e[key]<=0: raise WorkflowError(f"Invalid execution {key}")


def run_mummer(genomes: dict[str,str], output_dir: Path, cache_dir: Path | None = None,
               policy_path: Path | None = None, workers: int | None = None) -> dict:
    if "QUERY" not in genomes or len(genomes)<2 or not all(isinstance(v,str) for v in genomes.values()):
        raise WorkflowError("MUMmer needs a target assembly named QUERY and at least one candidate assembly")
    policy_path = policy_path or CONFIG_DIR/"mummer-comparison-policy.json"
    policy = load_json(policy_path)
    if workers is not None: policy["execution"]["workers"] = workers
    validate_policy(policy)
    output_dir.mkdir(parents=True,exist_ok=True)
    cache = (cache_dir or output_dir/"cache").resolve(); cache.mkdir(parents=True,exist_ok=True)
    target = Path(genomes["QUERY"]).resolve(); target_hash = sha256_file(target)
    lengths = fasta_lengths(target)
    tools = tool_identity(policy)
    started = time.monotonic()
    def run(item):
        sample,value = item
        try:
            result = compare_pair(target,Path(value).resolve(),cache,policy,tools,target_hash,lengths)
        except (WorkflowError,OSError,ValueError,subprocess.SubprocessError) as error:
            result = {"status":"ERROR","error":str(error),"cache_hit":False,"snp_distance":None}
        return {"sample":sample,**result}
    with ThreadPoolExecutor(max_workers=policy["execution"]["workers"]) as pool:
        rows = list(pool.map(run,sorted((k,v) for k,v in genomes.items() if k!="QUERY")))
    result = {"backend":"mummer","status":"PASS" if all(x["status"]=="COMPARED" for x in rows) else "PARTIAL",
              "comparisons":rows,"commands":[c for row in rows if not row.get("cache_hit") for c in row.get("commands",[])],
              "policy":policy,"policy_sha256":sha256_file(policy_path),"tools":tools,
              "implementation_sha256":sha256_file(Path(__file__)),"cache_directory":str(cache),
              "alignment_jobs":sum(x["status"]=="COMPARED" and not x["cache_hit"] for x in rows),
              "cache_hits":sum(x.get("cache_hit",False) for x in rows),
              "elapsed_seconds":time.monotonic()-started,"comparison_scope":"target_to_candidate"}
    write_json(output_dir/"run_mummer.json",result)
    return result


def main() -> int:
    p=argparse.ArgumentParser()
    p.add_argument("--target",required=True); p.add_argument("--candidate-manifest",required=True)
    p.add_argument("--output",required=True); p.add_argument("--cache-dir"); p.add_argument("--policy"); p.add_argument("--workers",type=int)
    a=p.parse_args()
    from fetch_reference_genomes import read_local_manifest
    refs=read_local_manifest(Path(a.candidate_manifest))
    genomes={"QUERY":str(Path(a.target).resolve())}
    for accession,row in refs.items():
        path=Path(row["path"])
        if sha256_file(path)!=row["sha256"]: raise WorkflowError(f"Assembly checksum mismatch: {accession}")
        if accession=="QUERY": raise WorkflowError("QUERY is reserved for the target")
        genomes[accession]=str(path)
    result=run_mummer(genomes,Path(a.output).resolve(),Path(a.cache_dir).resolve() if a.cache_dir else None,Path(a.policy).resolve() if a.policy else None,a.workers)
    return 0 if result["status"]=="PASS" else 2

if __name__=="__main__":
    raise SystemExit(main())
