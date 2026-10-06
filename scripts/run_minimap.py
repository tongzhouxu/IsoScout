#!/usr/bin/env python3
"""Cached minimap2 assembly alignment plus upstream paftools target-only calls."""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import fcntl
from bisect import bisect_right
import math
import os
from pathlib import Path
import re
import subprocess
import time
from common import CONFIG_DIR,WorkflowError,load_json,require_executable,sha256_file,write_json
from run_mummer import key_for,fasta_lengths,union_length,command_run
from shared_target_regions import save_evidence


def prepare_paf(path: Path, output: Path, tl: dict, ql: dict, policy: dict) -> dict:
    rows=[]; qi={k:[] for k in ql}; all_rows=0
    with path.open() as handle:
        for line in handle:
            f=line.rstrip("\n").split("\t");all_rows+=1
            if len(f)<12:raise WorkflowError("Malformed PAF output")
            q,t=f[0],f[5]
            if q not in ql or t not in tl or int(f[1])!=ql[q] or int(f[6])!=tl[t]:raise WorkflowError("PAF input identity/length mismatch")
            qs,qe,ts,te=map(int,(f[2],f[3],f[7],f[8]))
            if not 0<=qs<qe<=ql[q] or not 0<=ts<te<=tl[t] or f[4] not in ('+','-'):raise WorkflowError("PAF coordinate outside input")
            tags={x.split(':',2)[0]:x.split(':',2)[2] for x in f[12:] if len(x.split(':',2))==3}
            if tags.get('tp') in ('S','i') or ('s1' in tags and 's2' not in tags):continue
            if int(f[10])<policy['min_alignment_length'] or int(f[11])<policy['min_mapping_quality']:continue
            if 'cs' not in tags:raise WorkflowError("PAF alignment lacks required cs tag")
            # Short cs permits exact endpoint validation without loading long sequence tags.
            x,y=ts,qs
            parts=re.findall(r'(?::[0-9]+|\*[A-Za-z][A-Za-z]|[+-][A-Za-z]+)',tags['cs'])
            if ''.join(parts)!=tags['cs']:raise WorkflowError("Unsupported or malformed cs tag")
            for part in parts:
                if part[0]==':':x+=int(part[1:]);y+=int(part[1:])
                elif part[0]=='*':x+=1;y+=1
                elif part[0]=='-':x+=len(part)-1
                else:y+=len(part)-1
            if x!=te or y!=qe:raise WorkflowError("cs tag does not span PAF alignment")
            rows.append((t,ts,line));qi[q].append((qs+1,qe))
    output.write_text(''.join(x[2] for x in sorted(rows,key=lambda x:(x[0],x[1]))))
    return {'candidate_aligned_bases':sum(union_length(v) for v in qi.values()),'paf_rows':all_rows,'qualifying_paf_rows':len(rows)}


def parse_calls(path: Path, tl: dict) -> dict:
    intervals={k:[] for k in tl};variants=[];ambiguous=0
    for line in path.read_text().splitlines():
        f=line.split('\t')
        if f[0]=='R' and len(f)==4:
            name,start,end=f[1],int(f[2]),int(f[3])
            if name not in tl or not 0<=start<end<=tl[name]:raise WorkflowError("Invalid paftools region")
            intervals[name].append((start+1,end))
        elif f[0]=='V' and len(f)==12:
            name,start,end=f[1],int(f[2]),int(f[3])
            if name not in tl or not 0<=start<=end<=tl[name]:raise WorkflowError("Invalid paftools variant")
            if int(f[4])!=1:continue
            variants.append((name,start,end,f[6].upper(),f[7].upper()))
        else:raise WorkflowError("Unexpected paftools call output")
    # Merge coverage once; per-variant interval scans are costly on fragmented inputs.
    covered = {}
    for name, spans in intervals.items():
        merged = []
        for start, end in sorted(spans):
            if merged and start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        covered[name] = ([a for a, _ in merged], [b for _, b in merged])
    positions=set();indels=0
    for name,start,end,ref,alt in variants:
        if len(ref)==len(alt)==1 and ref in 'ACGT' and alt in 'ACGT':
            if end!=start+1 or ref==alt:raise WorkflowError("Invalid substitution")
            starts, ends = covered[name]
            region = bisect_right(starts, start + 1) - 1
            if region < 0 or start + 1 > ends[region]:continue
            site=(name,start+1)
            if site in positions:raise WorkflowError("Duplicate target substitution")
            positions.add(site)
        elif ref=='-' or alt=='-':indels+=len(alt) if ref=='-' else len(ref)
        else:ambiguous+=1
    aligned=sum(union_length(v) for v in intervals.values())
    return {'snp_distance':len(positions),'indel_bases':indels,'ambiguous_difference_rows':ambiguous,
            'target_aligned_bases':aligned,
            'target_unique_alignment_intervals':intervals, 'snp_positions':sorted(positions)}


def compare_pair(target,candidate,cache,policy,tools,target_hash,tl):
    started=time.monotonic();ch=sha256_file(candidate)
    ident={'schema':'minimap2-assembly-v1','target_sha256':target_hash,'candidate_sha256':ch,'alignment':policy['alignment'],'tools':tools}
    key=key_for(ident);directory=cache/key;directory.mkdir(parents=True,exist_ok=True)
    parser_hash=sha256_file(Path(__file__))
    call_key=key_for({'variant_calling':policy['variant_calling'],'parser':parser_hash})
    call_dir=directory/call_key;call_dir.mkdir(exist_ok=True)
    with (directory/'pair.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        done=directory/'alignment.json';paf=directory/'alignment.paf'
        hit=False;commands=[]
        if done.exists():
            record=load_json(done)
            hit=record.get('identity')==ident and paf.is_file() and sha256_file(paf)==record.get('paf_sha256')
        if not hit:
            done.unlink(missing_ok=True)
            cmd=[tools['minimap2']['path'],'-x',policy['alignment']['preset'],'-c','--cs','--secondary=no','-t','1',str(target),str(candidate)]
            command_run(cmd,directory,'alignment.paf',policy['execution']['command_timeout_seconds'])
            # Verify the inputs did not change while the aligner was reading them.
            if sha256_file(target)!=target_hash or sha256_file(candidate)!=ch:raise WorkflowError('Input changed during alignment')
            commands.append(cmd)
            write_json(done,{'identity':ident,'commands':[cmd],'paf_sha256':sha256_file(paf)})
        metric_file=call_dir/'metrics.json';calls=call_dir/'variants.tsv';filtered=call_dir/'eligible.paf'
        previous=load_json(metric_file) if metric_file.exists() else None
        metrics_hit=(hit and previous and all((call_dir/n).exists() and sha256_file(call_dir/n)==h for n,h in previous.get('artifacts',{}).items()) and set(previous.get('artifacts',{}))=={'variants.tsv','eligible.paf'})
        if metrics_hit:metrics=previous['metrics']
        else:
            metric_file.unlink(missing_ok=True)
            ql=fasta_lengths(candidate);stats=prepare_paf(paf,filtered,tl,ql,policy['variant_calling'])
            v=policy['variant_calling']
            cmd=[tools['k8']['path'],tools['paftools.js']['path'],'call','-l',str(v['min_alignment_length']),'-L',str(v['min_alignment_length']),'-q',str(v['min_mapping_quality']),str(filtered)]
            command_run(cmd,call_dir,'variants.tsv',policy['execution']['command_timeout_seconds']);commands.append(cmd)
            metrics={**stats,**parse_calls(calls,tl),'target_bases':sum(tl.values()),'candidate_bases':sum(ql.values())}
            ta,ca=metrics['target_aligned_bases'],metrics['candidate_aligned_bases']
            metrics.update(target_aligned_fraction=ta/sum(tl.values()),candidate_aligned_fraction=ca/sum(ql.values()),
                           snps_per_target_aligned_mb=metrics['snp_distance']*1e6/ta if ta else None,
                           rate_denominator='target bases covered uniquely by qualifying alignments; not a callable-site count')
            write_json(metric_file,{'metrics':metrics,'artifacts':{n:sha256_file(call_dir/n) for n in ('variants.tsv','eligible.paf')}})
        evidence=save_evidence(call_dir/'site_evidence.json',metrics,target_hash,ch,tl)
    summary = {k: v for k, v in metrics.items() if k not in {'snp_positions', 'target_unique_alignment_intervals'}}
    return {**summary,'site_evidence':evidence,'status':'COMPARED','cache_key':key,'cache_hit':hit,'interpretation_cache_hit':bool(metrics_hit),
            'target_sha256':target_hash,'candidate_sha256':ch,'artifact_directory':str(directory),'calls_directory':str(call_dir),
            'commands':commands,'elapsed_seconds':time.monotonic()-started}


def run_minimap(genomes,output_dir,cache_dir=None,policy_path=None,workers=None):
    if 'QUERY' not in genomes or len(genomes)<2 or not all(isinstance(v,str) for v in genomes.values()):raise WorkflowError('Assembly target QUERY and candidates are required')
    policy_path=policy_path or CONFIG_DIR/'assembly-comparison-policy.json';policy=load_json(policy_path)
    if workers is not None:policy['execution']['workers']=workers
    if policy['backend']!='minimap2' or policy['alignment']['preset']!='asm5':raise WorkflowError('Unsupported minimap2 policy')
    for value in policy['execution'].values():
        if type(value) is not int or value<=0:raise WorkflowError('Invalid execution budget')
    for name,value in policy['variant_calling'].items():
        if type(value) is not int or value<0 or (name=='min_alignment_length' and value==0):raise WorkflowError('Invalid variant-calling setting')
    if policy['variant_calling']['min_mapping_quality']>254:raise WorkflowError('Invalid mapping quality')
    for value in policy['comparability'].values():
        if not math.isfinite(value) or not 0<value<=1:raise WorkflowError('Invalid coverage criterion')
    tools={name:{'path':require_executable(name)} for name in ('minimap2','paftools.js','k8')}
    version=subprocess.run([tools['minimap2']['path'],'--version'],capture_output=True,text=True,check=True,timeout=30).stdout.strip()
    if version!=policy['minimap2_version']:raise WorkflowError('minimap2 version does not match pinned policy')
    for row in tools.values():row['sha256']=sha256_file(Path(row['path']))
    tools['minimap2']['version']=version
    output_dir.mkdir(parents=True,exist_ok=True);cache=(cache_dir or output_dir/'cache').resolve();cache.mkdir(parents=True,exist_ok=True)
    target=Path(genomes['QUERY']).resolve();th=sha256_file(target);tl=fasta_lengths(target);started=time.monotonic()
    def run(item):
        sample,path=item
        try:result=compare_pair(target,Path(path).resolve(),cache,policy,tools,th,tl)
        except (WorkflowError,OSError,ValueError,subprocess.SubprocessError) as error:result={'status':'ERROR','error':str(error),'cache_hit':False,'snp_distance':None}
        return {'sample':sample,**result}
    with ThreadPoolExecutor(max_workers=policy['execution']['workers']) as pool:rows=list(pool.map(run,sorted((k,v) for k,v in genomes.items() if k!='QUERY')))
    result={'backend':'minimap2','status':'PASS' if all(x['status']=='COMPARED' for x in rows) else 'PARTIAL','comparisons':rows,
            'commands':[cmd for row in rows for cmd in row.get('commands',[])],'policy':policy,'policy_sha256':sha256_file(policy_path),'tools':tools,
            'implementation_sha256':sha256_file(Path(__file__)),'cache_directory':str(cache),'alignment_jobs':sum(x['status']=='COMPARED' and not x['cache_hit'] for x in rows),
            'cache_hits':sum(x.get('cache_hit',False) for x in rows),'elapsed_seconds':time.monotonic()-started,'comparison_scope':'target_to_candidate'}
    write_json(output_dir/'run_minimap.json',result);return result


def main() -> int:
    """Compare a supplied candidate manifest without a Mashpit screen."""
    import argparse
    from fetch_reference_genomes import read_local_manifest
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', required=True)
    parser.add_argument('--candidate-manifest', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--cache-dir')
    parser.add_argument('--policy')
    parser.add_argument('--workers', type=int)
    args = parser.parse_args()
    refs = read_local_manifest(Path(args.candidate_manifest))
    genomes = {'QUERY': str(Path(args.target).resolve())}
    for accession, row in refs.items():
        path = Path(row['path'])
        if sha256_file(path) != row['sha256']:
            raise WorkflowError(f'Assembly checksum mismatch: {accession}')
        genomes[accession] = str(path)
    result = run_minimap(
        genomes, Path(args.output).resolve(),
        Path(args.cache_dir).resolve() if args.cache_dir else None,
        Path(args.policy).resolve() if args.policy else None, args.workers,
    )
    return 0 if result['status'] == 'PASS' else 2


if __name__ == '__main__':
    raise SystemExit(main())
