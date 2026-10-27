# Raspberry Pi full-analysis pilot

Use a 64-bit Linux Pi with an existing Python 3.11 IsoScout/Mashpit environment,
preassembled FASTA inputs, and one complete database package at a time. Raw-read
processing is not supported on Pi. No Docker or AI assistant is needed.

After updating IsoScout, install its native assembly-comparison tools:

```bash
sudo apt-get install -y curl unzip time
source ~/isoscout-pi-venv/bin/activate
bash scripts/install_native_alignment.sh "$VIRTUAL_ENV"
```

The installer uses an isolated native-tools environment managed by micromamba.
It installs minimap2, paftools, k8, and NCBI datasets; see [setup.md](setup.md).

For the prepared paper benchmark kit, copy `scripts/run_pi_pilot.sh` into
`~/pi-benchmark/`. The kit layout is `bundles/<organism>/inputs.json`,
`bundles/<organism>/queries/`, and `bundles/<organism>/database/`. Keep the complete
database folder, including its representative and membership metadata.

```bash
cd ~/pi-benchmark
source ~/isoscout-pi-venv/bin/activate
bash run_pi_pilot.sh cronobacter
```

This runs the first manifest query through QC, Mashpit, and **default minimap2
SNP comparison**, with one comparison worker. It excludes the query accession
from reference selection. Reference assemblies are downloaded as needed; network
time is included in the elapsed time. It does not enable extra cluster-member
expansion. Each invocation creates a fresh `results/<organism>-isoscout-full.*`
directory containing logs, GNU time measurements, endpoint temperature/throttling
flags, versions, and `analysis/{result.json,report.md,provenance.json}`.

Check `analysis_mode: full`, `SNP backend: minimap2`, and the overall and SNP
statuses. `INCOMPLETE` or `STOPPED` with exit code 2 is not a successful full run.
A completed comparison can be ambiguous; preserve that finding. Missing downloads
or failed comparisons must not be counted as SNP-supported results.

Start with **one full query per organism** to demonstrate feasibility. These five
pilots do not establish timing distributions or typing accuracy. The separate
Mashpit timing design uses ten queries per organism and three repetitions;
`run_benchmark.py` measures Mashpit alone. Do not label those timings as full
IsoScout runtime. Existing screening-only results remain screening-only results.

Copy results back to the server and verify they arrived before deleting the old
organism's database bundle on Pi. Full-analysis directories also contain downloaded
references and alignment artifacts; archive those directories before removing
them if space is needed. Then transfer the next complete bundle and repeat.

## Separate stage timings

CLI 1.5.0.dev2 records `analysis/stage_timings.json`, also embedded in
`analysis/result.json`. The pilot prints these values and saves `summary.txt`:

- `mashpit`: the Mashpit invocation and its wrapper processing.
- `reference_retrieval`: fetching/local lookup, retries and reference verification.
- `snp_comparison`: the backend comparison call, including minimap2, paftools
  variant calling and comparison processing; excludes reference retrieval.

All values are elapsed wall seconds measured with a monotonic clock, accumulated
across rounds. `calls` records attempts and `failed_calls` counts calls that raised an exception.
Partial results returned normally are described by the analysis status, not this
exception counter; unexecuted times are null. Partial timings are retained when
an ordinary workflow error is caught. Force-killing the process can prevent final
files being written. Cache reuse is included in the observed times.

These are not per-stage CPU or memory measurements. Their sum excludes QC,
selection, diagnostics, downstream interpretation/reporting, provenance and
startup. Use GNU time's `time.txt` for the whole CLI's elapsed time, CPU and peak
RSS. Older runs cannot be given exact missing stage timings retroactively.

## Full batch across the prepared queries

Keep `run_pi_batch.py` and `run_pi_pilot.sh` together in `~/pi-benchmark`.
The runner processes every query in the bundle, serially, in a saved randomized
order. It verifies query hashes, writes `runs.tsv` and `summary.json` after every
run, and stops on failure with partial artifacts retained. The summary identifies
planned, recorded and successful runs. It never overwrites an existing batch.
There is no automatic resume or dedicated warm-up. Each invocation uses a fresh
output/comparison cache; filesystem caches are not cleared. Retrieval time is
included in total time and recorded separately.

For limited Pi storage, run one repetition (ten distinct queries) at a time:

```bash
cd ~/pi-benchmark
source ~/isoscout-pi-venv/bin/activate
python run_pi_batch.py --organism cronobacter --repeats 1 \
  --seed 20261027 --output results/cronobacter-full-r1
```

Copy the entire completed batch back to the server and verify the copy before
removing that batch's results on Pi. Keep the database until all repetitions for
that organism finish. Repeat with fresh output names `cronobacter-full-r2` and
`cronobacter-full-r3`, using seeds `20261028` and `20261029`. Each separate batch
numbers its repetition locally from 1; preserve the r1/r2/r3 batch names.
Three batches give 30 full runs per organism; five organisms give 150 runs over
50 distinct queries. Repetitions are not independent biological observations.
If storage permits keeping all artifacts, `--repeats 3` runs them in one batch.
Downloaded genomes and alignment caches can occupy much more space than the
input database; archive results between batches on a small SD card.

After all three batches are safely copied, delete the old organism's database
bundle on Pi, transfer the next complete bundle, and repeat. For the figure,
use each query's median of its three repetitions, show total and separate stage
wall times, and retain warnings, thermal flags and failed attempts. Do not mix
these full-analysis measurements with `run_benchmark.py` Mashpit-only timings.
