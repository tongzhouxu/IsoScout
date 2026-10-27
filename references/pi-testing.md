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
