# IsoScout

[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](LICENSE)

IsoScout screens a bacterial isolate against local [Mashpit](https://github.com/tongzhouxu/mashpit) databases of NCBI Pathogen Detection SNP clusters, with optional higher-resolution comparisons against candidate assemblies.

Its **AI-agent skill** helps you set up an analysis, run it from a plain-language request, and understand the evidence and limitations in the report. The **command-line interface (CLI)** exposes the same workflow for direct use and automation. Both use the same versioned scientific code and policies; the assistant does not replace the analysis or decide its thresholds.

| How you want to work | Start here |
| --- | --- |
| Describe your task to an assistant and discuss the results | [Skill-assisted quick start](#skill-assisted-quick-start) |
| Install and run IsoScout yourself | [Command-line quick start](#command-line-quick-start) |
| Run without Docker using your own tool environment | [Native execution](#native-execution) |
| Screen an existing assembly on a Raspberry Pi | [Raspberry Pi: assemblies only](#raspberry-pi-assemblies-only) |

Supported organism keys are `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, and `cronobacter`. A candidate cluster is screening evidence, not confirmation of outbreak membership. See [limitations](references/limitations.md).

## Skill-assisted quick start

Use this path if you want an assistant to manage the commands and explain the results. Keep the IsoScout checkout: [SKILL.md](SKILL.md), `references/`, and `agents/` provide its instructions and supporting material. Installing the CLI alone does not register the skill with an assistant.

### 1. Prepare your computer and assistant

Use a terminal-capable assistant, such as Claude Code or Codex, that can read local files and execute commands. Point it to the checkout's `SKILL.md`; it can then follow IsoScout's workflow instructions. See your assistant's documentation for automatic skill discovery or installation.

For the default setup, you need:

- Python 3.9 or newer; **Python 3.11 is recommended**.
- A running Docker installation: [Docker Desktop](https://www.docker.com/products/docker-desktop/) on a desktop computer, or Docker Engine on Linux. Run the command-line steps inside WSL2 on Windows.
- Disk space for the container, extracted databases, and outputs. Download only the organism database(s) you need.
- One assembled genome (`.fa`, `.fasta`, or `.fna`), or paired Illumina FASTQ files on a workstation/server with sufficient memory for assembly. The configured SKESA memory budget is 16 GB, in addition to system overhead.

The current container targets `linux/amd64`; Apple Silicon uses emulation. **Raspberry Pi users should use the native, assembly-only path below. Raw-read processing is not supported on the Pi.** Setup and run times depend on downloads, database size, input and requested stages.

### 2. Ask the assistant to set up IsoScout

Replace the bracketed text in this prompt. If you already cloned the repository, give its path so the assistant can reuse it. Use `UNKNOWN` if you do not know the organism independently.

```text
Please set up IsoScout for me.

Repository: https://github.com/tongzhouxu/IsoScout
Existing checkout, if any: [PATH, or NONE]
Organism: [salmonella, ecoli_shigella, listeria, campylobacter,
           cronobacter, or UNKNOWN]

Reuse my existing checkout, or clone the repository if I have none.
Read SKILL.md and README.md. Check my computer's requirements, create
.venv in the checkout, and install the CLI there. For this workstation,
build the current Docker image as isoscout:local. Use the absolute path
to .venv/bin/isoscout in later commands.

Download and verify the matching databases-v2 package under
~/.isoscout/databases/. If the organism is UNKNOWN, prepare all five
supported packages for local MLST routing. Do not guess the organism.

Explain what you configured and tell me when it is ready.
```

This prompt uses the default workstation/container setup. If a native environment is already configured, explicitly tell the assistant to use `--executor native` and give its executable and database paths. Changing the interface does not remove the [Pi limitations](#raspberry-pi-assemblies-only).

### 3. Ask for a screen and an explanation

```text
Read the IsoScout SKILL.md in [CHECKOUT PATH] and screen:
[ABSOLUTE PATH TO ONE ASSEMBLY, or R1 AND R2 FASTQ PATHS]

Organism: [KNOWN ORGANISM, or UNKNOWN]
Output directory: [NEW DIRECTORY PATH]

Start with the screening workflow. Read report.md and result.json,
and explain the candidate clusters, alternatives or ties, QC stops,
and limitations. Distinguish completed analysis from failed or
unperformed stages. Keep the structured results and provenance.
```

The assistant should preserve QC stops and unresolved findings, explain which references were examined, and avoid presenting a screening result as outbreak confirmation. You can discuss [optional refinement](#optional-refinement) after reviewing the screen.

The analysis runs locally. The screening workflow does not upload the query genome, but your chosen assistant may send conversation or tool content to its provider. Reference downloads also require network access when assemblies are not supplied locally.

## Command-line quick start

An AI assistant is not required. The examples below use a Bash-compatible shell on Linux, macOS, or WSL2.

### 1. Install the CLI

From a checkout, with Python 3.9 or newer (3.11 recommended):

```bash
git clone https://github.com/tongzhouxu/IsoScout.git
cd IsoScout
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
isoscout --version
isoscout screen --help
```

If you already have a checkout, start from that directory. This CLI package is installed from source or a built wheel; it is not yet published to PyPI. `python -m isoscout` also works after installation. In a new shell, activate the environment again or use the absolute path to `.venv/bin/isoscout`.

The package includes the workflow scripts and policy files. **Install a runtime and a database as well:** `pip install .` does not install Docker, an image, Mashpit, other analysis tools, or databases.

### 2. Choose the runtime

| Runtime | Use it when | Setup |
| --- | --- | --- |
| Docker (default) | You want the reference tool environment on a workstation/server | Build the image below |
| Native (`--executor native`) | You manage the tools yourself, including on ARM64 | Follow [native execution](#native-execution), then continue to step 3 |

For Docker, start Docker and build from the IsoScout checkout:

```bash
docker build --platform linux/amd64 \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  -f container/Dockerfile -t isoscout:local .
```

The recipe extends a pinned baseline image and includes this checkout's workflow. The image is **AMD64-only**; the launcher requests that architecture explicitly. On Apple Silicon it runs under emulation. On a Pi, use native execution instead.

The container pins Mashpit to upstream commit `538d3421302fe6dd129780605b8ff5dedbf4c046`. Tool versions, policies, and the actual runtime image identity are recorded for reproducibility; see [setup details](references/setup.md).

### 3. Install an organism database

Both runtimes use the same [databases-v2 packages](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-v2). Set `org` to the organism you know your sample belongs to. This example downloads, verifies, and extracts one package; it stops if verification fails and leaves your shell in its original directory.

```bash
(
  set -euo pipefail
  mkdir -p ~/.isoscout/databases
  cd ~/.isoscout/databases
  org=salmonella
  curl -fLO "https://github.com/tongzhouxu/IsoScout/releases/download/databases-v2/${org}.tar.gz"
  curl -fLO https://github.com/tongzhouxu/IsoScout/releases/download/databases-v2/checksums.sha256.txt
  grep "  ${org}.tar.gz$" checksums.sha256.txt | shasum -a 256 -c -
  tar -xzf "${org}.tar.gz"
)
```

The layout is `~/.isoscout/databases/<organism>/`, containing `database.json`, the `.db` and `.sig`, and any bundled membership files. Transfer the **whole organism folder** if copying from another machine. The workflow checks database metadata and recorded checksums on each run.

Only the selected organism's package is required when you supply `--organism`. Automatic routing requires all databases it may select. Use `--database-root /path/to/databases` for a different location; point it at the parent of the organism folders, not at a `.db` file. Building a database is a separate [Mashpit](https://github.com/tongzhouxu/mashpit) task.

### 4. Run a screen

**From an existing assembly:**

```bash
isoscout screen /absolute/path/to/sample.fasta \
  --organism salmonella \
  --output /absolute/path/to/results/sample-screen
```

**From paired Illumina reads — workstation/server only, not Raspberry Pi:**

```bash
isoscout screen /absolute/path/to/sample_R1.fastq.gz \
  /absolute/path/to/sample_R2.fastq.gz \
  --organism salmonella \
  --output /absolute/path/to/results/sample-read-screen
```

The read path runs fastp → SKESA → QUAST before screening the resulting assembly. Use a machine with sufficient memory and the full tool environment for these stages.

These examples use Docker by default. Add `--executor native` only when the corresponding tools are installed locally. **The output directory must be new.** Omit `--organism` only when using local MLST routing with its tools and databases installed; routing is limited to the five supported groups and is not a general species-identification method.

Open `report.md` in the output directory and check `result.json` for stage status; see [reading the result](#reading-the-result). The existing `python3 scripts/launch_isoscout.py ...` entry point remains available to source-checkout users.

## Native execution

Native mode runs the bundled Python workflow using tools already installed on `PATH`. Use Python 3.11 for the documented Mashpit environment. Install Mashpit and the tools required for your input/stages in an appropriate environment, using the [setup reference](references/setup.md) and [version configuration](config/workflow.json). The container's `alignment-linux-64.lock` is an x86-64 lock, not an ARM64 installer.

| Requested stage | Tools needed in addition to the IsoScout CLI |
| --- | --- |
| Assembly screening with a known organism | Compatible Mashpit and its Python dependencies |
| Diagnostic plots | Matplotlib in the Python environment running IsoScout |
| Automatic organism routing | `mlst` and its bundled schemes |
| Paired-read processing, on a workstation/server | `fastp`, SKESA and QUAST |
| Assembly refinement | Pinned minimap2, `paftools.js`, `k8`; NCBI `datasets` if references are not supplied locally |
| Legacy paired-read refinement | SKA2 |

Check the environment before running an analysis:

```bash
command -v mashpit
mashpit --help
python -c 'import matplotlib'
isoscout screen --executor native --help
```

Then run an assembly screen:

```bash
isoscout screen --executor native /absolute/path/to/sample.fasta \
  --organism salmonella \
  --database-root "$HOME/.isoscout/databases" \
  --output /absolute/path/to/results/sample-native-screen
```

Native installations need their own verification. Neither installing the CLI nor successfully printing help establishes that every analysis dependency works. Container-only options such as `--image` and `--print-command` do not apply to native mode. Nextflow is not required.

### Raspberry Pi: assemblies only

**Raw-read processing is not supported on Raspberry Pi in the current IsoScout setup. Do not supply paired FASTQ files on the Pi.** The configured read workflow uses SKESA with a 16-GB memory budget, exceeding the RAM of the 8-GB Pi 4, and the complete ARM64 read-processing stack has not been validated. Assemble and QC reads on a workstation or server, then transfer the resulting FASTA assembly to the Pi.

The intended Pi use is **native screening of an existing assembly**, with the organism supplied explicitly and one database installed at a time:

- Use a 64-bit ARM OS and Python 3.11, with ARM-compatible Mashpit dependencies installed separately. Some dependencies may require source builds.
- Use `--executor native` as in the assembly example above. The current Docker image targets AMD64.
- Copy the complete organism database folder, including `database.json`; allow room for its extracted size and output files.
- Keep SNP refinement and member expansion on a workstation/server for this deployment. Their full Pi execution and performance have not been validated.

**Pi validation is still pending:** the CLI is pure Python and its Python requirement is compatible with this setup, but a complete assembly-screening run still needs verification on the device. The raw-read limitation is deployment guidance; the CLI does not currently detect a Pi and reject FASTQ automatically.

## Reading the result

Both interfaces and both runtimes produce the same primary files:

| File | What to use it for |
| --- | --- |
| `report.md` | Read the QC outcome, candidate clusters, scores, alternatives and limitations; ask the assistant to explain it if using the skill |
| `result.json` | Inspect the authoritative structured result, including `status`, `stop_reason`, stage results and warnings |
| `provenance.json` | Retain input/database checksums, configured and observed tool identities, commands, CLI version and executor |

Exit status zero means the workflow completed, possibly with warnings. A QC stop or required-stage failure returns nonzero. Optional refinement failure can leave a successful screen with warnings, so check the requested stage's result as well as the exit code. Docker runs additionally record the actual image identity. Keep the full output directory for review and reproducibility.

## Optional refinement


Add `--snp-resolve` to compare the target assembly directly with selected
candidate assemblies using minimap2 and its upstream paftools variant caller.
`--assembly-manifest` supplies local accession/path/sha256 records; only missing
candidates are downloaded. `--comparison-cache` reuses completed pairs during
expansion or a subsequent run. Reports retain SNPs, indels, coverage, ties and
unresolved comparisons. No candidate-to-candidate matrix or tree is calculated
for assemblies. MUMmer4 is available with `--snp-backend mummer`; paired reads
retain legacy SKA2. See [assembly comparison](references/assembly-comparison.md)
and [selection policy](references/snp-resolution.md).

Add `--snp-selection-mode all_returned` to compare every valid unique returned representative when the recorded 200-attempt budget fits. Supply `--query-accession` if the query itself appears in the database. Add `--snp-expand` to enable SNP refinement plus bounded exploration of additional members from the exact NCBI cluster release. The revised default prioritizes SNP-leading clusters after broad representative screening, reserves rotating slots for alternatives, and requires focused qualifying coverage before a stability stop. `--refinement-policy PATH` accepts an explicit versioned policy file with recorded budgets. Each round records allocation and stopping rationale; stable sampling is not exhaustive search or proof of nearest-genome recovery. For paired-read inputs, SKA2 uses cleaned reads directly, although Mashpit still needs the generated assembly. See [cluster expansion](references/cluster-expansion.md) and [software verification](references/verification.md). The implementation and comparability defaults require biological validation before strain-assignment claims.

The legacy `databases-v1` packages do not include full cluster membership.
All five [databases-v2](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-v2) packages were built and verified with the
[database package builder](references/setup.md#database-packages-with-local-cluster-membership),
then published with archive checksums and a release manifest. They bundle exact-release metadata
and cluster-membership tables so expansion can work after NCBI removes that
release. Existing packages continue to use exact-release retrieval while it is
available. Candidate assemblies must be supplied locally or downloaded separately.

## CDC REP designations

The [databases-rep-2026-10](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-rep-2026-10)
release adds a `REP_CLUSTER` table to each Mashpit `.db` (Salmonella
`PDG000000002.4263`, E. coli/Shigella `PDG000000004.6366`, Campylobacter
`PDG000000003.2927`). Columns: `PDS_acc` (primary key), `rep_designation`,
`organism`, `anchor_pnusa`, `anchor_pdt`, `evidence_status`, `source_url`. The
release's `rep_cluster_annotations.json` lists the same annotations. Install
these packages as for `databases-v2`, using the `databases-rep-2026-10` tag.

When the table exists, the screen reads it directly; you no longer need to
query it with `sqlite3`. `result.json` gains a `rep_designations` block, also
saved as `rep_designations.json`:

- `final_clusters` — the final called cluster, which is the SNP-resolved
  cluster or every unresolved SNP-supported cluster when refinement ran, and
  Mashpit's top candidate otherwise (`final_call_basis` records which).
- `candidate_clusters` — Mashpit's best and displayed alternative clusters.
- Each entry has `status` `LABELLED` (with `rep_designation`, `source_url`,
  anchors and `evidence_status`) or `NO_LABEL`.
- `status` is `AVAILABLE`, `TABLE_ABSENT` or `UNREADABLE`.

Matching uses the exact PDS accession, including its version. A missing table
or row means *no REP label in this database*, not *not a REP strain*. A REP
label means CDC has linked the called NCBI SNP cluster to a REP strain; it does
not confirm that the isolate belongs to the outbreak strain. `report.md` adds a
"CDC REP designations" section saying this. With `databases-v2` packages, the
screen's status, warnings and calls are unchanged; the block only records that
the table is absent.

## Similarity-distribution diagnostics

Each successfully parsed Mashpit query now adds a rank–similarity plot and numerical diagnostics to the report. These show score gaps at ranks 50/100/200, the number of representatives within Mashpit's sketch tolerance of the best score, and whether retrieval or the current SNP selection policy may cut through a near-tie. They preview policy 3.0.1: include the global top 50 unique nonself representatives, use up to 50 more slots for alternative-cluster coverage and global rank, and cap total reference attempts (including member expansion) at 200. Sketch near-ties describe scores but do not exclude SNP candidates. `--snp-selection-mode all_returned` attempts all returned representatives if they fit the recorded 200-attempt budget. Every inclusion, resource omission, download, and qualifying or excluded comparison is recorded; these computational defaults still require biological validation.

Full ranked scores, cluster composition, source checksums, and cutoff flags are saved in `similarity_distribution/summary.json`. See [the diagnostics guide](references/similarity-distribution.md) for interpretation and the standalone command for previous runs.

## Reference docs

- [SKILL.md](SKILL.md) — the agent-facing contract: when to run what, how to report results
- [references/setup.md](references/setup.md) — container/database contract in full
- [references/workflow.md](references/workflow.md) — the fixed pipeline, stage by stage
- [references/database-routing.md](references/database-routing.md) — organism routing rules
- [references/qc-policy.md](references/qc-policy.md) — QC thresholds and their evidence basis
- [references/mashpit-interpretation.md](references/mashpit-interpretation.md) — how a Mashpit result is labeled
- [references/snp-resolution.md](references/snp-resolution.md) — the optional SNP-refinement step
- [references/cluster-expansion.md](references/cluster-expansion.md) — focused member expansion and resource limits
- [references/verification.md](references/verification.md) — software test commands and interpretation limits
- [references/limitations.md](references/limitations.md) — scope and scientific limitations
- [CHANGELOG.md](CHANGELOG.md) — release history
- [CITATION.cff](CITATION.cff) — how to cite this skill and Mashpit itself
