# IsoScout

[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](LICENSE)

IsoScout screens a bacterial isolate (a genome you sequenced) against public outbreak-cluster data, to tell you which known cluster of related bacteria it's most similar to. You do not need to know how to code to use it — you need to install two ordinary applications, then just describe what you want in plain English to an AI coding assistant (Claude Code or Codex), and it does the rest. That's what this page walks you through first. If you're comfortable with the command line already, skip to [Manual setup](#manual-setup) below.

## Quick start (no coding experience needed)

### What you'll need

- A Mac, Windows (with WSL2), or Linux computer with at least **10GB of free disk space** and a decent internet connection (you'll be downloading a few GB).
- A genome file for the bacterial isolate you want to screen: either one assembled genome (a `.fasta`/`.fna` file), or two raw paired-end Illumina read files (commonly named something like `sample_R1.fastq.gz` and `sample_R2.fastq.gz`). If you don't have one yet and just want to try the process, ask the assistant in Step 3 below to fetch a public example genome for you from NCBI instead of using your own file.
- Roughly 15-40 minutes for one-time setup (mostly waiting for downloads), then a few minutes per genome you screen after that.

### Step 1 — Install Docker Desktop

This is the program that actually runs the analysis in an isolated, reproducible way. It's a normal application, installed like any other:

1. Go to [docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop/) and download it for your computer.
2. Install it the way you'd install any app (double-click, follow the prompts).
3. Open the Docker Desktop app once and leave it running in the background. You'll know it worked if you see a little whale icon in your menu bar / system tray. It needs to stay running any time you use this skill.

### Step 2 — Install Claude Code or Codex

This is the AI assistant that will read your instructions, download things, and run the analysis for you — so you don't have to type any of the technical commands yourself.

- **Claude Code**: see [claude.com/claude-code](https://claude.com/claude-code) for the current install instructions for your computer.
- **Codex**: see OpenAI's current Codex documentation for install instructions.

Either one works equally well with this skill. Once installed, open it (this usually means: open the Terminal app on your computer, type the assistant's name, e.g. `claude`, and press Enter — the install instructions above will confirm the exact command for your setup) so you have a chat window open with it.

### Step 3 — Ask it to set everything up

Copy the box below, paste it into your chat with Claude Code or Codex, replace `[ORGANISM]` with whichever of `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, or `cronobacter` matches your sample (if you're not sure, just pick `cronobacter` for a first try — its database is much smaller and faster to download), and send it:

```
Please set up IsoScout for me:

1. Clone https://github.com/tongzhouxu/IsoScout
   into a folder on my computer.
2. Read that repository's README.md and follow its "Manual setup"
   instructions to pull the published Docker image, and to download
   and checksum-verify the [ORGANISM] database from its GitHub
   Release into ~/.isoscout/databases/.

Explain each step briefly as you go, and tell me clearly once
everything is downloaded and ready.
```

It will likely pause partway through to ask your permission before running certain commands (downloading files, running Docker) — that's normal and expected; read what it says it's about to do, and approve it if it matches the steps above.

### Step 4 — Ask it to screen your genome

Once setup is done, copy this box, fill in the path to your genome file(s) and the organism, and send it:

```
Please screen this genome using IsoScout, which you
just set up: [PATH TO YOUR .fasta FILE, OR YOUR TWO .fastq.gz FILES]

Organism: [ORGANISM]

Once it finishes, open the report.md file it produces and explain
the result to me in plain English — in particular, tell me which
cluster it's closest to, how confident that is, and what you'd
recommend I do next.
```

A first run typically takes a few minutes. If you don't know the organism ahead of time, say so instead of naming one — it can auto-detect it from your genome, it'll just take slightly longer.

### If something goes wrong

- **"Cannot connect to the Docker daemon"** — Docker Desktop isn't running. Reopen the Docker Desktop app and wait for the whale icon to stop animating, then try again.
- **Running out of disk space** — each organism's database ranges from ~3MB (cronobacter) to ~530MB (ecoli_shigella) compressed; you only need the one(s) you're actually screening against, not all five.
- **Anything else** — paste the error message back to the assistant and ask it to diagnose and fix it; that's exactly the kind of thing it's good at.

## What this actually does

Screens a bacterial isolate — raw paired Illumina reads or an assembly — against [Mashpit](https://github.com/tongzhouxu/mashpit), a MinHash-sketch database of NCBI Pathogen Detection SNP clusters, and optionally confirms a candidate at true SNP resolution with [ska2](https://github.com/bacpop/ska.rust). It tells you which known cluster of bacteria your sample most resembles — useful for narrowing down a possible outbreak connection, though it is a screening tool, not proof of one (see [references/limitations.md](references/limitations.md)).

Every command, threshold, and parameter it runs is fixed in version-controlled config (`config/*.json`) — the AI assistant only runs one script and reports the result back to you; it never invents a bioinformatics command or a cutoff on its own. See [SKILL.md](SKILL.md) for the full technical contract it follows.

Supports `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, and `cronobacter`. Nothing is uploaded during a screen — your genome stays on your own computer. Mashpit is installed from a pinned upstream commit (`538d3421302fe6dd129780605b8ff5dedbf4c046c`), not the older published PyPI release.

## Manual setup

For anyone who'd rather run the commands themselves instead of asking an assistant to.

You need two things: the container image, and a Mashpit database for at least one organism.

The container image uses `ghcr.io/tongzhouxu/isoscout`; the examples tag it locally as `isoscout:local`.

**1. Get the image** — either pull the pre-built one:

```bash
# The registry package is currently private; authorized GHCR login is required.
docker pull --platform linux/amd64 ghcr.io/tongzhouxu/isoscout:latest
docker tag ghcr.io/tongzhouxu/isoscout:latest isoscout:local
```

or build it yourself:

```bash
docker build --platform linux/amd64 --tag isoscout:local --file container/Dockerfile .
```

Build check (2026-10-01): a clean dependency install currently fails with a conda dependency-resolution conflict. The local IsoScout image was built by retaining the existing image's installed dependencies and replacing its bundled project files; its entrypoint and all 36 unit tests passed. The renamed image was published to `ghcr.io/tongzhouxu/isoscout:latest` on 2026-10-01.

`--platform linux/amd64` is required everywhere here, not just on Apple Silicon: the image is only published/buildable for `linux/amd64` because `quast=5.3.0` has no native `linux/arm64` build compatible with the pinned Python 3.11. It runs fine under emulation on Apple Silicon; omitting the flag there pulls/builds nothing since Docker defaults to your host's native architecture.

**2. Get a database** — download the pre-built ones from [Releases](../../releases/tag/databases-v1):

```bash
mkdir -p ~/.isoscout/databases && cd ~/.isoscout/databases
for org in salmonella ecoli_shigella listeria campylobacter cronobacter; do
  curl -LO "https://github.com/tongzhouxu/IsoScout/releases/download/databases-v1/${org}.tar.gz"
done
curl -LO https://github.com/tongzhouxu/IsoScout/releases/download/databases-v1/checksums.sha256.txt
shasum -a 256 -c checksums.sha256.txt   # verify before extracting
for f in *.tar.gz; do tar -xzf "$f"; done
```

You only need the organism(s) you actually plan to screen against — each is independent. `screen_isolate.py` additionally verifies a per-file checksum recorded in each organism's `database.json` on every run, so a corrupted or tampered database is caught automatically, not just at download time.

Database creation and updating are out of scope for this skill; see [mashpit](https://github.com/tongzhouxu/mashpit) itself for that.

## Run a screen

```bash
docker run --rm \
  --volume "/absolute/path/to/data:/data:ro" \
  --volume "$HOME/.isoscout/databases:/databases:ro" \
  --volume "/absolute/path/to/results:/results" \
  isoscout:local \
  /data/sample.fasta \
  --organism salmonella \
  --database-root /databases \
  --output /results/sample-screen
```

- Input is one assembly (`.fa`/`.fasta`/`.fna`), or two paired-end FASTQ files — replace the assembly argument with both paths (e.g. `/data/sample_R1.fastq.gz /data/sample_R2.fastq.gz`) to run the fixed fastp → SKESA → QUAST assembly workflow first.
- Omit `--organism` to auto-detect it with local `mlst` against its bundled PubMLST schemes instead of asserting it.
- The output directory must not already exist — nothing gets silently overwritten.

Add `--snp-resolve` to also download the relevant representative genomes from NCBI and compute exact pairwise SNP distances with ska2 once a Mashpit candidate is found — a Neighbor-Joining tree, a per-cluster distance summary, and a confidence comparison between the nearest and next-nearest cluster. This is the only step that reaches out to the network (for public reference genomes; the query itself is never uploaded), so it's opt-in. See [references/snp-resolution.md](references/snp-resolution.md).

### Reading the result

Three files land in the output directory:

- **`report.md`** — a plain-language summary for a non-technical reader: organism determination, read QC when the input was raw reads, Mashpit's candidate clusters and scores plus its own Mash-based tree, and (with `--snp-resolve`) the ska2 SNP tables, confidence statement, and SNP tree — both trees rendered as PNGs with the query highlighted.
- **`result.json`** — the structured, authoritative result. Use `status`, `stop_reason`, and `user_summary`.
- **`provenance.json`** — checksums, pinned tool versions, and every command actually run, for reproducibility.

## Test without biological tools or databases

```bash
PYTHONPYCACHEPREFIX=/tmp/isoscout_pycache python3 -m unittest discover -s tests -v
```

These unit tests mock external bioinformatics execution and don't need Docker, Mashpit, or a real database. A real end-to-end run additionally needs the container and at least one downloaded database, per Manual setup above.

## Reference docs

- [SKILL.md](SKILL.md) — the agent-facing contract: when to run what, how to report results
- [references/setup.md](references/setup.md) — container/database contract in full
- [references/workflow.md](references/workflow.md) — the fixed pipeline, stage by stage
- [references/database-routing.md](references/database-routing.md) — organism routing rules
- [references/qc-policy.md](references/qc-policy.md) — QC thresholds and their evidence basis
- [references/mashpit-interpretation.md](references/mashpit-interpretation.md) — how a Mashpit result is labeled
- [references/snp-resolution.md](references/snp-resolution.md) — the optional ska2 SNP-resolution step
- [references/limitations.md](references/limitations.md) — scope and scientific limitations
- [CHANGELOG.md](CHANGELOG.md) — release history
- [CITATION.cff](CITATION.cff) — how to cite this skill and Mashpit itself
