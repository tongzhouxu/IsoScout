# Changelog

All notable changes to this project are documented here. Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

- Workflow 1.8.0 / skill 1.4.0 / refinement policy 2.0.0 replaces flat member round-robin with qualifying-SNP-led focused expansion and rotating reserved alternative capacity. It tracks per-cluster new qualifying coverage, blocks premature stability, separates isolated comparability failures from unresolved clusters, and records allocation and uncertainty in reports/provenance. Initial selection, Mashpit return limits, and default resource limits remain unchanged. Real pilot and independent confirmation remain unrun here.

- Removed research comparison scripts, historical policy copies, and study-specific fixtures from the skill repository. The study will be developed separately for the paper; operational regression tests and software verification remain here.

- Selection policy 3.0.1 keeps the same candidate choices and budgets while renaming the `all_returned` inclusion reason to `all_returned_mode`.

- Workflow 1.7.0 / selection policy 3.0.0 removes the sketch-tolerance cluster gate. The global first 50 nonself unique returned references are guaranteed when eligible and budgeted; adaptive mode uses up to 100 returned references with explicit cluster-leader coverage, and all-returned mode requires enough budget for every returned reference. Expansion now considers members of every selected cluster. Previews, actual coverage, reports, and provenance distinguish resource omissions, failed downloads, comparability exclusions, and nearest ties. No biological accuracy claim.

- Workflow 1.6.0 / skill 1.2.0 adds `package_database_release.py` to create new database bundles with exact-release metadata and cluster membership, verified against every representative and recorded with checksums. Expansion uses validated bundled tables offline, while existing `databases-v1` packages retain exact-release download support. Corrupt bundles fail closed.

- Built, verified, and published all five [`databases-v2`](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-v2) packages with the exact NCBI Pathogen Detection membership and metadata tables (Salmonella `PDG000000002.3864`, E. coli/Shigella `PDG000000004.5775`, Listeria `PDG000000001.4519`, Campylobacter `PDG000000003.2682`, Cronobacter `PDG000000043.415`). Archive checksums, database integrity, and offline member discovery passed. The updated container remains unpublished.

- Raised the refinement metadata download ceiling from 512 MiB to 1 GiB (policy 1.0.1): the published Salmonella release's available metadata TSV is 786,821,557 bytes and would otherwise be rejected. Confirmed exact-release metadata and membership URLs for all five published databases on 2026-10-02.

- Audited all 90 tests and consolidated four redundant cases into stronger existing checks, retaining their distinct assertions. Removed repeated assertions and unused test imports. All 86 remaining tests pass in the pinned-tool container; measured in-process production execution coverage is unchanged.

- Workflow 1.5.0 adds opt-in exact-release cluster-member expansion and bounded refinement feedback with per-round audits, explicit incomplete coverage, and preservation of successful comparisons after later failures. Cleaned paired reads now feed SKA2 directly. Comparability gates exclude low-overlap rankings; tied minima remain ambiguous and zero SNPs no longer imply identity. Added an offline real-SKA2 integration test using seeded random DNA. Biological validation remains required.

- Workflow 1.4.0 enables adaptive SNP selection (policy 2.0.0): initial target 50, boundary near-tie expansion, alternative-cluster coverage, and hard ceiling 200. Deduplicates repeated accessions, records every decision, and exposes resource/retrieval limits and download losses. Diagnostics and actual refinement use the same selector. Added 13 adaptive-selection tests, including preview/download agreement. These defaults are not biologically validated.

- Workflow 1.3.0 adds similarity-distribution diagnostics: ranked scores/gaps, sketch-tolerance counts, cluster composition, retrieval-limit warnings, and a preview of the existing SNP selector. Reports include a rank–similarity plot and JSON link. A standalone script can analyze prior Mashpit runs using their recorded query settings. Candidate selection is unchanged. Added 16 diagnostic and integration tests.

- Renamed the project and GitHub repository to IsoScout, with skill identifier `isoscout`. Updated the container image, build paths, environment name, and report title. The default database directory is now `~/.isoscout/databases`, and the environment override is `ISOSCOUT_DATABASE_ROOT`; existing databases can still be selected with `--database-root`.

- Added `LICENSE` (GPL v2, matching Mashpit's own license verbatim) — the repo had none before, which meant it was "all rights reserved" by default despite being public. `CITATION.cff` now declares `license: GPL-2.0-only` to match.
- Validated the raw paired-end Illumina reads path end-to-end for the first time, against real data (a real Cronobacter SRA run, `SRR1614321`, not a database representative): fastp → SKESA → QUAST → Mashpit correctly recovered the isolate's true NCBI Pathogen Detection cluster (`PDS000112223.1`, score 0.986) from 102x real coverage. Every other verification in this project up to this point used pre-built assemblies; the reads path itself was previously only unit-tested with tiny synthetic/mocked data.
- `report.md` now shows the read QC line (estimated coverage, Q30 fraction, read pairs) when the input was reads, and shows plain-language input-type labels ("raw paired-end Illumina reads" / "an existing genome assembly") instead of the raw internal enum value (`illumina_paired_fastq`) — both gaps found by that same real-data run.

## [1.0.0] - 2026-08-17

Initial release.

### Core screening pipeline

- Deterministic, LLM-agnostic entrypoint (`scripts/screen_isolate.py`) accepting one FASTA/FNA assembly or a paired-end Illumina FASTQ pair.
- Fixed read QC and assembly workflow: `fastp` → `SKESA` → `QUAST`, with a hard 20x estimated-coverage floor.
- Broad and organism-specific assembly QC against evidence-based thresholds (`config/qc-policy.json`; see `references/qc-policy.md`).
- Organism routing either by explicit `--organism` or local `mlst` auto-detection against its bundled PubMLST schemes (`references/database-routing.md`).
- Fixed-profile Mashpit query (`--number 200 --threshold 0.85 --tie-tolerance-hashes 2`) against a pinned Mashpit commit (`538d3421302fe6dd129780605b8ff5dedbf4c046c`).
- Structural, non-invented interpretation of Mashpit's result: top candidate, ambiguity (`near_top`), and a `below_threshold` flag catching the case where Mashpit's own top hit is noise-level (found via real testing: a cross-genus query still produced a nominal "candidate" at score 0.007).
- Five supported organisms: `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, `cronobacter`. (`ecoli_shigella`, not `ecoli`, since the underlying PubMLST scheme genuinely can't distinguish the two genera.)

### Optional SNP resolution (`--snp-resolve`)

- Selects representative genomes from Mashpit's top cluster plus any `near_top` clusters, capped at 100 representatives per cluster / 100 genomes total (round-robin across clusters) so a cluster with hundreds or thousands of representatives can't turn into an unbounded ska2 run.
- Downloads representatives from NCBI via one batched `datasets download` call, retried on failure.
- Runs `ska2` (`ska build` / `ska distance`) for exact pairwise SNP distances between the query and every downloaded representative, and between representatives themselves.
- Reports a per-cluster distance summary, a raw-ratio confidence comparison between the nearest and next-nearest cluster (no invented strength label), and a Neighbor-Joining tree (pure-Python implementation, verified against a hand-computed additive case) with tips labeled by cluster.
- Renders both the ska2 SNP tree and Mashpit's own Mash-similarity tree to PNG, query highlighted, embedded in the report.

### Reporting

- `report.md`: a plain-language summary (organism determination, Mashpit's candidate clusters and scores, SNP tables/tree/confidence statement when `--snp-resolve` ran) for a reader who isn't going to open JSON.
- `result.json` / `provenance.json`: structured result and full reproducibility record (checksums, pinned tool versions, every command run) as before.

### Distribution

- Pre-built, checksummed Mashpit databases for all five organisms published on the [`databases-v1`](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-v1) release.
- Container image published to `ghcr.io/tongzhouxu/isoscout` via GitHub Actions on version tags (`linux/amd64` only — `quast=5.3.0` has no `linux/arm64` build for the pinned Python 3.11).

### Verification

Validated end-to-end against real production Mashpit databases (Listeria, Salmonella, Cronobacter) with ground truth cross-referenced against NCBI Pathogen Detection's own isolate/cluster files (not just Mashpit self-consistency), including inside the actual published container. That process surfaced and fixed several real, pre-existing issues along the way: a `quast`/Python 3.11/`linux-arm64` build incompatibility, a missing C compiler blocking a pip dependency's build, a `sourmash`/`pkg_resources` break from an unpinned `setuptools`, and a provenance-collection bug that recorded a crashed tool's traceback as if it were a version string.

36 unit tests, all mocking external bioinformatics execution (no Docker/database required to run them).
