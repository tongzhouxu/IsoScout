# Limitations

IsoScout supports one bacterial isolate represented by paired-end Illumina FASTQ/FASTQ.GZ or an existing FASTA/FNA assembly. It does not support ONT, PacBio, hybrid, interleaved, single-end, or metagenomic input.

The workflow does not confirm outbreak membership, transmission, source attribution, clinical significance, or treatment decisions. Mashpit MinHash similarity can prioritize cluster representatives but does not replace validated SNP/cgMLST analysis and epidemiological investigation.

Local MLST auto-detection and the QC policy require validation with positive controls, unsupported near neighbors, contaminated isolates, and site-specific sequencing distributions. MLST scheme selection is not a general taxonomic classifier and can misassign close relatives or poor/mixed assemblies. A user-supplied organism bypasses this check. Cronobacter QC bounds especially require empirical calibration. The workflow retains warnings rather than masking these limitations.

The bundled `ecoli`/`ecoli_2` PubMLST schemes are shared by Escherichia and Shigella in the upstream scheme metadata, so automatic routing cannot reliably distinguish those genera. This is why the routed database key is `ecoli_shigella`, not `ecoli`: treat a selected `ecoli_shigella` database as a routing choice, not definitive species identification, whether reached via `--organism` or local MLST.

Building the underlying Mashpit `.db` and `.sig` is outside scope, while
packaging those files with exact-release membership tables is supported. Results
are only reproducible when the exact database archive and checksums, tool image,
code revision (especially for a bind-mounted checkout), downloaded reference
assemblies, and run provenance are retained.

Optional `--snp-resolve` (see [references/snp-resolution.md](snp-resolution.md)) refines a Mashpit candidate with ska2 pairwise SNP distances against representative genomes re-downloaded from NCBI. It is a screening refinement bounded by whatever representatives the underlying Mashpit database happened to select, not a validated outbreak-confirmation pipeline, and it introduces a network dependency the rest of this workflow does not otherwise have.

`--snp-expand` extends beyond representatives within selected returned clusters, using exact-release membership and bounded rounds. Metadata gaps, absent assemblies, failed downloads, truncated searches, and unexamined members remain explicit. It cannot assign universal strain identity or establish that no closer unexamined isolate exists. Direct-read SKA2 refinement is supported, but the initial Mashpit screen still requires an assembly. See [benchmarking.md](benchmarking.md) for the limits of software verification versus biological validation.

The legacy `databases-v1` packages do not bundle full membership tables;
their expansion still depends on NCBI retaining the matching release. The
published `databases-v2` packages include checked copies of
those tables. This preserves cluster membership but not the downloaded genome
assemblies used in each SNP run, which must also be kept for exact replay.
