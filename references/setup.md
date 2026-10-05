# Setup

Run analyses in the pinned container built from `container/Dockerfile`. Mashpit is installed from upstream commit `538d3421302fe6dd129780605b8ff5dedbf4c046c`, not from the older published PyPI artifact. Runtime package installation is not part of an analysis.

Pull the published image:

```bash
# The package is public; no GHCR login is needed to pull it.
docker pull --platform linux/amd64 ghcr.io/tongzhouxu/isoscout:latest
docker tag ghcr.io/tongzhouxu/isoscout:latest isoscout:local
```

The published image (2026-10-01) predates this checkout's adaptive selection,
cluster expansion, and bundled membership support. A clean build from
`container/Dockerfile` currently fails at conda dependency resolution. For
running the updated checkout with the existing pinned tools, run from the
repository root and add `--volume "$PWD:/opt/isoscout:ro"` to `docker run`.
An earlier checkout was verified with its then-current 89 tests, including the
real SKA2 integration test, on 2026-10-02. Retain the checkout commit with the results; the configured
container digest identifies the baseline image rather than the mounted code.

`--platform linux/amd64` is required on Apple Silicon (`quast=5.3.0` has no `linux/arm64` build for the pinned Python 3.11); it runs fine under emulation there.

Run an assembly with a known organism:

```bash
docker run --rm --platform linux/amd64 \
  --volume "/absolute/path/to/data:/data:ro" \
  --volume "$HOME/.isoscout/databases:/databases:ro" \
  --volume "/absolute/path/to/results:/results" \
  isoscout:local \
  /data/sample.fasta --organism salmonella \
  --database-root /databases --output /results/sample-screen
```

Omit `--organism` to use local MLST routing. The output directory must not already exist.

Set `ISOSCOUT_DATABASE_ROOT` or pass `--database-root`. The directory must contain:

```text
databases/
├── salmonella/
├── ecoli_shigella/
├── listeria/
├── campylobacter/
└── cronobacter/
```

Pre-built, checksummed databases for all five are published at
[databases-v2](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-v2). This example downloads one
organism; set `org` to the database key you need:

```bash
mkdir -p ~/.isoscout/databases && cd ~/.isoscout/databases
org=salmonella
curl -fLO "https://github.com/tongzhouxu/IsoScout/releases/download/databases-v2/${org}.tar.gz"
curl -fLO https://github.com/tongzhouxu/IsoScout/releases/download/databases-v2/checksums.sha256.txt
grep "  ${org}.tar.gz$" checksums.sha256.txt | shasum -a 256 -c -
tar -xzf "${org}.tar.gz"
```

Each organism directory must contain the Mashpit `<name>.db` and `<name>.sig` files plus `database.json`:

```json
{
  "name": "salmonella",
  "version": "PDG-build-identifier",
  "build_date": "2026-08-01",
  "checksum": "sha256-of-primary-artifact",
  "checksum_file": "salmonella.db",
  "source": "NCBI Pathogen Detection via Mashpit build"
}
```

The wrapper checks the name, required metadata, database and signature files,
and their recorded checksums. When a package includes membership tables, it
also checks their release, paths, sizes, and checksums. It never downloads,
builds, updates, or substitutes a database during a screen.

### Database packages with local cluster membership

The legacy [databases-v1](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-v1)
release predates member expansion. Its five packages can still screen isolates
and can expand clusters while their exact NCBI releases remain available. On
2026-10-02, all five `databases-v2` archives were built, verified, and
published. Their exact NCBI versions are:

| Organism | NCBI Pathogen Detection release |
| --- | --- |
| Salmonella | `PDG000000002.3864` |
| E. coli/Shigella | `PDG000000004.5775` |
| Listeria | `PDG000000001.4519` |
| Campylobacter | `PDG000000003.2682` |
| Cronobacter | `PDG000000043.415` |

To build another durable package, run
`scripts/package_database_release.py` against an extracted database. The script
fetches that database's exact-release metadata and membership TSVs, checks
every existing representative against both tables, compresses the files, and
writes a new package and optional tar archive. It leaves the input database
untouched:

```bash
python3 scripts/package_database_release.py \
  --database-dir /path/to/extracted/salmonella \
  --output-dir /path/to/new-release/salmonella \
  --archive /path/to/new-release/salmonella.tar.gz
```

The package's `database.json` gains a `membership_snapshot` object containing
the PDG version, NCBI species directory, relative filenames, source URLs,
original and compressed sizes and SHA-256 checksums. It also records the
Mashpit `.sig` checksum. The `.db` and `.sig`
remain the same. IsoScout verifies packaged membership files before screening
and reads them without contacting NCBI. Reference genome assemblies used for
SKA2 still need to be retrieved and retained with the run. To package already
saved exact-release TSVs, pass both `--metadata-tsv` and
`--all-isolates-tsv`; this mode records their provenance as supplied local files.
Choose a new output directory and keep the original `databases-v1` assets.
For a multi-organism release, create one archive per organism and publish a
SHA-256 list for the completed archives. Verify it before extracting:

```bash
cd /path/to/new-release
shasum -a 256 -c checksums.sha256.txt
```

Only the selected organism database is required for a run — each is independent, so download just the one(s) you need. Local `mlst` 2.35.0 and its bundled PubMLST schemes are installed in the image; no scheme download or sequence upload occurs at runtime.

The `container_digest` in `config/workflow.json` still identifies the published
baseline image. It does not identify the code in a bind-mounted checkout. For
exact replay, retain that checkout commit or source snapshot, the database
archive/checksums, downloaded reference assemblies, and the run output. Update
the digest and publish a new image when a clean image build is available.
