# Software verification

Build the CLI source distribution and wheel, then test the installed artifact:

```bash
python -m pip install build
python -m build
python -m venv /tmp/isoscout-wheel-check
/tmp/isoscout-wheel-check/bin/python -m pip install --no-deps dist/*.whl
/tmp/isoscout-wheel-check/bin/python tests/check_installed_cli.py -v
```

Use a fresh test environment. Three focused checks run outside the checkout:
installed package contents and startup, one native QC-stop run compared with the
source entry point, and Docker argument/provenance forwarding and exit status
using a fake Docker executable. Input-validation edge cases remain in the
existing workflow suite. The wheel checks do not require Docker or bioinformatics
dependencies and do not establish native biological-tool or ARM compatibility.
CI runs them on Python 3.9 and 3.11 before allowing container publication.

The wheel maps `scripts/` to the private `isoscout._workflow` package and
`config/` to `isoscout.config`; there is only one source for each implementation
and policy. The public CLI runs those scripts in a separate process, preserving
their existing imports and relative policy paths. Treat the private modules as
implementation details, not a stable Python API. Source script paths remain
supported for the skill and existing container entry point.

CLI package version `1.5.0.dev1` identifies the unreleased packaging work.
Scientific workflow and policy versions remain independently recorded in
`config/*.json`. Changing packaging does not validate a new tool environment.

For a quick host check, run the focused unit and mocked workflow suite from the repository root:

```bash
PYTHONPYCACHEPREFIX=/tmp/isoscout_pycache python3 -m unittest discover -s tests -v
```

The canonical release check builds the current image and runs the same suite plus
the three offline tool integrations on seeded DNA with known differences, indels,
and reverse complements:

```bash
docker build --platform linux/amd64 \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  -f container/Dockerfile -t isoscout:test .
docker run --rm --platform linux/amd64 \
  --entrypoint /usr/local/bin/_entrypoint.sh \
  --env ISOSCOUT_TOOL_TESTS=1 \
  isoscout:test python3 -m unittest discover -s /opt/isoscout/tests -v
```

The suite keeps one representative case for each important boundary: input and
database rejection, routing, fixed invocation, candidate/tie interpretation,
selection budgets, focused expansion, cache behavior, and assembly/read distance
parsing. It does not measure strain-assignment accuracy on real isolates. A
cluster label match or stable result among examined references cannot establish
recovery of the globally nearest genome.
