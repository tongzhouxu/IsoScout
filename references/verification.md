# Software verification

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
