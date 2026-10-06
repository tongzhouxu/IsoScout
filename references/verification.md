# Software verification

Run the unit and mocked workflow tests from the repository root:

```bash
PYTHONPYCACHEPREFIX=/tmp/isoscout_pycache python3 -m unittest discover -s tests -v
```

Optional offline integration tests exercise pinned SKA2, minimap2/paftools and
MUMmer executables on seeded DNA with known differences, indels and reverse
complements. Cache tests verify incremental work, invalidation and failures. Run it in the configured container:

```bash
ISOSCOUT_TOOL_TESTS=1 python3 -m unittest discover -s tests -v
```

These tests check software behavior, including selection budgets, focused
expansion, ties, exclusions, and assembly/read distance parsing. They do not
measure strain-assignment accuracy on real isolates. A cluster label match or
stable result among examined references cannot establish recovery of the
globally nearest genome.
