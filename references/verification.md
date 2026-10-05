# Software verification

Run the unit and mocked workflow tests from the repository root:

```bash
PYTHONPYCACHEPREFIX=/tmp/isoscout_pycache python3 -m unittest discover -s tests -v
```

The optional offline integration test uses the pinned SKA2 executable on seeded
DNA with known 1-SNP and 3-SNP differences. Run it in the configured container:

```bash
ISOSCOUT_TOOL_TESTS=1 python3 -m unittest discover -s tests -v
```

These tests check software behavior, including selection budgets, focused
expansion, ties, exclusions, and assembly/read distance parsing. They do not
measure strain-assignment accuracy on real isolates. A cluster label match or
stable result among examined references cannot establish recovery of the
globally nearest genome.
