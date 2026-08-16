Tests for the `strainmake` CLI, written with [pytest](https://docs.pytest.org).
They use temporary files only, and run in a few seconds.

To run them, change directory to the root of the git repository and run:

```sh
pip install -e ".[dev]"
pytest
```

| File | Covers |
| :--- | :--- |
| [`test_cli.py`](test_cli.py) | entry point: `--help`, `--version`, sub-commands |
| [`test_init.py`](test_init.py) | `strainmake init` |
| [`test_build.py`](test_build.py) | `strainmake build`, i.e. Snakefile generation |
| [`test_run.py`](test_run.py) | `strainmake run` () including `--offline` and `--create-envs-only`) |
| [`test_databases.py`](test_databases.py) | where reference databases live, and reading the config |
| [`test_prepare.py`](test_prepare.py) | `strainmake prepare` sub-commands |
| [`test_report.py`](test_report.py) | `strainmake report` |
| [`test_unit.py`](test_unit.py) | helper scripts, and Snakemake dry-runs of the workflow |

The dry-run tests need `snakemake` on the `PATH`. Without it they are skipped
rather than failed, so check the `skipped` count in the summary line if you
expect them to run.