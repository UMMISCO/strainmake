Several Python scripts in the pipeline have unit tests.
These tests use the [unittest](https://docs.python.org/3/library/unittest.html) library.
To run them, change directory to the root of the git repository and run:

```sh
python3 -m unittest discover -s workflow/scripts/test/
```

Needed software:

- Prodigal