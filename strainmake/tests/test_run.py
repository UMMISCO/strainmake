"""
Tests for `strainmake run` - Snakemake wrapper.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from strainmake import REPO_ROOT
from strainmake.cli.main import app
from strainmake.cli.run_cmd import _DEFAULT_SNAKEFILE, _PREPROCESSED_SNAKEFILE

runner = CliRunner()


def _invoke_run(args: list, returncode: int = 0):
    """Invoke the run command with subprocess.run mocked."""
    mock_result = MagicMock()
    mock_result.returncode = returncode
    with patch("strainmake.cli.run_cmd.subprocess.run", return_value=mock_result) as mock_sub:
        result = runner.invoke(app, ["run"] + args, catch_exceptions=False)
    return result, mock_sub


class TestDefaultSnakefile:
    def test_default_snakefile_points_to_workflow(self):
        """The default Snakefile should be the one in the workflow directory."""
        assert _DEFAULT_SNAKEFILE == REPO_ROOT / "workflow" / "Snakefile"

    def test_preprocessed_snakefile_exists(self):
        """The Snakefile for already preprocessed sequences should exist."""
        assert _PREPROCESSED_SNAKEFILE == REPO_ROOT / "workflow" / "Snakefile.already_preprocessed_seq.smk"


class TestRunCommand:
    def test_calls_snakemake(self):
        """Running `strainmake run` should call Snakemake."""
        result, mock_sub = _invoke_run(["--cores", "2"])
        assert mock_sub.called
        cmd = mock_sub.call_args[0][0]
        assert "snakemake" in cmd

    def test_use_conda_always_present(self):
        """The `--use-conda` flag should always be included when calling Snakemake, regardless of other options."""
        _, mock_sub = _invoke_run(["--cores", "4"])
        cmd = mock_sub.call_args[0][0]
        assert "--use-conda" in cmd

    def test_cores_forwarded(self):
        """The number of cores specified with `--cores` should be forwarded to Snakemake."""
        _, mock_sub = _invoke_run(["--cores", "16"])
        cmd = mock_sub.call_args[0][0]
        assert "--cores" in cmd
        assert "16" in cmd

    def test_dry_run_flag(self):
        """The `--dry-run` flag should be forwarded to Snakemake when specified."""
        _, mock_sub = _invoke_run(["--dry-run"])
        cmd = mock_sub.call_args[0][0]
        assert "--dry-run" in cmd

    def test_configfile_forwarded(self, tmp_path):
        """The `--configfile` option should be forwarded to Snakemake with the correct path."""
        cfg = tmp_path / "config.yaml"
        cfg.write_text("samples: data/metadata.tsv\n")
        _, mock_sub = _invoke_run(["--configfile", str(cfg)])
        cmd = mock_sub.call_args[0][0]
        assert "--configfile" in cmd
        assert str(cfg) in cmd

    def test_already_preprocessed_uses_correct_snakefile(self):
        """The `--already-preprocessed` flag should use the correct Snakefile."""
        _, mock_sub = _invoke_run(["--already-preprocessed"])
        cmd = mock_sub.call_args[0][0]
        assert "already_preprocessed" in " ".join(cmd)

    def test_extra_args_forwarded(self):
        """Extra arguments after `--` should be forwarded to Snakemake."""
        _, mock_sub = _invoke_run(["--cores", "4", "--", "--rerun-incomplete"])
        cmd = mock_sub.call_args[0][0]
        assert "--rerun-incomplete" in cmd

    def test_exit_code_propagated(self):
        """If the Snakemake subprocess returns a non-zero exit code, `strainmake run` should also exit with a non-zero code."""
        result, _ = _invoke_run([], returncode=1)
        assert result.exit_code == 1

    def test_default_cores_is_4(self):
        """If no `--cores` option is provided, the default should be 4 cores."""
        _, mock_sub = _invoke_run([])
        cmd = mock_sub.call_args[0][0]
        idx = cmd.index("--cores")
        assert cmd[idx + 1] == "4"
