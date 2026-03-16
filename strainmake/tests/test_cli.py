"""
Tests for the strainmake CLI entry point: --help, --version, banner, sub-groups.
"""

import re

from typer.testing import CliRunner

from strainmake import __version__
from strainmake.cli.main import app

runner = CliRunner()


def _invoke(args: list[str]):
    """Invoke CLI with deterministic terminal settings for CI-stable help output."""
    return runner.invoke(
        app,
        args,
        color=False,
        env={"NO_COLOR": "1", "TERM": "dumb", "COLUMNS": "200", "LINES": "40"},
    )


def _norm(text: str) -> str:
    """Normalize help output by stripping ANSI escapes and collapsing whitespace."""
    ansi = re.compile(r"\x1B\[[0-?]*[ -/]*[@-~]")
    text = ansi.sub("", text)
    return re.sub(r"\s+", " ", text).strip()


class TestVersion:
    def test_version_flag_short(self):
        """Test that the short version flag (-V) works and includes the version number in the output."""
        result = _invoke(["-V"])
        assert result.exit_code == 0
        assert __version__ in _norm(result.output)

    def test_version_flag_long(self):
        """Test that the long version flag (--version) works and includes the version number in the output."""
        result = _invoke(["--version"])
        out = _norm(result.output)
        assert result.exit_code == 0
        assert "StrainMake" in out
        assert __version__ in out


class TestHelp:
    def test_root_help(self):
        """Test that the root help command works and lists all top-level commands."""
        result = _invoke(["--help"])
        out = _norm(result.output)
        assert result.exit_code == 0
        # all top-level commands must be listed
        for cmd in ("init", "build", "run", "report", "prepare"):
            assert cmd in out

    def test_init_help(self):
        """Test that the 'init' command help works and lists its options."""
        result = _invoke(["init", "--help"])
        out = _norm(result.output)
        assert result.exit_code == 0
        assert "strainmake init" in out
        assert ("--samples" in out) or ("-s" in out)
        assert ("--output" in out) or ("-o" in out)

    def test_build_help(self):
        """Test that the 'build' command help works and lists its options."""
        result = _invoke(["build", "--help"])
        out = _norm(result.output)
        assert result.exit_code == 0
        assert "strainmake build" in out
        assert ("--config" in out) or ("-c" in out)
        assert ("--output" in out) or ("-o" in out)

    def test_run_help(self):
        """Test that the 'run' command help works and lists its options."""
        result = _invoke(["run", "--help"])
        out = _norm(result.output)
        assert result.exit_code == 0
        assert "strainmake run" in out
        assert ("--cores" in out) or ("-c" in out)
        assert ("--dry-run" in out) or ("-n" in out)

    def test_report_help(self):
        """Test that the 'report' command help works and lists its options."""
        result = _invoke(["report", "--help"])
        out = _norm(result.output)
        assert result.exit_code == 0
        assert "strainmake report" in out
        assert ("--results" in out) or ("-r" in out)
        assert ("--output" in out) or ("-o" in out)

    def test_prepare_help(self):
        """Test that the 'prepare' command help works and lists its sub-commands."""
        result = _invoke(["prepare", "--help"])
        out = _norm(result.output)
        assert result.exit_code == 0
        for sub in ("import-preprocessed", "split-samples", "gather-assemblies", "gather-binning"):
            assert sub in out

    def test_prepare_subcommand_help(self):
        """Test that each 'prepare' sub-command help works and lists its options."""
        for sub in ("import-preprocessed", "split-samples", "gather-assemblies", "gather-binning"):
            result = _invoke(["prepare", sub, "--help"])
            assert result.exit_code == 0, f"{sub} --help failed: {_norm(result.output)}"


class TestNoArgsShowsHelp:
    def test_no_args_no_error(self):
        """Test that running the CLI with no arguments does not raise an unhandled exception and shows help."""
        # no_args_is_help=True should print help, not raise an exception
        result = _invoke([])
        out = _norm(result.output).lower()
        # exit code 0 or 2 depending on typer version, but never an unhandled exception
        assert "strainmake" in out or result.exit_code in (0, 2)
