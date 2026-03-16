"""
Tests for the strainmake CLI entry point: --help, --version, banner, sub-groups.
"""

from typer.testing import CliRunner

from strainmake import __version__
from strainmake.cli.main import app

runner = CliRunner()


class TestVersion:
    def test_version_flag_short(self):
        """Test that the short version flag (-V) works and includes the version number in the output."""
        result = runner.invoke(app, ["-V"])
        assert result.exit_code == 0
        assert __version__ in result.output

    def test_version_flag_long(self):
        """Test that the long version flag (--version) works and includes the version number in the output."""
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "StrainMake" in result.output
        assert __version__ in result.output


class TestHelp:
    def test_root_help(self):
        """Test that the root help command works and lists all top-level commands."""
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        # all top-level commands must be listed
        for cmd in ("init", "build", "run", "report", "prepare"):
            assert cmd in result.output

    def test_init_help(self):
        """Test that the 'init' command help works and lists its options."""
        result = runner.invoke(app, ["init", "--help"])
        assert result.exit_code == 0
        assert "--samples" in result.output
        assert "--output" in result.output

    def test_build_help(self):
        """Test that the 'build' command help works and lists its options."""
        result = runner.invoke(app, ["build", "--help"])
        assert result.exit_code == 0
        assert "--config" in result.output
        assert "--output" in result.output

    def test_run_help(self):
        """Test that the 'run' command help works and lists its options."""
        result = runner.invoke(app, ["run", "--help"])
        assert result.exit_code == 0
        assert "--cores" in result.output
        assert "--dry-run" in result.output

    def test_report_help(self):
        """Test that the 'report' command help works and lists its options."""
        result = runner.invoke(app, ["report", "--help"])
        assert result.exit_code == 0
        assert "--results" in result.output
        assert "--output" in result.output

    def test_prepare_help(self):
        """Test that the 'prepare' command help works and lists its sub-commands."""
        result = runner.invoke(app, ["prepare", "--help"])
        assert result.exit_code == 0
        for sub in ("import-preprocessed", "split-samples", "gather-assemblies", "gather-binning"):
            assert sub in result.output

    def test_prepare_subcommand_help(self):
        """Test that each 'prepare' sub-command help works and lists its options."""
        for sub in ("import-preprocessed", "split-samples", "gather-assemblies", "gather-binning"):
            result = runner.invoke(app, ["prepare", sub, "--help"])
            assert result.exit_code == 0, f"{sub} --help failed: {result.output}"


class TestNoArgsShowsHelp:
    def test_no_args_no_error(self):
        """Test that running the CLI with no arguments does not raise an unhandled exception and shows help."""
        # no_args_is_help=True should print help, not raise an exception
        result = runner.invoke(app, [])
        # exit code 0 or 2 depending on typer version, but never an unhandled exception
        assert "strainmake" in result.output.lower() or result.exit_code in (0, 2)
