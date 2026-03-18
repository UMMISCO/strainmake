"""
Tests for `strainmake report` - MultiQC report generation.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from strainmake.cli.main import app
from strainmake.cli.report_cmd import _DEFAULT_MULTIQC_CONFIG, _load_handle_results

runner = CliRunner()

_PATCH = "strainmake.cli.report_cmd._load_handle_results"


def _make_dirs(tmp_path):
    results = tmp_path / "results"
    logs = tmp_path / "logs"
    output = tmp_path / "output"
    results.mkdir()
    logs.mkdir()
    return results, logs, output


class TestDefaultMultiqcConfig:
    def test_default_config_path_exists(self):
        """The default multiqc_config.yaml should exist at the expected path."""
        assert _DEFAULT_MULTIQC_CONFIG.exists(), (
            f"Default multiqc_config.yaml not found at {_DEFAULT_MULTIQC_CONFIG}"
        )

    def test_load_handle_results_raises_on_missing_module(self):
        """_load_handle_results should propagate ImportError when the module is absent."""
        import sys
        saved = sys.modules.pop(
            "workflow.scripts.multiqc_results.generate_multiqc_report", None
        )
        try:
            # Ensure the import actually fails by removing any cached module
            with patch.dict("sys.modules", {"workflow.scripts.multiqc_results.generate_multiqc_report": None}):
                with pytest.raises((ImportError, ModuleNotFoundError)):
                    _load_handle_results()
        finally:
            if saved is not None:
                sys.modules["workflow.scripts.multiqc_results.generate_multiqc_report"] = saved


class TestReportCommand:
    def _args(self, results, logs, output, extra=None):
        """Helper to construct the base CLI arguments for the report command, with optional extra flags."""
        base = [
            "report",
            "--results", str(results),
            "--logs", str(logs),
            "--output", str(output),
        ]
        return base + (extra or [])

    def test_calls_handle_results(self, tmp_path):
        """Running `strainmake report` should call the handle_results functions."""
        results, logs, output = _make_dirs(tmp_path)
        mock_handler = MagicMock()
        mock_class = MagicMock(return_value=mock_handler)

        with patch(_PATCH, return_value=mock_class):
            result = runner.invoke(app, self._args(results, logs, output))

        assert result.exit_code == 0, result.output
        mock_class.assert_called_once()
        mock_handler.collect_all_results.assert_called_once()
        mock_handler.prepare_results.assert_called_once()
        mock_handler.generate_report.assert_called_once_with(dry_run=False)

    def test_dry_run_forwarded(self, tmp_path):
        """The --dry-run flag should be forwarded to the report handler's generate_report method."""
        results, logs, output = _make_dirs(tmp_path)
        mock_handler = MagicMock()
        mock_class = MagicMock(return_value=mock_handler)

        with patch(_PATCH, return_value=mock_class):
            runner.invoke(app, self._args(results, logs, output, ["--dry-run"]))

        mock_handler.generate_report.assert_called_once_with(dry_run=True)

    def test_ani_forwarded_to_constructor(self, tmp_path):
        """The --ani option should be forwarded to the report handler's constructor."""
        results, logs, output = _make_dirs(tmp_path)
        mock_handler = MagicMock()
        mock_class = MagicMock(return_value=mock_handler)

        with patch(_PATCH, return_value=mock_class):
            runner.invoke(app, self._args(results, logs, output, ["--ani", "97"]))

        call_kwargs = mock_class.call_args[1]
        assert call_kwargs.get("ani") == 97

    def test_output_dir_created_if_missing(self, tmp_path):
        """The output directory should be created if it does not exist."""
        results, logs, _ = _make_dirs(tmp_path)
        output = tmp_path / "new" / "deep" / "dir"

        with patch(_PATCH, return_value=MagicMock(return_value=MagicMock())):
            runner.invoke(app, self._args(results, logs, output))

        assert output.exists()

    def test_missing_results_dir_fails(self, tmp_path):
        """If the results directory does not exist, the command should fail with a non-zero exit code."""
        logs = tmp_path / "logs"
        logs.mkdir()
        result = runner.invoke(
            app,
            [
                "report",
                "--results", str(tmp_path / "nonexistent"),
                "--logs", str(logs),
                "--output", str(tmp_path / "out"),
            ],
        )
        assert result.exit_code != 0

    def test_import_error_gives_friendly_message(self, tmp_path):
        """If an ImportError occurs, the command should fail with a friendly message."""
        results, logs, output = _make_dirs(tmp_path)

        with patch(_PATCH, side_effect=ImportError("multiqc not found")):
            result = runner.invoke(app, self._args(results, logs, output))

        assert result.exit_code != 0
        assert result.exit_code != 0
