"""
Tests for `strainmake init` - configuration YAML generation.
"""

from pathlib import Path
from unittest.mock import patch

import pytest
import yaml
from typer.testing import CliRunner

from strainmake.cli.main import app
from strainmake.cli.init_cmd import (
    SECTION_KEYS,
    _get_section,
    _prompt_for_section,
    _prompt_recursive,
    load_defaults,
)

runner = CliRunner()

# Simulate the user pressing Enter (accept defaults) for each section.
_ACCEPT_ALL_INPUT = "\n".join(["y"] * len(SECTION_KEYS)) + "\n"


class TestLoadDefaults:
    def test_returns_dict(self):
        """load_defaults should return a dictionary parsed from the template YAML."""
        defaults = load_defaults()
        assert isinstance(defaults, dict)

    def test_contains_expected_top_level_keys(self):
        """The loaded defaults should contain expected top-level keys like 'fastp', 'assembly', and 'binning'."""
        defaults = load_defaults()
        for key in ("fastp", "assembly", "binning"):
            assert key in defaults, f"Expected '{key}' in defaults"

    def test_raises_if_template_missing(self, tmp_path):
        """If the template YAML is missing, load_defaults should raise a FileNotFoundError."""
        from strainmake.cli import init_cmd as _mod
        original = _mod.TEMPLATE_CONFIG
        _mod.TEMPLATE_CONFIG = tmp_path / "nonexistent.yaml"
        try:
            with pytest.raises(FileNotFoundError):
                _mod.load_defaults()
        finally:
            _mod.TEMPLATE_CONFIG = original


class TestGetSection:
    def test_extracts_correct_keys(self):
        """_get_section should extract the specified keys from the defaults dictionary."""
        defaults = {"fastp": {"threads": 4}, "binning": {}, "other": 1}
        result = _get_section(defaults, ["fastp", "binning"])
        assert set(result.keys()) == {"fastp", "binning"}

    def test_ignores_missing_keys(self):
        """_get_section should ignore keys that are not present in the defaults dictionary."""
        result = _get_section({"a": 1}, ["a", "b"])
        assert "b" not in result


class TestPromptRecursive:
    def test_accepts_default_for_scalar(self):
        """When the user presses Enter without typing anything, the default scalar value should be returned."""
        with patch("builtins.input", return_value=""):
            result = _prompt_recursive(42, prefix="threads")
        assert result == 42

    def test_overrides_scalar(self):
        """When the user types a value, it should override the default scalar value."""
        with patch("builtins.input", return_value="8"):
            result = _prompt_recursive(4, prefix="threads")
        assert result == 8

    def test_accepts_default_for_list(self):
        """When the user presses Enter without typing anything, the default list value should be returned."""
        with patch("builtins.input", return_value=""):
            result = _prompt_recursive(["megahit"], prefix="assembler")
        assert result == ["megahit"]

    def test_overrides_list(self):
        """When the user types a comma-separated list, it should override the default list value."""
        with patch("builtins.input", return_value="megahit metaspades"):
            result = _prompt_recursive(["megahit"], prefix="assembler")
        assert result == ["megahit", "metaspades"]

    def test_nested_dict(self):
        """When the user presses Enter without typing anything, the default nested dictionary should be returned."""
        data = {"threads": 4, "mem": 16}
        # Accept defaults for all nested fields
        with patch("builtins.input", return_value=""):
            result = _prompt_recursive(data)
        assert result == data


class TestPromptForSection:
    def test_accept_defaults_with_y(self):
        """When the user inputs 'y', the default data for that section should be returned."""
        data = {"threads": 4}
        with patch("builtins.input", return_value="y"):
            result = _prompt_for_section("ASSEMBLY", data)
        assert result == data

    def test_accept_defaults_with_empty(self):
        """When the user presses Enter without typing anything, the default data for that section should be returned."""
        data = {"threads": 4}
        with patch("builtins.input", return_value=""):
            result = _prompt_for_section("ASSEMBLY", data)
        assert result == data

    def test_skip_with_x(self):
        """When the user inputs 'x', the section should be skipped and None returned."""
        data = {"threads": 4}
        with patch("builtins.input", return_value="x"):
            result = _prompt_for_section("ASSEMBLY", data)
        assert result is None


class TestInitCommand:
    def test_creates_output_yaml(self, tmp_path, sample_metadata_tsv):
        """Running `strainmake init` with valid inputs should create the output YAML file."""
        output = tmp_path / "out" / "config.yaml"
        result = runner.invoke(
            app,
            ["init", "--samples", str(sample_metadata_tsv), "--output", str(output)],
            input=_ACCEPT_ALL_INPUT,
        )
        assert result.exit_code == 0, result.output
        assert output.exists()

    def test_output_contains_samples_key(self, tmp_path, sample_metadata_tsv):
        """The output YAML should contain the 'samples' key with the correct value."""
        output = tmp_path / "config.yaml"
        runner.invoke(
            app,
            ["init", "--samples", str(sample_metadata_tsv), "--output", str(output)],
            input=_ACCEPT_ALL_INPUT,
        )
        cfg = yaml.safe_load(output.read_text())
        assert "samples" in cfg
        assert cfg["samples"] == str(sample_metadata_tsv)

    def test_output_contains_lr_seq_format(self, tmp_path, sample_metadata_tsv):
        """The output YAML should contain the 'lr_seq_format' key with the correct value."""
        output = tmp_path / "config.yaml"
        runner.invoke(
            app,
            [
                "init",
                "--samples", str(sample_metadata_tsv),
                "--lr-format", "fasta",
                "--output", str(output),
            ],
            input=_ACCEPT_ALL_INPUT,
        )
        cfg = yaml.safe_load(output.read_text())
        assert cfg["lr_seq_format"] == "fasta"

    def test_skipped_section_absent_from_output(self, tmp_path, sample_metadata_tsv):
        """When a section is skipped, its keys should be absent from the output YAML."""
        output = tmp_path / "config.yaml"
        # Sections: PREPROCESSING, ASSEMBLY, GENE_CATALOG, BINNING, TAXO_PROFILING, STRAINS_PROFILING
        # "x" on the 3rd prompt skips GENE_CATALOG (keys: mmseqs2, representative_genes)
        user_input = "y\ny\nx\ny\ny\ny\n"
        result = runner.invoke(
            app,
            ["init", "--samples", str(sample_metadata_tsv), "--output", str(output)],
            input=user_input,
        )
        assert result.exit_code == 0, result.output
        cfg = yaml.safe_load(output.read_text())
        # GENE_CATALOG keys (3rd section) should not appear
        for key in SECTION_KEYS["GENE_CATALOG"]:
            assert key not in cfg, f"Expected '{key}' to be absent (section was skipped)"

    def test_output_contains_lr_technology(self, tmp_path, sample_metadata_tsv):
        """lr_technology is part of PREPROCESSING and must appear in the output."""
        output = tmp_path / "config.yaml"
        runner.invoke(
            app,
            ["init", "--samples", str(sample_metadata_tsv), "--output", str(output)],
            input=_ACCEPT_ALL_INPUT,
        )
        cfg = yaml.safe_load(output.read_text())
        assert "lr_technology" in cfg

    def test_output_contains_strains_profiling(self, tmp_path, sample_metadata_tsv):
        """strains_profiling section must appear when STRAINS_PROFILING is accepted."""
        output = tmp_path / "config.yaml"
        runner.invoke(
            app,
            ["init", "--samples", str(sample_metadata_tsv), "--output", str(output)],
            input=_ACCEPT_ALL_INPUT,
        )
        cfg = yaml.safe_load(output.read_text())
        assert "strains_profiling" in cfg

    def test_missing_samples_file_fails(self, tmp_path):
        """If the provided samples TSV file does not exist, the command should fail with a non-zero exit code."""
        result = runner.invoke(
            app,
            ["init", "--samples", str(tmp_path / "nonexistent.tsv"), "--output", "config.yaml"],
        )
        assert result.exit_code != 0
