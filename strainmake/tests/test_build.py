"""
Tests for `strainmake build` - Snakefile generation from config YAML.
"""

import pytest
from typer.testing import CliRunner

from strainmake.cli.main import app
from strainmake.cli.build_cmd import _TEMPLATE_FILE

runner = CliRunner()


class TestTemplateFile:
    def test_template_exists(self):
        """The Snakefile.j2 template must exist at the expected path."""
        assert _TEMPLATE_FILE.exists(), (
            f"Snakefile.j2 template not found at {_TEMPLATE_FILE}"
        )


class TestBuildCommand:
    def test_creates_snakefile(self, tmp_path, sample_config_yaml):
        """Running `strainmake build` should create the specified Snakefile output."""
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(sample_config_yaml), "--output", str(output)],
        )
        assert result.exit_code == 0, result.output
        assert output.exists()

    def test_snakefile_contains_configfile_directive(self, tmp_path, sample_config_yaml):
        """The rendered Snakefile should contain a reference to the config YAML (for Snakemake's --configfile)."""
        output = tmp_path / "Snakefile"
        runner.invoke(
            app,
            ["build", "--config", str(sample_config_yaml), "--output", str(output)],
        )
        content = output.read_text()
        assert "configfile:" in content

    def test_snakefile_contains_config_path(self, tmp_path, sample_config_yaml):
        """The rendered Snakefile should contain the path to the config YAML (for debugging)."""
        output = tmp_path / "Snakefile"
        runner.invoke(
            app,
            ["build", "--config", str(sample_config_yaml), "--output", str(output)],
        )
        content = output.read_text()
        assert str(sample_config_yaml) in content

    def test_snakefile_contains_rule_all(self, tmp_path, sample_config_yaml):
        """The rendered Snakefile should contain a rule all."""
        output = tmp_path / "Snakefile"
        runner.invoke(
            app,
            ["build", "--config", str(sample_config_yaml), "--output", str(output)],
        )
        assert "rule all:" in output.read_text()

    def test_snakefile_contains_generation_timestamp(self, tmp_path, sample_config_yaml):
        """The rendered Snakefile should contain a generation timestamp."""
        output = tmp_path / "Snakefile"
        runner.invoke(
            app,
            ["build", "--config", str(sample_config_yaml), "--output", str(output)],
        )
        # Template renders '# generated on <ISO datetime>'
        assert "generated on" in output.read_text()

    def test_output_in_subdirectory_is_created(self, tmp_path, sample_config_yaml):
        """The output path can include subdirectories that don't exist yet - they should be created automatically."""
        output = tmp_path / "nested" / "dir" / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(sample_config_yaml), "--output", str(output)],
        )
        assert result.exit_code == 0, result.output
        assert output.exists()

    def test_missing_config_fails(self, tmp_path):
        """Running `strainmake build` with a missing config file should fail."""
        result = runner.invoke(
            app,
            ["build", "--config", str(tmp_path / "nonexistent.yaml"), "--output", "Snakefile"],
        )
        assert result.exit_code != 0

    def test_missing_template_fails(self, tmp_path, sample_config_yaml, monkeypatch):
        """Graceful error when template file is absent (e.g. wrong install)."""
        import strainmake.cli.build_cmd as _mod
        from pathlib import Path

        original = _mod._TEMPLATE_FILE
        _mod._TEMPLATE_FILE = tmp_path / "nonexistent.j2"
        try:
            result = runner.invoke(
                app,
                ["build", "--config", str(sample_config_yaml), "--output", str(tmp_path / "Snakefile")],
            )
            assert result.exit_code != 0
        finally:
            _mod._TEMPLATE_FILE = original


class TestPerformOption:
    """Tests for the --perform step-filtering option."""

    def _build(self, config_yaml, tmp_path, *extra_args):
        """Helper to run `strainmake build` with given config and extra args, returning the rendered Snakefile content for inspection."""
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(config_yaml), "--output", str(output)] + list(extra_args),
        )
        assert result.exit_code == 0, result.output
        return output.read_text()

    def test_no_perform_includes_assembly_when_in_config(self, tmp_path, sample_config_yaml):
        """
        When --perform is not given, the presence of assembly steps should be determined by 
        the config YAML (notably the presence of assembly samples).
        """
        content = self._build(sample_config_yaml, tmp_path)
        assert "03_assembly" in content

    def test_perform_assembly_includes_assembly(self, tmp_path, full_config_yaml):
        """When --perform assembly is given, assembly steps should be included regardless of config content."""
        content = self._build(full_config_yaml, tmp_path, "--perform", "assembly")
        assert "03_assembly" in content

    def test_perform_assembly_excludes_binning(self, tmp_path, full_config_yaml):
        """When --perform assembly is given, binning steps should be excluded even if present in the config."""
        content = self._build(full_config_yaml, tmp_path, "--perform", "assembly")
        assert "results/05_binning" not in content

    def test_perform_binning_excludes_assembly(self, tmp_path, full_config_yaml):
        """When --perform binning is given, assembly steps should be excluded even if present in the config."""
        content = self._build(full_config_yaml, tmp_path, "--perform", "binning")
        assert "results/03_assembly" not in content
        assert "results/05_binning" in content

    def test_perform_preprocessing_includes_qc(self, tmp_path, full_config_yaml):
        """ When --perform preprocessing is given, QC and preprocessing steps should be included regardless of config content. """
        content = self._build(full_config_yaml, tmp_path, "--perform", "preprocessing")
        assert "results/01_qc" in content
        assert "results/02_preprocess" in content

    def test_perform_assembly_excludes_qc(self, tmp_path, full_config_yaml):
        """When --perform is given without preprocessing, QC lines are absent."""
        content = self._build(full_config_yaml, tmp_path, "--perform", "assembly")
        assert "results/01_qc" not in content

    def test_perform_taxo_profiling(self, tmp_path, full_config_yaml):
        """When --perform taxo_profiling is given, taxonomic profiling steps should be included regardless of config content."""
        content = self._build(full_config_yaml, tmp_path, "--perform", "taxo_profiling")
        assert "results/09_taxonomic_profiling" in content
        assert "results/03_assembly" not in content

    def test_perform_strain_profiling(self, tmp_path, full_config_yaml):
        """When --perform strain_profiling is given, strain profiling steps should be included regardless of config content."""
        content = self._build(full_config_yaml, tmp_path, "--perform", "strain_profiling")
        assert "results/10_strain_profiling" in content
        assert "results/05_binning" not in content

    def test_perform_multiple_steps(self, tmp_path, full_config_yaml):
        """When multiple --perform flags are given, the union of specified steps should be included."""
        content = self._build(
            full_config_yaml, tmp_path,
            "--perform", "assembly", "--perform", "binning",
        )
        assert "results/03_assembly" in content
        assert "results/05_binning" in content
        assert "results/08_bins_postprocessing" not in content

    def test_strain_profiling_tool_implicitly_activates_step(self, tmp_path, full_config_yaml):
        """--strain-profiling without --perform strain_profiling should still render the block."""
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(full_config_yaml), "--output", str(output),
             "--perform", "assembly", "--perform", "binning", "--perform", "postprocessing",
             "--strain-profiling", "instrain"],
        )
        assert result.exit_code == 0, result.output
        content = output.read_text()
        assert "inStrain" in content
        assert "floria" not in content

    def test_post_processing_tool_implicitly_activates_step(self, tmp_path, full_config_yaml):
        """--post-processing without --perform postprocessing should still render the block."""
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(full_config_yaml), "--output", str(output),
             "--perform", "assembly", "--perform", "binning",
             "--post-processing", "gtdbtk"],
        )
        assert result.exit_code == 0, result.output
        content = output.read_text()
        assert "gtdb_tk" in content
        assert "results/05_binning" in content


class TestPostProcessingTools:
    """Tests for the --post-processing sub-tool filter."""

    def _build_pp(self, full_config_yaml, tmp_path, *tools):
        """Helper to run `strainmake build` with given config and post-processing tools, returning the rendered Snakefile content for inspection."""
        args = ["--perform", "postprocessing"]
        for t in tools:
            args += ["--post-processing", t]
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(full_config_yaml), "--output", str(output)] + args,
        )
        assert result.exit_code == 0, result.output
        return output.read_text()

    def test_no_post_processing_flag_renders_all_tools(self, tmp_path, full_config_yaml):
        """When --post-processing is not given, all post-processing tools should be included regardless of config content."""
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(full_config_yaml), "--output", str(output),
             "--perform", "postprocessing"],
        )
        assert result.exit_code == 0, result.output
        content = output.read_text()
        assert "gtdb_tk" in content
        assert "checkm1" in content
        assert "carveme" in content
        assert "bakta" in content

    def test_gtdbtk_only(self, tmp_path, full_config_yaml):
        """When --post-processing gtdbtk is given, only gtdbtk should be included."""
        content = self._build_pp(full_config_yaml, tmp_path, "gtdbtk")
        assert "gtdb_tk" in content
        assert "checkm1" not in content
        assert "carveme" not in content
        assert "bakta" not in content

    def test_coverage_maps_to_checkm1(self, tmp_path, full_config_yaml):
        """When --post-processing coverage is given, the checkm1 rule should be included (coverage is an alias for checkm1 within the CLI)."""
        content = self._build_pp(full_config_yaml, tmp_path, "coverage")
        assert "checkm1" in content
        assert "gtdb_tk" not in content

    def test_carveme_only(self, tmp_path, full_config_yaml):
        """When --post-processing carveme is given, only carveme should be included."""
        content = self._build_pp(full_config_yaml, tmp_path, "carveme")
        assert "carveme" in content
        assert "gtdb_tk" not in content
        assert "bakta" not in content

    def test_bakta_only(self, tmp_path, full_config_yaml):
        """When --post-processing bakta is given, only bakta should be included."""
        content = self._build_pp(full_config_yaml, tmp_path, "bakta")
        assert "bakta" in content
        assert "carveme" not in content

    def test_drep_always_present(self, tmp_path, full_config_yaml):
        """dRep is unconditional within postprocessing."""
        content = self._build_pp(full_config_yaml, tmp_path, "gtdbtk")
        assert "dRep" in content

    def test_multiple_pp_tools(self, tmp_path, full_config_yaml):
        """When multiple --post-processing flags are given, the union of specified tools should be included."""
        content = self._build_pp(full_config_yaml, tmp_path, "gtdbtk", "carveme")
        assert "gtdb_tk" in content
        assert "carveme" in content
        assert "checkm1" not in content
        assert "bakta" not in content


class TestStrainProfilingTools:
    """Tests for the --strain-profiling sub-tool filter."""

    def _build_sp(self, full_config_yaml, tmp_path, *tools):
        """Helper to run `strainmake build` with given config and strain profiling tools, returning the rendered Snakefile content for inspection."""
        args = ["--perform", "strain_profiling", "--perform", "postprocessing"]
        for t in tools:
            args += ["--strain-profiling", t]
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(full_config_yaml), "--output", str(output)] + args,
        )
        assert result.exit_code == 0, result.output
        return output.read_text()

    def test_no_strain_profiling_flag_renders_all_tools(self, tmp_path, full_config_yaml):
        """When --strain-profiling is not given, all strain profiling tools should be included regardless of config content."""
        output = tmp_path / "Snakefile"
        result = runner.invoke(
            app,
            ["build", "--config", str(full_config_yaml), "--output", str(output),
             "--perform", "strain_profiling", "--perform", "postprocessing"],
        )
        assert result.exit_code == 0, result.output
        content = output.read_text()
        assert "inStrain" in content
        assert "floria" in content

    def test_instrain_only(self, tmp_path, full_config_yaml):
        """When --strain-profiling instrain is given, only instrain should be included."""
        content = self._build_sp(full_config_yaml, tmp_path, "instrain")
        assert "inStrain" in content
        assert "floria" not in content

    def test_floria_only(self, tmp_path, full_config_yaml):
        """When --strain-profiling floria is given, only floria should be included."""
        content = self._build_sp(full_config_yaml, tmp_path, "floria")
        assert "floria" in content
        assert "inStrain" not in content
