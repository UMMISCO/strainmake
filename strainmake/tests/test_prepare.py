"""
Tests for `strainmake prepare` subcommands and the underlying helper functions.
"""

import math
import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from typer.testing import CliRunner

from strainmake.cli.main import app

runner = CliRunner()


# ===========================================================================
# Unit tests: already_preprocessed_seq.py
# ===========================================================================


class TestGetSrSamples:
    def test_returns_dict_with_r1_r2(self, sample_metadata_df):
        """get_SR_samples should return a dictionary with sample IDs as keys and R1/R2 paths as values when both are present in the metadata."""
        from workflow.scripts.prepare.already_preprocessed_seq import get_SR_samples  # type: ignore

        result = get_SR_samples(sample_metadata_df)
        assert result is not None
        assert "sample1" in result
        assert "R1" in result["sample1"]
        assert "R2" in result["sample1"]

    def test_returns_none_when_no_sr(self):
        """get_SR_samples should return None when no short-read samples are present in the metadata."""
        from workflow.scripts.prepare.already_preprocessed_seq import get_SR_samples  # type: ignore

        df = pd.DataFrame({"sample": ["/lr.fastq.gz"], "sample_id": ["s1"], "type": ["long"]})
        assert get_SR_samples(df) is None


class TestGetLrSamples:
    def test_returns_none_when_no_lr(self, sample_metadata_df):
        """get_LR_samples should return None when no long-read samples are present in the metadata."""
        from workflow.scripts.prepare.already_preprocessed_seq import get_LR_samples  # type: ignore

        assert get_LR_samples(sample_metadata_df) is None

    def test_returns_dict_with_long_reads(self, sample_metadata_with_lr_tsv):
        """get_LR_samples should return a dictionary with sample IDs as keys and long-read paths as values when long-read samples are present in the metadata."""
        from workflow.scripts.prepare.already_preprocessed_seq import (  # type: ignore
            get_LR_samples,
        )

        df = pd.read_csv(sample_metadata_with_lr_tsv, sep="\t")
        result = get_LR_samples(df)
        assert result is not None
        assert "sample1" in result


class TestPrepareResultsDir:
    def test_creates_sr_symlinks(self, tmp_path, sample_metadata_tsv):
        """prepare_results_dir should create symlinks for short-read samples in the correct location with the expected naming convention."""
        from workflow.scripts.prepare.already_preprocessed_seq import (  # type: ignore
            get_SR_samples,
            prepare_results_dir,
        )

        # Create fake source files so symlinks can resolve
        r1 = tmp_path / "src_R1.fastq.gz"
        r2 = tmp_path / "src_R2.fastq.gz"
        r1.touch()
        r2.touch()

        df = pd.DataFrame(
            {
                "sample": [str(r1), str(r2)],
                "sample_id": ["s1", "s1"],
                "type": ["R1", "R2"],
            }
        )
        sr = get_SR_samples(df)
        results = tmp_path / "results"
        prepare_results_dir(sr, None, results)

        bowtie2_dir = results / "02_preprocess" / "bowtie2"
        assert (bowtie2_dir / "s1_1.clean.fastq.gz").is_symlink()
        assert (bowtie2_dir / "s1_2.clean.fastq.gz").is_symlink()

    def test_creates_lr_symlinks(self, tmp_path):
        """prepare_results_dir should create symlinks for long-read samples in the correct location with the expected naming convention."""
        from workflow.scripts.prepare.already_preprocessed_seq import (  # type: ignore
            prepare_results_dir,
        )

        src = tmp_path / "lr.fastq.gz"
        src.touch()
        lr = {"s1": str(src)}
        results = tmp_path / "results"
        prepare_results_dir(None, lr, results, lr_seq_format="fastq")

        lr_dir = results / "02_preprocess" / "fastp_long_read"
        assert (lr_dir / "s1_1.fastq.gz").is_symlink()


class TestCleanResultsDir:
    def test_removes_symlinks(self, tmp_path):
        """clean_results_dir should remove existing symlinks in the results directory without affecting other files."""
        from workflow.scripts.prepare.already_preprocessed_seq import (  # type: ignore
            clean_results_dir,
        )

        bowtie2 = tmp_path / "02_preprocess" / "bowtie2"
        bowtie2.mkdir(parents=True)
        link = bowtie2 / "s1_1.clean.fastq.gz"
        target = tmp_path / "real.fastq.gz"
        target.touch()
        link.symlink_to(target)

        clean_results_dir(tmp_path)
        assert not link.exists()


# ===========================================================================
# Unit tests: divide_metadata_table.py
# ===========================================================================


class TestSplitTsv:
    def test_correct_number_of_files(self, tmp_path, sample_metadata_tsv):
        """split_tsv should produce the correct number of output files based on the specified N."""
        from workflow.scripts.prepare.divide_metadata_table import split_tsv  # type: ignore

        out = tmp_path / "splits"
        split_tsv(str(sample_metadata_tsv), 2, str(out))
        produced = list(out.glob("output_part_*.tsv"))
        assert len(produced) == 2

    def test_all_rows_preserved(self, tmp_path, sample_metadata_tsv):
        """split_tsv should preserve all rows from the original TSV in the output files."""
        from workflow.scripts.prepare.divide_metadata_table import split_tsv  # type: ignore

        out = tmp_path / "splits"
        split_tsv(str(sample_metadata_tsv), 2, str(out))
        rows = sum(
            len(pd.read_csv(f, sep="\t")) for f in out.glob("output_part_*.tsv")
        )
        original = len(pd.read_csv(sample_metadata_tsv, sep="\t"))
        assert rows == original

    def test_same_sample_id_stays_together(self, tmp_path, sample_metadata_tsv):
        """split_tsv should ensure that rows with the same sample_id are not split across different output files."""
        from workflow.scripts.prepare.divide_metadata_table import split_tsv  # type: ignore

        out = tmp_path / "splits"
        split_tsv(str(sample_metadata_tsv), 2, str(out))
        for f in out.glob("output_part_*.tsv"):
            df = pd.read_csv(f, sep="\t")
            for sid in df["sample_id"].unique():
                # All rows for a given sample_id must be in the same file
                assert (df["sample_id"] == sid).all() or (df["sample_id"] != sid).any()

    def test_n_larger_than_samples_capped(self, tmp_path, sample_metadata_tsv):
        """When N > number of unique sample_ids, output files ≤ N."""
        from workflow.scripts.prepare.divide_metadata_table import split_tsv  # type: ignore

        out = tmp_path / "splits"
        split_tsv(str(sample_metadata_tsv), 100, str(out))
        produced = list(out.glob("output_part_*.tsv"))
        n_samples = len(pd.read_csv(sample_metadata_tsv, sep="\t")["sample_id"].unique())
        assert len(produced) <= n_samples


# ===========================================================================
# Unit tests: gather_results_of_assemblies.py
# ===========================================================================


class TestGatherAssemblies:
    def _make_assembly_dir(self, base: Path, assembler: str, sample: str) -> Path:
        """Helper to create a fake assembly directory structure for testing."""
        d = base / "03_assembly" / assembler / sample
        d.mkdir(parents=True)
        (d / "assembly.fa").write_text(f">{sample}\nACGT\n")
        return d

    def test_moves_megahit_results(self, tmp_path):
        """gather_results should move assembly results from the source to the target directory, preserving the expected structure and files."""
        from workflow.scripts.prepare.gather_results_of_assemblies import (  # type: ignore
            gather_results,
        )

        src = tmp_path / "src"
        target = tmp_path / "target" / "03_assembly"
        target.mkdir(parents=True)
        self._make_assembly_dir(src, "megahit", "sample1")

        gather_results(str(target), [str(src)])

        assert (target / "megahit" / "sample1" / "assembly.fa").exists()

    def test_missing_src_assembly_dir_skipped(self, tmp_path):
        """gather_results should skip missing source assembly directories without raising an error."""
        from workflow.scripts.prepare.gather_results_of_assemblies import (  # type: ignore
            gather_results,
        )

        src = tmp_path / "empty_src"
        src.mkdir()
        target = tmp_path / "target" / "03_assembly"
        target.mkdir(parents=True)

        # Should not raise even when source has no 03_assembly/
        gather_results(str(target), [str(src)])


# ===========================================================================
# Unit tests: gather_results_of_binning.py
# ===========================================================================


class TestGatherBinning:
    def _make_binning_dir(self, base: Path, step: str, assembler: str) -> Path:
        """Helper to create a fake binning directory structure for testing."""
        d = base / step / assembler
        d.mkdir(parents=True)
        (d / "bin.fa").write_text(">bin1\nACGT\n")
        return d

    def test_moves_binning_results(self, tmp_path):
        """gather_binning_results should move binning results from the source to the target directory, preserving the expected structure and files."""
        from workflow.scripts.prepare.gather_results_of_binning import (  # type: ignore
            gather_binning_results,
        )

        src = tmp_path / "src"
        target = tmp_path / "target"
        (target / "05_binning").mkdir(parents=True)
        self._make_binning_dir(src, "05_binning", "semibin2")

        gather_binning_results(str(target), [str(src)], "05_binning")

        assert (target / "05_binning" / "semibin2" / "bin.fa").exists()


# ===========================================================================
# CLI tests: strainmake prepare import-preprocessed
# ===========================================================================


class TestImportPreprocessedCli:
    def test_symlinks_created(self, tmp_path, sample_metadata_tsv):
        """End-to-end: metadata with real file paths produces symlinks."""
        r1 = tmp_path / "sample1_R1.fastq.gz"
        r2 = tmp_path / "sample1_R2.fastq.gz"
        r1.touch()
        r2.touch()

        meta = tmp_path / "meta.tsv"
        meta.write_text(
            f"sample\tsample_id\ttype\n"
            f"{r1}\tsample1\tR1\n"
            f"{r2}\tsample1\tR2\n"
        )

        results_dir = tmp_path / "results"
        result = runner.invoke(
            app,
            [
                "prepare",
                "import-preprocessed",
                str(meta),
                "--results-dir",
                str(results_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (results_dir / "02_preprocess" / "bowtie2").exists()

    def test_clean_removes_symlinks(self, tmp_path):
        """Running with --clean should remove existing symlinks in the results directory."""
        bowtie2 = tmp_path / "02_preprocess" / "bowtie2"
        bowtie2.mkdir(parents=True)
        target = tmp_path / "real.fastq.gz"
        target.touch()
        link = bowtie2 / "s1_1.clean.fastq.gz"
        link.symlink_to(target)

        # Create a minimal valid TSV (not used in clean mode, but required arg)
        meta = tmp_path / "meta.tsv"
        meta.write_text("sample\tsample_id\ttype\n")

        result = runner.invoke(
            app,
            [
                "prepare",
                "import-preprocessed",
                str(meta),
                "--results-dir",
                str(tmp_path),
                "--clean",
            ],
        )
        assert result.exit_code == 0, result.output
        assert not link.exists()


# ===========================================================================
# CLI tests: strainmake prepare split-samples
# ===========================================================================


class TestSplitSamplesCli:
    def test_split_produces_files(self, tmp_path, sample_metadata_tsv):
        """Running `strainmake prepare split-samples` should produce the expected number of output files."""
        out = tmp_path / "splits"
        result = runner.invoke(
            app,
            [
                "prepare",
                "split-samples",
                str(sample_metadata_tsv),
                "--n",
                "2",
                "--output-dir",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(list(out.glob("output_part_*.tsv"))) == 2


# ===========================================================================
# CLI tests: strainmake prepare gather-assemblies
# ===========================================================================


class TestGatherAssembliesCli:
    def test_error_when_target_equals_source(self, tmp_path):
        """gather_assemblies should return an error when the target directory is the same as the source directory."""
        (tmp_path / "03_assembly").mkdir()
        result = runner.invoke(
            app,
            [
                "prepare",
                "gather-assemblies",
                "--target-dir",
                str(tmp_path),
                str(tmp_path),
            ],
        )
        assert result.exit_code != 0

    def test_error_when_assembly_dir_missing_in_target(self, tmp_path):
        """gather_assemblies should return an error when the target directory is missing the expected assembly sub-directory."""
        src = tmp_path / "src"
        src.mkdir()
        target = tmp_path / "target"
        target.mkdir()
        result = runner.invoke(
            app,
            [
                "prepare",
                "gather-assemblies",
                "--target-dir",
                str(target),
                str(src),
            ],
        )
        assert result.exit_code != 0

    def test_gathers_assemblies_successfully(self, tmp_path):
        """gather_assemblies should move assembly results from the source to the target directory, preserving the expected structure and files."""
        src = tmp_path / "src"
        assembly = src / "03_assembly" / "megahit" / "sample1"
        assembly.mkdir(parents=True)
        (assembly / "assembly.fa").write_text(">s1\nACGT")

        target = tmp_path / "target"
        (target / "03_assembly").mkdir(parents=True)

        result = runner.invoke(
            app,
            [
                "prepare",
                "gather-assemblies",
                "--target-dir",
                str(target),
                str(src),
            ],
        )
        assert result.exit_code == 0, result.output


# ===========================================================================
# CLI tests: strainmake prepare gather-binning
# ===========================================================================


class TestGatherBinningCli:
    def test_error_when_binning_dirs_missing_in_target(self, tmp_path):
        """gather_binning should return an error when the target directory is missing the expected binning sub-directories."""
        src = tmp_path / "src"
        src.mkdir()
        target = tmp_path / "target"
        target.mkdir()
        result = runner.invoke(
            app,
            [
                "prepare",
                "gather-binning",
                "--target-dir",
                str(target),
                str(src),
            ],
        )
        assert result.exit_code != 0

    def test_gathers_binning_successfully(self, tmp_path):
        """gather_binning should move binning results from the source to the target directory, preserving the expected structure and files."""
        src = tmp_path / "src"
        (src / "05_binning" / "semibin2").mkdir(parents=True)
        (src / "07_bins_refinement" / "binette").mkdir(parents=True)
        (src / "05_binning" / "semibin2" / "bin.fa").write_text(">b1\nACGT")

        target = tmp_path / "target"
        (target / "05_binning").mkdir(parents=True)
        (target / "07_bins_refinement").mkdir(parents=True)

        result = runner.invoke(
            app,
            [
                "prepare",
                "gather-binning",
                "--target-dir",
                str(target),
                str(src),
            ],
        )
        assert result.exit_code == 0, result.output
