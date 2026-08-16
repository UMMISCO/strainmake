"""Legacy `.test/unit` suite migrated to pytest."""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

from strainmake.tests.data import (
    generate_metadata_table,
    make_abundance_table,
    make_genome_list_table,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
PREPARE_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "prepare" / "already_preprocessed_seq.py"
DEDUP_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "deduplicate_contigs_name.py"
TEMPLATE_CONFIG = REPO_ROOT / "config" / "template_config.yaml"
FASTA_TEST = REPO_ROOT / "strainmake" / "tests" / "data" / "fasta" / "test.fa"

HAS_SNAKEMAKE = shutil.which("snakemake") is not None


# NOTE: DAG building used to need a placeholder GTDB-Tk directory, because the
# reference path came from `variables:` in workflow/envs/gtdb_tk.yaml and no rule
# could produce it. Every reference database is now declared in the `databases:`
# config section and produced by a rule in workflow/rules/00_databases.smk, so
# Snakemake resolves them on its own and no placeholder is required.


@pytest.fixture(name="simulated_reads_dir")
def _fixture_simulated_reads_dir(tmp_path: Path) -> Path:
    """Create minimal fake SR/LR files used to generate metadata tables."""
    reads = tmp_path / "reads"
    reads.mkdir(parents=True, exist_ok=True)

    files = [
        "fake_illumina_R1.SAMPLE1.fastq.gz",
        "fake_illumina_R2.SAMPLE1.fastq.gz",
        "fake_nanopore_sample0_aligned_reads.SAMPLE1.fastq.gz",
        "fake_nanopore_sample0_aligned_reads.SAMPLE1.fasta.gz",
    ]

    for name in files:
        with gzip.open(reads / name, "wt", encoding="utf-8") as fh:
            fh.write("@r1\nACGT\n+\nIIII\n" if name.endswith("fastq.gz") else ">c1\nACGT\n")

    return reads


@pytest.fixture(name="legacy_metadata_tables")
def _fixture_legacy_metadata_tables(tmp_path: Path, simulated_reads_dir: Path) -> dict[str, Path]:
    """Generate legacy metadata TSVs for SR/LR/ALL scenarios."""
    out = {}

    sr_df = generate_metadata_table.generate_metadata(str(simulated_reads_dir), "SR")
    sr_path = tmp_path / "metadata_sr.tsv"
    generate_metadata_table.write_metadata_tsv(sr_df, sr_path)
    out["sr"] = sr_path

    lr_df = generate_metadata_table.generate_metadata(str(simulated_reads_dir), "LR")
    lr_path = tmp_path / "metadata_lr.tsv"
    generate_metadata_table.write_metadata_tsv(lr_df, lr_path)
    out["lr"] = lr_path

    all_df = generate_metadata_table.generate_metadata(str(simulated_reads_dir), "all")
    all_path = tmp_path / "metadata_all.tsv"
    generate_metadata_table.write_metadata_tsv(all_df, all_path)
    out["all"] = all_path

    all_fasta_df = generate_metadata_table.generate_metadata(
        str(simulated_reads_dir), "all", lr_format="fasta"
    )
    all_fasta_path = tmp_path / "metadata_all_lr_fasta.tsv"
    generate_metadata_table.write_metadata_tsv(all_fasta_df, all_fasta_path)
    out["all_fasta"] = all_fasta_path

    sr_not_paired = pd.read_csv(all_path, sep="\t")
    sr_not_paired = sr_not_paired[sr_not_paired["type"] != "R2"]
    sr_not_paired_path = tmp_path / "metadata_sr_not_paired.tsv"
    sr_not_paired.to_csv(sr_not_paired_path, sep="\t", index=False)
    out["sr_not_paired"] = sr_not_paired_path

    return out


def _load_template_config() -> dict:
    with open(TEMPLATE_CONFIG, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _with_instrain_nested(config: dict) -> None:
    sp = config.setdefault("strains_profiling", {})
    instrain = sp.get("instrain", {}) or {}
    if "profile" not in instrain or "compare" not in instrain:
        threads = instrain.get("threads", 4)
        sp["instrain"] = {
            "profile": {"threads": threads},
            "compare": {"threads": threads},
        }


def _write_case_config(tmp_path: Path, case_name: str, samples_path: Path, lr_format: str = "fastq") -> Path:
    cfg = _load_template_config()
    cfg["samples"] = str(samples_path)
    cfg["lr_seq_format"] = lr_format

    if case_name == "1_SR_only":
        cfg["assembly"]["assembler"] = ["megahit", "metaspades"]
        cfg["assembly"]["long_read_assembler"] = None
        cfg["assembly"]["hybrid_assembler"] = None
    elif case_name == "2_LR_only":
        cfg["assembly"]["assembler"] = None
        cfg["assembly"]["long_read_assembler"] = ["metaflye"]
        cfg["assembly"]["hybrid_assembler"] = None
    elif case_name == "3_ALL_seq":
        cfg["assembly"]["assembler"] = ["megahit", "metaspades"]
        cfg["assembly"]["long_read_assembler"] = ["metaflye"]
        cfg["assembly"]["hybrid_assembler"] = ["hybridspades"]
    elif case_name == "4_SR_fail_if_LR_only":
        cfg["assembly"]["assembler"] = ["megahit", "metaspades"]
        cfg["assembly"]["long_read_assembler"] = None
        cfg["assembly"]["hybrid_assembler"] = None
    elif case_name == "5_LR_fail_if_SR_only":
        cfg["assembly"]["assembler"] = None
        cfg["assembly"]["long_read_assembler"] = ["metaflye"]
        cfg["assembly"]["hybrid_assembler"] = None
    elif case_name == "6_SR_not_paired":
        cfg["assembly"]["assembler"] = ["megahit", "metaspades"]
        cfg["assembly"]["long_read_assembler"] = ["metaflye"]
        cfg["assembly"]["hybrid_assembler"] = ["hybridspades"]
    elif case_name in {"7_ALL_already_preprocessed", "8_ALL_LR_FASTA"}:
        cfg["assembly"]["assembler"] = ["megahit", "metaspades"]
        cfg["assembly"]["long_read_assembler"] = ["metaflye"]
        cfg["assembly"]["hybrid_assembler"] = ["hybridspades"]

    _with_instrain_nested(cfg)

    out = tmp_path / f"{case_name}.yaml"
    with open(out, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return out


def _run_snakemake_dryrun(config_path: Path) -> subprocess.CompletedProcess[str]:
    """Run Snakemake in dry-run mode with the given config and capture output."""
    return subprocess.run(
        [
            "snakemake",
            "-c",
            "4",
            "--use-conda",
            "--configfile",
            str(config_path),
            "-s",
            str(SNAKEFILE),
            "--directory",
            str(REPO_ROOT),
            "--rerun-incomplete",
            "-n",
        ],
        check=False,
        capture_output=True,
        text=True,
    )


def parse_snakemake_dryrun_output(output: str) -> pd.DataFrame:
    """Parse Snakemake dry-run job summary table."""
    job_table_lines = []
    in_table = False

    for line in output.splitlines():
        if line.startswith("job") and "count" in line:
            in_table = True
            continue
        if in_table and line.strip().startswith("total "):
            job_table_lines.append(line.strip())
            break
        if in_table:
            job_table_lines.append(line.strip())

    job_data = []
    for job_line in job_table_lines[1:-1]:
        parts = re.split(r"\s{2,}", job_line)
        if len(parts) == 2:
            job_data.append(parts)

    if not job_data:
        return pd.DataFrame(columns=["job", "count"])

    df = pd.DataFrame(job_data, columns=["job", "count"])
    df["count"] = pd.to_numeric(df["count"])
    return df


def test_generate_metadata_table_sr(simulated_reads_dir: Path, tmp_path: Path):
    out = tmp_path / "metadata.tsv"
    rc = generate_metadata_table.main([str(simulated_reads_dir), "SR", "--output", str(out)])
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert set(df["type"].tolist()) == {"R1", "R2"}


def test_generate_metadata_table_lr_fasta(simulated_reads_dir: Path, tmp_path: Path):
    """Test that LR metadata can be generated with FASTA format and correct file extensions."""
    out = tmp_path / "metadata.tsv"
    rc = generate_metadata_table.main(
        [str(simulated_reads_dir), "LR", "--lr-format", "fasta", "--output", str(out)]
    )
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert set(df["type"].tolist()) == {"long"}
    assert df["sample"].str.endswith(".fasta.gz").all()


def test_make_abundance_table(tmp_path: Path):
    """Test that the abundance table can be generated from a set of genomes."""
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    for i in range(5):
        (genomes / f"g{i}.fna").write_text(">c1\nACGT\n", encoding="utf-8")

    out = tmp_path / "abundance.tsv"
    rc = make_abundance_table.main(["--genome-dir", str(genomes), "--output", str(out)])
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 5


def test_make_genome_list_table(tmp_path: Path):
    """Test that the genome list table and DNA type table can be generated from a set of genomes."""
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    for i in range(5):
        (genomes / f"g{i}.fna").write_text(">c1\nACGT\n", encoding="utf-8")

    genomes_out = tmp_path / "genomes_list_for_nanosim"
    dna_out = tmp_path / "dna_type_list_for_nanosim.tsv"

    rc = make_genome_list_table.main(
        [
            "--genome-dir",
            str(genomes),
            "--genomes-output",
            str(genomes_out),
            "--dna-type-output",
            str(dna_out),
        ]
    )
    assert rc == 0
    assert genomes_out.exists()
    assert dna_out.exists()


legacy_integration = pytest.mark.skipif(
    not HAS_SNAKEMAKE,
    reason="Snakemake is required to run legacy integration tests.",
)

@legacy_integration
@pytest.mark.parametrize(
    ("case_name", "metadata_key", "expected_returncode", "expected_job_counts"),
    [
        ("1_SR_only", "sr", 0, {"reads_mapping": 2}),
        ("2_LR_only", "lr", 0, {}),
        ("3_ALL_seq", "all", 0, {"reads_mapping": 2, "reads_mapping_LR": 1}),
        ("4_SR_fail_if_LR_only", "lr", 1, {}),
        ("5_LR_fail_if_SR_only", "sr", 1, {}),
        ("6_SR_not_paired", "sr_not_paired", 1, {}),
    ],
)
def test_legacy_dryrun_cases(
    tmp_path: Path,
    legacy_metadata_tables: dict[str, Path],
    case_name: str,
    metadata_key: str,
    expected_returncode: int,
    expected_job_counts: dict[str, int],
):
    """Test Snakemake dry-run for various legacy cases with expected job counts."""
    cfg = _write_case_config(tmp_path, case_name, legacy_metadata_tables[metadata_key])
    result = _run_snakemake_dryrun(cfg)

    combined_output = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == expected_returncode, combined_output

    if expected_returncode == 1:
        assert "ValueError" in combined_output
        return

    jobs = parse_snakemake_dryrun_output(combined_output)
    for job_name, expected_count in expected_job_counts.items():
        match = jobs[jobs["job"] == job_name]
        assert not match.empty, f"Job {job_name} not found in summary:\n{jobs}"
        assert int(match["count"].values[0]) == expected_count


@legacy_integration
def test_given_default_config_dryrun(tmp_path: Path, legacy_metadata_tables: dict[str, Path]):
    """Test Snakemake dry-run with the default configuration."""
    cfg = _write_case_config(tmp_path, "3_ALL_seq", legacy_metadata_tables["all"])
    result = _run_snakemake_dryrun(cfg)
    combined_output = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == 0, combined_output


@legacy_integration
def test_reads_already_preprocessed_implementation(tmp_path: Path, legacy_metadata_tables: dict[str, Path]):
    """Test that the prepare script creates expected symlinks for already preprocessed reads and can clean them up."""
    metadata = legacy_metadata_tables["all"]
    results_dir = tmp_path / "results"

    setup = subprocess.run(
        [
            "python3",
            str(PREPARE_SCRIPT),
            "--results_dir",
            str(results_dir),
            str(metadata),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert setup.returncode == 0, f"{setup.stdout}\n{setup.stderr}"

    expected_symlinks = [
        results_dir / "02_preprocess" / "bowtie2" / "SAMPLE1_1.clean.fastq.gz",
        results_dir / "02_preprocess" / "bowtie2" / "SAMPLE1_2.clean.fastq.gz",
        results_dir / "02_preprocess" / "fastp_long_read" / "SAMPLE1_1.fastq.gz",
    ]
    for symlink_path in expected_symlinks:
        assert symlink_path.exists()
        assert symlink_path.is_symlink()

    cleanup = subprocess.run(
        [
            "python3",
            str(PREPARE_SCRIPT),
            "--clean",
            "--results_dir",
            str(results_dir),
            str(metadata),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert cleanup.returncode == 0, f"{cleanup.stdout}\n{cleanup.stderr}"


@legacy_integration
def test_lr_reads_fasta_prepare_mode(tmp_path: Path, legacy_metadata_tables: dict[str, Path]):
    metadata = legacy_metadata_tables["all_fasta"]
    results_dir = tmp_path / "results"

    setup = subprocess.run(
        [
            "python3",
            str(PREPARE_SCRIPT),
            "--results_dir",
            str(results_dir),
            "--long-reads-seq-format",
            "fasta",
            str(metadata),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert setup.returncode == 0, f"{setup.stdout}\n{setup.stderr}"

    expected_lr = results_dir / "02_preprocess" / "fastp_long_read" / "SAMPLE1_1.fasta.gz"
    assert expected_lr.exists()
    assert expected_lr.is_symlink()


def test_sequences_renaming_in_fasta(tmp_path: Path):
    """Ensure sequence headers are renamed by `deduplicate_contigs_name.py`."""
    tmp_fasta_dir = tmp_path / "fasta"
    tmp_fasta_dir.mkdir(parents=True, exist_ok=True)
    fasta_test = tmp_fasta_dir / "test_tmp_for_test.fa"
    shutil.copy2(FASTA_TEST, fasta_test)

    before_headers = [
        line.strip()[1:].strip()
        for line in fasta_test.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    assert before_headers == ["contig_1", "contig_2", "contig_3"]

    result = subprocess.run(
        ["python3", str(DEDUP_SCRIPT), str(tmp_fasta_dir)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"

    after_headers = [
        line.strip()[1:].strip()
        for line in fasta_test.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    assert after_headers == [
        "test_tmp_for_test_contig_1",
        "test_tmp_for_test_contig_2",
        "test_tmp_for_test_contig_3",
    ]
