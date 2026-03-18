"""
Shared pytest fixtures for strainmake CLI tests.
"""

from pathlib import Path

import pytest
import yaml


# ---------------------------------------------------------------------------
# Sample data
# ---------------------------------------------------------------------------

SAMPLE_METADATA = (
    "sample\tsample_id\ttype\n"
    "/data/sample1_R1.fastq.gz\tsample1\tR1\n"
    "/data/sample1_R2.fastq.gz\tsample1\tR2\n"
    "/data/sample2_R1.fastq.gz\tsample2\tR1\n"
    "/data/sample2_R2.fastq.gz\tsample2\tR2\n"
)

SAMPLE_METADATA_WITH_LR = (
    "sample\tsample_id\ttype\n"
    "/data/sample1_R1.fastq.gz\tsample1\tR1\n"
    "/data/sample1_R2.fastq.gz\tsample1\tR2\n"
    "/data/sample1.fastq.gz\tsample1\tlong\n"
)

MINIMAL_CONFIG: dict = {
    "samples": "data/config_data.tsv",
    "lr_seq_format": "fastq",
    "assembly": {
        "assembler": ["megahit"],
        "long_read_assembler": None,
        "hybrid_assembler": None,
        "megahit": {
            "threads": 4,
            "min_contig_len": 500,
            "other_params": "",
        },
    },
    "fastp": {
        "compression": 2,
        "minimal_read_length": 50,
        "qualified_quality_phred": 15,
        "other_params": "",
    },
}

FULL_CONFIG: dict = {
    "samples": "data/config_data.tsv",
    "lr_seq_format": "fastq",
    "lr_technology": "pacbio",
    "assembly": {
        "assembler": ["megahit"],
        "long_read_assembler": None,
        "hybrid_assembler": None,
        "megahit": {"threads": 4, "min_contig_len": 500, "other_params": ""},
    },
    "fastp": {"compression": 2, "minimal_read_length": 50, "qualified_quality_phred": 15, "other_params": ""},
    "fastp_long_read": {"compression": 2, "minimal_read_length": 1000, "qualified_quality_phred": 12, "other_params": ""},
    "representative_genes": {"minimal_gene_length": 100},
    "mmseqs2": {"sequence_identity_threshold": 0.95, "alignment_coverage_shorter_sequence": 0.9, "threads": 4},
    "binning": {
        "binner": ["metabat2"],
        "long_read_binner": None,
        "minimap2": {"threads": 4},
        "samtools": {"threads": 4},
        "metabat2": {"min_contig_size": 1500, "minimum_mean_coverage": 1, "min_bin_size": 200000, "threads": 4},
    },
    "checkm2": {"threads": 4},
    "bins_refinement": {"binette": {"threads": 4, "low_mem": ""}},
    "bins_postprocessing": {
        "drep": {"ani": [95, 97], "comparison_algorithm": "ANImf", "other_args": "", "threads": 4},
        "genes_prediction": {"prodigal": {"threads": 4, "ani": 95}},
        "gtdbtk": {"threads": 4, "other_args": ""},
        "genomes_quality_filtration": {
            "checkm2": {"threads": 4},
            "filtration": {"min_completeness": 75, "max_contamination": 10},
        },
        "profiling": {"checkm1": {"threads": 4}},
        "carveme": {"threads": 4, "solver": "scip"},
        "bakta": {"threads": 4, "parallel_jobs": 1},
    },
    "taxonomic_profiling": {
        "metaphlan": {"threads": 4},
        "meteor": {"threads": 4, "downsize": [10000000]},
    },
    "strains_profiling": {
        "minimap2": {"threads": 4},
        "freebayes": {"min_alternate_count": 1, "min_alternate_fraction": 0.01},
        "instrain": {"threads": 4},
        "floria": {"threads": 4},
    },
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_metadata_tsv(tmp_path: Path) -> Path:
    """A minimal metadata TSV with two short-read samples."""
    p = tmp_path / "metadata.tsv"
    p.write_text(SAMPLE_METADATA)
    return p


@pytest.fixture
def sample_metadata_with_lr_tsv(tmp_path: Path) -> Path:
    """A metadata TSV with one short-read + one long-read sample."""
    p = tmp_path / "metadata_lr.tsv"
    p.write_text(SAMPLE_METADATA_WITH_LR)
    return p


@pytest.fixture
def sample_config_yaml(tmp_path: Path, sample_metadata_tsv: Path) -> Path:
    """A minimal config YAML suitable for testing `strainmake build`."""
    cfg = dict(MINIMAL_CONFIG)
    cfg["samples"] = str(sample_metadata_tsv)
    p = tmp_path / "config.yaml"
    p.write_text(yaml.dump(cfg, sort_keys=False))
    return p


@pytest.fixture
def full_config_yaml(tmp_path: Path, sample_metadata_tsv: Path) -> Path:
    """A full config YAML with all pipeline sections present."""
    import copy
    cfg = copy.deepcopy(FULL_CONFIG)
    cfg["samples"] = str(sample_metadata_tsv)
    p = tmp_path / "full_config.yaml"
    p.write_text(yaml.dump(cfg, sort_keys=False))
    return p


@pytest.fixture
def sample_metadata_df(sample_metadata_tsv: Path):
    """Loaded DataFrame of the sample metadata fixture."""
    import pandas as pd

    return pd.read_csv(sample_metadata_tsv, sep="\t")
