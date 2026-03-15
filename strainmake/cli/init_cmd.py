"""
strainmake init - interactively generate a configuration YAML for the pipeline.

This command replaces the standalone `workflow/scripts/config_generator/config_generator.py`
script. It reads defaults directly from `config/template_config.yaml` at runtime, removing
the fragile subprocess-based `generate_defaults.py` call in the original script.
"""

from pathlib import Path
from typing import Any, Optional

import typer
import yaml

from strainmake import REPO_ROOT

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TEMPLATE_CONFIG: Path = REPO_ROOT / "config" / "template_config.yaml"

# maps human-readable section names to the YAML top-level keys they cover.
SECTION_KEYS: dict[str, list[str]] = {
    "PREPROCESSING": ["fastp", "fastp_long_read", "bowtie2", "downsizing_for_hybrid", "lr_technology"],
    "ASSEMBLY": ["assembly", "quast"],
    "GENE_CATALOG": ["mmseqs2", "representative_genes"],
    "BINNING": ["binning", "checkm2", "bins_refinement", "bins_postprocessing"],
    "TAXO_PROFILING": ["taxonomic_profiling"],
    "STRAINS_PROFILING": ["strains_profiling"],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_defaults() -> dict:
    """Load default configuration values from the template YAML."""
    if not TEMPLATE_CONFIG.exists():
        raise FileNotFoundError(
            f"Template config not found: {TEMPLATE_CONFIG}\n"
            "Make sure you run `strainmake init` from the repo root."
        )
    with open(TEMPLATE_CONFIG) as fh:
        return yaml.safe_load(fh)


def _get_section(defaults: dict, keys: list[str]) -> dict:
    """Extract a subset of top-level keys from the defaults dict."""
    return {k: v for k, v in defaults.items() if k in keys}


def _prompt_recursive(data: Any, prefix: str = "") -> Any:
    """Recursively prompt the user to confirm or override each config value."""
    # for dicts, we recurse into each key and build a new dict with the same structure
    if isinstance(data, dict):
        result: dict = {}
        for key, value in data.items():
            full_key = f"{prefix}.{key}" if prefix else key
            result[key] = _prompt_recursive(value, prefix=full_key)
        return result
    # for lists, we can accept a comma-separated string or keep the default
    if isinstance(data, list):
        typer.echo(f"  {prefix} (default: {data})")
        val = input("> ").strip()
        return data if val == "" else val.split()
    # scalar
    typer.echo(f"  {prefix} (default: {data!r})")
    val = input("> ").strip()
    if val == "":
        return data
    try:
        return eval(val, {}, {})  # noqa: S307 – intentional for config values
    except Exception:
        return val


def _prompt_for_section(section_name: str, default_data: dict) -> Optional[dict]:
    """
    Ask the user how to handle one config section.

    Choices:
      Y / <Enter>  – keep defaults
      n            – edit interactively
      x            – skip (omit section from output)
    """
    typer.echo(
        f"\n  Would you like to use default parameters for "
        f"[bold]{section_name}[/bold]? [Y/n/x]"
    )
    choice = input("> ").strip().lower()

    # we return None for skipped sections, which signals the caller to omit them from the final config
    if choice == "x":
        typer.echo(f"  ❌  Section \"{section_name}\" skipped.")
        return None
    if choice in ("y", ""):
        return default_data

    typer.echo(f"  Editing section: {section_name}")
    return _prompt_recursive(default_data)


# ---------------------------------------------------------------------------
# Typer command
# ---------------------------------------------------------------------------


def init(
    samples: Path = typer.Option(
        ...,
        "--samples",
        "-s",
        help="Path to the sample metadata TSV (columns: sample, sample_id, type).",
        exists=True,
        readable=True,
    ),
    lr_seq_format: str = typer.Option(
        "fastq",
        "--lr-format",
        help="Format of long-read sequences: 'fastq' or 'fasta'.",
    ),
    output: Path = typer.Option(
        "config.yaml",
        "--output",
        "-o",
        help="Path to write the generated configuration YAML.",
    ),
) -> None:
    """
    Interactively generate a [bold]config.yaml[/bold] for the StrainMake pipeline.

    Walks through each pipeline section (preprocessing, assembly, …) and lets
    you accept defaults, edit values, or skip a section entirely.
    """
    # loading defaults from the template YAML
    try:
        defaults = load_defaults()
    except FileNotFoundError as exc:
        typer.echo(f"❌  {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo("\n📄  Generating StrainMake configuration …\n")

    config: dict = {
        "samples": str(samples),
        "lr_seq_format": lr_seq_format,
    }

    for section_name, keys in SECTION_KEYS.items():
        section_defaults = _get_section(defaults, keys)
        section_data = _prompt_for_section(section_name, section_defaults)
        if section_data is not None:
            config.update(section_data)

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as fh:
        yaml.dump(config, fh, sort_keys=False)

    typer.echo(f"\n✅  Configuration written to [bold]{output}[/bold]")
