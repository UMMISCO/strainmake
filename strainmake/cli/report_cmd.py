"""
strainmake report - generate a MultiQC report from StrainMake pipeline results.

Wraps the HandleResults class in
workflow/scripts/multiqc_results/generate_multiqc_report.py.
The import is intentionally lazy so that missing `multiqc` installations do
not break the entire CLI.
"""

from pathlib import Path

import typer

from strainmake import REPO_ROOT

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_MULTIQC_CONFIG: Path = (
    REPO_ROOT / "workflow" / "scripts" / "multiqc_results" / "multiqc_config.yaml"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_handle_results():
    """Lazy-load HandleResults; isolated into a function so tests can patch it."""
    import sys

    repo_root_str = str(REPO_ROOT)
    if repo_root_str not in sys.path:
        sys.path.insert(0, repo_root_str)

    from workflow.scripts.multiqc_results.generate_multiqc_report import (  # type: ignore[import]
        HandleResults,
    )

    return HandleResults


# ---------------------------------------------------------------------------
# Typer command
# ---------------------------------------------------------------------------


def report(
    results_dir: Path = typer.Option(
        ...,
        "--results",
        "-r",
        help="Directory containing StrainMake pipeline results (the `results/` folder).",
        exists=True,
        file_okay=False,
    ),
    log_dir: Path = typer.Option(
        ...,
        "--logs",
        "-l",
        help="Directory containing pipeline logs (the `logs/` folder).",
        exists=True,
        file_okay=False,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Directory where the MultiQC report(s) will be written.",
    ),
    ani: int = typer.Option(
        95,
        "--ani",
        help="ANI threshold used for MAG dereplication (must match the pipeline run).",
    ),
    multiqc_config: Path = typer.Option(
        _DEFAULT_MULTIQC_CONFIG,
        "--multiqc-config",
        help="Path to the MultiQC YAML configuration file.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        "-n",
        help="Print the MultiQC commands without executing them.",
    ),
) -> None:
    """
    Generate [bold]MultiQC report(s)[/bold] from StrainMake pipeline results.

    Collects QC artefacts from preprocessing, assembly, binning and annotation
    steps, then runs MultiQC.  One report is produced per assembler when
    multiple assemblers were used.
    """
    try:
        HandleResults = _load_handle_results()
    except ImportError as exc:
        typer.echo(
            f"❌  Could not import report dependencies: {exc}\n"
            "Ensure the strainmake conda environment is activated or install multiqc.",
            err=True,
        )
        raise typer.Exit(1) from exc

    output_dir.mkdir(parents=True, exist_ok=True)

    handler = HandleResults(
        results_dir=str(results_dir),
        log_dir=str(log_dir),
        output_dir=str(output_dir),
        ani=ani,
        multiqc_config=str(multiqc_config),
    )
    handler.collect_all_results()
    handler.prepare_results()
    handler.generate_report(dry_run=dry_run)

    if not dry_run:
        typer.echo("✅  MultiQC report generation complete.")
