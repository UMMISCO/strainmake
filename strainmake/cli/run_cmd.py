"""
strainmake run - run the StrainMake pipeline via Snakemake.

Wraps `snakemake --use-conda` with sensible defaults.
Any extra positional arguments are forwarded verbatim to Snakemake.
"""

import subprocess
from pathlib import Path
from typing import List, Optional

import typer

from strainmake import REPO_ROOT

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_SNAKEFILE: Path = REPO_ROOT / "workflow" / "Snakefile"
_PREPROCESSED_SNAKEFILE: Path = (
    REPO_ROOT / "workflow" / "Snakefile.already_preprocessed_seq.smk"
)


# ---------------------------------------------------------------------------
# Typer command
# ---------------------------------------------------------------------------


def run(
    ctx: typer.Context,
    snakefile: Path = typer.Option(
        _DEFAULT_SNAKEFILE,
        "--snakefile",
        "-s",
        help="Path to the Snakefile to execute.",
    ),
    cores: int = typer.Option(
        4,
        "--cores",
        "-c",
        min=1,
        help="Number of CPU cores to use.",
    ),
    configfile: Optional[Path] = typer.Option(
        None,
        "--configfile",
        help="Path to a config YAML (passed via --configfile to Snakemake).",
        exists=True,
        readable=True,
    ),
    already_preprocessed: bool = typer.Option(
        False,
        "--already-preprocessed",
        help=(
            "Start from already-preprocessed reads "
            "(uses Snakefile.already_preprocessed_seq.smk)."
        ),
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        "-n",
        help="Dry run: print rules that would be executed without running them.",
    ),
) -> None:
    """
    Run the [bold]StrainMake[/bold] pipeline using Snakemake.

    Wraps [cyan]snakemake --use-conda[/cyan] with sensible defaults.
    Any unknown options/flags are passed through to Snakemake unchanged.
    Use [cyan]--[/cyan] to separate StrainMake options from Snakemake ones:

    \b
        strainmake run --cores 16
        strainmake run --cores 8 --configfile my_config.yaml
        strainmake run --already-preprocessed --cores 8
        strainmake run --dry-run
        strainmake run --cores 16 -- --rerun-incomplete --keep-going
    """
    effective_snakefile = _PREPROCESSED_SNAKEFILE if already_preprocessed else snakefile

    cmd: List[str] = [
        "snakemake",
        "--snakefile", str(effective_snakefile),
        "--use-conda",
        "--cores", str(cores),
    ]

    if configfile:
        cmd += ["--configfile", str(configfile)]

    if dry_run:
        cmd.append("--dry-run")

    if ctx.args:
        cmd.extend(ctx.args)

    typer.echo(f"  Running: {' '.join(cmd)}\n")

    result = subprocess.run(cmd)
    raise typer.Exit(code=result.returncode)
