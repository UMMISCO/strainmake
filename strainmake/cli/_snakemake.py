"""
Shared helper for the commands that shell out to Snakemake.

Keeps the "is Snakemake installed?" answer in one place.
"""

import subprocess
from typing import List

import typer


def run_snakemake(cmd: List[str]) -> int:
    """
    Run a Snakemake command, reporting a missing executable clearly.

    Parameters:
    cmd (list): the full command line, starting with 'snakemake'

    Returns:
    int: Snakemake's exit code
    """

    typer.echo(f"  Running: {' '.join(cmd)}\n")

    try:
        return subprocess.run(cmd).returncode
    except FileNotFoundError:
        typer.echo(
            "❌  `snakemake` was not found on your PATH.\n"
            "    StrainMake drives Snakemake, it does not bundle it. See\n"
            "    https://snakemake.readthedocs.io/en/stable/getting_started/installation.html",
            err=True,
        )
        raise typer.Exit(127)
