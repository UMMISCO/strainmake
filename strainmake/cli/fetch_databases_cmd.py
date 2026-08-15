"""
strainmake fetch-databases - download the pipeline's reference databases.

Meant to be run once, from a machine that has internet access, so that the
pipeline itself can then run somewhere that does not, but has access to the downloaded databases. 
Each database lands at its own `databases: <name>:` path, and none of 
them are project-specific, so they can be shared across projects.
"""

from pathlib import Path
from typing import List, Optional

import typer

from strainmake import REPO_ROOT
from strainmake.cli._snakemake import run_snakemake
from strainmake.databases import (
    DATABASE_KEYS,
    conda_prefix,
    inspect_databases,
    load_config,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_SNAKEFILE: Path = REPO_ROOT / "workflow" / "Snakefile"


# ---------------------------------------------------------------------------
# Typer command
# ---------------------------------------------------------------------------


def fetch_databases(
    config: Path = typer.Option(
        ...,
        "--config",
        "-c",
        help="Path to the StrainMake configuration YAML.",
        exists=True,
        readable=True,
    ),
    only: Optional[List[str]] = typer.Option(
        None,
        "--only",
        help=(
            "Fetch only these databases. Pass once per database, e.g. "
            "`--only checkm2 --only bakta`. Defaults to everything the "
            "configuration needs."
        ),
    ),
    all_databases: bool = typer.Option(
        False,
        "--all",
        help="Fetch every known database, even ones this configuration does not use.",
    ),
    cores: int = typer.Option(
        1,
        "--cores",
        min=1,
        help="Number of CPU cores to use (a few databases are unpacked in parallel).",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        "-n",
        help="Show what would be downloaded, without downloading anything.",
    ),
) -> None:
    """
    Download the reference databases the pipeline needs.

    Run this on a machine with internet access. Each database lands at its own
    [cyan]databases:[/cyan] path, which other projects can point at too.

    \b
    Examples:
        # everything the configuration calls for
        strainmake fetch-databases --config config.yaml

        # just one database
        strainmake fetch-databases --config config.yaml --only checkm2

        # see what is missing without downloading
        strainmake fetch-databases --config config.yaml --dry-run
    """
    config_data = load_config(config)
    statuses = inspect_databases(config_data)

    if only:
        unknown = sorted(set(only) - set(DATABASE_KEYS))
        if unknown:
            typer.echo(
                f"❌  Unknown database(s): {', '.join(unknown)}\n"
                f"    Known databases: {', '.join(DATABASE_KEYS)}",
                err=True,
            )
            raise typer.Exit(1)
        selected = [statuses[key] for key in only]
    elif all_databases:
        selected = list(statuses.values())
    else:
        selected = [status for status in statuses.values() if status.required]

    # databases StrainMake cannot fetch on its own are reported, never attempted
    manual = [status for status in selected if not status.downloadable]
    fetchable = [status for status in selected if status.downloadable]

    typer.echo("StrainMake reference databases\n")
    for status in selected:
        if status.present:
            mark, note = "✅", "already present"
        elif status.downloadable:
            mark, note = "⬇️ ", "will be downloaded"
        else:
            mark, note = "⚠️ ", "must be installed manually"
        typer.echo(f"  {mark} {status.label:<20} {status.path}  ({note})")

    if manual:
        typer.echo("\nDatabases StrainMake cannot download for you:")
        for status in manual:
            if status.present:
                continue
            typer.echo(f"  • {status.label}: place it at {status.path}")
            if status.manual_source:
                typer.echo(f"    see {status.manual_source}")

    todo = [status for status in fetchable if not status.present]
    if not todo:
        typer.echo("\n✅  Nothing left to download.")
        # a manual database that is still missing is worth a non-zero exit, so
        # that a setup script notices rather than proceeding to a doomed run
        raise typer.Exit(1 if any(not s.present for s in manual) else 0)

    targets = [status.target for status in todo]

    cmd: List[str] = [
        "snakemake",
        "--snakefile", str(_DEFAULT_SNAKEFILE),
        "--configfile", str(config),
        "--use-conda",
        "--cores", str(cores),
        # databases are shared state: never let a half-finished download stand in
        "--rerun-incomplete",
        *targets,
    ]

    # reuse the environments the pipeline itself will use, rather than building a
    # second copy of them somewhere else
    prefix = conda_prefix(config_data)
    if prefix:
        cmd += ["--conda-prefix", prefix]

    if dry_run:
        cmd.append("--dry-run")

    returncode = run_snakemake(cmd)

    if returncode == 0 and not dry_run:
        typer.echo(
            "\n✅  Databases ready. Once the conda environments are built too\n"
            "    (`strainmake run --config <config.yaml> --create-envs-only`),\n"
            "    the pipeline can run without internet access:\n"
            "    strainmake run --config <config.yaml> --offline --cores <N>"
        )

    raise typer.Exit(code=returncode)
