"""
strainmake - unified CLI entry point for the StrainMake metagenomic pipeline.

Usage:
    strainmake [--version] <command> [args...]

Commands:
    init     Interactively generate a configuration YAML
    build    Generate a Snakefile from a configuration YAML
    run      Run the pipeline via Snakemake
    report   Generate a MultiQC report from pipeline results
    prepare  Prepare input data and results (subgroup)

    fetch-databases  Download the reference databases (run where there is internet)
"""

from typing import Optional

import typer

from strainmake import __version__
from strainmake.cli import (
    init_cmd,
    build_cmd,
    run_cmd,
    report_cmd,
    prepare_cmd,
    fetch_databases_cmd,
)

# ---------------------------------------------------------------------------
# Root application
# ---------------------------------------------------------------------------

app = typer.Typer(
    name="strainmake",
    help="[bold green]StrainMake[/bold green]: reproducible hybrid metagenomics with MAG recovery and strain-level resolution.",
    add_completion=True,
    no_args_is_help=True,
    rich_markup_mode="rich",
)

# ---------------------------------------------------------------------------
# Prepare sub-group
# ---------------------------------------------------------------------------

prepare_app = typer.Typer(
    name="prepare",
    help="Prepare input data and results for the pipeline.",
    no_args_is_help=True,
    rich_markup_mode="rich",
)

app.add_typer(prepare_app, name="prepare")

prepare_app.command("import-preprocessed")(prepare_cmd.import_preprocessed)
prepare_app.command("split-samples")(prepare_cmd.split_samples)
prepare_app.command("gather-assemblies")(prepare_cmd.gather_assemblies)
prepare_app.command("gather-binning")(prepare_cmd.gather_binning)

# ---------------------------------------------------------------------------
# Top-level commands
# ---------------------------------------------------------------------------

app.command("init")(init_cmd.init)
app.command("build")(build_cmd.build)
app.command(
    "run",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)(run_cmd.run)
app.command("report")(report_cmd.report)
app.command("fetch-databases")(fetch_databases_cmd.fetch_databases)

# ---------------------------------------------------------------------------
# --version flag
# ---------------------------------------------------------------------------

BANNER = """
╔══════════════════════════════════════════╗
║           [bold cyan]StrainMake[/bold cyan] v{version:<20}  ║
╚══════════════════════════════════════════╝
"""


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"StrainMake {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        "-V",
        help="Show version and exit.",
        callback=_version_callback,
        is_eager=True,
    ),
) -> None:
    """
    [bold green]StrainMake[/bold green]: reproducible hybrid metagenomics with MAG recovery and strain-level resolution.

    Run [cyan]strainmake COMMAND --help[/cyan] for help on a specific command.
    """


if __name__ == "__main__":
    app()
