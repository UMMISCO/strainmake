"""
strainmake prepare - subcommands for preparing input data and combining results.

Each command wraps a function from the corresponding script in
workflow/scripts/prepare/, keeping the business logic in one place while
providing a modern, discoverable CLI interface.

Commands
--------
import-preprocessed   Create symlinks for already-preprocessed FASTQ files.
split-samples         Split a metadata TSV across N sub-tables.
gather-assemblies     Merge assembly results from multiple split runs.
gather-binning        Merge binning results from multiple split runs.
"""

import sys
from pathlib import Path
from typing import List

import typer

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _import_prepare(module_name: str):
    """
    Lazily import a module from workflow.scripts.prepare.

    The import is deferred so that a missing optional dependency (e.g. pandas)
    only fails when the relevant subcommand is actually invoked.
    """
    import importlib

    from strainmake import REPO_ROOT

    repo_root_str = str(REPO_ROOT)
    if repo_root_str not in sys.path:
        sys.path.insert(0, repo_root_str)

    try:
        return importlib.import_module(f"workflow.scripts.prepare.{module_name}")
    except ModuleNotFoundError as exc:
        typer.echo(
            f"❌  Could not import prepare module '{module_name}': {exc}\n"
            "Make sure pandas is installed: pip install pandas",
            err=True,
        )
        raise typer.Exit(1) from exc


# ---------------------------------------------------------------------------
# import-preprocessed
# ---------------------------------------------------------------------------


def import_preprocessed(
    tsv_file: Path = typer.Argument(
        ...,
        help="Path to the sample metadata TSV (columns: sample, sample_id, type).",
        exists=True,
        readable=True,
    ),
    results_dir: Path = typer.Option(
        Path("results"),
        "--results-dir",
        help="Root directory where the `results/` tree will be created.",
    ),
    lr_seq_format: str = typer.Option(
        "fastq",
        "--lr-format",
        help="Format of long-read sequences: 'fastq' or 'fasta'.",
    ),
    clean: bool = typer.Option(
        False,
        "--clean",
        help="Remove previously created symlinks instead of creating new ones.",
    ),
) -> None:
    """
    Set up [bold]results/[/bold] symlinks for already-preprocessed FASTQ files.

    Reads the metadata TSV and creates symbolic links inside
    [dim]results/02_preprocess/bowtie2/[/dim] (short reads) and
    [dim]results/02_preprocess/fastp_long_read/[/dim] (long reads) so that the
    pipeline can start from the assembly step.

    After running this command, execute:

    \\b
        strainmake run --already-preprocessed --cores <N>
    """
    import pandas as pd

    mod = _import_prepare("already_preprocessed_seq")
    results_dir = results_dir.resolve()

    if clean:
        mod.clean_results_dir(results_dir)
        typer.echo("✅  Symlinks removed.")
        return

    metadata = pd.read_csv(tsv_file, sep="\t")
    sr_samples = mod.get_SR_samples(metadata)
    lr_samples = mod.get_LR_samples(metadata)
    mod.prepare_results_dir(sr_samples, lr_samples, results_dir, lr_seq_format)
    typer.echo("✅  Symlinks created in results/02_preprocess/")


# ---------------------------------------------------------------------------
# split-samples
# ---------------------------------------------------------------------------


def split_samples(
    input_file: Path = typer.Argument(
        ...,
        help="Path to the input metadata TSV to split.",
        exists=True,
        readable=True,
    ),
    n: int = typer.Option(
        ...,
        "--n",
        "-n",
        min=1,
        help="Number of output sub-tables to produce.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        help="Directory where the split TSV files will be written.",
    ),
) -> None:
    """
    Split a metadata TSV into [bold]N[/bold] sub-tables.

    Rows with the same [cyan]sample_id[/cyan] are always kept in the same
    output file.  Useful before launching parallel pipeline runs on a cluster.

    Example:

    \\b
        strainmake prepare split-samples metadata.tsv --n 4 --output-dir splits/
    """
    mod = _import_prepare("divide_metadata_table")
    mod.split_tsv(str(input_file), n, str(output_dir))
    typer.echo(f"✅  Split into {n} sub-tables under {output_dir}")


# ---------------------------------------------------------------------------
# gather-assemblies
# ---------------------------------------------------------------------------


def gather_assemblies(
    target_dir: Path = typer.Option(
        ...,
        "--target-dir",
        "-t",
        help="Root results directory where assemblies will be gathered into.",
        exists=True,
        file_okay=False,
    ),
    res_dirs: List[Path] = typer.Argument(
        ...,
        help="Source results directories from individual split runs.",
    ),
) -> None:
    """
    Merge [bold]assembly results[/bold] from multiple split runs into one directory.

    Use this after running `strainmake prepare split-samples` and launching
    separate pipeline runs per sub-table.  Assemblies from the [dim]03_assembly/[/dim]
    sub-directories of each source are moved into the target's [dim]03_assembly/[/dim].

    Example:

    \\b
        strainmake prepare gather-assemblies \\
            --target-dir results/ \\
            results_part1/ results_part2/ results_part3/
    """
    import os

    mod = _import_prepare("gather_results_of_assemblies")

    target = os.path.abspath(str(target_dir))
    sources = [os.path.abspath(str(d)) for d in res_dirs]

    # sanity checks to avoid accidentally deleting data: the target directory must not be one of the sources, 
    # and it must already contain a 03_assembly/ folder since the script will move files into it.
    if target in sources:
        typer.echo("❌  Target directory must not be one of the source directories.", err=True)
        raise typer.Exit(1)

    assembly_target = os.path.join(target, "03_assembly")
    if not os.path.exists(assembly_target):
        typer.echo(f"❌  {assembly_target} does not exist in the target directory.", err=True)
        raise typer.Exit(1)

    typer.echo(f"📂  Gathering assemblies from {len(sources)} source(s) → {assembly_target}")
    mod.gather_results(assembly_target, sources)
    typer.echo("✅  Assembly results gathered.")


# ---------------------------------------------------------------------------
# gather-binning
# ---------------------------------------------------------------------------


def gather_binning(
    target_dir: Path = typer.Option(
        ...,
        "--target-dir",
        "-t",
        help="Root results directory where binning results will be gathered into.",
        exists=True,
        file_okay=False,
    ),
    res_dirs: List[Path] = typer.Argument(
        ...,
        help="Source results directories from individual split runs.",
    ),
) -> None:
    """
    Merge [bold]binning results[/bold] from multiple split runs into one directory.

    Handles both [dim]05_binning/[/dim] and [dim]07_bins_refinement/[/dim].

    Example:

    \\b
        strainmake prepare gather-binning \\
            --target-dir results/ \\
            results_part1/ results_part2/
    """
    import os

    mod = _import_prepare("gather_results_of_binning")

    target = os.path.abspath(str(target_dir))
    sources = [os.path.abspath(str(d)) for d in res_dirs]

    # sanity checks to avoid accidentally deleting data: the target directory must not be one of the sources
    if target in sources:
        typer.echo("❌  Target directory must not be one of the source directories.", err=True)
        raise typer.Exit(1)

    for folder in ("05_binning", "07_bins_refinement"):
        full = os.path.join(target, folder)
        if not os.path.exists(full):
            typer.echo(f"❌  {full} does not exist in the target directory.", err=True)
            raise typer.Exit(1)

    typer.echo(f"📂  Gathering binning results from {len(sources)} source(s) → {target}")
    mod.gather_binning_results(target, sources, "05_binning")
    mod.gather_binning_results(target, sources, "07_bins_refinement")
    typer.echo("✅  Binning results gathered.")
