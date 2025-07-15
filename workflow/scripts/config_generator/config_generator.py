""" 
CLI for interactively generating the configuration YAML file for the pipeline.
"""

import typer
import yaml
from pathlib import Path
from typing import Any
from defaults import *
import subprocess

app = typer.Typer()

ALL_SECTIONS = {
    "PREPROCESSING": PREPROCESSING,
    "ASSEMBLY": ASSEMBLY,
    "GENE_CATALOG": GENE_CATALOG,
    "BINNING": BINNING,
    "TAXO_PROFILING": TAXO_PROFILING,
}

def prompt_for_section(section_name: str, default_data: dict) -> dict:
    """
    Prompts the user to decide how to handle configuration parameters for a given section.

    Prompts:
        - "Y" or Enter: Use default parameters.
        - "n": Edit parameters interactively.
        - "x": Skip this section.
    """
    typer.echo(f"\n Would you like to use default parameters for \"{section_name}\"? [Y/n/x]")
    choice = input("> ").strip().lower()

    if choice == "x":
        typer.echo(f"❌ Section \"{section_name}\" will be skipped.")
        return None
    if choice in ("y", ""):
        return default_data

    typer.echo(f"Editing section: {section_name}")
    return prompt_recursive(default_data)

def prompt_recursive(data: Any, prefix: str = "") -> Any:
    """
    Recursively prompts the user to input values for each field in a nested data structure.
    """
    if isinstance(data, dict):
        result = {}
        for key, value in data.items():
            full_key = f"{prefix}.{key}" if prefix else key
            result[key] = prompt_recursive(value, prefix=full_key)
        return result
    elif isinstance(data, list):
        typer.echo(f"{prefix} (default: {data})")
        val = input("> ")
        return data if val.strip() == "" else val.split()
    else:
        typer.echo(f"{prefix} (default: {data}): ")
        val = input("> ")
        if val.strip() == "":
            return data
        try:
            return eval(val, {}, {})
        except Exception:
            return val
        
@app.command()
def generate(
    samples: Path = typer.Option(..., help="Path to sample metadata file needed by the pipeline (TSV)"),
    lr_seq_format: str = typer.Option("fastq", help="Format of long reads: 'fastq' or 'fasta'"),
    output: Path = typer.Option("config.yaml", help="Path to write the final YAML")
):
    """
    Generate a configuration YAML file for the pipeline.
    """
    # Run generate_defaults.py to create defaults.yaml
    generate_script = Path(__file__).parent / "generate_defaults.py"
    typer.echo("🔄 Generating default YAML settings...")
    subprocess.run(["python", str(generate_script)], check=True)

    typer.echo("\n📄 Config file initialized.\n")

    config = {
        "samples": str(samples),
        "lr_seq_format": lr_seq_format,
    }

    for section_name, default_data in ALL_SECTIONS.items():
        section_data = prompt_for_section(section_name, default_data)
        if section_data is not None:
            config.update(section_data)

    with open(output, "w") as f:
        yaml.dump(config, f, sort_keys=False)

    typer.echo(f"\n✅ {output} written to disk.")

if __name__ == "__main__":
    app()
