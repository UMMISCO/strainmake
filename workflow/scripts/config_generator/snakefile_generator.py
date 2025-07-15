""" 
Script to generate the suited Snakefile based on the generated configuration by `config_generator.py`.

Update this script and the Snakefile template when the main Snakifle changes.
"""

import typer
import yaml
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

BASE_SNAKEFILE = Path(__file__).parent.parent.parent / "Snakefile"

app = typer.Typer()

@app.command()
def generate(
    config: Path = typer.Option(..., help="Path to the generated YAML configuration file"),
    output: Path = typer.Option("Snakefile", help="Path to write the Snakefile")
):
    """
    Generate a Snakefile with relevant final parts based on the given configuration (YAML).
    """

    # loading the YAML configuration file
    with open(config) as f:
        config_data = yaml.safe_load(f)

    # defining the path to the Jinja2 Snakefile template
    template = Path("Snakefile.j2")

    # setting up the Jinja2 environment to load templates from the template's directory
    env = Environment(loader=FileSystemLoader(template.parent))
    template_obj = env.get_template(template.name)

    # getting the current date and time for metadata
    generation_datetime = datetime.now().isoformat()

    # rendering the template with the configuration data and metadata
    rendered = template_obj.render(
        config=config_data,
        config_path=config,
        generation_datetime=generation_datetime
    )

    # writing the rendered Snakefile to the specified output path
    with open(output, "w") as f:
        f.write(rendered)

    # printing a success message to the user
    typer.echo(f"✅ Snakefile generated at {output}")

if __name__ == "__main__":
    app()