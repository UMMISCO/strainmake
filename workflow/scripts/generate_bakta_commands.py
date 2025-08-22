# a CLI that, using GTDB-Tk bacterial summary, generates the Bakta commands to 
# annotate/plot a set of MAGs
# these commands are written to a .txt file whose path is given as an argument

import argparse
import os
import pandas as pd
import shlex

def generate_bakta_commands_annot(gtdb_summary_path: str, genomes_dir: str, out_dir: str, extension: str = ".fa", threads: int = 1, bakta_database: str = "") -> list:
    """ 
    Reads the GTDB-Tk summary file and generates Bakta commands for each MAG.
    """

    gtdb_summary_df = pd.read_csv(gtdb_summary_path, sep="\t")
    bakta_commands = []    

    for index, row in gtdb_summary_df.iterrows():
        mag_name = row['user_genome']
        classification = row['classification']
        taxonomy = classification.split(";")

        # constructing the path to the genome file
        path_to_genome = os.path.join(genomes_dir, f"{mag_name}{extension}")

        # constructing the output directory for the Bakta annotation
        out_dir_genome = os.path.join(out_dir, mag_name)

        # if there is an assigned genus in the annotation we use it
        genus = None
        for taxon in taxonomy:
            if taxon.startswith("g__") and len(taxon) > 3:
                genus = shlex.quote(taxon[3:])
                genus_arg = f"--genus {genus}"
                break
            # if there is no genus, it should be set to empty
            else:
                genus_arg = ""

        # same for species
        species = None
        for taxon in taxonomy:
            if taxon.startswith("s__") and len(taxon) > 3:
                # using shlex.quote to safely escape the species name for shell commands
                species_full = taxon[3:]
                species_parts = species_full.split()
                if species_parts and genus and species_parts[0] == genus:
                    species = shlex.quote(" ".join(species_parts[1:]))
                else:
                    species = shlex.quote(species_full)
                species_arg = f"--species {species}"
                break
            # if there is no species, it should be set to empty
            else:
                species_arg = ""

        # construct the command, in case there is no genus or species, nothing is added since it would be empty
        command = f"bakta {genus_arg} {species_arg} --skip-plot --output {out_dir_genome} --prefix {mag_name} --threads {threads} --verbose --db {bakta_database} {path_to_genome}"
        bakta_commands.append(command)


    return bakta_commands

def generate_bakta_commands_plot(gtdb_summary_path: str, bakta_annot_dir: str, out_dir: str) -> list:
    """
    Reads the GTDB-Tk summary file and generates Bakta commands for plotting genomes.
    """

    gtdb_summary_df = pd.read_csv(gtdb_summary_path, sep="\t")
    bakta_plot_commands = []

    for index, row in gtdb_summary_df.iterrows():
        mag_name = row['user_genome']

        # constructing the output directory for the Bakta annotation
        out_dir_genome = os.path.join(out_dir, mag_name)
    
        # we should have the following Bakta annotation output (JSON)
        results_json = os.path.join(bakta_annot_dir, mag_name, f"{mag_name}.json")

        if os.path.exists(results_json):
            # constructing the command for plotting
            command = f"bakta_plot --output {out_dir_genome} --type cog --verbose {results_json}"
            bakta_plot_commands.append(command)
        else:
            raise FileNotFoundError(f"Bakta annotation results for {mag_name} not found at {results_json}")
        
    return bakta_plot_commands

def main():
    """
    CLI logic
    """
    parser = argparse.ArgumentParser(description="Generate Bakta commands from GTDB-Tk summary.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # subcommand for generating Bakta commands for annotation
    bakta_parser = subparsers.add_parser("bakta_annot", help="Generate Bakta commands from GTDB-Tk summary")
    bakta_parser.add_argument("--gtdb_tk", required=True, help="Path to GTDB-Tk summary file")
    bakta_parser.add_argument("--genomes_dir", required=True, help="Directory containing genome files")
    bakta_parser.add_argument("--output_commands", required=True, help="Output file to write Bakta commands")
    bakta_parser.add_argument("--extension", required=True, help="Genome file extension (e.g., .fa)")
    bakta_parser.add_argument("--bakta_database", required=True, help="Path to the Bakta database directory")
    bakta_parser.add_argument("--threads", type=int, default=1, help="Number of threads for Bakta (default: 1)")
    bakta_parser.add_argument("--output_dir", default=".", help="Output directory for Bakta annotations (default: current directory)")

    # subcommand for generating Bakta commands for plotting genomes
    plot_parser = subparsers.add_parser("bakta_plot", help="Generate Bakta commands for plotting genomes from GTDB-Tk summary")
    plot_parser.add_argument("--gtdb_tk", required=True, help="Path to GTDB-Tk summary file")
    plot_parser.add_argument("--bakta_annot_dir", required=True, help="Directory containing Bakta annotation results (command `bakta_annot` output)")
    plot_parser.add_argument("--output_commands", required=True, help="Output file to write Bakta plot commands")
    plot_parser.add_argument("--output_dir", default=".", help="Output directory for Bakta plots (default: current directory)")

    args = parser.parse_args()

    if args.command == "bakta_annot":
        commands = generate_bakta_commands_annot(
            gtdb_summary_path=args.gtdb_tk,
            genomes_dir=args.genomes_dir,
            out_dir=args.output_dir,
            extension=args.extension,
            threads=args.threads,
            bakta_database=args.bakta_database
        )

        with open(args.output_commands, "w") as f:
            for cmd in commands:
                f.write(cmd.strip() + "\n")
    elif args.command == "bakta_plot":
        commands = generate_bakta_commands_plot(
            gtdb_summary_path=args.gtdb_tk,
            bakta_annot_dir=args.bakta_annot_dir,
            out_dir=args.output_dir
        )
        
        with open(args.output_commands, "w") as f:
            for cmd in commands:
                f.write(cmd.strip() + "\n")


if __name__ == "__main__":
    main()