# a CLI that will find all supported results by MultiQC, renaming them if needed using symlinks, generating one or several MultiQC reports in case where several assembler were used

import os
import typer
import subprocess
import shutil
import pandas as pd


class HandleResults:
    """
    This class contains one method by supported tool, i.e.:

    - Bowtie2
    - FastQC
    - fastp (fastplong is not supported by MultiQC for now: https://github.com/MultiQC/MultiQC/issues/3121)
    - MEGAHIT
    - QUAST
    - CheckM2
    - GTDB-Tk
    - Bakta

    (check the README for the summary of what the tools are doing)

    Each one of these methods collect the paths of the suited results files, and store them.
    To collect everything automatically, use the `collect_all_results` method.
    Once it is done, use the `prepare_results` method to prepare the results for MultiQC.

    It also contains information about found assembly methods (in order to establish the number of MultiQC reports to generate)

    The `generate_report` method will generate the MultiQC report(s) based on the collected results.
    """

    def __init__(self, results_dir, log_dir, output_dir, ani: int = 95, multiqc_config: str = "multiqc_config.yaml"):
        self.results_dir = results_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.ani = ani
        self.multiqc_config = multiqc_config

        # ensuring results_dir/results exists
        if not os.path.isdir(self.results_dir):
            raise FileNotFoundError(f"Required directory not found: {self.results_dir}")

        # recovering the used assemblers
        self.assemblers = []
        assembly_dir = os.path.join(self.results_dir, "03_assembly")
        if os.path.isdir(assembly_dir):
            valid_assemblers = {
                "hybridspades",
                "megahit",
                "metaflye",
                "metaspades",
                "hylight",
            }
            self.assemblers = [
                name
                for name in os.listdir(assembly_dir)
                if os.path.isdir(os.path.join(assembly_dir, name))
                and name in valid_assemblers
            ]
        else:
            self.assemblers = []

    def collect_all_results(self):
        """
        Collect all results from the specified directories.
        This method calls each specific result collection method.
        """

        self.bowtie2()
        self.fastqc()
        self.fastp()
        self.quast()
        self.checkm2()
        self.gtdbtk()
        self.bakta()

    def prepare_results(self):
        """
        Prepare the results for MultiQC.
        This method creates symlinks, with different name if needed, to the results files in a common directory structure.
        """

        # creating a directory for MultiQC results
        multiqc_dir = os.path.join(self.output_dir, "multiqc_results")
        os.makedirs(multiqc_dir, exist_ok=True)

        # processing results
        # Bowtie2 results
        if self.bowtie2_results:
            if len(self.bowtie2_results) > 0:
                bowtie2_dir = os.path.join(multiqc_dir, "bowtie2")
                os.makedirs(bowtie2_dir, exist_ok=True)

                for sample_name, file_path in self.bowtie2_results.items():
                    symlink_path = os.path.join(bowtie2_dir, f"{sample_name}.txt")
                    if not os.path.exists(symlink_path):
                        os.symlink(file_path, symlink_path)

        # FastQC results
        if self.fastqc_results:
            if len(self.fastqc_results) > 0:
                fastqc_dir = os.path.join(multiqc_dir, "fastqc")
                os.makedirs(fastqc_dir, exist_ok=True)

                for sample_name, file_path in self.fastqc_results.items():
                    symlink_path = os.path.join(fastqc_dir, f"{sample_name}_fastqc.zip")
                    if not os.path.exists(symlink_path):
                        os.symlink(file_path, symlink_path)

        # fastp results
        if self.fastp_results:
            if len(self.fastp_results) > 0:
                fastp_dir = os.path.join(multiqc_dir, "fastp")
                os.makedirs(fastp_dir, exist_ok=True)

                for sample_name, file_path in self.fastp_results.items():
                    symlink_path = os.path.join(fastp_dir, f"{sample_name}.json")
                    if not os.path.exists(symlink_path):
                        os.symlink(file_path, symlink_path)

        # QUAST results
        if self.quast_results:
            if len(self.quast_results) > 0:
                quast_dir = os.path.join(multiqc_dir, "quast")
                os.makedirs(quast_dir, exist_ok=True)

                for assembler, samples in self.quast_results.items():
                    assembler_dir = os.path.join(quast_dir, assembler)
                    os.makedirs(assembler_dir, exist_ok=True)

                    for sample_name, file_path in samples.items():
                        sample_dir = os.path.join(assembler_dir, sample_name)
                        os.makedirs(sample_dir, exist_ok=True)
                        dest_path = os.path.join(sample_dir, "report.tsv")
                        if not os.path.exists(dest_path):
                            # we copy the file because we need to modify the header
                            shutil.copy2(file_path, dest_path)

                            # replacing the second column header 'assembly' with the sample name in the copied report.tsv
                            df = pd.read_csv(dest_path, sep="\t")
                            if df.columns.size > 1:
                                cols = list(df.columns)
                                cols[1] = sample_name
                                df.columns = cols
                                df.to_csv(dest_path, sep="\t", index=False)

        # CheckM2 results
        if self.checkm2_results:
            if len(self.checkm2_results) > 0:
                checkm2_dir = os.path.join(multiqc_dir, "checkm2")
                os.makedirs(checkm2_dir, exist_ok=True)

                for assembler, file_path in self.checkm2_results.items():
                    symlink_path = os.path.join(checkm2_dir, f"{assembler}.tsv")
                    if not os.path.exists(symlink_path):
                        os.symlink(file_path, symlink_path)

        # GTDB-Tk results
        if self.gtdbtk_results:
            if len(self.gtdbtk_results) > 0:
                gtdbtk_dir = os.path.join(multiqc_dir, "gtdbtk")
                os.makedirs(gtdbtk_dir, exist_ok=True)

                for assembler, reports in self.gtdbtk_results.items():
                    if isinstance(reports, dict):
                        for key, file_path in reports.items():
                            if key == "bact":
                                symlink_path = os.path.join(
                                    gtdbtk_dir, f"{assembler}_bact.tsv"
                                )
                            elif key == "ar":
                                symlink_path = os.path.join(
                                    gtdbtk_dir, f"{assembler}_ar.tsv"
                                )
                            if not os.path.exists(symlink_path):
                                os.symlink(file_path, symlink_path)

        # Bakta results
        if self.bakta_results:
            if len(self.bakta_results) > 0:
                bakta_dir = os.path.join(multiqc_dir, "bakta")
                os.makedirs(bakta_dir, exist_ok=True)

                for assembler, mags in self.bakta_results.items():
                    assembler_dir = os.path.join(bakta_dir, assembler)
                    os.makedirs(assembler_dir, exist_ok=True)

                    for mag_name, file_path in mags.items():
                        symlink_path = os.path.join(assembler_dir, f"{mag_name}.txt")
                        if not os.path.exists(symlink_path):
                            os.symlink(file_path, symlink_path)

    def bowtie2(self):
        """
        Bowtie2 mapping stats on reference during the decontamination step
        """
        directory = os.path.join(self.log_dir, "02_preprocess/bowtie2/")

        if os.path.exists(directory):
            results = {}
            for file in os.listdir(directory):
                if file.endswith(".stderr") and not "index" in file:
                    sample_name = file.replace(".stderr", "")
                    results[sample_name] = os.path.join(directory, file)
        else:
            results = {}

        self.bowtie2_results = results

    def fastqc(self):
        directory = os.path.join(self.results_dir, "02_preprocess/fastqc")

        if os.path.exists(directory):
            results = {}
            for file in os.listdir(directory):
                if file.endswith(".clean_fastqc.zip"):
                    sample_name = file.replace(".clean_fastqc.zip", "")
                    results[sample_name] = os.path.join(directory, file)
        else:
            results = {}

        self.fastqc_results = results

    def fastp(self):
        directory = os.path.join(self.results_dir, "02_preprocess/fastp")

        if os.path.exists(directory):
            results = {}
            for file in os.listdir(directory):
                if file.endswith("_report.json"):
                    sample_name = file.replace("_report.json", "")
                    results[sample_name] = os.path.join(directory, file)
        else:
            results = {}

        self.fastp_results = results

    def quast(self):
        directory = os.path.join(self.results_dir, "04_assembly_qc/quast")

        if os.path.exists(directory):
            results = {}
            for assembler in self.assemblers:
                results[assembler] = {}
                # searching reports for this assembler
                assembler_dir = os.path.join(directory, assembler)
                if os.path.exists(assembler_dir):
                    for sample in os.listdir(assembler_dir):
                        sample_dir = os.path.join(assembler_dir, sample)
                        if os.path.isdir(sample_dir):
                            report_file = os.path.join(sample_dir, "report.tsv")
                            if os.path.exists(report_file):
                                results[assembler][sample] = report_file
        else:
            results = {}

        self.quast_results = results

    def checkm2(self):
        """
        CheckM2 is used several times in the pipeline, here we only collect results on MAGs produced from all samples
        """
        directory = os.path.join(
            self.results_dir,
            "08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/",
            str(self.ani),
        )

        if os.path.exists(directory):
            results = {}
            for assembler in self.assemblers:
                results[assembler] = {}
                # searching reports for this assembler
                assembler_dir = os.path.join(directory, assembler, "checkm2")
                if os.path.exists(assembler_dir):
                    report_file = os.path.join(assembler_dir, "quality_report.tsv")
                    if os.path.exists(report_file):
                        results[assembler] = report_file
        else:
            results = {}

        self.checkm2_results = results

    def gtdbtk(self):
        directory = os.path.join(
            self.results_dir, "08_bins_postprocessing/gtdb_tk/", str(self.ani)
        )

        if os.path.exists(directory):
            results = {}
            for assembler in self.assemblers:
                results[assembler] = {}
                # searching reports for this assembler
                assembler_dir = os.path.join(directory, assembler)
                if os.path.exists(assembler_dir):
                    report_file1 = os.path.join(
                        assembler_dir, "gtdbtk.bac120.summary.tsv"
                    )
                    if os.path.exists(report_file1):
                        results[assembler]["bact"] = report_file1

                    report_file2 = os.path.join(
                        assembler_dir, "gtdbtk.ar53.summary.tsv"
                    )
                    if os.path.exists(report_file2):
                        results[assembler]["ar"] = report_file2
        else:
            results = {}

        self.gtdbtk_results = results

    def bakta(self):
        directory = os.path.join(
            self.results_dir, "08_bins_postprocessing/bakta", str(self.ani)
        )

        if os.path.exists(directory):
            results = {}
            for assembler in self.assemblers:
                results[assembler] = {}
                # searching reports for this assembler
                assembler_dir = os.path.join(directory, assembler, "annotation")
                if os.path.exists(assembler_dir):
                    for mag in os.listdir(assembler_dir):
                        mag_dir = os.path.join(assembler_dir, mag)
                        if os.path.isdir(mag_dir):
                            report_file = os.path.join(mag_dir, f"{mag}.txt")
                            if os.path.exists(report_file):
                                if assembler not in results:
                                    results[assembler] = {}
                                results[assembler][mag] = report_file
        else:
            results = {}

        self.bakta_results = results

    def generate_report(self, dry_run: bool = False):
        """
        Generates a MultiQC report based on the collected results. 
        It may generates multiple reports if several assemblers were used in parallel
        """

        # there are two possible cases:
        # 1. only one assembler was used, we generate one report
        if len(self.assemblers) == 1:
            # ensuring we have the MultiQC configuration file
            if not os.path.exists(self.multiqc_config):
                raise FileNotFoundError(f"MultiQC config file not found: {self.multiqc_config}")

            multiqc_command = [
                "multiqc",
                "--config",
                self.multiqc_config,
                self.output_dir
            ]

            if dry_run:
                print("Dry run: MultiQC command would be:", " ".join(multiqc_command))
            else:
                print("Running MultiQC command:", " ".join(multiqc_command))
                subprocess.run(multiqc_command)
        else:
            # preparing one command for each assembler
            # it will also collect results assembler-dependent to build the command (eg. the MAGs produced from the given assembler only) 
            for assembler in self.assemblers:
                multiqc_dir = os.path.join(self.output_dir, f"multiqc_report_{assembler}")
                os.makedirs(multiqc_dir, exist_ok=True)

                # folders or files to include for MultiQC report
                folders_to_include = []

                if len(self.bowtie2_results) > 0:
                    folders_to_include.append(os.path.join(self.output_dir, "multiqc_results", "bowtie2"))
                if len(self.fastqc_results) > 0:
                    folders_to_include.append(os.path.join(self.output_dir, "multiqc_results", "fastqc"))
                if len(self.fastp_results) > 0:
                    folders_to_include.append(os.path.join(self.output_dir, "multiqc_results", "fastp"))
                if assembler in self.quast_results and len(self.quast_results[assembler]) > 0:
                    folders_to_include.append(os.path.join(self.output_dir, "multiqc_results", "quast", assembler))
                if assembler in self.checkm2_results:
                    folders_to_include.append(os.path.join(self.output_dir, "multiqc_results", "checkm2", f"{assembler}.tsv"))
                if assembler in self.gtdbtk_results:
                    folder = os.path.join(self.output_dir, "multiqc_results", "gtdbtk")
                    for file in os.listdir(folder):
                        if file.startswith(assembler) and file.endswith(".tsv"):
                            folders_to_include.append(os.path.join(folder, file))
                if assembler in self.bakta_results:
                    folders_to_include.append(os.path.join(self.output_dir, "multiqc_results", "bakta", assembler))
                
                multiqc_command = [
                    "multiqc",
                    "--config",
                    self.multiqc_config,
                    "-o", multiqc_dir
                ] + folders_to_include

                if dry_run:
                    print("Dry run: MultiQC command would be:", " ".join(multiqc_command))
                else:
                    print("Running MultiQC command:", " ".join(multiqc_command))
                    subprocess.run(multiqc_command)

app = typer.Typer()

@app.command()
def generate_report(
    results_dir: str = typer.Option(..., help="Directory containing the results of the pipeline."),
    log_dir: str = typer.Option(..., help="Directory containing the logs of the pipeline."),
    output_dir: str = typer.Option(..., help="Directory where the MultiQC report will be generated."),
    ani: int = typer.Option(95, help="ANI threshold used for dereplicating MAGs."),
    multiqc_config: str = typer.Option("multiqc_config.yaml", help="Path to the MultiQC configuration file."),
    dry_run: bool = typer.Option(False, "--dry-run", "-d", help="If set, only prints the MultiQC command without executing it.")
):
    """
    Generates MultiQC report(s) based on the pipeline results collected from the specified directories.
    """
    handler = HandleResults(results_dir, log_dir, output_dir, ani, multiqc_config)
    handler.collect_all_results()
    handler.prepare_results()
    handler.generate_report(dry_run)

if __name__ == "__main__":
    app()
