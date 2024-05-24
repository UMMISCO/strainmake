# this rule is to remove: one should have already download the data before running the script, since that is 
# very heavy

# rule gtdb_tk_download_ref_data:
#     output: directory("results/08_bins_postprocessing/gtdb_tk/reference/data/release220")
#     log:
#         stdout = "logs/08_bins_postprocessing/gtdb_tk/data_download.stdout",
#         stderr = "logs/08_bins_postprocessing/gtdb_tk/data_download.stderr"
#     params:
#         ref_data_online = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz",
#         save_location = "results/08_bins_postprocessing/gtdb_tk/reference/data/gtdbtk_data"
#     threads: 2
#     shell:
#         """
#         curl -C - -L --retry 10 \
#             --output {params.save_location}.tar.gz {params.ref_data_online} \
#             > {log.stdout} 2> {log.stderr} \
#         && \
#         tar xvzf {params.save_location}.tar.gz -C results/08_bins_postprocessing/gtdb_tk/reference/data/ \
#         """

rule gtdb_tk_taxonomic_annotation:
    input:
        # folder with refined bins
        refined_bins = "results/07_bins_refinement/binette/{assembler}/{sample}",
        ref_data = "data/08_bins_postprocessing/gtdb_tk/reference/data/release220"
    output: directory("results/08_bins_postprocessing/gtdb_tk/{assembler}/{sample}")
    conda:
        "../envs/gtdb_tk.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/gtdb_tk/{assembler}/{sample}/classify.stdout",
        stderr = "logs/08_bins_postprocessing/gtdb_tk/{assembler}/{sample}/classify.stderr"
    params:
        threads = config['bins_postprocessing']['gtdbtk']['threads']
    threads: 4
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.refined_bins}/final_bins --cpus {params.threads} --out_dir {output} \
            --extension ".fa" \
            --skip_ani_screen \
            > {log.stdout} 2> {log.stderr}
        """
