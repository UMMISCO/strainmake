rule gtdb_tk_taxonomic_annotation:
    input:
        # folder with refined bins
        refined_bins = "results/07_bins_refinement/binette/{assembler}/{sample}",
        ref_data = "data/gtdb_tk/release220"
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
