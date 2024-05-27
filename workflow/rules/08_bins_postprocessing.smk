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
    threads: config['bins_postprocessing']['gtdbtk']['threads']
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.refined_bins}/final_bins --cpus {threads} --out_dir {output} \
            --extension ".fa" \
            --skip_ani_screen --pplacer_cpus 1 \
            > {log.stdout} 2> {log.stderr}
        """
