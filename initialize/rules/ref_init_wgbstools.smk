### Genome prepatation for wgbstools use ###

# ref_init_wgbstools: Initialize reference genome for wgbstools
# update this comment to match config needs: change ref.fasta and ref.wgbstools_ref_name in config.yaml to the name of the desired reference genome fasta file and name for wgbstools to use for the reference genome
rule ref_init_wgbstools:
    input:
        fasta_path = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]),

    output:
        directory(expand("{root}/{wgbstools_ref_dir}/{fasta}", root = config["root"], wgbstools_ref_dir = config["ref"]["wgbstools_idx_dir"], fasta = config["ref"]["fasta"]))
    
    log:
        expand("logs/initialize_rules/ref_init_wgbstools--{fasta}.log", fasta=config["ref"]["fasta"])

    conda:
        "../../environment_files/wgbstools.yaml"

    threads: 2
    
    resources:
        mem_mb=2000

    params:
        genome_name = config["ref"]["fasta"],

    shell:
        """
        wgbstools init_genome {params.genome_name} --fasta_path {input.fasta_path} > {log} 2>&1
        """
