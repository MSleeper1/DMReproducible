### Genome preparation for bismark ###

# ref_index_bismark: index the reference genome for bismark
# change ref.dir and ref.fasta in config.yaml to use a different reference genome
rule ref_index_bismark:
    input:
        fasta_path = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]),

    output:
        directory(expand("{root}/{genomes_dir}/{genome}/bismark/Bisulfite_Genome/", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"])),
        expand("{root}/{genomes_dir}/{genome}/bismark/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]),
        expand("{root}/{genomes_dir}/{genome}/bismark/genomic_nucleotide_frequencies.txt", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),

    log:
        stdout = expand("logs/initialize_rules/ref_index_bismark--{fasta}.out", fasta = config["ref"]["fasta"]),
        stderr = expand("logs/initialize_rules/ref_index_bismark--{fasta}.err", fasta = config["ref"]["fasta"])
    
    conda:
        "../../environment_files/bismark.yaml"
    
    # shadow:
    #     "shallow"

    params:
        bismk_args = config["prep_args"]["bismark_genome_prep"],
        bismk_dir = expand("{root}/{genomes_dir}/{genome}/bismark", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        bismk_fasta_path = expand("{root}/{genomes_dir}/{genome}/bismark/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"])

    shell:
        """
        mkdir -p {params.bismk_dir} 
        cp {input.fasta_path} {params.bismk_fasta_path} 
        echo "indexing {params.bismk_fasta_path} for bismark"
        bismark_genome_preparation {params.bismk_args} {params.bismk_dir} > {log.stdout} 2> {log.stderr}
        echo "prepping genomic nucleotide frequencies"
        bam2nuc --genomic_composition_only --genome_folder {params.bismk_dir}
        echo "done"
        """