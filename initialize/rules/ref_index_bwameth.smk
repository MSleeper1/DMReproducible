### Genome preparation for bwameth ###

# ref_index_bwameth: index the reference genome for bwameth
# change ref.dir and ref.fasta in config.yaml to use a different reference genome
rule ref_index_bwameth:
    input:
        fasta_path = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]),

    output:
        bwa_out_files = expand("{root}/{genomes_dir}/{genome}/bwameth/{fasta}.fa.gz.bwameth.{suf}", suf=["c2t","c2t.bwt","c2t.pac","c2t.ann","c2t.amb","c2t.sa"], root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]),
        bwa_fasta = expand("{root}/{genomes_dir}/{genome}/bwameth/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"])

    log:
        expand("logs/initialize_rules/ref_index_bwameth--{fasta}.log", fasta = config["ref"]["fasta"])

    conda:
        "../../environment_files/bwameth.yaml"

    threads: 4

    resources:
        mem_mb=8000

    params:
        bwa_dir = expand("{root}/{genomes_dir}/{genome}/bwameth/", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"])

    shell:
        """
        mkdir -p {params.bwa_dir}
        echo "copying {input.fasta_path} to {output.bwa_fasta}"
        cp {input.fasta_path} {output.bwa_fasta}
        echo "indexing {output.bwa_fasta} for bwameth mapping"
        bwameth.py index {output.bwa_fasta} > {log} 2>&1
        echo "done"
        """