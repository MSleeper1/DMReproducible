### Genome reference file download rules ###

# get_ref: download the reference genome from the ref goldenPath specified in the config file
# change ref.goldenPath and ref.fasta in the config file to the desired goldenPath and fasta file name expected
rule get_ref:
    output:
        zipped_fasta = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]),
        unzipped_fasta = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"])
    
    log:
        expand("logs/initialize_rules/get_ref--{fasta}.log", fasta = config["ref"]["fasta"])

    conda:
        "../../environment_files/rsync.yaml"
    
    params:
        gold_path = config["ref"]["goldenPath"], 
        out_dir = expand("{root}/{genomes_dir}/{genome}/", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        rsync_args = config["prep_args"]["rsync_get_ref"]

    shell:
        """
        mkdir -p {params.out_dir}
        echo "downloading {params.gold_path} to {params.out_dir}"
        rsync {params.rsync_args} rsync:{params.gold_path} {params.out_dir} > {log} 2>&1
        echo "unzipping {output.zipped_fasta} to {output.unzipped_fasta}"
        if [ -f {output.zipped_fasta} ]; then
            gunzip -c {output.zipped_fasta} > {output.unzipped_fasta}
        else
            echo "ref file already unzipped"
        fi 
        echo "done"
        """

# get_ref_gtf: download the reference gtf from the gtf goldenPath specified in the config file
# change ref.gtf_goldenPath and ref.gtf in the config file to the desired goldenPath and gtf file name expected
rule get_ref_gtf:
    output:
        zipped_gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        unzipped_gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])
    
    log:
        expand("logs/initialize_rules/get_ref_gtf--{gtf}.log", gtf = config["ref"]["gtf"])

    conda:
        "../../environment_files/rsync.yaml"
    
    params:
        gold_path = config["ref"]["gtf_goldenPath"], 
        out_dir = expand("{root}/{genomes_dir}/{genome}/", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        rsync_args = config["prep_args"]["rsync_get_ref"]

    shell:
        """
        mkdir -p {params.out_dir}
        echo "downloading {params.gold_path} to {params.out_dir}"
        rsync {params.rsync_args} rsync:{params.gold_path} {params.out_dir} > {log} 2>&1
        echo "unzipping {output.zipped_gtf} to {output.unzipped_gtf}"
        gunzip -c {output.zipped_gtf} > {output.unzipped_gtf}
        echo "done"
        """