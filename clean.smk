rule clean_all_genome_preps:
    params:
        genome_preps = expand("{genomes_dir}", fasta = config["genomes_dir"])
    shell: "rm -rf {params.genome_preps}"

rule clean_bwa_index:
    params:
        bwa_index = expand("{idx_dir}", idx_dir = config["ref"]["bwa_idx_dir"])
    shell: "rm -rf {params.bwa_index}"

rule clean_bismark_index:
    params:
        bismark_index = expand("{idx_dir}", idx_dir = config["ref"]["bismark_idx_dir"])
    shell: "rm -rf {params.bismark_index}"

rule clean_raw_sequence_files:
    params: 
        raw_seq_files = expand("{data_dir}/raw_sequence_files/", data_dir = config["data"]["dir"])
    shell: "rm -rf {params.raw_seq_files}"

rule clean_logs:
    params:
        rule_logs_dir = "../logs/" # note: cluster logs are in a different location (cluster logs can be found in the same directory the snakefile was run)
    shell: "rm -rf {params.rule_logs_dir}"
