### SRA download rules ###
# This file contains the rules to download the raw sequence files from the SRA database using the fasterq-dump tool.

# get_samples_se: rule to download single-end raw sequence files from the SRA database.
rule get_samples_se:
    output:
        expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", root = config["root"], data_dir = config["data_dir"])
    
    log:
       "logs/initialize_rules/get_samples_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/sra-download.yaml"

    params:
        temp_dir = expand("{root}/{data_dir}/temp/get_samples_se-{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}-{{accession}}", root = config["root"], data_dir = config["data_dir"]),
        output_dir = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir = config["data_dir"])

    wildcard_constraints:
        layout = "se"

    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {params.temp_dir} 
        echo "downloading {wildcards.accession} to {params.output_dir}"
        fasterq-dump --temp {params.temp_dir} -O {params.output_dir} {wildcards.accession} > {log} 2>&1
        echo "done"
        echo "removing temp directory"
        rm -rf {params.temp_dir}
        echo "done"
        """

# get_samples_pe: rule to download paired-end raw sequence files from the SRA database.
rule get_samples_pe:
    output:
        r1 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", root = config["root"], data_dir = config["data_dir"]),
        r2 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", root = config["root"], data_dir = config["data_dir"])
    
    log:
        "logs/initialize_rules/get_samples_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/sra-download.yaml"

    params:
        temp_dir = expand("{root}/{data_dir}/temp/get_samples_se-{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}-{{accession}}", root = config["root"], data_dir = config["data_dir"]),
        output_dir = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir = config["data_dir"])

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        mkdir -p {params.output_dir} 
        mkdir -p {params.temp_dir} 
        echo "downloading {wildcards.accession} to {params.output_dir}" > {log}
        fasterq-dump --temp {params.temp_dir} -O {params.output_dir} {wildcards.accession} >> {log} 2>&1
        echo "done"
        echo "removing temp directory" >> {log}
        rm -rf {params.temp_dir}
        echo "done" >> {log}
        """