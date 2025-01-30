### BWA-METH alignment rules ###
# bwa-meth is a tool for aligning bisulfite-treated reads to a reference genome. It is a wrapper around BWA and is designed to work with WGBS data.
# samtools is used to convert the sam file to a bam file and the sam file is removed to save space
# input: trimmed fastq files and reference genome
# output: aligned bam files and alignment report (captured standard error report of bwameth)
# LOGGING STANDARD OUT OF BWA-METH WILL CAPTURE ALL OUTPUTS BECAUSE BWAMETH.PY SENDS OUTPUT FILE TO STDOUT

# bwameth_mapping_se_pipe: align single-end reads to reference genome using bwameth
rule bwameth_mapping_se:
    input:
        read_se = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", root = config["root"], data_dir=config["data_dir"]),
        index = expand("{root}/{genomes_dir}/{genome}/bwameth/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"])

    output:
        bam = expand("{root}/{data_dir}/03_aligned_bwameth/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", root = config["root"], data_dir=config["data_dir"]),
        bwa_report = expand("{root}/{rep_dir}/03_bwameth/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bwameth_report.txt", root = config["root"], rep_dir=config["reports_dir"])
         
    log: 
        "logs/secondary_rules/03_bwameth_mapping_se/03_bwameth_mapping_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
        
    conda: 
        "../../environment_files/bwameth.yaml"
    
    # shadow: 
    #     "shallow"

    threads: 3

    params: 
        accession = "{accession}",
        sam = expand("{root}/{data_dir}/03_aligned_bwameth/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.sam", root = config["root"], data_dir=config["data_dir"])

    shell:
        """
        echo "aligning {input.read_se} to {input.index}" > {log}
        bwameth.py --threads {threads} --reference {input.index} {input.read_se} > {params.sam} 2> {output.bwa_report}
        echo "done with alignment" >> {log}
        echo "converting {params.sam} to {output.bam}" >> {log}
        samtools view -S -b {params.sam} > {output.bam}
        echo "alignment and conversion complete for {params.accession}" >> {log} 
        echo "removing {params.sam} file to save space" >> {log}
        rm {params.sam} >> {log} 2>> {log}
        echo "done" >> {log}
        """

# bwameth_mapping_pe: align paired-end reads to reference genome using bwameth
rule bwameth_mapping_pe:
    input:
        read_1 = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_val_1.fq", root = config["root"], data_dir=config["data_dir"]),
        read_2 = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_val_2.fq", root = config["root"], data_dir=config["data_dir"]),
        index = expand("{root}/{genomes_dir}/{genome}/bwameth/{fasta}.fa.gz", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"])

    output:
        bam = expand("{root}/{data_dir}/03_aligned_bwameth/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", root = config["root"], data_dir=config["data_dir"]),
        bwa_report = expand("{root}/{rep_dir}/03_bwameth/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bwameth_report.txt", root = config["root"], rep_dir=config["reports_dir"])
         
    log:
        "logs/secondary_rules/03_bwameth_mapping_pe/03_bwameth_mapping_pe-{ref}--{patient_id}-{group}-{layout}-{srx_id}-{accession}.log"

    conda:
        "../../environment_files/bwameth.yaml"

    # shadow: 
    #     "shallow"

    threads: 6

    params: 
        accession = "{accession}",
        sam = expand("{root}/{data_dir}/03_aligned_bwameth/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.sam", root = config["root"], data_dir=config["data_dir"])

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        echo "aligning {input.read_1} and {input.read_2} to {input.index}" > {log}
        bwameth.py --threads {threads} --reference {input.index} {input.read_1} {input.read_2} > {params.sam} 2> {output.bwa_report}
        echo "done with alignment" >> {log}
        echo "converting {params.sam} to {output.bam}" >> {log}
        samtools view -S -b {params.sam} > {output.bam}
        echo "alignment and conversion complete for {params.accession}" >> {log} 
        echo "removing {params.sam} file to save space" >> {log}
        rm {params.sam} >> {log} 2>> {log}
        echo "done" >> {log}
        """
