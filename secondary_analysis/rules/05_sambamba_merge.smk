
### sambamba merge rule for deduped bam files ###
# sambamba is a tool for working with SAM/BAM files. It is faster than samtools and can be used to merge bam files.

# sambamba_merge: determines if merging is needed based on SRX ID. 
# If only one file is associated with the SRX ID, a symlink is created. 
# If more than one file is associated with the SRX ID, the files are merged using sambamba merge. 
# The merged file is then used for the next step in the pipeline.
# input: deduped bam files (04)
# output: merged bam file (05)

# # has been tested with files that do not need to be merged
# # still need to test rule for sambamba merge with files that need to be merged
# symlinks were experiencing latency wait time errors. To avoid specifying latency wait time for all jobs, I added a sleep 10 seconds command to the shell script. This will allow the symlink to be created before the next rule_all is executed.

rule sambamba_merge_bwameth:
    input:
       bams = lambda wildcards: expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{accession}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"], accession = sample_info[sample_info["srx_id"] == wildcards.srx_id]["accession"].tolist())
    
    output:
        merged_bam = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam.bai", root = config["root"], data_dir=config["data_dir"])

    log:
        "logs/secondary_rules/05_sambamba_merge_bwameth/05_sambamba_merge_bwameth-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/sambamba.yaml"

    # shadow:
    #     "shallow"

    threads: 3
    
    wildcard_constraints:
        srx_id = "|".join(sample_info["srx_id"].tolist()),
        accession = "|".join(sample_info["accession"].tolist())

    params:
        outdir = expand("{root}/{data_dir}/05_merged_sambamba_bwa", root = config["root"], data_dir=config["data_dir"])

    shell:
        """
        echo "Adding files to array..." > {log}
        FILES=()
        for i in {input.bams}; do
            FILES+=($i)
        done
        echo "Checking number of files in array..." >> {log}
        input_len=${{#FILES[@]}}
        echo "Number of files in array: $input_len" >> {log}
        echo "Merge files if more than one file is present..." >> {log}
        if [ $input_len -gt 1 ]; then
            echo "Merging files: {input.bams} into {output.merged_bam} with sambamba merge..." >> {log}
            sambamba merge -t {threads} {output.merged_bam} {input.bams} >> {log} 2>> {log}
            sambamba index -t {threads} {output.merged_bam} {output.bai} >> {log} 2>> {log}
        else
            echo "Only one file, no need to merge. Creating symlink for file at expected output location: {output.merged_bam}" >> {log}
            cp -f {input.bams} {output.merged_bam}
            cp -f {input.bams}.bai {output.bai}
        fi
        echo "Done" >> {log}
        sleep 10
        """

rule sambamba_merge_bismark:
    input:
        bams = lambda wildcards: expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{accession}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"], accession = sample_info[sample_info["srx_id"] == wildcards.srx_id]["accession"].tolist())
    
    output:
        merged_bam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])

    log:
        "logs/secondary_rules/05_sambamba_merge_bismark/05_sambamba_merge_bismark-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/sambamba.yaml"

    # shadow:
    #     "shallow"

    threads: 3
    
    wildcard_constraints:
        srx_id = "|".join(sample_info["srx_id"].tolist()),
        accession = "|".join(sample_info["accession"].tolist())

    params:
        outdir = expand("{root}/{data_dir}/05_merged_sambamba_bis", root = config["root"], data_dir=config["data_dir"])

    shell:
        """
        echo "Adding files to array..." > {log}
        FILES=()
        for i in {input.bams}; do
            FILES+=($i)
        done
        echo "Checking number of files in array..." >> {log}
        input_len=${{#FILES[@]}}
        echo "Number of files in array: $input_len" >> {log}
        echo "Merge files if more than one file is present..." >> {log}
        if [ $input_len -gt 1 ]; then
            echo "Merging files: {input.bams} into {output.merged_bam} with sambamba merge..." >> {log}
            sambamba merge -t {threads} {output.merged_bam} {input.bams} >> {log} 2>> {log}
            sambamba index -t {threads} {output.merged_bam} {output.merged_bam}.bai >> {log} 2>> {log}
        else
            echo "Only one file, no need to merge. Creating symlink for file at expected output location: {output.merged_bam}" >> {log}
            cp -f {input.bams} {output.merged_bam}
        fi
        echo "Done" >> {log}
        sleep 10
        """