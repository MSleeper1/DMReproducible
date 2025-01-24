# Cleans up large intermediate files produced by the secondary analysis
# workflow. This is a separate rule to avoid deleting files that are still being used
# Removes all intermediate data files not including initial fastq files or final methylation data files
# only removes files that are specified by sample_info file

# Note that `snakemake clean` for secondary workflow will remove all output files created by the workflow
# This rule can be used to only remove intermediate files and directories

#### import modules ####
import pandas as pd
import secondary_helper_functions as hf

#### assign config ####
configfile: "../config.yaml"

#### sample info ####
sample_info = hf.get_sample_info_df(config["root"] + "/" + config["samples_tsv"])
sample_info_se = sample_info[sample_info['layout'] == 'se']
sample_info_pe = sample_info[sample_info['layout'] == 'pe']

# clean up intermediate files
secondary_data_files = [
    expand("{root}/{data_dir}/02_trimmed_trim_galore/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed.fq", root = config["root"], data_dir=config["data_dir"], se=sample_info_se.itertuples()), # trim_galore se fastq output
    expand("{root}/{data_dir}/02_trimmed_trim_galore/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed.fq", root = config["root"], data_dir=config["data_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # trim_galore pe reports
    expand("{root}/{data_dir}/03_aligned_bismark_bwt2/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2.bam", root = config["root"], data_dir=config["data_dir"], se=sample_info_se.itertuples()), # mapped bam output
    expand("{root}/{data_dir}/03_aligned_bismark_bwt2/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2.bam", root = config["root"], data_dir=config["data_dir"], pe=sample_info_pe.itertuples()), # mapped bam output
    expand("{root}/{data_dir}/03_aligned_bwameth/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed.bam", root = config["root"], data_dir=config["data_dir"], se=sample_info_se.itertuples()), # bwameth_mapping se output
    expand("{root}/{data_dir}/03_aligned_bwameth/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed.bam", root = config["root"], data_dir=config["data_dir"], pe=sample_info_pe.itertuples()), # bwameth_mapping pe output
    expand("{root}/{data_dir}/04_bismark_deduped/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()), # bismark deduped bam output
    expand("{root}/{data_dir}/04_deduped_sambamba/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()), # sambamba deduped bam output
    expand("{root}/{data_dir}/04_deduped_sambamba/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.bai", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()), # sambamba deduped bam index output
    expand("{root}/{data_dir}/05_merged_sambamba_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/05_merged_sambamba_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.bai", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/05_merged_sambamba_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/06_merged_deduped_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.bai", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/06_merged_deduped_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/06_merged_deduped_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.bai", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples())
]

rule clean_secondary_int_files:
    input:
        secondary_data_files
    output:
        touch("{root}/{data_dir}/clean_secondary_int_files.done", root = config["root"], data_dir=config["data_dir"])
    log:
        "logs/clean/clean_secondary_int_files.log"
    shell:
        """
        echo -e "Cleaning up intermediate files produced by the secondary analysis workflow\n" > {log}
        echo "Processed on:" >> {log}
        echo $(date) >> {log}
        echo -e "\nFiles to be removed:" >> {log}
        echo "{input}" >> {log}
        rm -f {input}
        echo -e "\nFinished" >> {log}
        cp {log} {output}
        """

rule all:
    input:
        expand("{root}/{data_dir}/clean_secondary_int_files.done", root = config["root"], data_dir=config["data_dir"])
