#!/usr/bin/env snakemake --cluster-config ../cluster.yaml --cluster "sbatch --parsable --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-per-task={cluster.cpus-per-task} --output={cluster.output} --error={cluster.error} --job-name={cluster.name}" --jobs 20 --use-conda --rerun-incomplete --printshellcmds
# Garden usage: snakemake --profile garden --use-conda --rerun-incomplete --printshellcmds --rerun-triggers mtime --jobs 2 -s sub_sample_fq.smk

#### snakemake subsample subworkflow pipeline ####
# author: Meghan M. Sleeper
# run on farm cluster with:  snakemake --cluster-config ../cluster.yaml --cluster "sbatch --parsable --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-per-task={cluster.cpus-per-task} --output={cluster.output} --error={cluster.error} --job-name={cluster.name}" --jobs 20 --use-conda --rerun-incomplete --printshellcmds --rerun-triggers mtime --dry-run
# run on garden biolab with: snakemake --profile garden --use-conda --rerun-incomplete --printshellcmds --rerun-triggers mtime --dry-run
# run as a subworkflow for the DMReproducible pipeline 
# run before secondary if you need to downsample for a more managable file size
# purpose:
#    - to downsample fq files downloaded by sra in initialization step
# input: 
#    - tsv file specifying data accession numbers and associated metadata for wgbs data on SRA
#    - config file specifying reference genome and file paths
# output: 
#    - downsampled fastq files
#
# when using the output downsampled files with secondary workflow, 
# you must put 'sub100000' trailing ref in tsv info sheet

#### import modules ####
import pandas as pd
import initialize_helper_functions as hf

#### assign config ####
configfile: "../config.garden.yaml"
# configfile: "../config.farm.yaml"

#### sample info ####
sample_info = hf.get_sample_info_df(config["root"] + "/" + config["samples_tsv"])
sample_info_se = sample_info[sample_info['layout'] == 'se']
sample_info_pe = sample_info[sample_info['layout'] == 'pe']

#### set up conditional files for rule all input ####
rule_all_input_list = [
    expand("{root}/{data_dir}/01_raw_sequence_files/{se.ref}sub100000--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq", root = config["root"], data_dir=config["data_dir"], se=sample_info_se.itertuples()), # sra_get_data se output
    expand("{root}/{data_dir}/01_raw_sequence_files/{pe.ref}sub100000--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", root = config["root"], data_dir=config["data_dir"], pe=sample_info_pe.itertuples(), suf={"_1.fastq", "_2.fastq"}), # sra_get_data pe R1 and R2 output
    ]

#### default rule ####
rule all:
     input:
          rule_all_input_list

# downsampling rules
rule downsample_fq_se:
    input:
        fq_file = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", root = config["root"], data_dir = config["data_dir"])
    output:
        sub_fq_file = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}sub100000--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", root = config["root"], data_dir = config["data_dir"])
    log:
        "logs/additional_analysis/subset100000-{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.log"
    conda:
        "../environment_files/seqtk.yaml"
    wildcard_constraints:
        layout = "se"
    params:
        outdir = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}sub100000--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/", root = config["root"], data_dir = config["data_dir"]),
        subset_num = 100000
    shell:
        """
        mkdir -p {params.outdir} 2>>{log}
        seqtk sample {input.fq_file} {params.subset_num} > {output.sub_fq_file}
        echo "done subsampling" >>{log}
        """

rule downsample_fq_pe:
    input:
        r1 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", root = config["root"], data_dir = config["data_dir"]),
        r2 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", root = config["root"], data_dir = config["data_dir"])
    output:
        sub_r1 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}sub100000--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", root = config["root"], data_dir = config["data_dir"]),
        sub_r2 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}sub100000--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", root = config["root"], data_dir = config["data_dir"])
    log:
        "logs/additional_analysis/subset100000-{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.log"
    conda:
        "../environment_files/seqtk.yaml"
    wildcard_constraints:
        layout = "pe"
    params:
        outdir = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}sub100000--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/", root = config["root"], data_dir = config["data_dir"]),
        subset_num = 100000
    shell:
        """
        mkdir -p {params.outdir} 2>>{log}
        seqtk sample -s100 {input.r1} {params.subset_num} > {output.sub_r1}
        seqtk sample -s100 {input.r2} {params.subset_num} > {output.sub_r2}
        echo "done subsampling" >>{log}
        """
        
onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")
