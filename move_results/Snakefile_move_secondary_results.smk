#!/usr/bin/env snakemake --cluster-config ../cluster.yaml --cluster "sbatch --parsable --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-per-task={cluster.cpus-per-task} --output={cluster.output} --error={cluster.error} --job-name={cluster.name}" --jobs 20 --use-conda --rerun-incomplete --printshellcmds

#### snakemake move data subworkflow pipeline ####
# Author: Meghan M. Sleeper
# run on farm cluster with:  snakemake -s Snakefile_move_secondary_results.smk --cluster-config ../cluster.yaml --cluster "sbatch --parsable --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-per-task={cluster.cpus-per-task} --output={cluster.output} --error={cluster.error} --job-name={cluster.name}" --jobs 20 --use-conda --rerun-incomplete --printshellcmds --rerun-triggers mtime --dry-run
# subworkflow for the DMReproducible pipeline that writes scripts for scp of files from remote to local
# Note: `ls -d -1 "$PWD/"215* | tr '\n' ','` can be used to get a comma separated list of files with absolute path for use iwth scp (in this example, 215 is the ref)
# The comma separated list can be used in the scp command to move multiple files at once as used in the shell commands below
# Output: bash scripts for moving files from the farm to the local drive

#### import modules ####
import pandas as pd
import secondary_helper_functions as hf

#### assign config ####
configfile: "config_move_results.yaml"

#### sample info ####
sample_info = hf.get_sample_info_df(config["root"] + "/" + config["samples_tsv"])
sample_info_se = sample_info[sample_info['layout'] == 'se']
sample_info_pe = sample_info[sample_info['layout'] == 'pe']

#### Set up input file lists ####

# File inputs 
pre_merge_files = [
    # expand("{root}/{rep_dir}/01_fastqc/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], se=sample_info_se.itertuples()),
    # expand("{root}/{rep_dir}/01_fastqc/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_1_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], pe=sample_info_pe.itertuples()),
    # expand("{root}/{rep_dir}/01_fastqc/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_2_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], pe=sample_info_pe.itertuples()),
    # expand("{root}/{rep_dir}/01_fastq_screen/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["txt", "html"], se=sample_info_se.itertuples()),
    # expand("{root}/{rep_dir}/01_fastq_screen/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_1_screen.{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["txt", "html"], pe=sample_info_pe.itertuples()),
    # expand("{root}/{rep_dir}/01_fastq_screen/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_2_screen.{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["txt", "html"], pe=sample_info_pe.itertuples()),      
    # expand("{root}/{rep_dir}/02_trim_galore/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se trimming report
    # expand("{root}/{rep_dir}/02_fastqc_post_trim/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc report post-trim
    # expand("{root}/{rep_dir}/02_fastqc_post_trim/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_val_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # trimmed fastq outputs
    # expand("{root}/{rep_dir}/02_trim_galore/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # fastqc reports post-trim
    # expand("{root}/{rep_dir}/03_bwameth/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bwameth_report.txt", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bwameth_mapping se report from stderr
    # expand("{root}/{rep_dir}/03_bwameth/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bwameth_report.txt", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bwameth_mapping pe report from stderr 
    # expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # sambamba deduplication report
    # expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.flagstat", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    # expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    # expand("{root}/{rep_dir}/04_feature_counts_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()),
    # expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    # expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # produced unless --no-per-base specified
    # expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz.csi", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # produced unless --no-per-base specified
    # expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # this named output is required for prefix parsing
    expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # sambamba deduplication report
]
post_merge_files = [
    expand("{root}/{rep_dir}/06_merged_deduped_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.log", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/06_feature_counts_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()), # featureCounts bwa
    expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], sample = sample_info.itertuples()), # fastqc report bwa
    expand("{root}/{rep_dir}/06_samtools_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.stats", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bwa
    expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.flagstat.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
    expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
    expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
]
# Directory inputs  
pre_merge_dirs = [
    expand("{root}/{rep_dir}/01_fastqc/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir = config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/01_fastq_screen/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/02_trim_galore/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/02_fastqc_post_trim/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/03_bwameth/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"],  sample = sample_info.itertuples()), 
    expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/04_feature_counts_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/04_qualimap_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
]
post_merge_dirs = [
    expand("{root}/{rep_dir}/06_qualimap_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), 
    expand("{root}/{data_dir}/07_wgbstools_betas_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mbias", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples())
]
### Data file inputs including bams and methylation data from wgbstools and methyldackel ###
bam_files = [
    expand("{root}/{data_dir}/06_merged_deduped_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/06_merged_deduped_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.bai", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples())
]
wgbstools_betas = [
    expand("{root}/{data_dir}/07_wgbstools_betas_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["pat.gz", "pat.gz.csi", "beta"], sample = sample_info.itertuples()), # wgbstools beta output
]
methyldackel = [
    expand("{root}/{data_dir}/07_methyldackel_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_CpG.bedGraph", root = config["root"], sample = sample_info.itertuples(), data_dir=config["data_dir"]),
    expand("{root}/{data_dir}/07_methyldackel_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_mergedContext_CpG.bedGraph", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"]),
    expand("{root}/{data_dir}/07_methyldackel_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_CpG.methylKit", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"]),
    expand("{root}/{data_dir}/07_methyldackel_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.cytosine_report.txt", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"])
]

dirs_to_tree = [
    expand("{root}/{rep_dir}/01_fastqc/", root = config["root"], rep_dir = config["reports_dir"]),
    expand("{root}/{rep_dir}/01_fastq_screen/", root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{rep_dir}/02_trim_galore/", root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{rep_dir}/02_fastqc_post_trim/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/03_bwameth/", root = config["root"], rep_dir=config["reports_dir"],  sample = sample_info.itertuples()), 
    expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # sambamba deduplication report
    expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/04_feature_counts_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/04_qualimap_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/06_qualimap_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), 
    expand("{root}/{rep_dir}/06_feature_counts_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # featureCounts bwa
    expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa/", root = config["root"], rep_dir = config["reports_dir"], sample = sample_info.itertuples()), # fastqc report bwa
    expand("{root}/{rep_dir}/06_samtools_post_merge_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bwa
    expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bwa/", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{data_dir}/06_merged_deduped_bwa/", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/07_wgbstools_betas_bwa/", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()), # wgbstools beta output
    expand("{root}/{data_dir}/07_methyldackel_bwa/", root = config["root"], sample = sample_info.itertuples(), data_dir=config["data_dir"]),
]

file_dirs_to_tree = [
    expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # sambamba deduplication report
    expand("{root}/{rep_dir}/04_feature_counts_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/06_merged_deduped_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()),
    expand("{root}/{rep_dir}/06_feature_counts_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # featureCounts bwa
    expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa/", root = config["root"], rep_dir = config["reports_dir"], sample = sample_info.itertuples()), # fastqc report bwa
    expand("{root}/{rep_dir}/06_samtools_post_merge_bwa/", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bwa
    expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bwa/", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    expand("{root}/{data_dir}/06_merged_deduped_bwa/", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
    expand("{root}/{data_dir}/07_wgbstools_betas_bwa/", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()), # wgbstools beta output
    expand("{root}/{data_dir}/07_methyldackel_bwa/", root = config["root"], sample = sample_info.itertuples(), data_dir=config["data_dir"]),
]

#### default rule ####
rule all:
    input:
        expand("{root}/scratch/move_files/{sample.ref}_scp_pre_merge_files.sh", root = config["root"], sample = sample_info.itertuples()), # scp commands for files to be moved to the local drive
        expand("{root}/scratch/move_files/{sample.ref}_scp_post_merge_files.sh", root = config["root"], sample = sample_info.itertuples()), # scp commands for dirs to be moved to the local drive
        expand("{root}/scratch/move_files/{sample.ref}_scp_pre_merge_dirs.sh", root = config["root"], sample = sample_info.itertuples()), # scp commands for files to be moved to the local drive
        expand("{root}/scratch/move_files/{sample.ref}_scp_post_merge_dirs.sh", root = config["root"], sample = sample_info.itertuples()), # scp commands for dirs to be moved to the local drive
        expand("{root}/scratch/move_files/{sample.ref}_scp_bams_commands.sh", root = config["root"], sample = sample_info.itertuples()),
        expand("{root}/scratch/move_files/{sample.ref}_scp_betas_commands.sh", root = config["root"], sample = sample_info.itertuples()),
        expand("{root}/scratch/move_files/{sample.ref}_scp_methyldackel_commands.sh", root = config["root"], sample = sample_info.itertuples()),
        expand("{root}/scratch/move_files/{sample.ref}_move_all.sh", root = config["root"], sample = sample_info.itertuples()),
        expand("{root}/scratch/move_files/trees/{sample.ref}_dirs_to_tree.sh", root = config["root"], sample = sample_info.itertuples()),
        expand("{root}/scratch/move_files/file_info/{sample.ref}_to_list.sh", root = config["root"], sample = sample_info.itertuples())

# Create a file with the directories to tree
rule dirs_to_tree:
    input:
        dirs_to_tree
    output:
        tree = expand("{root}/scratch/move_files/trees/{{ref}}_files_to_list.txt", root = config["root"]),
        script = expand("{root}/scratch/move_files/trees/{{ref}}_dirs_to_tree.sh", root = config["root"])
    params:
        ref = lambda wildcards, input: wildcards.ref,
        files=lambda wildcards, input: ' '.join(list(set(input))),
        tree_dir = expand("{root}/scratch/move_files/trees/", root = config["root"])
    shell:
        """
        echo "# Ref: {params.ref} script to build tree" > {output.script}
        echo "mkdir -p {params.tree_dir}" >> {output.script}
        echo "echo 'Dataset {params.ref}: farm structure of directories and files moved to local via scp\n' > {output.tree}" >> {output.script}
        for file in {params.files}; do
            echo "tree -h -n --noreport" $file"/{params.ref}--* >> {output.tree}" >> {output.script}
            echo "\n >> {output.tree}" >> {output.script}
        done
        sleep 5
        bash {output.script}
        """

rule files_to_ls:
    input:
        pre_merge_files,
        post_merge_files,
        pre_merge_dirs,
        post_merge_dirs,
        bam_files,
        wgbstools_betas,
        methyldackel
    output:
        expand("{root}/scratch/move_files/file_info/{{ref}}_to_list.sh", root = config["root"])
    params:
        ref = lambda wildcards, input: wildcards.ref,
        files=lambda wildcards, input: ' '.join(list(set(input))),
        list_out = expand("{root}/scratch/move_files/file_info/{{ref}}_dirs_to_tree.txt", root = config["root"]),
        list_dir = expand("{root}/scratch/move_files/file_info/", root = config["root"])
    shell:
        """
        mkdir -p {params.list_dir}
        echo "# Ref: {params.ref} script to list files" > {output}
        echo "mkdir -p {params.list_dir}" >> {output}
        echo "echo 'Dataset {params.ref}: farm directories and files moved to local via scp\n' > {params.list_out}" >> {output}
        echo "ls -lh {params.files} >> {params.list_out}" >> {output}
        sleep 5
        bash {output}
        """


rule pre_merge_files_command:
    input:
        pre_merge_files
    output:
        expand("{root}/scratch/move_files/{{ref}}_scp_pre_merge_files.sh", root = config["root"])
    params:
        files=lambda wildcards, input: ','.join(list(set(input))),
        dir="{ref}/pre-merge",
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        echo "cd {params.local}" > {output}
        echo "mkdir -p {params.dir}" >> {output}
        echo "cd {params.dir}" >> {output}
        echo "scp 'farm:{params.files}' ." >> {output}
        sleep 5
        """

rule pre_merge_dirs_command:
    input:
        pre_merge_dirs
    output:
        expand("{root}/scratch/move_files/{{ref}}_scp_pre_merge_dirs.sh", root = config["root"])
    params:
        files=lambda wildcards, input: ','.join(list(set(input))),
        dir="{ref}/pre-merge",
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        echo "cd {params.local}" > {output}
        echo "mkdir -p {params.dir}" >> {output}
        echo "cd {params.dir}" >> {output}
        echo "scp -r 'farm:{params.files}' ." >> {output}
        sleep 5
        """

rule post_merge_files_command:
    input:
        post_merge_files
    output:
        expand("{root}/scratch/move_files/{{ref}}_scp_post_merge_files.sh", root = config["root"])
    params:
        files=lambda wildcards, input: ','.join(list(set(input))),
        dir="{ref}/post-merge",
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        echo "cd {params.local}" > {output}
        echo "mkdir -p {params.dir}" >> {output}
        echo "cd {params.dir}" >> {output}
        echo "scp 'farm:{params.files}' ." >> {output}
        sleep 5
        """

rule post_merge_dirs_command:
    input:
        post_merge_dirs
    output:
        expand("{root}/scratch/move_files/{{ref}}_scp_post_merge_dirs.sh", root = config["root"])
    params:
        files=lambda wildcards, input: ','.join(list(set(input))),
        dir="{ref}/post-merge",
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        echo "cd {params.local}" > {output}
        echo "mkdir -p {params.dir}" >> {output}
        echo "cd {params.dir}" >> {output}
        echo "scp -r 'farm:{params.files}' ." >> {output}
        sleep 5
        """

rule move_bams_command:
    input:
        bam_files
    output:
        expand("{root}/scratch/move_files/{{ref}}_scp_bams_commands.sh", root = config["root"])
    params:
        files=lambda wildcards, input: ','.join(list(set(input))),
        dir="06_merged_deduped_bwa",
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        echo "cd {params.local}" > {output}
        echo "mkdir -p {params.dir}" >> {output}
        echo "cd {params.dir}" >> {output}
        echo "scp 'farm:{params.files}' ." >> {output}
        sleep 5
        """

rule move_betas_command:
    input:
        wgbstools_betas
    output:
        expand("{root}/scratch/move_files/{{ref}}_scp_betas_commands.sh", root = config["root"])
    params:
        files=lambda wildcards, input: ','.join(list(set(input))),
        dir="07_wgbstools_betas",
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        echo "cd {params.local}" > {output}
        echo "mkdir -p {params.dir}" >> {output}
        echo "cd {params.dir}" >> {output}
        echo "scp 'farm:{params.files}' ." >> {output}
        sleep 5
        """

rule move_methyldackel_command:
    input:
        methyldackel
    output:
        expand("{root}/scratch/move_files/{{ref}}_scp_methyldackel_commands.sh", root = config["root"])
    params:
        files=lambda wildcards, input: ','.join(list(set(input))),
        dir="07_methyldackel_bwa",
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        echo "cd {params.local}" > {output}
        echo "mkdir -p {params.dir}" >> {output}
        echo "cd {params.dir}" >> {output}
        echo "scp 'farm:{params.files}' ." >> {output}
        sleep 5
        """

rule join_scripts:
    input:
        pre_files = expand("{root}/scratch/move_files/{{ref}}_scp_pre_merge_files.sh", root = config["root"]), 
        post_files = expand("{root}/scratch/move_files/{{ref}}_scp_post_merge_files.sh", root = config["root"]),
        pre_dirs = expand("{root}/scratch/move_files/{{ref}}_scp_pre_merge_dirs.sh", root = config["root"]), 
        post_dirs = expand("{root}/scratch/move_files/{{ref}}_scp_post_merge_dirs.sh", root = config["root"]),
        bams = expand("{root}/scratch/move_files/{{ref}}_scp_bams_commands.sh", root = config["root"]),
        betas = expand("{root}/scratch/move_files/{{ref}}_scp_betas_commands.sh", root = config["root"]),
        methyld = expand("{root}/scratch/move_files/{{ref}}_scp_methyldackel_commands.sh", root = config["root"])
    output:
        expand("{root}/scratch/move_files/{{ref}}_move_all.sh", root = config["root"])
    params:
        local="/Volumes/Research_Backup/DMReproducible_Data/"
    shell:
        """
        cat {input.pre_files} > {output}
        cat {input.pre_dirs} >> {output}
        cat {input.post_files} >> {output}
        cat {input.post_dirs} >> {output}
        cat {input.betas} >> {output}
        cat {input.methyld} >> {output}
        cat {input.bams} >> {output}
        sleep 8
        """

onsuccess:
    print("Workflow finished, no error")
onerror:
    print("An error occurred")