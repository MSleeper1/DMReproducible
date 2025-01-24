# Clean up intermediate directories produced by the secondary analysis workflow
# This rule will remove all intermediate directories produced by the secondary analysis workflow
# Note that `snakemake clean` for secondary workflow will remove all output files created by the workflow
# This rule can be used to only remove intermediate directories and files

#### import modules ####
import pandas as pd
import secondary_helper_functions as hf

#### assign config ####
configfile: "../config.yaml"

#### sample info ####
sample_info = hf.get_sample_info_df(config["root"] + "/" + config["samples_tsv"])
sample_info_se = sample_info[sample_info['layout'] == 'se']
sample_info_pe = sample_info[sample_info['layout'] == 'pe']

#### clean up intermediate directories ####
secondary_data_directories = [
    directory(expand("{root}/{data_dir}/02_trimmed_trim_galore/", root = config["root"], data_dir=config["data_dir"])), # trim_galore data directory
    directory(expand("{root}/{data_dir}/03_aligned_bismark_bwt2/" root = config["root"], data_dir=config["data_dir"])), # mapped data directory
    directory(expand("{root}/{data_dir}/03_aligned_bwameth/" root = config["root"], data_dir=config["data_dir"])), # bwameth_mapping data directory
    directory(expand("{root}/{data_dir}/04_bismark_deduped/" root = config["root"], data_dir=config["data_dir"])), # bismark deduped data directory
    directory(expand("{root}/{data_dir}/04_deduped_sambamba/", root = config["root"], data_dir=config["data_dir"])), # sambamba deduped data directory
    directory(expand("{root}/{data_dir}/05_merged_sambamba_bwa/", root = config["root"], data_dir=config["data_dir"])),
    directory(expand("{root}/{data_dir}/05_merged_sambamba_bis/", root = config["root"], data_dir=config["data_dir"])),
    directory(expand("{root}/{data_dir}/06_merged_deduped_bis/", root = config["root"], data_dir=config["data_dir"])),
    directory(expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/", root = config["root"], data_dir=config["data_dir"])),
    directory(expand("{root}/{data_dir}/06_merged_deduped_bwa/", root = config["root"], data_dir=config["data_dir"]))
]

rule clean_secondary_int_dirs:
    input:
        secondary_data_directories
    output:
        touch("{root}/{data_dir}/clean_secondary_int_dirs.done", root = config["root"], data_dir=config["data_dir"])
    log:
        "logs/clean/clean_secondary_int_dirs.log"
    shell:
        """
        echo -e "Cleaning up intermediate directories produced by the secondary analysis workflow\n" > {log}
        echo "Processed on:" >> {log}
        echo $(date) >> {log}
        echo -e "\nDirectories to be removed:" >> {log}
        echo "{input}" >> {log}
        rm -rf {input}
        echo -e "\nFinished" >> {log}
        cp {log} {output}
        """

rule all:
    input:
        expand("{root}/{data_dir}/clean_secondary_int_dirs.done", root = config["root"], data_dir=config["data_dir"])