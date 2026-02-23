# ### sambamba rules for sorting, deduplexing, and indexing ###

# # sambamba is a high performance, robust, and fast tool for working with SAM and BAM files. It is a faster alternative to samtools and picard tools.
# # sambamba sort: Sorts the input bam file by coordinates. The output is a sorted bam file.
# # sambamba markdup: Marks duplicates in the input bam file. It also removes duplicates if the --remove-duplicates option is used.
# # sambamba index: Indexes the input bam file. The output is a .bai file.
# # sambamba does not produce a log file by default. Therefore, the log file is created using the shell command by capturing standard output and standard error.

# # sambamba_sort_index_markdups rule: sorts, deduplexes, and indexes the input bam file using sambamba.
# # rule input: bwameth aligned bam file.
# # rule output: sorted/deduplexed bam file, .bai file, and a log file.
# rule sambamba_sort_index_markdups:
#     input: 
#         bam = expand("{root}/{data_dir}/03_aligned_bwameth/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", root = config["root"], data_dir=config["data_dir"])

#     output:
#         bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"]),
#         bai = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.bai", root = config["root"], data_dir=config["data_dir"]),
#         report = expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}-{{accession}}.log", root = config["root"], rep_dir=config["reports_dir"])

#     log:
#         expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}-{{accession}}.log", root = config["root"], rep_dir=config["reports_dir"])

#     conda:
#         "../../environment_files/sambamba.yaml"
    
#     # shadow:
#     #     "shallow"

#     params:
#         temp_dir = expand("{root}/{data_dir}/temp/sambamba/{{ref}}--{{accession}}", root = config["root"], data_dir=config["data_dir"]),
#         sorted_bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted.bam", root = config["root"], data_dir=config["data_dir"]),

#     threads: 3

#     shell:
#         """
#         mkdir -p {params.temp_dir}
#         echo "Sorting bam file..." > {log}
#         sambamba sort -t {threads} -o {params.sorted_bam} --tmpdir {params.temp_dir} {input.bam} >> {log} 2>> {log}
#         echo "Sorting complete. Now marking duplicates..." >> {log} 
#         sambamba markdup -t {threads} --remove-duplicates --tmpdir {params.temp_dir} {params.sorted_bam} {output.bam} >> {log} 2>> {log}
#         echo "Marking duplicates complete. Now indexing..." >> {log}
#         sambamba index -t {threads} {output.bam} {output.bai} >> {log} 2>> {log}
#         echo "Indexing complete. Now removing intermediate bam files and tem directory.." >> {log}
#         rm -f {params.sorted_bam} >> {log} 2>> {log}
#         rm -f {params.sorted_bam}.bai >> {log} 2>> {log}
#         rm -rf {params.temp_dir} >> {log} 2>> {log}
#         echo "Intermediate files: {params.sorted_bam}, {params.sorted_bam}.bai, and temporary directory: {params.temp_dir} have been removed. Done." >> {log}
#         """
