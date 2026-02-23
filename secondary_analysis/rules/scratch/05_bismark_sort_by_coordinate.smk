
# ### bismark_sort_by_coordinate rule ###
# rule bismark_sort_by_coordinate_post_merge:
#     input:
#         bam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"])

#     output:
#         bam = expand("{root}/{data_dir}/05_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_bismark.deduplicated_sorted_merged.bam", root = config["root"], data_dir=config["data_dir"]),
#         bai = expand("{root}/{data_dir}/05_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_bismark.deduplicated_sorted_merged.bam.bai", root = config["root"], data_dir=config["data_dir"])

#     log:
#         "logs/secondary_rules/05_bismark_sort_by_coordinate/05_bismark_sort_by_coordinate-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
#     conda:
#         "../../environment_files/samtools.yaml"
    
#     params:
#         temp_dir = expand("{root}/{data_dir}/temp/samtools/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir=config["data_dir"]),
#         sorted_bam = expand("{root}/{data_dir}/05_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_bismark.deduplicated_sorted_merged.bam", root = config["root"], data_dir=config["data_dir"])

#     threads: 3

#     shell:
#         """
#         mkdir -p {params.temp_dir}
#         echo "Sorting bam file..." > {log}
#         samtools sort -@ {threads} -o {output.bam} -T {params.temp_dir} {input.bam} >> {log} 2>> {log}
#         echo "Sorting complete. Now indexing..." >> {log}
#         samtools index {output.bam} >> {log} 2>> {log}
#         echo "Indexing complete. Now removing intermediate bam files and tem directory.." >> {log}
#         rm -rf {params.temp_dir} >> {log} 2>> {log}
#         echo "Intermediate files: {output.bam}, {output.bam}.bai, and temporary directory: {params.temp_dir} have been removed. Done." >> {log}
#         """