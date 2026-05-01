### 05 deduplication, sorting, and indexing rules ###
#--------------------------------------------------------------#
# Rules to deduplicate bismark mapped bam files after merging  #
#--------------------------------------------------------------#

# deduplicate_bismark rule: deduplicates the input bam file using deduplicate_bismark.
rule bismark_deduplicate_post_merge:
    input: 
        expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])
        
    output:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        report = expand("{root}/{rep_dir}/06_merged_bismark_deduplication/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.deduplication_report.txt", root = config["root"], rep_dir=config["reports_dir"])

    log:
        "logs/secondary_rules/06_bismark_deduplicate_post_merge/06_bismark_deduplicate-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/bismark.yaml"

    threads: 2

    resources:
        mem_mb=4000

    params:
        output_dir = expand("{root}/{data_dir}/06_merged_deduped_bis/", root = config["root"], data_dir=config["data_dir"]),
        base_name = "{ref}--{patient_id}-{group}-{srx_id}-{layout}",
        report_dir = expand("{root}/{rep_dir}/06_merged_bismark_deduplication/", root = config["root"], rep_dir=config["reports_dir"]),
        report = "{ref}--{patient_id}-{group}-{srx_id}-{layout}.deduplication_report.txt"

    shell:
        """
        echo "making output directory {params.output_dir}" > {log}
        mkdir -p {params.output_dir} 2>>{log}
        echo "running bismark deduplication on {input}"
        deduplicate_bismark --bam {input} --output_dir {params.output_dir} --outfile {params.base_name} >> {log} 2>> {log}
        echo "done with deduplication"
        echo "renaming deduplicated bam file to {output.bam}" >> {log}
        mv -f -v {params.output_dir}{params.base_name}.deduplicated.bam {output.bam} >> {log} 2>> {log}
        echo "making report directory: {params.report_dir}" >> {log}
        mkdir -p {params.report_dir} 2>>{log}
        echo "moving {params.output_dir}{params.base_name}.deduplication_report.txt to {params.report_dir}" >> {log}
        mv -f -v --target-directory={params.report_dir} {params.output_dir}{params.report} >> {log} 2>> {log}
        echo "done"
        """
### bismark_sort_by_coordinate post merge rule ###
rule bismark_sort_by_coordinate_post_merge:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])  
        
    output:
        bam = expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam.bai", root = config["root"], data_dir=config["data_dir"])

    log:
        "logs/secondary_rules/06_bismark_sort_by_coordinate/06_bismark_sort_by_coordinate-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/samtools.yaml"

    threads: 2

    resources:
        mem_mb=4000
        
    params:
        temp_dir = expand("{root}/{data_dir}/temp/samtools/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir=config["data_dir"]),
        sorted_bam = expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])

    threads: 3

    shell:
        """
        mkdir -p {params.temp_dir}
        echo "Sorting bam file..." > {log}
        samtools sort -@ {threads} -o {output.bam} -T {params.temp_dir} {input.bam} >> {log} 2>> {log}
        echo "Sorting complete. Now indexing..." >> {log}
        samtools index {output.bam} >> {log} 2>> {log}
        echo "Indexing complete. Now removing intermediate bam files and tem directory.." >> {log}
        rm -rf {params.temp_dir} >> {log} 2>> {log}
        echo "Intermediate files: {output.bam}, {output.bam}.bai, and temporary directory: {params.temp_dir} have been removed. Done." >> {log}
        """


#------------------------------------------------#
# Rules to deduplicate bwameth aligned bam files #
#------------------------------------------------#
### sambamba rules for sorting, deduplexing, and indexing ###

# sambamba is a high performance, robust, and fast tool for working with SAM and BAM files. It is a faster alternative to samtools and picard tools.
# sambamba sort: Sorts the input bam file by coordinates. The output is a sorted bam file.
# sambamba markdup: Marks duplicates in the input bam file. It also removes duplicates if the --remove-duplicates option is used.
# sambamba index: Indexes the input bam file. The output is a .bai file.
# sambamba does not produce a log file by default. Therefore, the log file is created using the shell command by capturing standard output and standard error.

# sambamba_sort_index_markdups rule: sorts, deduplexes, and indexes the input bam file using sambamba.
# rule input: bwameth aligned bam file.
# rule output: sorted/deduplexed bam file, .bai file, and a log file.
rule sambamba_sort_index_markdups_post_merge_bwa:
    input: 
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])
        
    output:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam.bai", root = config["root"], data_dir=config["data_dir"]),
        report = expand("{root}/{rep_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.log", root = config["root"], rep_dir=config["reports_dir"])

    log:
        "logs/secondary_rules/06_sambamba_sort_index_markdups_post_merge_bwa/06_sambamba_sort_index_markdups_post_merge_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/sambamba.yaml"

    params:
        out_dir = expand("{root}/{data_dir}/06_merged_deduped_bwa/", root = config["root"], data_dir=config["data_dir"]),
        temp_dir = expand("{root}/{data_dir}/temp/sambamba/{{ref}}--{{srx_id}}", root = config["root"], data_dir=config["data_dir"]),
        sorted_bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_sorted.bam", root = config["root"], data_dir=config["data_dir"]),

    threads: 2

    resources:
        mem_mb=4000

    shell:
        """
        mkdir -p {params.temp_dir}
        mkdir -p {params.out_dir}
        touch {log}
        echo "Sorting bam file..." > {log}
        sambamba sort -t {threads} -o {params.sorted_bam} --tmpdir {params.temp_dir} {input.bam} >> {log} 2>> {log}
        echo "Sorting complete. Now marking duplicates..." >> {log} 
        sambamba markdup -t {threads} --remove-duplicates --tmpdir {params.temp_dir} {params.sorted_bam} {output.bam} >> {log} 2>> {log}
        echo "Marking duplicates complete. Now indexing..." >> {log}
        sambamba index -t {threads} {output.bam} {output.bai} >> {log} 2>> {log}
        echo "Indexing complete. Now removing intermediate bam files and tem directory.." >> {log}
        rm -f {params.sorted_bam} >> {log} 2>> {log}
        rm -f {params.sorted_bam}.bai >> {log} 2>> {log}
        rm -rf {params.temp_dir} >> {log} 2>> {log}
        echo "Intermediate files: {params.sorted_bam}, {params.sorted_bam}.bai, and temporary directory: {params.temp_dir} have been removed. Done." >> {log}
        cp {log} {output.report}
        """
