# rule to segment, index, and create table tsv summarizing beta values by segment
rule segment_betas:
    input:
        betas = pre_processing_workflow(expand("{data_dir}/06_wgbstools_betas/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.{suf}", data_dir=config["data_dir"], suf=["beta"], sample=sample_info.itertuples()))

    output:
        blocks = expand("{data_dir}/06_wgbstools/blocks.bed", data_dir=config["data_dir"]),
        index = expand("{data_dir}/06_wgbstools/blocks.bed.{suf}", data_dir=config["data_dir"], suf=["gz", "gz.tbi"]),
        table = expand("{data_dir}/06_wgbstools/segments.tsv", data_dir=config["data_dir"])

    log:
        "../pre-processing/logs/rule-logs/06_wgbstools_segment.log"

    conda:
        "../environment_files/wgbstools.yaml"

    # shadow:
    #     "shallow"

    shell: 
        """
        wgbstools segment --betas {input.betas} --min_cpg 3 --max_bp 2000 -o {output.blocks}
        wgbstools index {output.blocks}
        wgbstools beta_to_table {output.blocks} --betas {input.betas} | column -t >> {output.table}
        """   


