# rule to segment, index, and create table tsv summarizing beta values by segment
rule segment_betas_bwa:
    input:
        betas = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["beta"], sample = sample_info.itertuples())

    output:
        blocks = expand("{root}/{data_dir}/07_wgbstools/blocks.{{ref}}.bed", data_dir=config["data_dir"], root = config["root"]),
        index = expand("{root}/{data_dir}/07_wgbstools/blocks.{{ref}}.bed.{suf}", data_dir=config["data_dir"], suf=["gz", "gz.tbi"], root = config["root"]),
        table = expand("{root}/{data_dir}/07_wgbstools/avg_meth.{{ref}}.tsv", data_dir=config["data_dir"], root = config["root"])

    log:
        "logs/tertiary_rules/{ref}--segment.log"

    conda:
        "../../environment_files/wgbstools.yaml"
    
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())

    params:
        beta_dir = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
        beta_files = expand("{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.{suf}", suf=["beta"], sample = sample_info.itertuples()),
        block_file = "blocks.{ref}.bed",
        meth_table = "avg_meth.{ref}.tsv"

    shell: 
        """
        echo "segmenting betas: {input.betas}" > {log}
        wgbstools segment --betas {input.betas} --min_cpg 3 --max_bp 2000 -o {params.block_file}
        echo " indexing blocks output" >> {log}
        wgbstools index {params.block_file}
        echo "extracting average methylation per block with beta to table" >> {log}
        wgbstools beta_to_table {params.block_file} --betas {input.betas} | column -t > {params.meth_table}
        echo "moving block file and meth table to {params.beta_dir}" >> {log} 2>> {log}
        mv -f -v --target-directory={params.beta_dir} {params.block_file} >> {log} 2>> {log}
        mv -f -v --target-directory={params.beta_dir} {params.meth_table} >> {log} 2>> {log}
        echo "done" >> {log} 2>> {log}
        pwd >> {log}
        """   