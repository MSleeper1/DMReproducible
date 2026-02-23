rule make_group_file:
    input:
        betas = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["beta"], sample = sample_info.itertuples())
    
    output:
        csv_file = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/{{ref}}--groups.csv", root = config["root"], data_dir=config["data_dir"])
    
    log:
        "logs/tertiary_rules/{ref}--groups_file.log"

    conda:
        "../../environment_files/wgbstools.yaml"
    
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())

    params:
        betas_file_names = expand("{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", sample = sample_info.itertuples())

    run:
        sample_info[['custom_id', 'group']].to_csv(output.csv_file, index=False)
        



    

