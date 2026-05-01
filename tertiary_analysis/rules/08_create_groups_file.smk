rule make_group_file:
    input:
        # listed as dependency to ensure beta files exist before group file is created
        betas = expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            sample = sample_info.itertuples()
        )
    output:
        csv_file = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{{ref}}--groups.csv",
            root = config["root"],
            data_dir = config["data_dir"],
        )

    log:
        "logs/tertiary_rules/{ref}--groups_file.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())

    resources:
        mem_mb=4000
    run:
        import os
        os.makedirs(os.path.dirname(output.csv_file[0]), exist_ok=True)
        sample_info[['custom_id', 'group']].to_csv(output.csv_file[0], index=False)

rule make_group_file_tVb:
    input:
        # listed as dependency to ensure beta files exist before group file is created
        betas = expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            sample = _tVb_samples.itertuples()
        )
    output:
        csv_file = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{{ref}}--groups.{target}.v.{background}.csv",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )

    log:
        "logs/tertiary_rules/{ref}--groups_file.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(_tVb_samples["ref"].unique().tolist())
    resources:
        mem_mb=4000
    run:
        import os
        os.makedirs(os.path.dirname(output.csv_file[0]), exist_ok=True)
        _tVb_samples[['custom_id', 'group']].to_csv(output.csv_file[0], index=False)

# rule make_group_file:
#     input:
#         betas = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["beta"], sample = sample_info.itertuples())
    
#     output:
#         csv_file = expand("{root}/{data_dir}/08_wgbstools_bwa/groups/{{ref}}--groups.csv", root = config["root"], data_dir=config["data_dir"])
    
#     log:
#         "logs/tertiary_rules/{ref}--groups_file.log"

#     conda:
#         "../../environment_files/wgbstools.yaml"
    
#     wildcard_constraints:
#         ref = "|".join(sample_info["ref"].tolist())

#     params:
#         out_dir = expand("{root}/{data_dir}/08_wgbstools_bwa/groups/", root = config["root"], data_dir=config["data_dir"])
    
#     run:
#         import os
#         os.makedirs(os.path.dirname(params.out_dir), exist_ok=True)
#         sample_info[['custom_id', 'group']].to_csv(output.csv_file, index=False)


    

