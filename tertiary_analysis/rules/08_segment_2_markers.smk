# rule to segment, index, and create table tsv summarizing beta values by segment
rule segment_2_markers_target_specified_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = _tVb_samples[_tVb_samples["ref"] == wildcards.ref].itertuples()
        ),
        pats = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.pat.gz",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = _tVb_samples[_tVb_samples["ref"] == wildcards.ref].itertuples()
        ),
        group_file = lambda wildcards: expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{{ref}}--groups.{target}.v.{background}.csv",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )

    output:
        index = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/blocks.{{ref}}.{target}.v.{background}.bed.{suf}",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"],
            suf = ["gz", "gz.tbi"]
        ),
        table = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/avg_meth.{{ref}}.{target}.v.{background}.tsv",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        ),
        table_grouped = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/avg_meth.grouped.{{ref}}.{target}.v.{background}.tsv",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        ),
        homog_out = directory(expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/homog.{{ref}}.{target}.v.{background}",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )),
        markers = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/markers.{{ref}}.{target}.v.{background}.log",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )

    log:
        "logs/tertiary_rules/{ref}--segment_2_markers_target_specified_bwa.log"

    conda:
        "../../environment_files/wgbstools.yaml"

    wildcard_constraints:
        ref = "|".join(_tVb_samples["ref"].unique().tolist())

    threads: 4
    
    resources:
        mem_mb=8000

    params:
        out_dir = lambda wildcards: f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{wildcards.ref}-{config['tert_args']['wgbstools_target_group']}_v_{config['tert_args']['wgbstools_background_group']}",
        target_group = config["tert_args"]["wgbstools_target_group"],
        background_group = config["tert_args"]["wgbstools_background_group"],
        segment_args = config["tert_args"]["wgbstools_segment_betas"],
        find_marker_args = config["tert_args"]["wgbstools_find_markers"],
        blocks = lambda wildcards: f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{wildcards.ref}-{config['tert_args']['wgbstools_target_group']}_v_{config['tert_args']['wgbstools_background_group']}/blocks.{wildcards.ref}.{config['tert_args']['wgbstools_target_group']}.v.{config['tert_args']['wgbstools_background_group']}.bed"

    shell:
        """
        echo "creating output directory if it does not exist" > {log}
        mkdir -p {params.out_dir}
        echo "done creating output directory" >> {log}

        echo "segmenting betas: {input.betas}" >> {log}
        wgbstools segment --betas {input.betas} -o {params.blocks} {params.segment_args} >> {log} 2>> {log}
        echo "Done with segmentation step" >> {log}

        echo "Counting homogenously methylated fragments with wgbstools homog" >> {log}
        wgbstools homog {input.pats} -b {params.blocks} -o {output.homog_out} --thresholds 0.25,0.75
        echo "homog step done" >> {log}

        echo "indexing blocks.bed output" >> {log}
        wgbstools index {params.blocks}
        echo "done indexing blocks file" >> {log}

        echo "extracting average methylation per block with beta_to_table" >> {log}
        wgbstools beta_to_table {params.blocks}.gz --betas {input.betas} --output {output.table}
        echo "done running beta_to_table" >> {log}

        echo "extracting average methylation per group per block with beta_to_table" >> {log}
        wgbstools beta_to_table {params.blocks}.gz --betas {input.betas} -g {input.group_file} --output {output.table_grouped}
        echo "done running beta_to_table merged by group" >> {log}

        echo "finding markers" >> {log}
        wgbstools find_markers --out_dir {params.out_dir} --blocks_path {params.blocks}.gz --betas {input.betas} --targets {params.target_group} --background {params.background_group} --groups_file {input.group_file} {params.find_marker_args} > {output.markers} 2>> {output.markers}
        echo "done finding markers" >> {log}
        """


        # echo "moving output files to {params.out_dir}" >> {log}
        # mv -f -v --target-directory={params.out_dir} {output.blocks}.gz >> {log} 2>> {log}
        # mv -f -v --target-directory={params.out_dir} {output.blocks}.gz.tbi >> {log} 2>> {log}
        # mv -f -v --target-directory={params.out_dir} {params.meth_table} >> {log} 2>> {log}
        # mv -f -v --target-directory={params.out_dir} {params.meth_table_grouped} >> {log} 2>> {log}
        # mv -f -v --target-directory={params.out_dir} {params.homog_out_dir} >> {log} 2>> {log}
        # echo "done" >> {log} 2>> {log}

# rule to segment, index, and create table tsv summarizing beta values by segment
rule segment_2_markers_one_v_all_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = sample_info[sample_info["ref"] == wildcards.ref].itertuples()
        ),
        pats = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.pat.gz",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = sample_info[sample_info["ref"] == wildcards.ref].itertuples()
        ),
        group_file = lambda wildcards: expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{ref}--groups.csv",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref
        )

    output:
        index = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed.{suf}",
            root = config["root"],
            data_dir = config["data_dir"],
            suf = ["gz", "gz.tbi"]
        ),
        table = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/avg_meth.{{ref}}.tsv",
            root = config["root"],
            data_dir = config["data_dir"]
        ),
        table_grouped = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/avg_meth.grouped.{{ref}}.tsv",
            root = config["root"],
            data_dir = config["data_dir"]
        ),
        homog_out = directory(expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/homog.{{ref}}",
            root = config["root"],
            data_dir = config["data_dir"]
        )),
        markers = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/markers.{{ref}}.log",
            root = config["root"],
            data_dir = config["data_dir"]
        )

    log:
        "logs/tertiary_rules/{ref}--segment_2_markers_one_v_all_bwa.log"

    conda:
        "../../environment_files/wgbstools.yaml"

    threads: 4
    
    resources:
        mem_mb=8000

    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())

    params:
        out_dir = lambda wildcards: f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{wildcards.ref}-all-groups",
        segment_args = config["tert_args"]["wgbstools_segment_betas"],
        find_marker_args = config["tert_args"]["wgbstools_find_markers"],
        blocks = lambda wildcards: f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{wildcards.ref}-all-groups/blocks.{wildcards.ref}.bed"

    shell:
        """
        echo "creating output directory if it does not exist" > {log}
        mkdir -p {params.out_dir}
        echo "done creating output directory" >> {log}

        echo "segmenting betas: {input.betas}" >> {log}
        wgbstools segment --betas {input.betas} -o {params.blocks} {params.segment_args} >> {log} 2>> {log}
        echo "Done with segmentation step" >> {log}

        echo "Counting homogenously methylated fragments with wgbstools homog" >> {log}
        wgbstools homog {input.pats} -b {params.blocks} -o {output.homog_out} --thresholds 0.25,0.75
        echo "homog step done" >> {log}

        echo "indexing blocks.bed output" >> {log}
        wgbstools index {params.blocks}
        echo "done indexing blocks file" >> {log}

        echo "extracting average methylation per block with beta_to_table" >> {log}
        wgbstools beta_to_table {params.blocks}.gz --betas {input.betas} --output {output.table}
        echo "done running beta_to_table" >> {log}

        echo "extracting average methylation per group per block with beta_to_table" >> {log}
        wgbstools beta_to_table {params.blocks}.gz --betas {input.betas} -g {input.group_file} --output {output.table_grouped}
        echo "done running beta_to_table merged by group" >> {log}

        echo "finding markers" >> {log}
        wgbstools find_markers --out_dir {params.out_dir} --blocks_path {params.blocks}.gz --betas {input.betas} --groups_file {input.group_file} {params.find_marker_args} > {output.markers} 2>> {output.markers}
        echo "done finding markers" >> {log}
        """
        # It is also possible to calculate the average methylation of each block for groups of beta files with the -g group_file.csv flag.


# # this rule could be run if more than 1 ref is present
# rule segment_2_markers_combined_refs_bwa:
#     input:
#         betas = expand(
#             "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             sample = sample_info.itertuples()
#         ),
#         pats = expand(
#             "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.pat.gz",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             sample = sample_info.itertuples()
#         ),
#         group_file = expand(
#             "{root}/{data_dir}/08_wgbstools_bwa/groups/{ref}--groups.csv",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             ref = sample_info["ref"].unique().tolist()
#         )

#     output:
#         blocks = expand(
#             "{root}/{data_dir}/08_wgbstools_bwa/{combined_ref}-all-groups/blocks.{combined_ref}.bed",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             combined_ref = _combined_ref
#         ),
#         index = expand(
#             "{root}/{data_dir}/08_wgbstools_bwa/{combined_ref}-all-groups/blocks.{combined_ref}.bed.{suf}",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             combined_ref = _combined_ref,
#             suf = ["gz", "gz.tbi"]
#         ),
#         table = expand(
#             "{root}/{data_dir}/08_wgbstools_bwa/{combined_ref}-all-groups/avg_meth.{combined_ref}.tsv",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             combined_ref = _combined_ref
#         ),
#         table_grouped = expand(
#             "{root}/{data_dir}/08_wgbstools_bwa/{combined_ref}-all-groups/avg_meth.grouped.{combined_ref}.tsv",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             combined_ref = _combined_ref
#         ),
#         homog_out = directory(expand(
#             "{root}/{data_dir}/08_wgbstools_bwa/{combined_ref}-all-groups/homog.{combined_ref}",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             combined_ref = _combined_ref
#         )),
#         markers = expand(
#             "{root}/{data_dir}/08_wgbstools_bwa/{combined_ref}-all-groups/markers.{combined_ref}.tsv",
#             root = config["root"],
#             data_dir = config["data_dir"],
#             combined_ref = _combined_ref
#         )

#     log:
#         f"logs/tertiary_rules/{_combined_ref}--segment_2_markers_combined_refs_bwa.log"

#     conda:
#         "../../environment_files/wgbstools.yaml"

#     params:
#         out_dir = f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{_combined_ref}-all-groups",
#         block_file = f"blocks.{_combined_ref}.bed",
#         meth_table = f"avg_meth.{_combined_ref}.tsv",
#         meth_table_grouped = f"avg_meth.grouped.{_combined_ref}.tsv",
#         homog_out_dir = f"homog.{_combined_ref}",
#         segment_args = config["tert_args"]["wgbstools_segment_betas"],
#         find_marker_args = config["tert_args"]["wgbstools_find_markers"]

#     shell:
#         """
#         echo "creating output directory if it does not exist" > {log}
#         mkdir -p {params.out_dir}
#         echo "done creating output directory" >> {log}

#         echo "segmenting betas: {input.betas}" >> {log}
#         wgbstools segment --betas {input.betas} -o {output.blocks} {params.segment_args}
#         echo "Done with segmentation step" >> {log}

#         echo "Counting homogenously methylated fragments with wgbstools homog" >> {log}
#         wgbstools homog {input.pats} -b {output.blocks} -o {params.homog_out_dir} --thresholds 0.25,0.75
#         echo "homog step done" >> {log}

#         echo "indexing blocks.bed output" >> {log}
#         wgbstools index {output.blocks}
#         echo "done indexing blocks file" >> {log}

#         echo "extracting average methylation per block with beta_to_table" >> {log}
#         wgbstools beta_to_table {output.blocks}.gz --betas {input.betas} --output {params.meth_table}
#         echo "done running beta_to_table" >> {log}

#         echo "extracting average methylation per group per block with beta_to_table" >> {log}
#         wgbstools beta_to_table {output.blocks}.gz --betas {input.betas} -g {input.group_file} --output {params.meth_table_grouped}
#         echo "done running beta_to_table merged by group" >> {log}

#         echo "finding markers" >> {log}
#         wgbstools find_markers --out_dir {params.out_dir} --blocks_path {output.blocks} --betas {input.betas} --groups_file {input.group_file} {params.find_marker_args}
#         echo "done finding markers" >> {log}

#         echo "moving output files to {params.out_dir}" >> {log}
#         mv -f -v --target-directory={params.out_dir} {output.blocks}.gz >> {log} 2>> {log}
#         mv -f -v --target-directory={params.out_dir} {output.blocks}.gz.tbi >> {log} 2>> {log}
#         mv -f -v --target-directory={params.out_dir} {params.meth_table} >> {log} 2>> {log}
#         mv -f -v --target-directory={params.out_dir} {params.meth_table_grouped} >> {log} 2>> {log}
#         mv -f -v --target-directory={params.out_dir} {params.homog_out_dir} >> {log} 2>> {log}
#         echo "done" >> {log} 2>> {log}
#         """