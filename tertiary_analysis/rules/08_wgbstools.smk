# ================================================================
# segment_2_markers_target_specified_bwa rules
# ================================================================

rule segment_target_specified_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = _tVb_samples[_tVb_samples["ref"] == wildcards.ref].itertuples()
        )
    output:
        blocks = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/blocks.{{ref}}.{target}.v.{background}.bed",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    log:
        "logs/tertiary_rules/{ref}--segment_target_specified_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(_tVb_samples["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    params:
        out_dir = lambda wildcards: f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{wildcards.ref}-{config['tert_args']['wgbstools_target_group']}_v_{config['tert_args']['wgbstools_background_group']}",
        segment_args = config["tert_args"]["wgbstools_segment_betas"]
    shell:
        """
        echo "creating output directory if it does not exist" > {log}
        mkdir -p {params.out_dir}
        echo "done creating output directory" >> {log}

        echo "segmenting betas: {input.betas}" >> {log}
        wgbstools segment --betas {input.betas} -o {output.blocks} {params.segment_args} >> {log} 2>> {log}
        echo "done with segmentation step" >> {log}
        """


rule homog_target_specified_bwa:
    input:
        pats = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.pat.gz",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = _tVb_samples[_tVb_samples["ref"] == wildcards.ref].itertuples()
        ),
        blocks = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/blocks.{{ref}}.{target}.v.{background}.bed",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    output:
        homog_out = directory(expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/homog.{{ref}}.{target}.v.{background}",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        ))
    log:
        "logs/tertiary_rules/{ref}--homog_target_specified_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(_tVb_samples["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "counting homogeneously methylated fragments with wgbstools homog" > {log}
        wgbstools homog {input.pats} -b {input.blocks} -o {output.homog_out} --thresholds 0.25,0.75 >> {log} 2>> {log}
        echo "homog step done" >> {log}
        """


rule index_target_specified_bwa:
    input:
        blocks = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/blocks.{{ref}}.{target}.v.{background}.bed",
            root = config["root"],
            data_dir = config["data_dir"],
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
        )
    log:
        "logs/tertiary_rules/{ref}--index_target_specified_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(_tVb_samples["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "indexing blocks.bed output" > {log}
        wgbstools index {input.blocks} >> {log} 2>> {log}
        echo "done indexing blocks file" >> {log}
        """


rule beta_to_table_target_specified_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = _tVb_samples[_tVb_samples["ref"] == wildcards.ref].itertuples()
        ),
        blocks_gz = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/blocks.{{ref}}.{target}.v.{background}.bed.gz",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    output:
        table = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/avg_meth.{{ref}}.{target}.v.{background}.tsv",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    log:
        "logs/tertiary_rules/{ref}--beta_to_table_target_specified_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(_tVb_samples["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "extracting average methylation per block with beta_to_table" > {log}
        wgbstools beta_to_table {input.blocks_gz} --betas {input.betas} --output {output.table} >> {log} 2>> {log}
        echo "done running beta_to_table" >> {log}
        """


rule beta_to_table_grouped_target_specified_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = _tVb_samples[_tVb_samples["ref"] == wildcards.ref].itertuples()
        ),
        blocks_gz = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/blocks.{{ref}}.{target}.v.{background}.bed.gz",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        ),
        group_file = lambda wildcards: expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{ref}--groups.{target}.v.{background}.csv",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    output:
        table_grouped = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/avg_meth.grouped.{{ref}}.{target}.v.{background}.tsv",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    log:
        "logs/tertiary_rules/{ref}--beta_to_table_grouped_target_specified_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(_tVb_samples["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "extracting average methylation per group per block with beta_to_table" > {log}
        wgbstools beta_to_table {input.blocks_gz} --betas {input.betas} -g {input.group_file} --output {output.table_grouped} >> {log} 2>> {log}
        echo "done running beta_to_table merged by group" >> {log}
        """


rule find_markers_target_specified_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = _tVb_samples[_tVb_samples["ref"] == wildcards.ref].itertuples()
        ),
        blocks_gz = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/blocks.{{ref}}.{target}.v.{background}.bed.gz",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        ),
        group_file = lambda wildcards: expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{ref}--groups.{target}.v.{background}.csv",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    output:
        markers = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-{target}_v_{background}/markers.{{ref}}.{target}.v.{background}.log",
            root = config["root"],
            data_dir = config["data_dir"],
            target = config["tert_args"]["wgbstools_target_group"],
            background = config["tert_args"]["wgbstools_background_group"]
        )
    log:
        "logs/tertiary_rules/{ref}--find_markers_target_specified_bwa.log"
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
        find_marker_args = config["tert_args"]["wgbstools_find_markers"]
    shell:
        """
        echo "finding markers" > {log}
        wgbstools find_markers --out_dir {params.out_dir} --blocks_path {input.blocks_gz} --betas {input.betas} --targets {params.target_group} --background {params.background_group} --groups_file {input.group_file} {params.find_marker_args} > {output.markers} 2>> {output.markers}
        echo "done finding markers" >> {log}
        """


# ================================================================
# segment_2_markers_one_v_all_bwa rules
# ================================================================

rule segment_one_v_all_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = sample_info[sample_info["ref"] == wildcards.ref].itertuples()
        )
    output:
        blocks = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed",
            root = config["root"],
            data_dir = config["data_dir"]
        )
    log:
        "logs/tertiary_rules/{ref}--segment_one_v_all_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    params:
        out_dir = lambda wildcards: f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{wildcards.ref}-all-groups",
        segment_args = config["tert_args"]["wgbstools_segment_betas"]
    shell:
        """
        echo "creating output directory if it does not exist" > {log}
        mkdir -p {params.out_dir}
        echo "done creating output directory" >> {log}

        echo "segmenting betas: {input.betas}" >> {log}
        wgbstools segment --betas {input.betas} -o {output.blocks} {params.segment_args} >> {log} 2>> {log}
        echo "done with segmentation step" >> {log}
        """


rule homog_one_v_all_bwa:
    input:
        pats = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.pat.gz",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = sample_info[sample_info["ref"] == wildcards.ref].itertuples()
        ),
        blocks = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed",
            root = config["root"],
            data_dir = config["data_dir"]
        )
    output:
        homog_out = directory(expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/homog.{{ref}}",
            root = config["root"],
            data_dir = config["data_dir"]
        ))
    log:
        "logs/tertiary_rules/{ref}--homog_one_v_all_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "counting homogeneously methylated fragments with wgbstools homog" > {log}
        wgbstools homog {input.pats} -b {input.blocks} -o {output.homog_out} --thresholds 0.25,0.75 >> {log} 2>> {log}
        echo "homog step done" >> {log}
        """


rule index_one_v_all_bwa:
    input:
        blocks = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed",
            root = config["root"],
            data_dir = config["data_dir"]
        )
    output:
        index = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed.{suf}",
            root = config["root"],
            data_dir = config["data_dir"],
            suf = ["gz", "gz.tbi"]
        )
    log:
        "logs/tertiary_rules/{ref}--index_one_v_all_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "indexing blocks.bed output" > {log}
        wgbstools index {input.blocks} >> {log} 2>> {log}
        echo "done indexing blocks file" >> {log}
        """


rule beta_to_table_one_v_all_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = sample_info[sample_info["ref"] == wildcards.ref].itertuples()
        ),
        blocks_gz = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed.gz",
            root = config["root"],
            data_dir = config["data_dir"]
        )
    output:
        table = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/avg_meth.{{ref}}.tsv",
            root = config["root"],
            data_dir = config["data_dir"]
        )
    log:
        "logs/tertiary_rules/{ref}--beta_to_table_one_v_all_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "extracting average methylation per block with beta_to_table" > {log}
        wgbstools beta_to_table {input.blocks_gz} --betas {input.betas} --output {output.table} >> {log} 2>> {log}
        echo "done running beta_to_table" >> {log}
        """


rule beta_to_table_grouped_one_v_all_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = sample_info[sample_info["ref"] == wildcards.ref].itertuples()
        ),
        blocks_gz = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed.gz",
            root = config["root"],
            data_dir = config["data_dir"]
        ),
        group_file = lambda wildcards: expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{ref}--groups.csv",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref
        )
    output:
        table_grouped = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/avg_meth.grouped.{{ref}}.tsv",
            root = config["root"],
            data_dir = config["data_dir"]
        )
    log:
        "logs/tertiary_rules/{ref}--beta_to_table_grouped_one_v_all_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        echo "extracting average methylation per group per block with beta_to_table" > {log}
        wgbstools beta_to_table {input.blocks_gz} --betas {input.betas} -g {input.group_file} --output {output.table_grouped} >> {log} 2>> {log}
        echo "done running beta_to_table merged by group" >> {log}
        """


rule find_markers_one_v_all_bwa:
    input:
        betas = lambda wildcards: expand(
            "{root}/{data_dir}/07_wgbstools_betas_bwa/{ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.beta",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref,
            sample = sample_info[sample_info["ref"] == wildcards.ref].itertuples()
        ),
        blocks_gz = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/blocks.{{ref}}.bed.gz",
            root = config["root"],
            data_dir = config["data_dir"]
        ),
        group_file = lambda wildcards: expand(
            "{root}/{data_dir}/08_wgbstools_bwa/groups/{ref}--groups.csv",
            root = config["root"],
            data_dir = config["data_dir"],
            ref = wildcards.ref
        )
    output:
        markers = expand(
            "{root}/{data_dir}/08_wgbstools_bwa/{{ref}}-all-groups/markers.{{ref}}.log",
            root = config["root"],
            data_dir = config["data_dir"]
        )
    log:
        "logs/tertiary_rules/{ref}--find_markers_one_v_all_bwa.log"
    conda:
        "../../environment_files/wgbstools.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].unique().tolist())
    threads: 4
    resources:
        mem_mb=8000
    params:
        out_dir = lambda wildcards: f"{config['root']}/{config['data_dir']}/08_wgbstools_bwa/{wildcards.ref}-all-groups",
        find_marker_args = config["tert_args"]["wgbstools_find_markers"]
    shell:
        """
        echo "finding markers" > {log}
        wgbstools find_markers --out_dir {params.out_dir} --blocks_path {input.blocks_gz} --betas {input.betas} --groups_file {input.group_file} {params.find_marker_args} > {output.markers} 2>> {output.markers}
        echo "done finding markers" >> {log}
        """