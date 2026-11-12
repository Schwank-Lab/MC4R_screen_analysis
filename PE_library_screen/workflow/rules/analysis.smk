# Analysis rule
rule analyze_library:
    input:
        target = config["output_folder"] + "{sample}_target_5trim.fastq.gz",
        spacer = config["output_folder"] + "{sample}_Spacer.fastq.gz",
        rtt = config["output_folder"] + "{sample}_RTTtrim.fastq.gz",
        template = config["library_design_file"],
    output:
        analysis = config["analysis_folder"] + "{date_prefix}_{library_name}_{sample}_analysisdf_focused.csv"
    log:
        config["log_folder"] + "{date_prefix}_{library_name}_python_analysis_log_{sample}.txt"
    conda:
        "../envs/self_target_PE_analysis.yaml"
    shell:
        """
        python -u scripts/self_targeting_analysis.py \
        --shortname {wildcards.sample} \
        --path {config[output_folder]} \
        --template_file {input.template} \
        --output_dir {config[analysis_folder]} \
        --date_prefix {wildcards.date_prefix} \
        --library_name {config[library_name]} \
        --barcode_start {config[barcode_start]} \
        --barcode_length {config[barcode_length]} \
        --min_quality {config[min_quality]} \
        --basesbeforenick {config[basesbeforenick]} \
        --baseafterflap {config[baseafterflap]} \
        --target_start_to_nick_distance {config[target_start_to_nick_distance]} \
        --nick_pos_in_wide_initial_target {config[nick_pos_in_wide_initial_target]} \
        --rtt_lookup_length {config[rtt_lookup_length]} \
        --read_length {config[read_length]} \
        > {log} 2>&1
        """