# Cutadapt rules

# Trim R1 files to target - this is the key file where barcodes will be extracted from
rule trim_r1_to_target:
    input:
        r1 = lambda wildcards: SAMPLE_FILES[wildcards.sample]['r1']
    output:
        target = config["output_folder"] + "{sample}_target_5trim.fastq.gz"
    log:
        config["log_folder"] + "{sample}_target_trim.log"
    conda:
        "../envs/self_target_PE_analysis.yaml"
    shell:
        """
        cutadapt -j {config[cutadapt_threads]} -g {config[target_adapter]} \
        --discard-untrimmed -o {output.target} {input.r1} > {log} 2>&1
        """

# Trim R2 files to spacer sequence (first step)
rule trim_r2_to_spacer_5:
    input:
        r2 = lambda wildcards: SAMPLE_FILES[wildcards.sample]['r2']
    output:
        spacer_5trim = config["output_folder"] + "{sample}_spacer_5trim.fastq.gz"
    log:
        config["log_folder"] + "{sample}_spacer_5trim.log"
    conda:
        "../envs/self_target_PE_analysis.yaml"
    shell:
        """
        cutadapt -j {config[cutadapt_threads]} -g {config[spacer_adapter]} \
        --discard-untrimmed -o {output.spacer_5trim} {input.r2} > {log} 2>&1
        """

# Trim spacer_5trim to get spacer
rule trim_spacer_to_final_length:
    input:
        spacer_5trim = config["output_folder"] + "{sample}_spacer_5trim.fastq.gz"
    output:
        spacer = config["output_folder"] + "{sample}_Spacer.fastq.gz"
    log:
        config["log_folder"] + "{sample}_spacer.log"
    conda:
        "../envs/self_target_PE_analysis.yaml"
    shell:
        """
        cutadapt -j {config[cutadapt_threads]} -a {config[spacer_adapter_end]} \
        -m {config[spacer_length]} -M {config[spacer_length]} --discard-untrimmed -o {output.spacer} {input.spacer_5trim} > {log} 2>&1
        """

# Trim R2 files to RTT
rule trim_r2_to_rtt:
    input:
        r2 = lambda wildcards: SAMPLE_FILES[wildcards.sample]['r2']
    output:
        rtt = config["output_folder"] + "{sample}_RTTtrim.fastq.gz"
    log:
        config["log_folder"] + "{sample}_rtt_trim.log"
    conda:
        "../envs/self_target_PE_analysis.yaml"
    shell:
        """
        cutadapt -j {config[cutadapt_threads]} -g {config[rtt_adapter]} \
        --discard-untrimmed -o {output.rtt} {input.r2} > {log} 2>&1
        """