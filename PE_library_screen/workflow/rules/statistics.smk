# Comprehensive statistics analysis rule

rule analyze_statistics:
    input:
        spacer_files = expand(config["output_folder"] + "{sample}_Spacer.fastq.gz", sample=SAMPLES),
        target_files = expand(config["output_folder"] + "{sample}_target_5trim.fastq.gz", sample=SAMPLES),
        analysis_logs = expand(config["log_folder"] + "{date_prefix}_{library_name}_python_analysis_log_{sample}.txt", 
                              sample=SAMPLES, date_prefix=DATE_PREFIX, library_name=config["library_name"]),
        cutadapt_logs = expand([
            config["log_folder"] + "{sample}_target_trim.log",
            config["log_folder"] + "{sample}_spacer.log",
            config["log_folder"] + "{sample}_spacer_5trim.log",
            config["log_folder"] + "{sample}_rtt_trim.log"
        ], sample=SAMPLES)
    output:
        main_report = config["analysis_folder"] + "{date_prefix}_{library_name}_statistics.csv",
        cutadapt_report = config["analysis_folder"] + "{date_prefix}_{library_name}_cutadapt_stats.csv"
    conda:
        "../envs/self_target_PE_analysis.yaml"
    params:
        library_size = config["library_size"],
        min_coverage = config["min_coverage"],
        sample_list = ",".join(SAMPLES),
        library_name = config["library_name"]
    shell:
        """
        python scripts/analyze_statistics.py \
        --fastq-dir {config[output_folder]} \
        --log-dir {config[log_folder]} \
        --output-dir {config[analysis_folder]} \
        --output-csv {output.main_report} \
        --sample-list "{params.sample_list}" \
        --date-prefix {wildcards.date_prefix} \
        --library-size {params.library_size} \
        --min-coverage {params.min_coverage} \
        --library-name {params.library_name}
        """