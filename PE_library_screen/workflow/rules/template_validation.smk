# Template validation rule

# Validate template file
rule validate_template:
    input:
        template = config["library_design_file"]
    output:
        validation_result = config["analysis_folder"] + "{date_prefix}_{library_name}_template_validation.txt"
    conda:
        "../envs/self_target_PE_analysis.yaml"
    params:
        required_columns = ",".join(config["required_columns"]) if "required_columns" in config else "uniquename,Spacer-Sequence,barcode,RTrevcomp,PBSrevcomp"
    shell:
        """
        python scripts/validate_template.py \
        --template {input.template} \
        --required-columns "{params.required_columns}" \
        --output {output.validation_result}
        """