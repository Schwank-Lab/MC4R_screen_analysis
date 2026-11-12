# Self-Targeting Library Analysis Pipeline

A Snakemake pipeline for processing and analyzing PRIDICT2.0-designed prime editing self-targeting library data from paired-end sequencing.

## Overview

This pipeline automates the process of:
1. Trimming FASTQ.gz files to extract UMI-barcode, target, spacer, and RTT sequences
    - UMI-barcode currently not used for downstream analysis (can be used to get "unique" library integrations into cell lines); library barcode is extracted from "target" sequence.
2. Validating template files to ensure they contain the required columns
3. Checking coverage to ensure sufficient sequencing depth
4. Analyzing editing efficiency and patterns

## Getting Started

### Prerequisites

- Anaconda or Miniconda
- Snakemake
- FASTQ files from paired-end sequencing of self-targeting libraries

### Setup

1. Clone this repository
2. Place your FASTQ files in the `fastqs/` directory
3. Configure your analysis by editing `config/self_targeting_config.csv`
4. Make sure your library design file is available at the location specified in the config

### Configuration

Edit `config/self_targeting_config.csv` to set parameters for your analysis. The configuration file is organized into sections:

#### Essential Parameters

| Parameter | Description |
|-----------|-------------|
| library_design_file | Path to your library design CSV |
| library_size | Number of members in your library |
| library_name | Name of the current library |
| date_prefix | Date prefix for output files (YYYYMMDD) |
| target_adapter, spacer_adapter, etc. | Adapter sequences for trimming |
| target_start_to_nick_distance | Distance in bases from the start of the target sequence to the nick site |
| required_columns | Comma-separated list of required columns in the library design file |

#### Path Parameters
| Parameter | Description |
|-----------|-------------|
| fastq_folder | Directory with input FASTQ.gz files |
| output_folder | Directory for intermediate files |
| analysis_folder | Directory for analysis results |
| log_folder | Directory for log files |

#### Performance and Sequence Parameters
Additional parameters control the performance and detailed sequence analysis settings. See the config file for the complete list of options.

The config file uses a clear, Excel-friendly format that can be easily edited in a spreadsheet application.

### Running the Pipeline

```bash
# Run locally with 4 cores
./run.sh --cores 4

# Run on SLURM cluster
./run.sh --cluster

# Perform a dry run to preview commands
./run.sh --dry-run
```

## Output Files

The pipeline generates several output files:

- **Trimmed FASTQ files**: `Output/{sample}_*.fastq.gz`
- **Analysis results**: `Analysis/{date}_{sample}_analysisdf_focused.csv`
- **Template validation**: `validation/{date}_template_validation.txt`
- **Coverage reports**: `coverage/{date}_{sample}_coverage.txt`

## Adding a New Library

To analyze a new library:

1. Create a library design CSV file with the required columns
2. Update `config/self_targeting_config.csv` with:
   - Path to the new library design file
   - Updated library size
   - Library-specific adapter sequences (if different)
3. Run the pipeline as usual

## Required Columns in Library Design File

The library design file must contain these columns:
- uniquename
- Spacer-Sequence
- barcode
- RTrevcomp
- PBSrevcomp
- wide_initial_target
- wide_mutated_target
- Mutation_Type
- Correction_Type
- Correction_Length
- RTlength
- Editing_Position

## Troubleshooting

- **Missing columns in template file**: Check the validation report and ensure your library design file contains all required columns
- **Low coverage warning**: Consider increasing sequencing depth or reducing library complexity
- **Incorrect sequence extraction**: If the analysis is not correctly identifying edited sequences, check the `target_start_to_nick_distance` parameter. This value represents the distance in bases from the start of the trimmed target sequence to the nick site and is crucial for proper sequence alignment.