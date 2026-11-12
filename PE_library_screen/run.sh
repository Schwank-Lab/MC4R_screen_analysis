#!/bin/bash
# Script to run the self-targeting library analysis pipeline
# Can be run locally or on a cluster

# Parse command-line arguments
mode="local"  # Default to local mode
cores=1
dry_run=false

# Help function
function show_help {
    echo "Usage: $0 [OPTIONS]"
    echo "Run the PE self-targeting analysis pipeline."
    echo ""
    echo "Options:"
    echo "  --cluster             Run on cluster using SLURM"
    echo "  --cores INT           Number of cores to use for local execution (default: 1)"
    echo "  --dry-run             Perform a dry run (don't execute commands)"
    echo "  --help                Show this help message and exit"
    echo ""
    echo "Example:"
    echo "  $0 --cores 4          # Run locally with 4 cores"
    echo "  $0 --cluster          # Run on SLURM cluster"
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --cluster)
            mode="cluster"
            shift
            ;;
        --cores)
            cores=$2
            shift 2
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        --help)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
done

# Construct Snakemake command
cmd="snakemake"

# Add --profile for cluster mode
if [ "$mode" = "cluster" ]; then
    echo "Running in cluster mode..."
    cmd="$cmd --profile profiles/slurm"
else
    echo "Running in local mode with $cores cores..."
    cmd="$cmd --cores $cores"
fi

# Add --use-conda to use conda environments
cmd="$cmd --use-conda"

# Add --dry-run if specified
if [ "$dry_run" = true ]; then
    cmd="$cmd --dry-run"
fi

# Run Snakemake
echo "Executing: $cmd"
eval $cmd

# Check for errors
if [ $? -ne 0 ]; then
    echo "Pipeline execution failed."
    exit 1
else
    echo "Pipeline execution completed successfully."
fi