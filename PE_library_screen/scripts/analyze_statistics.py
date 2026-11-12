#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Comprehensive statistics analysis for self-targeting library data.
Combines coverage analysis, summary statistics from logs, and cutadapt metrics.
"""

import pandas as pd
import gzip
import re
from Bio import SeqIO
import argparse
import os


def count_fastq_records(fastq_file):
    """Count the number of records in a FASTQ file."""
    try:
        with gzip.open(fastq_file, 'rt') as handle:
            return sum(1 for _ in SeqIO.parse(handle, 'fastq'))
    except Exception as e:
        print(f"Error counting records in {fastq_file}: {e}")
        raise


def check_fastq_stats(spacer_file, target_file, library_size):
    """
    Get basic statistics from FASTQ files.
    """
    try:
        # Count reads in spacer and target files
        spacer_count = count_fastq_records(spacer_file)
        target_count = count_fastq_records(target_file)

        return {
            "spacer_read_count": spacer_count,
            "target_read_count": target_count
        }

    except Exception as e:
        print(f"Error checking FASTQ stats: {e}")
        raise


def extract_summary_stats_from_log(log_file):
    """Extract summary statistics from an analysis log file."""
    stats = {
        'unique_matched_variants': 0,
        'total_reads': 0,
        'successfully_assigned_reads': 0,
        'recombined_reads_extension': 0,
        'percentage_recombined_extension': 0,
        'recombined_reads_scaffold': 0,
        'percentage_recombined_scaffold': 0
    }

    try:
        with open(log_file, 'r') as f:
            log_content = f.read()

            # Extract statistics using regex
            unique_match = re.search(r'Unique matched variants: (\d+)',
                                     log_content)
            if unique_match:
                stats['unique_matched_variants'] = int(unique_match.group(1))

            total_match = re.search(r'Total reads: (\d+)', log_content)
            if total_match:
                stats['total_reads'] = int(total_match.group(1))

            assigned_match = re.search(r'Successfully assigned reads: (\d+)',
                                       log_content)
            if assigned_match:
                stats['successfully_assigned_reads'] = int(
                    assigned_match.group(1))

            recombined_match_extension = re.search(r'Recombined reads at extension: (\d+)',
                                                   log_content)
            if recombined_match_extension:
                stats['recombined_reads_extension'] = int(recombined_match_extension.group(1))

            percentage_match_extension = re.search(r'Percentage extension recombined: ([\d\.]+)%',
                                                   log_content)
            if percentage_match_extension:
                stats['percentage_recombined_extension'] = float(
                    percentage_match_extension.group(1))

            recombined_match_scaffold = re.search(r'Recombined reads at scaffold: (\d+)',
                                                  log_content)
            if recombined_match_scaffold:
                stats['recombined_reads_scaffold'] = int(recombined_match_scaffold.group(1))

            percentage_match_scaffold = re.search(r'Percentage scaffold recombined: ([\d\.]+)%',
                                                  log_content)
            if percentage_match_scaffold:
                stats['percentage_recombined_scaffold'] = float(
                    percentage_match_scaffold.group(1))
        return stats

    except Exception as e:
        print(f"Error extracting statistics from log file {log_file}: {e}")
        raise


def extract_cutadapt_stats(log_file):  # noqa: C901
    """Extract cutadapt statistics from a log file."""
    stats = {
        'command': '',
        'total_reads': 0,
        'reads_with_adapters': 0,
        'adapter_percentage': 0,
        'reads_too_short': 0,
        'reads_too_long': 0,
        'reads_written': 0,
        'written_percentage': 0,
        'total_bp_processed': 0,
        'total_bp_written': 0,
        'bp_written_percentage': 0,
        'adapter_sequence': '',
        'processing_time': ''
    }

    try:
        with open(log_file, 'r') as f:
            content = f.read()

            # Extract command
            cmd_match = re.search(r'Command line parameters: (.*)', content)
            if cmd_match:
                stats['command'] = cmd_match.group(1)

            # Extract adapter sequence
            seq_match = re.search(r'Sequence: ([^;]+);', content)
            if seq_match:
                stats['adapter_sequence'] = seq_match.group(1)

            # Extract summary stats
            total_match = re.search(
                r'Total reads processed:\s+([0-9,]+)', content)
            if total_match:
                stats['total_reads'] = int(
                    total_match.group(1).replace(',', ''))

            adapter_match = re.search(
                r'Reads with adapters:\s+([0-9,]+) \(([0-9.]+)%\)', content)
            if adapter_match:
                stats['reads_with_adapters'] = int(
                    adapter_match.group(1).replace(',', ''))
                stats['adapter_percentage'] = float(adapter_match.group(2))

            short_match = re.search(
                r'Reads that were too short:\s+([0-9,]+) \(([0-9.]+)%\)',
                content)
            if short_match:
                stats['reads_too_short'] = int(
                    short_match.group(1).replace(',', ''))

            long_match = re.search(
                r'Reads that were too long:\s+([0-9,]+) \(([0-9.]+)%\)',
                content)
            if long_match:
                stats['reads_too_long'] = int(
                    long_match.group(1).replace(',', ''))

            written_match = re.search(
                (r'Reads written \(passing filters\):\s+([0-9,]+) '
                 r'\(([0-9.]+)%\)'),
                content)
            if written_match:
                stats['reads_written'] = int(
                    written_match.group(1).replace(',', ''))
                stats['written_percentage'] = float(written_match.group(2))

            bp_proc_match = re.search(
                r'Total basepairs processed:\s+([0-9,]+)', content)
            if bp_proc_match:
                stats['total_bp_processed'] = int(
                    bp_proc_match.group(1).replace(',', ''))

            bp_written_match = re.search(
                r'Total written \(filtered\):\s+([0-9,]+) bp \(([0-9.]+)%\)',
                content)
            if bp_written_match:
                stats['total_bp_written'] = int(
                    bp_written_match.group(1).replace(',', ''))
                stats['bp_written_percentage'] = float(
                    bp_written_match.group(2))
        return stats

    except Exception as e:
        print(f"Error extracting cutadapt stats from {log_file}: {e}")
        raise


def get_cutadapt_logs_for_sample(sample, log_dir):
    """Get all cutadapt log files for a sample."""
    log_types = {
        'target': f"{sample}_target_trim.log",
        'spacer': f"{sample}_spacer.log",
        'spacer_5trim': f"{sample}_spacer_5trim.log",
        'rtt': f"{sample}_rtt_trim.log"
    }

    logs = {}

    for log_type, filename in log_types.items():
        log_path = os.path.join(log_dir, filename)
        if os.path.exists(log_path):
            logs[log_type] = log_path

    return logs


def collect_cutadapt_stats(sample, log_dir):
    """Collect statistics from all cutadapt logs for a sample."""
    logs = get_cutadapt_logs_for_sample(sample, log_dir)
    stats = {}

    for log_type, log_path in logs.items():
        log_stats = extract_cutadapt_stats(log_path)
        stats[log_type] = log_stats

    return stats


def process_sample(sample, fastq_dir, log_dir,
                   date_prefix, library_size, min_coverage, library_name):
    """Process a single sample and return its statistics."""
    spacer_file = os.path.join(fastq_dir, f"{sample}_Spacer.fastq.gz")
    target_file = os.path.join(fastq_dir, f"{sample}_target_5trim.fastq.gz")
    log_file = os.path.join(
        log_dir, f"{date_prefix}_{library_name}_python_analysis_log_{sample}.txt")

    # Check that fastq files exist
    if not os.path.exists(spacer_file) or not os.path.exists(target_file):
        print(f"Warning: FASTQ files for sample {sample} not found")
        fastq_stats = {
            "spacer_read_count": 0,
            "target_read_count": 0
        }
    else:
        # Get FASTQ stats
        fastq_stats = check_fastq_stats(spacer_file, target_file, library_size)

    # Extract summary stats from log file
    if os.path.exists(log_file):
        summary_stats = extract_summary_stats_from_log(log_file)
    else:
        print(f"Warning: Log file {log_file} not found")
        summary_stats = {
            'unique_matched_variants': 0,
            'total_reads': 0,
            'successfully_assigned_reads': 0,
            'recombined_reads_extension': 0,
            'percentage_recombined_extension': 0,
            'recombined_reads_scaffold': 0,
            'percentage_recombined_scaffold': 0
        }

    # Collect cutadapt statistics
    cutadapt_stats = collect_cutadapt_stats(sample, log_dir)

    # Calculate additional metrics
    assigned_percentage = 0
    if summary_stats['total_reads'] > 0:
        assigned_percentage = (summary_stats['successfully_assigned_reads'] /
                               summary_stats['total_reads'] * 100)

    useful_coverage = (summary_stats['successfully_assigned_reads'] /
                       library_size) if library_size > 0 else 0

    # Determine coverage status based on useful coverage
    if useful_coverage >= min_coverage:
        coverage_status = "PASSED"
    else:
        coverage_status = "WARNING"

    # Combine all statistics
    stats = {
        "sample": sample,
        # Library info
        "library_size": library_size,
        "library_name": library_name,
        "min_coverage_required": min_coverage,

        # FASTQ file stats
        "spacer_read_count": fastq_stats["spacer_read_count"],
        "target_read_count": fastq_stats["target_read_count"],

        # Variant assignment stats
        "unique_matched_variants": summary_stats["unique_matched_variants"],
        "total_reads_analyzed": summary_stats["total_reads"],
        "successfully_assigned_reads": summary_stats[
            "successfully_assigned_reads"],
        "recombined_reads_extension": summary_stats["recombined_reads_extension"],
        "recombined_reads_scaffold": summary_stats["recombined_reads_scaffold"],
        "assigned_percentage": round(assigned_percentage, 2),
        "recombined_percentage_extension": summary_stats["percentage_recombined_extension"],
        "recombined_percentage_scaffold": summary_stats["percentage_recombined_scaffold"],
        "useful_coverage": round(useful_coverage, 2),
        "coverage_status": coverage_status
    }

    # Add cutadapt stats for each step
    for log_type, log_stats in cutadapt_stats.items():
        stats[f"{log_type}_adapter"] = log_stats["adapter_sequence"]
        stats[f"{log_type}_reads_total"] = log_stats["total_reads"]
        stats[f"{log_type}_reads_written"] = log_stats["reads_written"]
        stats[f"{log_type}_written_percentage"] = log_stats[
            "written_percentage"]

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Create comprehensive statistics report")
    parser.add_argument(
        '--fastq-dir', required=True, help="Directory containing FASTQ files")
    parser.add_argument(
        '--log-dir', required=True, help="Directory containing log files")
    parser.add_argument(
        '--output-dir', required=True, help="Directory for output files")
    parser.add_argument(
        '--output-csv', required=True, help="Output CSV file path")
    parser.add_argument(
        '--sample-list', required=True,
        help="Comma-separated list of sample names")
    parser.add_argument(
        '--date-prefix', required=True, help="Date prefix used in log files")
    parser.add_argument(
        '--library-size', type=int, required=True,
        help="Number of library members")
    parser.add_argument(
        '--min-coverage', type=int, required=True,
        help="Minimum desired coverage")
    parser.add_argument(
        '--library-name', required=True, help="Name of the library")
    args = parser.parse_args()

    # Parse sample list
    samples = [s.strip() for s in args.sample_list.split(',')]

    # Process each sample
    all_data = []
    for sample in samples:
        sample_data = process_sample(
            sample,
            args.fastq_dir,
            args.log_dir,
            args.date_prefix,
            args.library_size,
            args.min_coverage,
            args.library_name
        )
        all_data.append(sample_data)

    # Create DataFrame
    df = pd.DataFrame(all_data)

    # Sort by sample name
    df = df.sort_values('sample')

    # Generate separate cutadapt stats file
    cutadapt_cols = [col for col in df.columns if any(x in col for x in [
        'adapter', 'reads_total', 'reads_written', 'written_percentage'])]
    if cutadapt_cols:
        cutadapt_df = df[['sample'] + cutadapt_cols]
        cutadapt_output = os.path.join(
            args.output_dir,
            f"{args.date_prefix}_{args.library_name}_cutadapt_stats.csv"
        )
        cutadapt_df.to_csv(cutadapt_output, index=False)
        print(f"Cutadapt statistics report created: {cutadapt_output}")

    # Keep only the main stats in the main report
    main_cols = [col for col in df.columns if col not in cutadapt_cols]
    main_df = df[main_cols]

    # Save to CSV
    main_df.to_csv(args.output_csv, index=False)
    print(f"Comprehensive statistics report created: {args.output_csv}")
    print(f"Processed {len(all_data)} samples")


if __name__ == "__main__":
    main()
