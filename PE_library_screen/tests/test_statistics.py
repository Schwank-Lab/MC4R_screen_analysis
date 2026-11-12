import os
import sys
import pytest
import pandas as pd

# Import the module to test
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.analyze_statistics import count_fastq_records, \
    check_fastq_stats, extract_cutadapt_stats, extract_summary_stats_from_log, \
    get_cutadapt_logs_for_sample, collect_cutadapt_stats, \
    process_sample, main    # noqa: E402


@pytest.fixture
def processed_test_files(test_data_path):
    """Paths to processed test FASTQ files."""
    test_data_path += "/fastq/trimmed"
    return {
        "spacer": os.path.join(test_data_path, "test_MC4R-lib-control1_Spacer.fastq.gz"),
        "target": os.path.join(test_data_path, "test_MC4R-lib-control1_target_5trim.fastq.gz")
    }


def test_count_fastq_records(test_fastq_files):
    """Test counting records in a valid FASTQ file."""
    # Use the test FASTQ file
    fastq_file = test_fastq_files["R1"]

    expected_count = 25000  # Known count for this test file

    # Count records
    count = count_fastq_records(fastq_file)

    # Verify exact count
    assert count == expected_count, f"Should count exactly {expected_count} records"


def test_count_fastq_records_nonexistent_file(capsys):
    """Test that counting records in a non-existent file prints error and raises exception."""
    with pytest.raises(FileNotFoundError):
        count_fastq_records("nonexistent_file.fastq.gz")

    # Check that the error was printed
    captured = capsys.readouterr()
    assert "Error counting records in nonexistent_file.fastq.gz" in captured.out


def test_check_fastq_stats(processed_test_files):
    """Test checking stats from processed FASTQ files."""
    spacer_file = processed_test_files["spacer"]
    target_file = processed_test_files["target"]
    library_size = 500

    stats = check_fastq_stats(spacer_file, target_file, library_size)

    assert "spacer_read_count" in stats
    assert "target_read_count" in stats
    assert stats["spacer_read_count"] == 23998  # Known count for this test file
    assert stats["target_read_count"] == 24959  # Known count for this test file


def test_check_fastq_stats_missing_file(capsys, processed_test_files):
    """Test handling when one file is missing."""
    spacer_file = processed_test_files["spacer"]
    target_file = "nonexistent_file.fastq.gz"
    library_size = 500

    with pytest.raises(FileNotFoundError):
        check_fastq_stats(spacer_file, target_file, library_size)

    # Check that the error was printed
    captured = capsys.readouterr()
    assert "Error checking FASTQ stats: " in captured.out


def test_extract_cutadapt_stats_from_real_log(test_data_path):
    """Test extract_cutadapt_stats with a real cutadapt log file."""
    test_data_path += "/logs"
    # Path to the actual log file
    log_file = os.path.join(test_data_path, "test_MC4R-lib-control1_spacer.log")

    # Extract stats
    stats = extract_cutadapt_stats(log_file)

    # Verify key values from the log are extracted correctly
    assert stats['command'] == (
        "-j 0 -a GTTTCAGA -m 19 -M 19 --discard-untrimmed "
        "-o test_MC4R-lib-control1_Spacer.fastq.gz "
        "test_MC4R-lib-control1_spacer_5trim.fastq.gz"
    )
    assert stats['total_reads'] == 24727
    assert stats['reads_with_adapters'] == 24321
    assert stats['adapter_percentage'] == 98.4
    assert stats['reads_too_short'] == 282
    assert stats['reads_too_long'] == 447
    assert stats['reads_written'] == 23998
    assert stats['written_percentage'] == 97.1
    assert stats['total_bp_processed'] == 3082568
    assert stats['total_bp_written'] == 455962
    assert stats['bp_written_percentage'] == 14.8
    assert stats['adapter_sequence'] == "GTTTCAGA"


def test_extract_cutadapt_stats_from_inexistent_log(capsys):
    """Test extract_cutadapt_stats with a inexsitent cutadapt log file."""

    log_file = "non_existent_file.txt"

    # Test with a non-existent file
    with pytest.raises(FileNotFoundError):
        extract_cutadapt_stats(log_file)

    # Check that the error was printed
    captured = capsys.readouterr()
    assert f"Error extracting cutadapt stats from {log_file}" in captured.out


def test_extract_summary_stats_from_real_log(test_data_path):
    """Test extract_summary_stats_from_log with a real analysis log file."""
    test_data_path += "/logs"
    # Path to the actual log file
    log_file = os.path.join(test_data_path, "20250319_fulllib_python_analysis_log_test_MC4R-lib-control1.txt")

    # Extract stats
    stats = extract_summary_stats_from_log(log_file)

    # Verify key values from the log are extracted correctly
    assert stats['unique_matched_variants'] == 488
    assert stats['total_reads'] == 322265
    assert stats['successfully_assigned_reads'] == 241124
    assert stats['recombined_reads_extension'] == 7738
    assert stats['percentage_recombined_extension'] == 2.4
    assert stats['recombined_reads_scaffold'] == 73403
    assert stats['percentage_recombined_scaffold'] == 23.34


def test_extract_summary_stats_nonexistent_file(capsys):
    """Test extract_summary_stats_from_log with a non-existent file."""
    # Test with a non-existent file, should raise an exception
    with pytest.raises(Exception):
        extract_summary_stats_from_log("non_existent_file.txt")

    # Check that the error was printed
    captured = capsys.readouterr()
    assert "Error extracting statistics from log file" in captured.out


def test_get_cutadapt_logs_for_sample(test_data_path):
    """Test getting cutadapt log files for a sample."""
    test_data_path += "/logs"
    # Use the sample name matching the test files
    sample = "test_MC4R-lib-control1"

    # Get logs using the test_data directory
    logs = get_cutadapt_logs_for_sample(sample, test_data_path)

    # Verify the correct log files were found
    assert len(logs) == 4
    assert "target" in logs
    assert "spacer" in logs
    assert "spacer_5trim" in logs
    assert "rtt" in logs

    # Verify the correct file paths
    assert os.path.basename(logs["target"]) == "test_MC4R-lib-control1_target_trim.log"
    assert os.path.basename(logs["spacer"]) == "test_MC4R-lib-control1_spacer.log"
    assert os.path.basename(logs["spacer_5trim"]) == "test_MC4R-lib-control1_spacer_5trim.log"
    assert os.path.basename(logs["rtt"]) == "test_MC4R-lib-control1_rtt_trim.log"

    # Test with a nonexistent sample
    nonexistent_logs = get_cutadapt_logs_for_sample("nonexistent-sample", test_data_path)
    assert len(nonexistent_logs) == 0


def test_collect_cutadapt_stats(test_data_path):
    """Test collecting cutadapt stats for a sample."""
    test_data_path += "/logs"

    # Use the sample name matching the test files
    sample = "test_MC4R-lib-control1"

    # Collect stats using the test_data directory
    stats = collect_cutadapt_stats(sample, test_data_path)

    # Verify we collected stats for all log types
    assert len(stats) == 4
    assert "target" in stats
    assert "spacer" in stats
    assert "spacer_5trim" in stats
    assert "rtt" in stats

    # Verify some log's stats
    assert stats["spacer"]["total_reads"] == 24727
    assert stats["target"]["total_reads"] == 25000
    assert stats["rtt"]["total_reads"] == 25000
    assert stats["spacer_5trim"]["total_reads"] == 25000
    assert stats["spacer"]["adapter_sequence"] == "GTTTCAGA"
    assert stats["target"]["reads_with_adapters"] == 24959
    assert stats["spacer"]["reads_written"] == 23998

    # Check a few key details from each log
    assert stats["target"]["adapter_sequence"] != stats["spacer"]["adapter_sequence"]
    assert stats["spacer_5trim"]["adapter_sequence"] != stats["spacer"]["adapter_sequence"]
    assert stats["rtt"]["adapter_sequence"] != stats["spacer"]["adapter_sequence"]

    # Test with a nonexistent sample
    nonexistent_stats = collect_cutadapt_stats("nonexistent-sample", test_data_path)
    assert len(nonexistent_stats) == 0


def test_process_sample(test_data_path, sample_config):
    """Test processing a single sample with real test files and compare logs with expected values."""
    # Setup test directories and parameters
    fastq_dir = os.path.join(test_data_path, "fastq/trimmed")
    log_dir = os.path.join(test_data_path, "logs")
    print(fastq_dir, log_dir)

    sample = "test_MC4R-lib-control1"
    date_prefix = "20250319"  # From the actual log filename
    library_size = sample_config["library_size"]  # 500 from config
    library_name = sample_config["library_name"]  # From config
    min_coverage = sample_config["min_coverage"]  # From config

    # Process the sample
    stats = process_sample(
        sample=sample,
        fastq_dir=fastq_dir,
        log_dir=log_dir,
        date_prefix=date_prefix,
        library_size=library_size,
        min_coverage=min_coverage,
        library_name=library_name
    )

    # Verify the basic stats are present and correct
    assert stats["sample"] == sample
    assert stats["library_size"] == library_size
    assert stats["library_name"] == library_name
    assert stats["min_coverage_required"] == min_coverage

    # Verify FASTQ stats
    assert stats["spacer_read_count"] == 23998
    assert stats["target_read_count"] == 24959

    # Verify analysis stats from the log file
    assert stats['unique_matched_variants'] == 488
    assert stats['total_reads_analyzed'] == 322265
    assert stats['successfully_assigned_reads'] == 241124
    assert stats['recombined_reads_extension'] == 7738
    assert stats['recombined_percentage_extension'] == 2.4
    assert stats['recombined_reads_scaffold'] == 73403
    assert stats['recombined_percentage_scaffold'] == 23.34

    # Verify coverage status
    assert stats["useful_coverage"] > 0
    assert stats["coverage_status"] in ["PASSED", "WARNING"]

    # Verify cutadapt stats are present
    assert "target_adapter" in stats
    assert "spacer_adapter" in stats
    assert "spacer_5trim_adapter" in stats
    assert "rtt_adapter" in stats

    # Verify some key adapter stats
    assert stats["spacer_reads_total"] == 24727
    assert stats["spacer_reads_written"] == 23998
    assert stats["target_reads_total"] == 25000
    assert stats["target_reads_written"] == 24959

    # Verify adapter sequences
    assert stats["spacer_adapter"] == "GTTTCAGA"
    assert "GTCATAGGAGTCAC".upper() in stats["target_adapter"]
    assert "gaaacaccG".upper() in stats["spacer_5trim_adapter"]
    assert "AGTCGGTGC".upper() in stats["rtt_adapter"]

    # Verify written percentages
    assert stats["spacer_written_percentage"] == 97.1
    assert stats["target_written_percentage"] > 90


def test_process_sample_passed_coverage(test_data_path, sample_config):
    """Test that coverage_status = PASSED branch is reached."""
    fastq_dir = os.path.join(test_data_path, "fastq/trimmed")
    log_dir = os.path.join(test_data_path, "logs")

    sample = "test_MC4R-lib-control1"
    date_prefix = "20250319"
    library_size = sample_config["library_size"]
    library_name = sample_config["library_name"]
    # Set min_coverage lower than the actual coverage to trigger PASSED
    min_coverage = 400  # Lower than the ~477 from the test data

    stats = process_sample(
        sample=sample,
        fastq_dir=fastq_dir,
        log_dir=log_dir,
        date_prefix=date_prefix,
        library_size=library_size,
        min_coverage=min_coverage,
        library_name=library_name
    )

    # Now specifically check for PASSED status
    assert stats["coverage_status"] == "PASSED"


def test_main_function(monkeypatch, test_data_path, tmpdir):
    """Test the main function with real test data and command-line arguments."""
    # Setup test directories
    fastq_dir = os.path.join(test_data_path, "fastq/trimmed")
    log_dir = os.path.join(test_data_path, "logs")
    output_dir = str(tmpdir)

    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define sample and output parameters
    sample = "test_MC4R-lib-control1"  # From existing tests
    date_prefix = "20250319"           # From existing tests
    library_size = 500                 # From config
    min_coverage = 500                 # From config
    library_name = "fulllib"              # From config
    output_csv = os.path.join(output_dir, f"{date_prefix}_{library_name}_statistics.csv")

    # Mock command-line arguments
    args = [
        'analyze_statistics.py',
        '--fastq-dir', fastq_dir,
        '--log-dir', log_dir,
        '--output-dir', output_dir,
        '--output-csv', output_csv,
        '--sample-list', sample,
        '--date-prefix', date_prefix,
        '--library-size', str(library_size),
        '--min-coverage', str(min_coverage),
        '--library-name', library_name
    ]
    monkeypatch.setattr('sys.argv', args)

    # Call the main function
    main()

    # Verify output files were created
    assert os.path.exists(output_csv), f"Main output file not created: {output_csv}"
    cutadapt_output = os.path.join(output_dir, f"{date_prefix}_{library_name}_cutadapt_stats.csv")
    assert os.path.exists(cutadapt_output), f"Cutadapt stats file not created: {cutadapt_output}"

    # Read and verify contents of main stats file
    main_df = pd.read_csv(output_csv)
    assert len(main_df) == 1, "Main stats file should contain one row"
    assert main_df.iloc[0]['sample'] == sample, f"Expected sample name {sample}, got {main_df.iloc[0]['sample']}"
    assert main_df.iloc[0][
        'library_size'] == library_size, f"Expected library size {library_size}, got {main_df.iloc[0]['library_size']}"
    assert main_df.iloc[0]['unique_matched_variants'] > 0, "Expected positive number of unique matched variants"

    # Read and verify contents of cutadapt stats file
    cutadapt_df = pd.read_csv(cutadapt_output)
    assert len(cutadapt_df) == 1, "Cutadapt stats file should contain one row"
    assert cutadapt_df.iloc[0][
        'sample'] == sample, f"Expected sample name {sample}, got {cutadapt_df.iloc[0]['sample']}"
    assert any('adapter' in col for col in cutadapt_df.columns), "Expected adapter columns in cutadapt stats"
    assert any('reads_total' in col for col in cutadapt_df.columns), "Expected reads_total columns in cutadapt stats"
    assert any(
        'reads_written' in col for col in cutadapt_df.columns), "Expected reads_written columns in cutadapt stats"