import os
import sys
import subprocess
import gzip
from Bio import SeqIO

# Import script to test cutadapt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.analyze_statistics import extract_cutadapt_stats  # noqa: E402


def test_cutadapt_install():
    """Test that cutadapt is installed and accessible."""
    result = subprocess.run(
        ["cutadapt", "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert result.returncode == 0, "Cutadapt command failed"
    # Just check that we got some non-empty output
    assert result.stdout.strip(), "No output from cutadapt version command"


def test_trim_r1_to_target(test_output_dirs, sample_config):
    """Test trimming R1 files to target using cutadapt."""
    subprocess.run(["cutadapt", "--version"], check=False)

    # Set up paths
    r1_file = os.path.join(test_output_dirs["fastq"], "MC4R-lib-control1_R1.fastq.gz")
    output_file = os.path.join(
        test_output_dirs["output"], "MC4R-lib-control1_target_5trim.fastq.gz"
    )
    log_file = os.path.join(
        test_output_dirs["log"], "MC4R-lib-control1_target_trim.log"
    )

    # Define the cutadapt command based on the Snakemake rule
    cmd = [
        "cutadapt",
        "-j",
        str(sample_config["cutadapt_threads"]),
        "-g",
        sample_config["target_adapter"],
        "--discard-untrimmed",
        "-o",
        output_file,
        r1_file,
    ]

    # Run cutadapt
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    # Check that command executed successfully
    assert result.returncode == 0, "Cutadapt command failed"

    # Check that output file was created
    assert os.path.exists(output_file), "Output file not created"

    # Check that log file was created and has content
    assert os.path.exists(log_file), "Log file not created"
    assert os.path.getsize(log_file) > 0, "Log file is empty"

    # Extract and check cutadapt stats
    stats = extract_cutadapt_stats(log_file)
    assert stats["total_reads"] > 0, "No reads processed"
    assert stats["reads_written"] > 0, "No reads written"
    assert (
        stats["adapter_sequence"] == sample_config["target_adapter"]
    ), "Incorrect adapter sequence in log"

    # Check that output file contains FASTQ data
    with gzip.open(output_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
        assert len(records) > 0, "No records in output file"


def test_trim_r2_to_spacer(test_output_dirs, sample_config):
    """Test trimming R2 files to spacer using cutadapt."""
    subprocess.run(["cutadapt", "--version"], check=False)

    # Set up paths
    r2_file = os.path.join(test_output_dirs["fastq"], "MC4R-lib-control1_R2.fastq.gz")
    output_file = os.path.join(
        test_output_dirs["output"], "MC4R-lib-control1_spacer_5trim.fastq.gz"
    )
    log_file = os.path.join(
        test_output_dirs["log"], "MC4R-lib-control1_spacer_5trim.log"
    )

    # Define the cutadapt command based on the Snakemake rule
    cmd = [
        "cutadapt",
        "-j",
        str(sample_config["cutadapt_threads"]),
        "-g",
        sample_config["spacer_adapter"],
        "--discard-untrimmed",
        "-o",
        output_file,
        r2_file,
    ]

    # Run cutadapt
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    # Check that command executed successfully
    assert result.returncode == 0, "Cutadapt command failed"

    # Check that output file was created and has content
    assert os.path.exists(output_file), "Output file not created"

    # Check that log file was created and has content
    assert os.path.exists(log_file), "Log file not created"
    assert os.path.getsize(log_file) > 0, "Log file is empty"

    # Check that output file contains FASTQ data
    with gzip.open(output_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
        assert len(records) > 0, "No records in output file"


def test_trim_spacer_to_final_length(test_output_dirs, sample_config):
    """Test trimming spacer to final length using cutadapt."""
    subprocess.run(["cutadapt", "--version"], check=False)

    # First, run the prerequisite trim_r2_to_spacer
    test_trim_r2_to_spacer(test_output_dirs, sample_config)

    # Set up paths for the second trimming
    input_file = os.path.join(
        test_output_dirs["output"], "MC4R-lib-control1_spacer_5trim.fastq.gz"
    )
    output_file = os.path.join(
        test_output_dirs["output"], "MC4R-lib-control1_Spacer.fastq.gz"
    )
    log_file = os.path.join(test_output_dirs["log"], "MC4R-lib-control1_spacer.log")

    # Define the cutadapt command based on the Snakemake rule
    cmd = [
        "cutadapt",
        "-j",
        str(sample_config["cutadapt_threads"]),
        "-a",
        sample_config["spacer_adapter_end"],
        "-m",
        str(sample_config["spacer_length"]),
        "-M",
        str(sample_config["spacer_length"]),
        "--discard-untrimmed",
        "-o",
        output_file,
        input_file,
    ]

    # Run cutadapt
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    # Check that command executed successfully
    assert result.returncode == 0, "Cutadapt command failed"

    # Check that output file was created and has content
    assert os.path.exists(output_file), "Output file not created"

    # Check that log file was created and has content
    assert os.path.exists(log_file), "Log file not created"
    assert os.path.getsize(log_file) > 0, "Log file is empty"

    # Check that output file contains FASTQ data with correct read length
    with gzip.open(output_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
        assert len(records) > 0, "No records in output file"

        # Check that reads have the expected length
        for record in records:
            assert (
                len(record.seq) == sample_config["spacer_length"]
            ), f"Read length {len(record.seq)} does not match expected length {sample_config['spacer_length']}"


def test_trim_r2_to_rtt(test_output_dirs, sample_config):
    """Test trimming R2 files to RTT using cutadapt."""
    subprocess.run(["cutadapt", "--version"], check=False)

    # Set up paths
    r2_file = os.path.join(test_output_dirs["fastq"], "MC4R-lib-control1_R2.fastq.gz")
    output_file = os.path.join(
        test_output_dirs["output"], "MC4R-lib-control1_RTTtrim.fastq.gz"
    )
    log_file = os.path.join(test_output_dirs["log"], "MC4R-lib-control1_rtt_trim.log")

    # Define the cutadapt command based on the Snakemake rule
    cmd = [
        "cutadapt",
        "-j",
        str(sample_config["cutadapt_threads"]),
        "-g",
        sample_config["rtt_adapter"],
        "--discard-untrimmed",
        "-o",
        output_file,
        r2_file,
    ]

    # Run cutadapt
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    # Check that command executed successfully
    assert result.returncode == 0, "Cutadapt command failed"

    # Check that output file was created and has content
    assert os.path.exists(output_file), "Output file not created"

    # Check that log file was created and has content
    assert os.path.exists(log_file), "Log file not created"
    assert os.path.getsize(log_file) > 0, "Log file is empty"

    # Check that output file contains FASTQ data
    with gzip.open(output_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
        assert len(records) > 0, "No records in output file"


def test_trim_with_nonmatching_adapter(test_output_dirs, sample_config):
    """Test trimming with an adapter that doesn't match any sequences."""
    subprocess.run(["cutadapt", "--version"], check=False)

    # Set up paths
    r1_file = os.path.join(test_output_dirs["fastq"], "MC4R-lib-control1_R1.fastq.gz")
    output_file = os.path.join(
        test_output_dirs["output"], "MC4R-lib-control1_nonmatching.fastq.gz"
    )
    log_file = os.path.join(
        test_output_dirs["log"], "MC4R-lib-control1_nonmatching.log"
    )

    # Define the cutadapt command with a nonmatching adapter
    nonmatching_adapter = "AAAAAAAAAAAAAAAAAAAA"  # Long poly-A sequence unlikely to match
    cmd = [
        "cutadapt",
        "-j", str(sample_config["cutadapt_threads"]),
        "-g", nonmatching_adapter,
        "-O", "20",  # Set overlap to 20 to require longer match
        "-e", "0",  # Set error rate to 0 to require exact match
        "--discard-untrimmed",
        "-o", output_file,
        r1_file,
    ]

    # Run cutadapt
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    # Command should execute successfully (return code 0) even if no sequences match
    assert result.returncode == 0, "Cutadapt command failed"

    # Check that log file was created and has content
    assert os.path.exists(log_file), "Log file not created"
    assert os.path.getsize(log_file) > 0, "Log file is empty"

    # Extract cutadapt stats
    stats = extract_cutadapt_stats(log_file)

    # Verify that no reads were written (or very few if any happen to match by chance)
    assert stats["reads_written"] == 0 or stats["written_percentage"] < 1.0, \
        "Unexpectedly high number of reads matched the nonmatching adapter"

    # Verify the output file exists but should be empty or nearly empty
    assert os.path.exists(output_file), "Output file not created"

    # Read the output file and verify it has no or very few records
    with gzip.open(output_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
        assert len(records) == 0 or len(records) < 5, \
            f"Expected 0 or very few records, but found {len(records)}"
