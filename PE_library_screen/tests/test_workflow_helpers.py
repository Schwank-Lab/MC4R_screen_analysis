import os
import sys
import pytest
import tempfile
from unittest.mock import patch

# Import the function to test
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.workflow_helpers import get_samples_from_sheet  # type: ignore # noqa: E402


def test_get_samples_from_sheet_underscore():
    """Test that sample sheet is read correctly."""
    # Create a temporary sample sheet
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("sample,r1_file,r2_file\n")
        f.write("sample1_rep1,sample1_R1.fastq.gz,sample1_R2.fastq.gz\n")
        f.write("sample2_rep2,sample2_R1.fastq.gz,sample2_R2.fastq.gz\n")
        temp_file = f.name

    try:
        # Mock os.path.exists to make it think files exist
        with patch('os.path.exists', return_value=True):
            fastq_folder = "/path/to/fastqs"
            with pytest.raises(ValueError, match="Sample names must not contain underscores"):
                get_samples_from_sheet(temp_file, fastq_folder)

    finally:
        # Clean up
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_get_samples_from_sheet_error_underscore():
    """Test that sample sheet is read correctly."""
    # Create a temporary sample sheet
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("sample,r1_file,r2_file\n")
        f.write("sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz\n")
        f.write("sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz\n")
        temp_file = f.name

    try:
        # Mock os.path.exists to make it think files exist
        with patch('os.path.exists', return_value=True):
            fastq_folder = "/path/to/fastqs"
            samples, sample_dict = get_samples_from_sheet(temp_file, fastq_folder)

            # Check samples list
            assert samples == ["sample1", "sample2"]

            # Check file mappings
            assert sample_dict["sample1"]["r1"] == os.path.join(fastq_folder, "sample1_R1.fastq.gz")
            assert sample_dict["sample1"]["r2"] == os.path.join(fastq_folder, "sample1_R2.fastq.gz")
            assert sample_dict["sample2"]["r1"] == os.path.join(fastq_folder, "sample2_R1.fastq.gz")
            assert sample_dict["sample2"]["r2"] == os.path.join(fastq_folder, "sample2_R2.fastq.gz")

    finally:
        # Clean up
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_get_samples_from_sheet_missing_file():
    """Test that function raises ValueError when sample sheet doesn't exist."""
    with pytest.raises(ValueError):
        get_samples_from_sheet("nonexistent_file.csv", "/path/to/fastqs")


def test_get_samples_from_sheet_missing_files():
    """Test that function raises ValueError when FASTQ files don't exist."""
    # Create a temporary sample sheet with non-existent files
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("sample,r1_file,r2_file\n")
        f.write("sample1,nonexistent_R1.fastq.gz,nonexistent_R2.fastq.gz\n")
        temp_file = f.name

    try:
        # Should raise ValueError because files don't exist
        with pytest.raises(ValueError):
            get_samples_from_sheet(temp_file, "/path/to/fastqs")
    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_get_samples_from_sheet_missing_columns():
    """Test that function raises ValueError when required columns are missing."""
    # Create a temporary sample sheet with missing columns
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("sample,r1_file\n")  # Missing r2_file column
        f.write("sample1,sample1_R1.fastq.gz\n")
        temp_file = f.name

    try:
        # Should raise ValueError because a required column is missing
        with pytest.raises(ValueError):
            get_samples_from_sheet(temp_file, "/path/to/fastqs")
    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_get_samples_from_sheet_duplicate_samples():
    """Test that function raises ValueError when duplicate sample names are found."""
    # Create a temporary sample sheet with duplicate sample names
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("sample,r1_file,r2_file\n")
        f.write("sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz\n")
        f.write("sample1,sample1_other_R1.fastq.gz,sample1_other_R2.fastq.gz\n")  # Duplicate
        temp_file = f.name

    try:
        # Should raise ValueError because of duplicate sample names
        with pytest.raises(ValueError):
            get_samples_from_sheet(temp_file, "/path/to/fastqs")
    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_get_samples_from_sheet_empty():
    """Test that function raises ValueError when no samples are found."""
    # Create an empty sample sheet (just headers)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("sample,r1_file,r2_file\n")
        temp_file = f.name

    try:
        # Should raise ValueError because no samples were found
        with pytest.raises(ValueError):
            get_samples_from_sheet(temp_file, "/path/to/fastqs")
    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)
