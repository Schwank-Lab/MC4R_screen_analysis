import os
import sys
import pytest
import tempfile
import pandas as pd

# Import the script to test
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.validate_template import validate_template, main  # noqa: E402


def test_validate_template(template_file, sample_config):
    """Test that validation passes with a valid template file."""
    required_columns = sample_config["required_columns"]

    # With updated validate_template, it should return a df and empty list of missing columns
    template_df, missing_columns = validate_template(template_file, required_columns)
    assert len(missing_columns) == 0, "Expected no missing columns"
    assert len(template_df) > 0, "Library size should be greater than 0"


def test_validate_template_missing_columns():
    """Test that validation fails when required columns are missing."""
    # Create a temporary file with missing columns
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("uniquename,Original_Sequence\n")
        f.write("test1,ACTGACTGACTGACTG\n")
        temp_file = f.name

    try:
        required_columns = [
            "uniquename",
            "Original_Sequence",
            "Edited-Sequences",
            "non_existent_column",
        ]

        # With updated validate_template, it should return missing columns
        template_df, missing_columns = validate_template(temp_file, required_columns)
        assert len(missing_columns) > 0, "Expected some missing columns"
        assert "Edited-Sequences" in missing_columns
        assert "non_existent_column" in missing_columns
    finally:
        # Clean up
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_validate_template_with_extra_columns():
    """Test validation with template containing extra columns (should pass)."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("uniquename,Original_Sequence,Edited-Sequences,extra_column1,extra_column2\n")
        f.write("test1,ACTGACTG,ACTAACTG,extra1,extra2\n")
        temp_file = f.name

    try:
        # We only require these columns
        required_columns = ["uniquename", "Original_Sequence", "Edited-Sequences"]

        # With updated validate_template, it should return no missing columns
        template_df, missing_columns = validate_template(temp_file, required_columns)
        assert len(missing_columns) == 0, "Expected no missing columns"
        assert len(template_df) == 1, "File size should be 1"
    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_validate_template_nonexistent_file():
    """Test that validation handles nonexistent files appropriately."""
    required_columns = ["uniquename"]

    with pytest.raises((FileNotFoundError, pd.errors.EmptyDataError)):
        validate_template("nonexistent_file.csv", required_columns)


def test_validate_template_main_success(monkeypatch, test_data_path, tmpdir):
    """Test the main function of validate_template.py with a valid template file."""
    # Set up test parameters
    template_file = os.path.join(test_data_path, "20240624_mc4r_library_designed_final.csv")
    output_file = os.path.join(tmpdir, "validation_result.txt")
    required_columns = "uniquename,Spacer-Sequence,barcode,RTrevcomp,PBSrevcomp,wide_initial_target"

    # Mock command-line arguments for a successful validation
    args = [
        'validate_template.py',
        '--template', template_file,
        '--required-columns', required_columns,
        '--output', output_file
    ]
    monkeypatch.setattr('sys.argv', args)

    # Call the main function (should pass validation and exit normally)
    main()

    # Verify output file was created with expected content
    assert os.path.exists(output_file), f"Output file not created: {output_file}"
    with open(output_file, 'r') as f:
        content = f.read()
    assert "PASSED" in content, f"Expected 'PASSED' in output, got: {content}"
    assert "All required columns present" in content, "Expected positive confirmation about columns"
    assert "Template Validation Report" in content, "Expected report header"


def test_validate_template_main_failure(monkeypatch, test_data_path, tmpdir, capsys):
    """Test the main function of validate_template.py with missing required columns."""
    # Set up test parameters
    template_file = os.path.join(test_data_path, "20240624_mc4r_library_designed_final.csv")
    output_file = os.path.join(tmpdir, "validation_failure.txt")

    # Specify some columns that definitely don't exist in the template
    non_existent_columns = "uniquename,non_existent_column1,non_existent_column2"

    # Mock command-line arguments for a failed validation
    args = [
        'validate_template.py',
        '--template', template_file,
        '--required-columns', non_existent_columns,
        '--output', output_file
    ]
    monkeypatch.setattr('sys.argv', args)

    # Call the main function with missing columns - should raise ValueError
    with pytest.raises(ValueError) as excinfo:
        main()
    assert "Template validation failed" in str(
        excinfo.value), "Expected ValueError with message about validation failure"

    # Verify output file was created with expected content
    assert os.path.exists(output_file), f"Output file not created: {output_file}"
    with open(output_file, 'r') as f:
        content = f.read()

    # Check for detailed failure information in the report
    assert "FAILED" in content, f"Expected 'FAILED' in output, got: {content}"
    assert "Missing required columns" in content, "Expected details about missing columns"
    assert "non_existent_column1" in content, "Expected specific missing column name"
    assert "non_existent_column2" in content, "Expected specific missing column name"

    # Check that error messages were printed
    captured = capsys.readouterr()
    assert "required columns are missing" in captured.out, "Expected error message about missing columns"


def test_validate_template_main_file_not_found(monkeypatch, tmpdir, capsys):
    """Test the main function of validate_template.py with a non-existent template file."""
    # Set up test parameters
    non_existent_template = os.path.join(tmpdir, "non_existent_file.csv")
    output_file = os.path.join(tmpdir, "validation_non_existent.txt")
    required_columns = "uniquename,Spacer-Sequence"

    # Mock command-line arguments for a non-existent template
    args = [
        'validate_template.py',
        '--template', non_existent_template,
        '--required-columns', required_columns,
        '--output', output_file
    ]
    monkeypatch.setattr('sys.argv', args)

    # Call the main function with a non-existent file - should raise ValueError
    with pytest.raises(ValueError) as excinfo:
        main()
    assert "validation failed" in str(excinfo.value).lower(), "Expected error message about validation failure"

    # Verify output file was created with expected content
    assert os.path.exists(output_file), f"Output file not created: {output_file}"
    with open(output_file, 'r') as f:
        content = f.read()

    # Check for detailed error information
    assert "FAILED" in content, f"Expected 'FAILED' in output, got: {content}"
    assert "Error reading template file" in content, "Expected details about file error"
