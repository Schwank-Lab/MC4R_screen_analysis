import os
import sys
import pytest
import yaml

# Import the script to test
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.read_config import read_csv_config, process_value, write_yaml_config  # noqa: E402


def test_process_value():
    """Test the process_value function handles different input types correctly."""
    # Test integer conversion
    assert process_value("123") == 123

    # Test float conversion
    assert process_value("123.45") == 123.45

    # Test string that should remain a string
    assert process_value("abc123") == "abc123"

    # Test empty string
    assert process_value("") == ""

    # Test None
    assert process_value(None) is None


def test_read_csv_config(test_data_path):
    """Test that the config reader correctly parses a CSV config file."""

    # Construct path to the config file in test_data directory
    config_file = os.path.join(test_data_path, "self_targeting_config.csv")

    # Read the config file
    config = read_csv_config(config_file)

    # Test basic parameter reading
    assert "library_name" in config
    assert config["library_name"] == "fulllib"
    assert "fastq_folder" in config
    assert config["fastq_folder"] == "fastqs/"

    # Test value type conversion
    assert isinstance(config["library_size"], int)
    assert config["library_size"] == 500
    assert isinstance(config["min_quality"], float)
    assert config["min_quality"] == 27.5
    assert isinstance(config["barcode_length"], int)
    assert config["barcode_length"] == 6
    assert isinstance(config["barcode_start"], int)
    assert config["barcode_start"] == 0

    # Test special handling of comma-separated lists
    assert "required_columns" in config
    assert isinstance(config["required_columns"], list)
    assert len(config["required_columns"]) > 0
    assert "uniquename" in config["required_columns"]
    assert "barcode" in config["required_columns"]

    # Test that adapter sequences remain strings
    assert isinstance(config["target_adapter"], str)
    assert config["target_adapter"] == "GTCATAGGAGTCAC"


def test_read_csv_config_nonexistent_file():
    """Test that the config reader handles nonexistent files appropriately."""
    with pytest.raises(Exception):
        read_csv_config("nonexistent_file.csv")


def test_write_yaml_config(tmpdir):
    """Test writing config to YAML file."""
    # Create a simple config dict
    config = {"key1": "value1", "key2": 123}
    output_file = os.path.join(tmpdir, "test_config.yaml")

    # Test writing
    write_yaml_config(config, output_file)

    # Verify file was created and contains expected content
    assert os.path.exists(output_file)
    with open(output_file, 'r') as f:
        loaded_config = yaml.safe_load(f)
    assert loaded_config == config


def test_write_yaml_config_invalid_path():
    """Test writing to an invalid path."""
    # Try to write to a location that can't be written to
    invalid_path = "/nonexistent/directory/file.yaml"  # Path that doesn't exist

    with pytest.raises(Exception):  # This will catch IOError, OSError, etc.
        write_yaml_config({"test": "value"}, invalid_path)


def test_read_config_main_function(monkeypatch, test_data_path, tmpdir):
    """Test the main function with command-line arguments."""
    csv_file = os.path.join(test_data_path, "self_targeting_config.csv")
    yaml_file = os.path.join(tmpdir, "output.yaml")
    monkeypatch.setattr('sys.argv', ['read_config.py', '--csv', csv_file, '--yaml', yaml_file])

    # Call the main function
    from scripts.read_config import main
    main()

    # Verify YAML was created
    assert os.path.exists(yaml_file)

    # Verify the content matches what read_csv_config produces
    with open(yaml_file, 'r') as f:
        loaded_config = yaml.safe_load(f)

    expected_config = read_csv_config(csv_file)
    assert loaded_config == expected_config
