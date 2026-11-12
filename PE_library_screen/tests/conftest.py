import os
import sys
import pytest
import shutil
import tempfile

# Add the project root directory to Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


@pytest.fixture
def test_data_path():
    """Return path to test data directory."""
    return os.path.join(os.path.dirname(__file__), "test_data")


@pytest.fixture
def template_file(test_data_path):
    """Path to the template file."""
    return os.path.join(test_data_path, "20240624_mc4r_library_designed_final.csv")


@pytest.fixture
def test_fastq_files(test_data_path):
    test_data_path += "/fastq"
    """Dict containing paths to test FASTQ files."""
    return {
        "R1": os.path.join(test_data_path, "test_MC4R-lib-control1_R1.fastq.gz"),
        "R2": os.path.join(test_data_path, "test_MC4R-lib-control1_R2.fastq.gz"),
    }


@pytest.fixture
def sample_config(test_data_path):
    """Load the actual config file for testing."""
    import os
    from scripts.read_config import read_csv_config

    config_file = os.path.join(test_data_path, "self_targeting_config.csv")
    config = read_csv_config(config_file)
    config["library_design_file"] = os.path.join(test_data_path, os.path.basename(config["library_design_file"]))
    return config


@pytest.fixture
def test_output_dirs():
    """Create temporary directories for test outputs."""
    temp_dir = tempfile.mkdtemp()

    # Create subdirectories
    fastq_dir = os.path.join(temp_dir, "fastqs")
    output_dir = os.path.join(temp_dir, "temp_output")
    log_dir = os.path.join(temp_dir, "logs")
    analysis_dir = os.path.join(temp_dir, "analysis")

    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(analysis_dir, exist_ok=True)

    # Copy test FASTQ files to the temporary fastq directory
    test_data = os.path.join(os.path.dirname(__file__), "test_data/fastq")
    for f in os.listdir(test_data):
        if f.endswith(".fastq.gz"):
            shutil.copy(
                os.path.join(test_data, f),
                os.path.join(fastq_dir, f.replace("test_", "")),
            )

    # Return dict with paths
    paths = {
        "root": temp_dir,
        "fastq": fastq_dir,
        "output": output_dir,
        "log": log_dir,
        "analysis": analysis_dir,
    }

    yield paths

    # Cleanup
    shutil.rmtree(temp_dir)
