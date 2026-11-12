import os
import sys
import pandas as pd
import pytest
import numpy as np

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Import the module to test
from scripts.jupyter_analysis_helper import (  # noqa: E402
    list_files, create_scatter_plot, correlation_heatmap, create_histogram)


def test_list_files(tmp_path):
    """Test the list_files function with a temporary directory."""
    # Create temporary files
    (tmp_path / "file1_combined_test.txt").write_text("content1")
    (tmp_path / "file2_test.txt").write_text("content2")
    (tmp_path / "file3_combined_data.csv").write_text("content3")

    # Test with default pattern ('combined')
    files = list_files(str(tmp_path))
    assert len(files) == 2
    assert "file1_combined_test.txt" in files
    assert "file3_combined_data.csv" in files
    assert "file2_test.txt" not in files

    # Test with different pattern
    files = list_files(str(tmp_path), pattern='test')
    assert len(files) == 2
    assert "file1_combined_test.txt" in files
    assert "file2_test.txt" in files
    assert "file3_combined_data.csv" not in files

    # Test with pattern that matches no files
    files = list_files(str(tmp_path), pattern='nonexistent')
    assert len(files) == 0


def test_list_files_empty_directory(tmp_path):
    """Test list_files with an empty directory."""
    # Test with empty directory
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()

    files = list_files(str(empty_dir))
    assert len(files) == 0


def test_create_scatter_plot(monkeypatch, capsys, tmp_path):
    """Test the create_scatter_plot function."""
    # Mock plt.show to prevent plots from displaying during tests
    monkeypatch.setattr('matplotlib.pyplot.show', lambda: None)

    # Mock titlesize global variable if it's not defined in the test environment
    if 'titlesize' not in globals():
        current_module = sys.modules[__name__]
        setattr(current_module, 'titlesize', 16)

    # Create a test DataFrame with some noise to avoid perfectly linear data
    # which causes gaussian_kde to fail with singular covariance matrix
    np.random.seed(42)  # For reproducibility
    x_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    y_values = [2.1, 3.9, 6.2, 7.8, 10.1, 12.2, 13.9]  # Not perfectly linear

    # Add some scatter points to avoid the KDE issue
    more_x = np.random.uniform(1, 7, 10)
    more_y = 2 * more_x + np.random.normal(0, 0.5, 10)

    x_values.extend(more_x)
    y_values.extend(more_y)

    data = {
        'x_values': x_values,
        'y_values': y_values
    }
    df = pd.DataFrame(data)

    # Mock gaussian_kde to avoid the actual computation
    def mock_kde(xy):
        def inner(points):
            # Just return some random density values
            return np.random.uniform(0, 1, len(points[0]))
        return inner

    monkeypatch.setattr('scipy.stats.gaussian_kde', mock_kde)

    # Test without saving
    create_scatter_plot(df, 'x_values', 'y_values', 'Test Plot', 'X', 'Y')

    # Capture print output
    captured = capsys.readouterr()
    assert "Full length of df:" in captured.out
    assert "Plotted number after dropping NaNs:" in captured.out

    # Verify correlation values are calculated and printed
    assert "Spearman correlation:" in captured.out
    assert "Pearson correlation:" in captured.out

    # Test with savepath
    save_path = tmp_path / "test_plot.png"

    # Replace savefig with a version that just creates an empty file
    def mock_savefig(path, *args, **kwargs):
        with open(path, 'w'):
            pass

    monkeypatch.setattr('matplotlib.pyplot.savefig', mock_savefig)

    create_scatter_plot(df, 'x_values', 'y_values', 'Test Plot', 'X', 'Y', savepath=str(save_path))
    assert os.path.exists(save_path)

    # Test with NaN values
    df_with_nan = df.copy()
    df_with_nan.loc[0, 'y_values'] = np.nan

    create_scatter_plot(df_with_nan, 'x_values', 'y_values', 'Test Plot with NaN', 'X', 'Y')

    # Capture print output again
    captured = capsys.readouterr()
    assert "Full length of df:" in captured.out
    assert "Plotted number after dropping NaNs:" in captured.out

    # Test with empty DataFrame after dropping NaNs
    df_empty = pd.DataFrame({'x_values': [np.nan, np.nan], 'y_values': [1, 2]})

    # Use pytest.raises as a context manager to catch the expected exception
    with pytest.raises(ValueError):
        create_scatter_plot(df_empty, 'x_values', 'y_values', 'Empty Test Plot', 'X', 'Y')

    # Capture print output for empty DataFrame
    captured = capsys.readouterr()
    assert "Full length of df: 2" in captured.out
    assert "Plotted number after dropping NaNs: 0" in captured.out


def test_correlation_heatmap(monkeypatch, capsys, tmp_path):
    """Test the replicate_correlation_heatmap function."""
    # Mock plt.show to prevent plots from displaying during tests
    monkeypatch.setattr('matplotlib.pyplot.show', lambda: None)

    # Mock seaborn functions
    def mock_set_theme(*args, **kwargs):
        pass

    def mock_heatmap(*args, **kwargs):
        class DummyAxes:
            def set(self, **kwargs):
                return self

            def tick_params(self, **kwargs):
                pass

            def set_facecolor(self, color):
                pass
        return DummyAxes()

    monkeypatch.setattr('seaborn.set_theme', mock_set_theme)
    monkeypatch.setattr('seaborn.heatmap', mock_heatmap)

    # Create test data with replicate columns (containing numbers in the names)
    data = {
        'HEK293T_PE2_D3_1_edited_percentage': [20.0, 30.0, 40.0, 50.0, 25.0],
        'HEK293T_PE2_D3_2_edited_percentage': [19.6, 29.4, 39.2, 49.0, 24.5],
        'HEK293T_PE2_D3_3_edited_percentage': [20.4, 30.6, 40.8, 51.0, 25.5]
    }
    df = pd.DataFrame(data)
    cols = list(data.keys())

    # Create labels for the axes
    labels = ['Rep1', 'Rep2', 'Rep3']

    # Test without saving
    correlation_heatmap(df, cols, 'Test Heatmap', x_labels=labels, y_labels=labels)

    # Test with savepath
    save_path = tmp_path / "test_heatmap.png"

    # Replace savefig with a version that just creates an empty file
    def mock_savefig(path, *args, **kwargs):
        with open(path, 'w'):
            pass

    monkeypatch.setattr('matplotlib.pyplot.savefig', mock_savefig)

    correlation_heatmap(df, cols, 'Test Heatmap', x_labels=labels, y_labels=labels, savepath=str(save_path))
    assert os.path.exists(save_path)

    # Test with all NaN values in one column
    df_with_nan = df.copy()
    df_with_nan['HEK293T_PE2_D3_1_edited_percentage'] = np.nan

    # Use pytest.raises as a context manager to catch the expected exception
    with pytest.raises(ValueError):
        correlation_heatmap(df_with_nan, cols, 'NaN Test Heatmap', x_labels=labels, y_labels=labels)

    # Test with empty DataFrame
    empty_df = pd.DataFrame(columns=cols)
    with pytest.raises(ValueError):
        correlation_heatmap(empty_df, cols, 'Empty Test Heatmap', x_labels=labels, y_labels=labels)


def test_create_histogram(monkeypatch, capsys, tmp_path):
    """Test the create_histogram function."""
    # Mock plt.show to prevent plots from displaying during tests
    monkeypatch.setattr('matplotlib.pyplot.show', lambda: None)

    # Mock titlesize global variable if it's not defined in the test environment
    if 'titlesize' not in globals():
        current_module = sys.modules[__name__]
        setattr(current_module, 'titlesize', 16)

    # Create a test DataFrame with some data
    data = {
        'edited_percentage': [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0],
        'unintended_editing_percentage': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    }
    df = pd.DataFrame(data)

    # Test without saving
    create_histogram(df, 'edited_percentage', 'Test Histogram')

    # Capture print output
    captured = capsys.readouterr()
    assert "Full length of df:" in captured.out
    assert "Plotted number after dropping NaNs:" in captured.out

    # Test with NaN values
    df_with_nan = df.copy()
    df_with_nan.loc[0, 'edited_percentage'] = np.nan

    create_histogram(df_with_nan, 'edited_percentage', 'Test Histogram with NaN')

    # Capture print output again
    captured = capsys.readouterr()
    assert "Full length of df:" in captured.out
    assert "Plotted number after dropping NaNs:" in captured.out
