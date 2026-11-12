import os
import pandas as pd


def get_samples_from_sheet(sample_sheet_path, fastq_folder):
    """Read sample sheet and return sample names and mapping dictionary."""
    if not os.path.exists(sample_sheet_path):
        raise ValueError(f"Sample sheet {sample_sheet_path} not found.")

    sample_df = pd.read_csv(sample_sheet_path)

    # Check required columns
    required_columns = ['sample', 'r1_file', 'r2_file']
    missing_columns = [col for col in required_columns if col not in sample_df.columns]
    if missing_columns:
        raise ValueError(f"Sample sheet is missing required columns: {', '.join(missing_columns)}")

    # Check for duplicate sample names
    if sample_df['sample'].duplicated().any():
        duplicate_samples = sample_df['sample'][sample_df['sample'].duplicated()].unique().tolist()
        raise ValueError(f"Duplicate sample names found in sample sheet: {', '.join(duplicate_samples)}")

    samples = sample_df['sample'].tolist()

    # Check for underscores in sample names (to avoid Snakemake wildcard issues)
    underscored_samples = [s for s in samples if "_" in s]
    if underscored_samples:
        raise ValueError(
            "Sample names must not contain underscores ('_') because Snakemake wildcards may parse them incorrectly.\n"
            f"The following sample names are problematic:\n{', '.join(underscored_samples)}\n"
            "Consider renaming your samples using dashes (-) or another character instead."
        )

    # Create dictionary mapping sample names to file paths
    sample_dict = {}
    missing_files = []

    for _, row in sample_df.iterrows():
        r1_path = os.path.join(fastq_folder, row['r1_file'])
        r2_path = os.path.join(fastq_folder, row['r2_file'])

        # Check if files exist
        if not os.path.exists(r1_path):
            missing_files.append(r1_path)
        if not os.path.exists(r2_path):
            missing_files.append(r2_path)

        sample_dict[row['sample']] = {
            'r1': r1_path,
            'r2': r2_path
        }

    # Raise error if any files are missing
    if missing_files:
        if len(missing_files) <= 5:
            raise ValueError("The following FASTQ files were not found:\n" +
                             "\n".join(missing_files))
        else:
            raise ValueError("{len(missing_files)} FASTQ files were not found, including:\n" +
                             "\n".join(missing_files[:5]) +
                             f"\n... and {len(missing_files) - 5} more")

    # Raise error if no samples found
    if not samples:
        raise ValueError("No samples found in the sample sheet.")

    return samples, sample_dict
