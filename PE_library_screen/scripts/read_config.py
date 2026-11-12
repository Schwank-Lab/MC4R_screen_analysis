#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read configuration from CSV file
and convert to a format suitable for Snakemake.
"""

import pandas as pd
import yaml
import argparse


def read_csv_config(config_file):
    """Read configuration from CSV file and convert to dictionary."""
    try:
        config_df = pd.read_csv(config_file)
        config_dict = {}

        # Process each parameter
        for _, row in config_df.iterrows():
            # Get parameter and value, ignoring the section column
            param = row['parameter']
            value = process_value(row['value'])

            # Add to config dictionary (flat structure)
            config_dict[param] = value

        # Special handling for comma-separated lists
        if 'required_columns' in config_dict:
            config_dict['required_columns'] = [
                col.strip() for col in config_dict['required_columns'].split(',')]

        return config_dict

    except Exception as e:
        print(f"Error reading config file: {e}")
        raise


def process_value(value):
    """Process a value from the CSV to convert to appropriate type."""
    if isinstance(value, str):
        # Convert numeric strings to numbers
        try:
            if '.' in value:
                return float(value)
            else:
                return int(value)
        except ValueError:
            # Keep as string if not a number
            return value
    return value


def write_yaml_config(config_dict, output_file):
    """Write configuration dictionary to YAML file."""
    try:
        with open(output_file, 'w') as yaml_file:
            yaml.dump(config_dict, yaml_file, default_flow_style=False)
    except Exception as e:
        raise IOError(f"Failed to write YAML config to {output_file}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert CSV config to YAML for Snakemake")
    parser.add_argument(
        '--csv', required=True, help="Input CSV config file")
    parser.add_argument(
        '--yaml', required=True, help="Output YAML config file")
    args = parser.parse_args()

    config_dict = read_csv_config(args.csv)
    write_yaml_config(config_dict, args.yaml)
    print(f"Configuration converted from {args.csv} to {args.yaml}")


if __name__ == "__main__":
    main()
