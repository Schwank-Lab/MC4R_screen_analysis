#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Validate the template file to ensure it has all required columns
and expected structure.
"""

import pandas as pd
import argparse
import datetime


def validate_template(template_file, required_columns):
    """
    Validate that the template file has all required columns
    and expected structure.

    Returns (template_df, missing_columns) if file can be read,
    otherwise raises an exception.
    """
    # Read the template file
    template_df = pd.read_csv(template_file)

    # Check if all required columns exist
    missing_columns = [
        col for col in required_columns if col not in template_df.columns]

    if missing_columns:
        print("Error: The following required columns are missing from "
              f"the template file: {', '.join(missing_columns)}")

    # Return the dataframe and any missing columns
    return template_df, missing_columns


def main():
    parser = argparse.ArgumentParser(
        description="Validate template file for self-targeting analysis")
    parser.add_argument(
        '--template', required=True,
        help="Template file path")
    parser.add_argument(
        '--required-columns', required=True,
        help="Comma-separated list of required columns")
    parser.add_argument(
        '--output', required=True,
        help="Output file to write validation results")
    args = parser.parse_args()

    # Convert comma-separated string to list
    required_columns = [
        col.strip() for col in args.required_columns.split(',')]

    # Prepare output content
    output_content = ["Template Validation Report"]
    output_content.append("=" * 30)
    output_content.append(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    output_content.append(f"Template file: {args.template}")
    output_content.append(f"Required columns: {', '.join(required_columns)}")
    output_content.append("-" * 30)

    is_valid = False
    library_size = 0

    try:
        # Validate template
        template_df, missing_columns = validate_template(args.template, required_columns)
        library_size = len(template_df)

        # Add validation results to report
        output_content.append(f"File loaded successfully. Found {library_size} rows.")
        output_content.append(f"Columns in file: {', '.join(template_df.columns)}")

        if missing_columns:
            output_content.append(f"FAILED: Missing required columns: {', '.join(missing_columns)}")
            is_valid = False
        else:
            output_content.append("PASSED: All required columns present.")
            is_valid = True

    except Exception as e:
        output_content.append(f"FAILED: Error reading template file: {e}")
        is_valid = False

    # Add summary
    output_content.append("-" * 30)
    output_content.append(f"Template file validation: {'PASSED' if is_valid else 'FAILED'}")
    if library_size > 0:
        output_content.append(f"Library size: {library_size}")

    # Write the output file with detailed info
    with open(args.output, 'w') as f:
        f.write('\n'.join(output_content))

    # Raise exception if validation failed
    if not is_valid:
        raise ValueError("Template validation failed. See report for details.")


if __name__ == "__main__":
    main()
