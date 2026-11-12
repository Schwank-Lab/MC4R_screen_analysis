#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analysis of self-targeting library.

@author: nimath
"""
import pandas as pd
from Bio import SeqIO
import gzip
import numpy as np
import os
import time
from Bio.Seq import Seq
import argparse
from itertools import combinations
from collections import Counter


def lookup(spacertemplate, exttemplate, barcodetemplate, rtt_lookup_length=6):
    """Generate lookup dictionaries for Spacer, RTTstart and barcodes."""
    # Create dictionaries mapping sequences to their indices
    spacerlookup = {}
    for i, z in enumerate(spacertemplate):
        spacerlookup.setdefault(z, []).append(i)

    rttlookup = {}
    for i, z in enumerate(exttemplate):
        # Use configurable length for RTT lookup
        rttlookup.setdefault(z[:rtt_lookup_length], []).append(i)

    barcodelookup = {}
    for i, z in enumerate(barcodetemplate):
        barcodelookup.setdefault(z, []).append(i)

    return spacerlookup, rttlookup, barcodelookup


def process_fastq_file(path, filename, diseasedict, match_func,
                       report_interval=100000):
    """Process a single FASTQ file with standardized reporting."""
    full_path = os.path.join(path, filename)
    count = 0
    counttemp = 0

    print(f'Processing {filename}...')
    with gzip.open(full_path, "rt") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            count += 1
            counttemp += 1
            match_func(seq_record, diseasedict)

            if counttemp == report_interval:
                print(f"Processed {count} reads")
                counttemp = 0

    print(f'Finished processing {filename}')
    return count


def sample_read_lengths(path, filename, max_reads=50):
    """Sample read lengths from a FASTQ file."""
    read_lengths = []
    full_path = os.path.join(path, filename)

    try:
        with gzip.open(full_path, "rt") as fasta_file:
            for i, seq_record in enumerate(SeqIO.parse(fasta_file, 'fastq')):
                if i >= max_reads:
                    break
                read_lengths.append(len(seq_record.seq))
    except Exception as e:
        print(f"Warning: Could not sample read lengths from {filename}: {e}")

    return read_lengths


def importfiles(shortname, path, spacerlookup, rttlookup, barcodelookup,
                filtered='', barcode_start=0, barcode_length=6, rtt_lookup_length=6):
    """Import FASTQ files and match reads against templates."""
    diseasedict = {}

    # Check read lengths to validate config (sample first 50 reads)
    read_lengths = sample_read_lengths(path, f"{shortname}_spacer_5trim{filtered}.fastq.gz", max_reads=50)

    # Define matching functions for each file type
    def match_spacer(seq_record, diseasedict):
        spacer = seq_record.seq
        spacermatch = spacerlookup.get("G"+str(spacer))
        identifier = seq_record.id
        if identifier in diseasedict:
            diseasedict[identifier]["spacermatch"] = spacermatch
        else:
            diseasedict[identifier] = {'spacermatch': spacermatch}

    def match_barcode(seq_record, diseasedict):
        barcode = seq_record.seq[barcode_start:barcode_start+barcode_length]
        barcodematch = barcodelookup.get(barcode)
        identifier = seq_record.id
        if identifier in diseasedict:
            diseasedict[identifier]["barcodematch"] = barcodematch

    def match_extension(seq_record, diseasedict):
        rtt = seq_record.seq
        exttemplate = str(rtt[:rtt_lookup_length])
        exttemplatematch = rttlookup.get(exttemplate)
        identifier = seq_record.id
        if identifier in diseasedict:
            diseasedict[identifier]["exttemplatematch"] = exttemplatematch
            # check if extension contains tevopreq sequence,
            # otherwise it might be recombination
            diseasedict[identifier]["tevopreq_in_extension"] = ('CGCGGTTCTATCTA' in str(rtt))  # gives True or False

    # Expose the functions for testing
    importfiles.match_spacer = match_spacer
    importfiles.match_barcode = match_barcode
    importfiles.match_extension = match_extension

    # Process each file type
    process_fastq_file(path, f"{shortname}_Spacer{filtered}.fastq.gz",
                       diseasedict, match_spacer)

    process_fastq_file(path, f"{shortname}_target_5trim{filtered}.fastq.gz",
                       diseasedict, match_barcode)

    process_fastq_file(path, f"{shortname}_RTTtrim{filtered}.fastq.gz",
                       diseasedict, match_extension)

    return diseasedict, read_lengths


def get_intermediate_strings(str1, str2):
    """
    Find all intermediate strings between a mutated and wildtype sequence.
    Used to quantify incomplete editing.
    """
    intermediate_strings = []
    changed_positions = []
    # Find positions that differ between the two strings
    diff_indexes = [i for i in range(len(str1)) if str1[i] != str2[i]]

    # Generate all possible intermediate combinations
    for i in range(len(diff_indexes)):
        for comb in combinations(diff_indexes, i + 1):
            intermediate = list(str1)
            changes = []
            # Apply changes one by one
            for index in comb:
                intermediate[index] = str2[index]
                changes.append(index)
            intermediate = "".join(intermediate)
            if intermediate != str1 and intermediate != str2:
                intermediate_strings.append(intermediate)
                changed_positions.append(changes)

    return intermediate_strings, changed_positions


def sorted_counter(lst):
    """Convert a list to a sorted Counter dictionary."""
    counter = Counter(lst)
    return dict(sorted(counter.items(), key=lambda x: x[1], reverse=True))


def preprocess_template(template_file):
    """Load and preprocess the template file."""
    template_df = pd.read_csv(template_file)

    # Process barcode and sequences
    template_df['barcode_rev'] = template_df['barcode'].apply(
        lambda x: str(Seq(x).reverse_complement())
    )
    template_df['RTTrevcomp_original'] = template_df['RTrevcomp']
    template_df['RTTrevcomp'] = template_df['RTrevcomp'].apply(
        lambda x: x.upper())  # edits are in lowercase
    template_df['extension'] = template_df['RTTrevcomp'] + \
        template_df['PBSrevcomp']
    template_df['extension_15bp'] = template_df['extension'].apply(
        lambda x: x[:15])
    template_df['Name'] = template_df['uniquename']

    # Process intermediate states for MultibpReplacement
    process_multibp_replacements(template_df)

    return template_df


def process_multibp_replacements(template_df):
    """Process intermediate states for MultibpReplacement variants."""
    # Calculate intermediate sequences
    template_df['wide_mutated_target_intermediate'] = template_df.apply(
        lambda x: get_intermediate_strings(x['wide_initial_target'],
                                           x['wide_mutated_target'])[0]
        if x['Mutation_Type'] == 'MultibpReplacement' else None, axis=1)

    # Calculate intermediate positions
    template_df['deepeditposition_intermediate'] = template_df.apply(
        lambda x: ['_'.join(str(item) for item in inner_list) for inner_list in
                   get_intermediate_strings(x['wide_initial_target'],
                                            x['wide_mutated_target'])[1]]
        if x['Mutation_Type'] == 'MultibpReplacement' else None, axis=1)


def initialize_template_counts(template_df):
    """Initialize count columns in the template dataframe."""
    template_df["editedcount"] = 0
    template_df["uneditedcount"] = 0
    template_df["indelcount"] = 0
    template_df["nickindelcount"] = 0
    template_df["intermediatecount"] = 0
    template_df["intermediatepositions"] = np.empty(
        (len(template_df), 0)).tolist()
    template_df["beforeflapindelcount"] = 0
    template_df["totalreads"] = 0
    template_df['barcodes'] = np.empty((len(template_df), 0)).tolist()
    template_df["skippedcount"] = 0  # for reads that were skipped due to short sequences
    return template_df


def get_column_indices(template_df):
    """Get indices of important columns in the template dataframe."""
    return {
        'uneditedcount': template_df.columns.get_loc("uneditedcount"),
        'editedcount': template_df.columns.get_loc("editedcount"),
        'indelcount': template_df.columns.get_loc("indelcount"),
        'nickindelcount': template_df.columns.get_loc("nickindelcount"),
        'intermediatecount': template_df.columns.get_loc("intermediatecount"),
        'intermediatepos': template_df.columns.get_loc(
            "intermediatepositions"),
        'beforeflapindel': template_df.columns.get_loc("beforeflapindelcount"),
        'wide_initial': template_df.columns.get_loc("wide_initial_target"),
        'wide_mutated': template_df.columns.get_loc("wide_mutated_target"),
        'wide_intermediate': template_df.columns.get_loc(
            "wide_mutated_target_intermediate"),
        'deepeditpos_inter': template_df.columns.get_loc(
            "deepeditposition_intermediate"),
        'RTTlength': template_df.columns.get_loc("RTlength"),
        'correction_type': template_df.columns.get_loc("Correction_Type"),
        'mutation_type': template_df.columns.get_loc("Mutation_Type"),
        'correction_length': template_df.columns.get_loc("Correction_Length"),
        'skippedcount': template_df.columns.get_loc("skippedcount")
    }


def process_matches(diseasedict, read_length=300):
    """Process matches to identify unique matches and recombinations."""
    # Convert to dataframe and process matches
    diseasedf = pd.DataFrame.from_dict(diseasedict, orient='index')

    # Filter and process matches
    diseasedf.dropna(
        subset=['spacermatch', 'exttemplatematch', 'barcodematch'],
        inplace=True)
    diseasedf['match'] = [list(set(a).intersection(set(b))) for a, b in zip(
        diseasedf.spacermatch, diseasedf.barcodematch)]
    diseasedf['matchwithextension'] = [list(set(a).intersection(
                                       set(b), set(c))) for a, b, c in
                                       zip(diseasedf.spacermatch,
                                           diseasedf.exttemplatematch,
                                           diseasedf.barcodematch)]
    diseasedf['recombination'] = diseasedf.apply(
        lambda x: len(x.matchwithextension) == 0, axis=1)

    # For longer reads (≥ 250bp), filter for tevopreq in extension
    if read_length >= 250:
        print("Filtering reads for tevopreq sequence in extension (long read mode)...")
        # only keep reads where we still have tevopreq in extension
        # otherwise it might be recombination of extension with target
        diseasedf_tevo = diseasedf[diseasedf['tevopreq_in_extension']].copy()
        print(f"Reads after tevopreq filtering: {len(diseasedf_tevo)} of {len(diseasedf)}",
              f"{len(diseasedf_tevo)/len(diseasedf)*100:.1f}%)")
    else:
        # For shorter reads, skip tevopreq filtering but add a warning
        print(f"WARNING: Skipping tevopreq sequence filtering (short read mode, {read_length}bp)")
        print("Recombination detection at the extension level will be less accurate.")
        diseasedf_tevo = diseasedf.copy()

    # Create dataframe of unique matches
    final_diseasedf = diseasedf_tevo[diseasedf_tevo['match'].map(
        lambda d: len(d)) > 0][['match', 'matchwithextension']]

    return diseasedf, diseasedf_tevo, final_diseasedf


def extract_sequences(targetend, variantindex, templatenumpy, correction_type,
                      RTTlength, correction_length,
                      target_start_to_nick_distance,
                      nick_pos_in_wide_initial_target, basesbeforenick,
                      baseafterflap, wide_initial_idx, wide_mutated_idx):
    """
    Extract sequence regions for analysis based on correction type.
    Uses a parameterized approach to handle different correction types
    with minimal code duplication.
    """
    # Determine offsets based on correction type
    if correction_type == 'Replacement':
        # No additional offset for replacement
        wt_length_offset = 0
        mut_length_offset = 0
    elif correction_type == 'Deletion':
        # Include the deletion length for WT, nothing for MUT
        wt_length_offset = correction_length
        mut_length_offset = 0
    elif correction_type == 'Insertion':
        # Subtract insertion length for WT, nothing for MUT
        wt_length_offset = -correction_length
        mut_length_offset = 0
    else:
        raise ValueError(f"Unknown correction type: {correction_type}")

    # Common starting positions - simplified calculation
    seq_start = target_start_to_nick_distance - basesbeforenick
    template_start = nick_pos_in_wide_initial_target - basesbeforenick

    # Calculate required lengths for each sequence
    target_required_length = seq_start + basesbeforenick + RTTlength + \
        max(wt_length_offset, mut_length_offset) + baseafterflap
    wt_template_required_length = template_start + basesbeforenick + RTTlength + wt_length_offset + baseafterflap
    mut_template_required_length = template_start + basesbeforenick + RTTlength + mut_length_offset + baseafterflap

    # Validate sequence lengths
    if len(targetend) < target_required_length:
        raise ValueError(f"Target sequence is too short. Length: {len(targetend)}, required: {target_required_length}")

    wt_template = templatenumpy[variantindex, wide_initial_idx]
    mut_template = templatenumpy[variantindex, wide_mutated_idx]

    if len(wt_template) < wt_template_required_length:
        raise ValueError(
            f"WT template sequence too short. Length: {len(wt_template)}, required: {wt_template_required_length}")

    if len(mut_template) < mut_template_required_length:
        raise ValueError(
            f"MUT template sequence too short. Length: {len(mut_template)}, required: {mut_template_required_length}")

    # Extract all sequences with appropriate offsets
    sequences = {}

    # Extract sequences from target read
    sequences['sequenceWT'] = targetend[seq_start:seq_start + basesbeforenick + RTTlength +
                                        wt_length_offset + baseafterflap]
    sequences['sequenceWT_beforeflap'] = targetend[seq_start:seq_start + basesbeforenick +
                                                   RTTlength +
                                                   wt_length_offset - 2]
    sequences['sequenceWT_nick'] = targetend[seq_start:seq_start + basesbeforenick + 2]

    sequences['sequenceMUT'] = targetend[seq_start:seq_start + basesbeforenick + RTTlength +
                                         mut_length_offset + baseafterflap]
    sequences['sequenceMUT_beforeflap'] = targetend[seq_start:seq_start + basesbeforenick +
                                                    RTTlength +
                                                    mut_length_offset - 2]
    sequences['sequenceMUT_nick'] = targetend[seq_start:seq_start + basesbeforenick + 2]

    # Extract control sequences from templates
    sequences['controlWT'] = templatenumpy[variantindex, wide_initial_idx][
        template_start:template_start + basesbeforenick + RTTlength + wt_length_offset +
        baseafterflap].upper()
    sequences['controlWT_beforeflap'] = templatenumpy[
        variantindex, wide_initial_idx][
        template_start:template_start + basesbeforenick + RTTlength +
        wt_length_offset - 2].upper()
    sequences['controlWT_nick'] = templatenumpy[
        variantindex, wide_initial_idx][
            template_start:template_start + basesbeforenick + 2].upper()

    sequences['controlMUT'] = templatenumpy[
        variantindex, wide_mutated_idx][
            template_start:template_start + basesbeforenick + RTTlength +
            mut_length_offset + baseafterflap].upper()
    sequences['controlMUT_beforeflap'] = templatenumpy[
        variantindex, wide_mutated_idx][
            template_start:template_start + basesbeforenick + RTTlength +
            mut_length_offset - 2].upper()
    sequences['controlMUT_nick'] = templatenumpy[
        variantindex, wide_mutated_idx][
            template_start:template_start + basesbeforenick + 2].upper()

    return sequences


def analyze_sequence(sequences, templatenumpy, variantindex, col_indices,
                     mutation_type, intermediate_data=None):
    """
    Analyze a sequence to determine edit status.
    Centralizes the analysis logic to avoid code duplication.
    """
    # Check for exact matches against controls
    if sequences['sequenceWT'] == sequences['controlWT']:
        templatenumpy[variantindex, col_indices['uneditedcount']] += 1
        return
    elif sequences['sequenceMUT'] == sequences['controlMUT']:
        templatenumpy[variantindex, col_indices['editedcount']] += 1
        return

    # If we get here, it's an indel
    templatenumpy[variantindex, col_indices['indelcount']] += 1

    # For multi-bp replacements, check for intermediate states
    if mutation_type == 'MultibpReplacement' and intermediate_data:
        intermediate_editlist = intermediate_data['editlist']
        deepeditposition_intermediatelist = intermediate_data['positions']

        # Check if sequence matches any intermediate state
        for i, elem in enumerate(intermediate_editlist):
            if (sequences['sequenceWT'] in elem or
                    sequences['sequenceMUT'] in elem):
                templatenumpy[variantindex,
                              col_indices['intermediatecount']] += 1
                templatenumpy[variantindex,
                              col_indices['intermediatepos']].append(
                    deepeditposition_intermediatelist[i])
                break

    # Check for indels around nick site
    if (sequences['sequenceWT_nick'] != sequences['controlWT_nick'] and
            sequences['sequenceMUT_nick'] != sequences['controlMUT_nick']):
        templatenumpy[variantindex, col_indices['nickindelcount']] += 1

    # Check for indels before flap
    if (sequences['sequenceWT_beforeflap'] != sequences['controlWT_beforeflap']
        and sequences['sequenceMUT_beforeflap'] !=
            sequences['controlMUT_beforeflap']):
        templatenumpy[variantindex, col_indices['beforeflapindel']] += 1


def calculate_statistics(templatedf):
    """Calculate summary statistics from analysis results."""
    # Calculate total reads
    totalreads = templatedf["uneditedcount"] + templatedf["editedcount"] +\
        templatedf["indelcount"]
    templatedf['totalreads'] = totalreads
    templatedf['totalreads'].replace(0, np.nan, inplace=True)

    # Calculate percentages
    templatedf['percentageedited'] = (templatedf["editedcount"] /
                                      templatedf['totalreads']) * 100
    templatedf['percentageunedited'] = (templatedf["uneditedcount"] /
                                        templatedf['totalreads']) * 100
    templatedf['percentageindel'] = (templatedf["indelcount"] /
                                     templatedf['totalreads']) * 100

    # Process intermediate positions
    templatedf['intermediatepositions'] = templatedf[
        'intermediatepositions'].apply(
        sorted_counter
    )

    return templatedf


def fix_column_names(templatedf):
    """Standardize column names to match expected output format."""
    column_mappings = {
        'RToverhanglength': 'RTToverhanglength',
        'RTlength': 'RTTlength',
        'RTseqoverhangrevcomp': 'RTTseqoverhangrevcomp',
        'RTTrevcomp_original': 'RTTrevcomp_lowercase',
        'RTrevcomp': 'RTTrevcomp',
        'RTmt': 'RTTmt',
        'RToverhangmt': 'RTToverhangmt',
        'RToverhangmatches': 'RTToverhangmatches',
        'Spacer-Sequence': 'Spacer',
        'protospacermt': 'spacermt'
    }

    for src, dst in column_mappings.items():
        templatedf[dst] = templatedf[src]

    return templatedf


def get_output_columns():
    """Return standard list of output columns in the correct order."""
    return [
        'Name', 'Original_Sequence', 'Edited-Sequences', 'Target-Strand',
        'Mutation_Type', 'Correction_Type', 'Correction_Length',
        'Editing_Position', 'PBSlength', 'RTToverhanglength', 'RTTlength',
        'EditedAllele', 'OriginalAllele', 'Spacer', 'PBSrevcomp',
        'RTTseqoverhangrevcomp', 'RTTrevcomp', 'RTTrevcomp_lowercase',
        'pegRNA', 'Editor_Variant', 'spacermt', 'extensionmt', 'RTTmt',
        'RTToverhangmt', 'PBSmt', 'original_base_mt', 'edited_base_mt',
        'original_base_mt_nan', 'edited_base_mt_nan', 'RTToverhangmatches',
        'wide_initial_target', 'wide_mutated_target',
        'wide_mutated_target_intermediate', 'protospacerlocation_only_initial',
        'PBSlocation', 'RT_initial_location', 'RT_mutated_location',
        'deepeditposition', 'deepeditposition_intermediate',
        'barcode', 'barcode_rev', 'extension', 'extension_15bp',
        'editedcount', 'uneditedcount', 'indelcount', 'intermediatecount',
        'intermediatepositions', 'nickindelcount', 'beforeflapindelcount',
        'skippedcount', 'totalreads', 'percentageedited', 'percentageunedited',
        'percentageindel', 'PRIDICT2_0_editing_Score_deep_HEK',
        'PRIDICT2_0_editing_Score_deep_K562'
    ]


def process_target_file(path, shortname, diseasedf, final_diseasedf,
                        templatenumpy, col_indices, filtered='', barcode_start=0,
                        barcode_length=6, min_quality=27.5, basesbeforenick=2,
                        baseafterflap=5, target_start_to_nick_distance=27,
                        nick_pos_in_wide_initial_target=26):
    """Process the target file to analyze editing outcomes."""
    readcounter = 0
    readcounttemp = 0
    # if reads are too short to quantify editing properly, skip them
    skipped_short_seq_count = 0
    skipped_short_ids = []
    start = time.time()
    target_file = f"{shortname}_target_5trim{filtered}.fastq.gz"

    print(f"Analyzing target file: {target_file}")
    with gzip.open(os.path.join(path, target_file), "rt") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            # Only analyze high quality reads
            qual_scores = seq_record.letter_annotations["phred_quality"]
            if qual_scores and (sum(qual_scores) /
                                len(qual_scores)) > min_quality:
                targetend = str(seq_record.seq)
                identifier = seq_record.id

                readcounter += 1
                readcounttemp += 1
                if readcounttemp == 100000:
                    print(f"Processed {readcounter} reads")
                    readcounttemp = 0

                # Skip if no match found
                if identifier not in final_diseasedf.index:
                    continue

                # Get variant index and parameters
                variantindex = diseasedf.at[identifier, 'match'][0]

                # Extract parameters for this variant
                RTTlength = int(templatenumpy[variantindex,
                                              col_indices['RTTlength']])
                correction_type = templatenumpy[variantindex,
                                                col_indices['correction_type']]
                mutation_type = templatenumpy[variantindex,
                                              col_indices['mutation_type']]
                correction_length = int(templatenumpy[
                    variantindex, col_indices['correction_length']])

                try:
                    # Extract sequences for analysis
                    sequences = extract_sequences(
                        targetend, variantindex, templatenumpy,
                        correction_type, RTTlength, correction_length,
                        target_start_to_nick_distance,
                        nick_pos_in_wide_initial_target,
                        basesbeforenick, baseafterflap,
                        col_indices['wide_initial'], col_indices['wide_mutated']
                    )

                    # Set up intermediate data if relevant
                    intermediate_data = None
                    if mutation_type == 'MultibpReplacement':
                        intermediate_data = {
                            'editlist': templatenumpy[
                                variantindex, col_indices['wide_intermediate']],
                            'positions': templatenumpy[
                                variantindex, col_indices['deepeditpos_inter']
                            ]
                        }

                    # Analyze the sequence
                    analyze_sequence(
                        sequences, templatenumpy, variantindex,
                        col_indices, mutation_type, intermediate_data
                    )
                except ValueError as ve:
                    skipped_short_seq_count += 1
                    skipped_short_ids.append((identifier, variantindex, str(ve)))
                    print(f"WARNING: Skipping read {identifier} (variant {variantindex}) - {ve}")
                    # Count skipped reads per variant
                    templatenumpy[variantindex, col_indices['skippedcount']] += 1
                    continue

    end = time.time()
    print(f'Analysis loop completed in {end-start:.2f} seconds')
    if skipped_short_seq_count:
        print(f"WARNING: Skipped {skipped_short_seq_count} reads due to short/unquantifiable extracted sequences.")
    return readcounter


def print_summary_stats(final_diseasedf, diseasedf, diseasedf_tevo, read_length=300):
    """Print summary statistics about the analysis."""

    # Validate input data - check for empty datasets
    if len(diseasedf) == 0:
        raise ValueError("No reads found for analysis")

    if len(diseasedf_tevo) == 0:
        raise ValueError("No reads found after extension filtering")

    print(f"Unique matched variants: {len(final_diseasedf['match'].unique())}")
    print(f"Total reads: {len(diseasedf)}")
    print(f"Successfully assigned reads: {len(final_diseasedf)}")
    if read_length >= 250:
        # For long reads, calculate extension recombination
        recombined_extension = len(diseasedf) - len(diseasedf_tevo)
        print(f"Recombined reads at extension: {recombined_extension}")
        percentage_extension_recombined = round((recombined_extension /
                                                len(diseasedf)) * 100, 2)
        print(f"Percentage extension recombined: {percentage_extension_recombined}%")
    else:
        # For short reads, note that this statistic is not available
        print("NOTE: Extension recombination statistics not available (short read mode)")

    # Scaffold recombination is always calculated
    recombined_scaffold = len(diseasedf_tevo) - len(final_diseasedf)
    print(f"Recombined reads at scaffold: {recombined_scaffold}")
    percentage_scaffold_recombined = round((recombined_scaffold /
                                           len(diseasedf_tevo)) * 100, 2)
    print(f"Percentage scaffold recombined: {percentage_scaffold_recombined}%")


def save_summary_stats(final_diseasedf,
                       diseasedf,
                       diseasedf_tevo,
                       output_dir,
                       date_prefix,
                       shortname,
                       read_length=300):
    """Save summary statistics to a CSV file."""

    # Validate input data - check for empty datasets
    if len(diseasedf) == 0:
        raise ValueError("No reads found for analysis")

    if len(diseasedf_tevo) == 0:
        raise ValueError("No reads found after extension filtering")

    unique_variants = len(final_diseasedf['match'].unique())
    total_reads = len(diseasedf)
    assigned_reads = len(final_diseasedf)

    # Create a dictionary of stats
    stats = {
        'sample': shortname,
        'unique_matched_variants': unique_variants,
        'total_reads': total_reads,
        'successfully_assigned_reads': assigned_reads,
        'read_length': read_length
    }

    # Add extension recombination statistics for long reads
    if read_length >= 250:
        recombined_extension = len(diseasedf) - len(diseasedf_tevo)
        percentage_extension_recombined = round((recombined_extension /
                                                 len(diseasedf)) * 100, 2)
        stats['recombined_reads_extension'] = recombined_extension
        stats['percentage_extension_recombined'] = percentage_extension_recombined

    # Scaffold recombination is always calculated
    recombined_scaffold = len(diseasedf_tevo) - len(final_diseasedf)
    percentage_scaffold_recombined = round((recombined_scaffold /
                                            len(diseasedf_tevo)) * 100, 2)

    stats['recombined_reads_scaffold'] = recombined_scaffold
    stats['percentage_scaffold_recombined'] = percentage_scaffold_recombined

    # Save as CSV
    stats_file = f"{output_dir}{date_prefix}_{shortname}_summary_stats.csv"
    pd.DataFrame([stats]).to_csv(stats_file, index=False)
    print(f"Summary statistics saved to {stats_file}")

    return stats


def main():
    """Main function to run the analysis."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Process self-targeting library data.')
    parser.add_argument(
        '--shortname', required=True, help='Short name of the sample')
    parser.add_argument(
        '--path', default='./Output/', help='Path to input files')
    parser.add_argument(
        '--template_file', required=True,
        help='Template file with library design')
    parser.add_argument(
        '--output_dir', default='Analysis/', help='Output directory')
    parser.add_argument(
        '--date_prefix', required=True, help='Date prefix for output files')
    parser.add_argument(
        '--library_name', required=True, help='Name of the library')
    parser.add_argument(
        '--barcode_start', type=int, default=0,
        help='Start position of the barcode in the read')
    parser.add_argument(
        '--barcode_length', type=int, required=True,
        help='Length of barcode in bases')
    parser.add_argument(
        '--min_quality', type=float, default=27.5,
        help='Minimum average Phred quality score for reads')
    parser.add_argument(
        '--basesbeforenick', type=int, default=2,
        help='Number of bases before the nick site to extract')
    parser.add_argument(
        '--baseafterflap', type=int, default=5,
        help='Number of bases after the flap to extract')
    parser.add_argument(
        '--target_start_to_nick_distance', type=int, default=27,
        help='Distance from target start to nick site')
    parser.add_argument(
        '--nick_pos_in_wide_initial_target', type=int, default=26,
        help='Position of nick within the wide_initial_target template '
             'sequence.')
    parser.add_argument(
        '--rtt_lookup_length', type=int, default=6,
        help='Number of leading bases from RTT sequence to use for lookup')
    parser.add_argument(
        '--read_length', type=int, default=300,
        help='Original read length from sequencer (150, 250, or 300bp)')
    args = parser.parse_args()

    # Set up paths and parameters
    path = args.path if args.path.endswith('/') else args.path + '/'
    output_dir = args.output_dir if args.output_dir.endswith('/') else \
        args.output_dir + '/'

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Use parameters from command-line arguments
    filtered = ''
    basesbeforenick = args.basesbeforenick
    baseafterflap = args.baseafterflap
    target_start_to_nick_distance = args.target_start_to_nick_distance
    nick_pos_in_wide_initial_target = \
        args.nick_pos_in_wide_initial_target

    # Process template file
    print(f"Loading template file: {args.template_file}")
    templatedf = preprocess_template(args.template_file)

    shortname = args.shortname

    print(f"\nProcessing sample: {shortname}")

    # Copy template to avoid modifications to original
    templatedf_sample = templatedf.copy()

    # Extract template sequences
    spacertemplate = templatedf_sample['Spacer-Sequence'].values.flatten()
    exttemplate = templatedf_sample['extension_15bp'].values.flatten()
    barcodetemplate = templatedf_sample['barcode_rev'].values.flatten()

    # Build lookup dictionaries with configurable RTT lookup length
    spacerlookup, rttlookup, barcodelookup = lookup(
        spacertemplate, exttemplate, barcodetemplate, args.rtt_lookup_length)
    del spacertemplate, exttemplate, barcodetemplate

    # Import and match files
    diseasedict, read_lengths = importfiles(shortname, path, spacerlookup, rttlookup,
                                            barcodelookup, filtered, args.barcode_start, args.barcode_length,
                                            args.rtt_lookup_length)

    # Validate read length configuration
    if read_lengths:
        median_length = np.median(read_lengths)
        if args.read_length >= 250 and median_length < 200:
            print("WARNING: Read length configuration mismatch!")
            print(f"Config specifies {args.read_length}bp, but median read length is {median_length}bp")
            raise ValueError(
                "Read length configuration mismatch. Please check your config and update read_length.")
        elif args.read_length <= 150 and median_length > 200:
            print("WARNING: Read length configuration mismatch!")
            print(f"Config specifies {args.read_length}bp, but median read length is {median_length}bp")
            print("Tevopreq filtering is disabled, but could be enabled. Consider updating read_length in config.")

    # Process matches
    diseasedf, diseasedf_tevo, final_diseasedf = process_matches(diseasedict, read_length=args.read_length)
    del diseasedict  # Free memory

    # Initialize template for analysis
    templatedf_sample = initialize_template_counts(templatedf_sample)

    # Get column indices
    col_indices = get_column_indices(templatedf_sample)

    # Convert to numpy for faster processing
    templatenumpy = templatedf_sample.to_numpy()

    # Process target file using parameters from config
    process_target_file(
        path, shortname, diseasedf, final_diseasedf, templatenumpy,
        col_indices, filtered, args.barcode_start, args.barcode_length,
        min_quality=args.min_quality, basesbeforenick=basesbeforenick,
        baseafterflap=baseafterflap,
        target_start_to_nick_distance=target_start_to_nick_distance,
        nick_pos_in_wide_initial_target=nick_pos_in_wide_initial_target
    )

    # Convert back to dataframe
    templatedf_sample = pd.DataFrame(
        data=templatenumpy,
        index=templatedf_sample.index.tolist(),
        columns=templatedf_sample.columns.tolist()
    )

    # Calculate statistics and fix column names
    templatedf_sample = calculate_statistics(templatedf_sample)
    templatedf_sample = fix_column_names(templatedf_sample)

    # Final processing
    final_diseasedf['match'] = final_diseasedf['match'].apply(lambda x: x[0])

    # Print summary statistics
    print_summary_stats(final_diseasedf, diseasedf, diseasedf_tevo,
                        read_length=args.read_length)

    # Save summary statistics
    save_summary_stats(final_diseasedf,
                       diseasedf,
                       diseasedf_tevo,
                       output_dir,
                       args.date_prefix,
                       shortname,
                       read_length=args.read_length)

    # Select output columns and save
    output_columns = get_output_columns()
    output_df = templatedf_sample[output_columns]

    # Save results
    output_file = f"{output_dir}{args.date_prefix}_{args.library_name}_{shortname}_analysisdf_focused.csv"
    output_df.to_csv(output_file)
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    main()
