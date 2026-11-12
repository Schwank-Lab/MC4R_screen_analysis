import os
import sys
import pandas as pd
from unittest.mock import patch
import io
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest
import tempfile
import gzip

# Import the script to test
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.self_targeting_analysis import (  # noqa: E402
    lookup, get_intermediate_strings, sorted_counter,
    get_column_indices, process_multibp_replacements, initialize_template_counts,
    preprocess_template, extract_sequences, importfiles, analyze_sequence,
    process_matches, calculate_statistics, fix_column_names, process_target_file,
    print_summary_stats, save_summary_stats, sample_read_lengths, main)


def test_lookup():
    """Test the lookup function creates correct lookup dictionaries."""
    # Create test input data
    spacertemplate = ["GACGTACGT", "ATCGATCG", "GACGTACGT"]  # Note: first and third are identical
    exttemplate = ["ACGTACGT", "TACGTACG", "GTACGTAC"]
    barcodetemplate = ["ACGT", "TACG", "CGTA"]

    # Call the function with default lookup length
    spacerlookup, rttlookup, barcodelookup = lookup(
        spacertemplate, exttemplate, barcodetemplate)

    # Check spacer lookup (with duplicate entry)
    assert "GACGTACGT" in spacerlookup
    assert len(spacerlookup["GACGTACGT"]) == 2
    assert 0 in spacerlookup["GACGTACGT"]
    assert 2 in spacerlookup["GACGTACGT"]

    # Check barcode lookup
    assert "ACGT" in barcodelookup
    assert barcodelookup["ACGT"] == [0]

    # Check RTT lookup (default length=6)
    assert "ACGTAC" in rttlookup  # First 6 characters of exttemplate[0]
    assert rttlookup["ACGTAC"] == [0]

    # Test with custom RTT lookup length
    _, rttlookup_custom, _ = lookup(
        spacertemplate, exttemplate, barcodetemplate, rtt_lookup_length=4)

    # Check custom length RTT lookup
    assert "ACGT" in rttlookup_custom  # First 4 characters
    assert rttlookup_custom["ACGT"] == [0]
    assert "TACG" in rttlookup_custom
    assert rttlookup_custom["TACG"] == [1]


def test_get_intermediate_strings():
    """Test the get_intermediate_strings function for finding intermediate sequences."""
    # Simple test case with single base change
    str1 = "ACGT"
    str2 = "ACTT"

    intermediates, positions = get_intermediate_strings(str1, str2)
    assert len(intermediates) == 0  # No intermediates for a single change
    assert len(positions) == 0

    # Test case with multiple changes
    str1 = "ACGT"
    str2 = "TCGA"

    intermediates, positions = get_intermediate_strings(str1, str2)
    assert len(intermediates) == 2  # 2^2 -2 = 2 intermediates

    # Check some expected intermediates are present
    assert "TCGT" in intermediates  # First position changed
    assert "ACGA" in intermediates  # Last position changed

    # Check positions
    for pos_list in positions:
        assert all(p in [0, 1, 2, 3] for p in pos_list)  # All positions should be 0-3

    # Check a more realistic sequence with fewer changes
    str1 = "ACGTACGTACGT"
    str2 = "ACGTTCGTACGT"  # Change at position 4

    intermediates, positions = get_intermediate_strings(str1, str2)
    assert len(intermediates) == 0  # Single change, no intermediates

    # Three changes
    str1 = "ACGTACGTACGT"
    str2 = "ACGTTCGCACGA"  # Changes at positions 4, 8, 11

    intermediates, positions = get_intermediate_strings(str1, str2)
    assert len(intermediates) == 6  # 2^(3) - 2 = 6 intermediates
    assert "ACGTTCGTACGT" in intermediates  # First change only
    assert "ACGTACGCACGT" in intermediates  # Second change only
    assert "ACGTACGTACGA" in intermediates  # Third change only

    # Check positions for one of the intermediates
    idx = intermediates.index("ACGTTCGTACGT")
    assert positions[idx] == [4]  # Position of the first change


def test_sorted_counter():
    """Test the sorted_counter function for converting lists to sorted counters."""
    # Test with a simple list
    test_list = ["A", "B", "A", "C", "B", "A", "D"]
    result = sorted_counter(test_list)

    # Check it's a dictionary
    assert isinstance(result, dict)

    # Check counts are correct
    assert result["A"] == 3
    assert result["B"] == 2
    assert result["C"] == 1
    assert result["D"] == 1

    # Check order is by frequency (descending)
    keys = list(result.keys())
    assert keys[0] == "A"  # Most frequent
    assert keys[1] == "B"
    assert keys[2] in ["C", "D"]  # Equal frequency, could be either order
    assert keys[3] in ["C", "D"]

    # Test with empty list
    assert sorted_counter([]) == {}

    # Test with list of numbers
    num_list = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]
    result = sorted_counter(num_list)
    assert result[5] == 3  # Most frequent
    assert result[1] == 2
    assert result[3] == 2

    # Check first item is the most frequent
    assert list(result.keys())[0] == 5


def test_get_column_indices():
    """Test the get_column_indices function."""
    # Create a sample template dataframe
    template_df = pd.DataFrame({
        "uneditedcount": [0, 0],
        "editedcount": [0, 0],
        "indelcount": [0, 0],
        "nickindelcount": [0, 0],
        "intermediatecount": [0, 0],
        "intermediatepositions": [[], []],
        "beforeflapindelcount": [0, 0],
        "skippedcount": [0, 0],
        "wide_initial_target": ["ACGT", "TGCA"],
        "wide_mutated_target": ["ACGa", "TGCt"],
        "wide_mutated_target_intermediate": [[], []],
        "deepeditposition_intermediate": [[], []],
        "RTlength": [10, 12],
        "Correction_Type": ["Replacement", "Replacement"],
        "Mutation_Type": ["Replacement", "Replacement"],
        "Correction_Length": [1, 1]
    })

    # Get column indices
    indices = get_column_indices(template_df)

    # Check that all expected keys are present
    expected_keys = [
        'uneditedcount', 'editedcount', 'indelcount', 'nickindelcount',
        'intermediatecount', 'intermediatepos', 'beforeflapindel',
        'wide_initial', 'wide_mutated', 'wide_intermediate', 'deepeditpos_inter',
        'RTTlength', 'correction_type', 'mutation_type', 'correction_length',
        'skippedcount'
    ]
    for key in expected_keys:
        assert key in indices

    # Check some specific indices
    assert indices['uneditedcount'] == 0
    assert indices['editedcount'] == 1
    assert indices['RTTlength'] == 12
    assert indices['correction_type'] == 13
    assert indices['mutation_type'] == 14
    assert indices['correction_length'] == 15


def test_process_multibp_replacements():
    """Test the process_multibp_replacements function."""
    # Create a sample template dataframe with TRUE multi-bp replacements
    template_df = pd.DataFrame({
        "wide_initial_target": ["ACGTACGT", "TGCATGCA", "AAGGCCTT"],
        "wide_mutated_target": ["tCGaACGT", "TGCgTGaA", "AAtagCTT"],  # Multiple changes in first two rows
        "Mutation_Type": ["MultibpReplacement", "MultibpReplacement", "SinglebpDeletion"],
    })

    # Process multibp replacements
    process_multibp_replacements(template_df)

    # Check that new columns were added
    assert "wide_mutated_target_intermediate" in template_df.columns
    assert "deepeditposition_intermediate" in template_df.columns

    # Check values for MultibpReplacement entries - should have intermediate states
    assert template_df.at[0, "wide_mutated_target_intermediate"] is not None
    assert len(template_df.at[0, "wide_mutated_target_intermediate"]) > 0  # Multiple changes = intermediates
    # Expect intermediates like "TCGTACGT" (first position changed only) and "ACGAACGT" (second position only)
    assert "tCGTACGT" in template_df.at[0, "wide_mutated_target_intermediate"]
    assert "ACGaACGT" in template_df.at[0, "wide_mutated_target_intermediate"]

    assert template_df.at[1, "wide_mutated_target_intermediate"] is not None
    assert len(template_df.at[1, "wide_mutated_target_intermediate"]) > 0  # Multiple changes = intermediates
    # Expect intermediates like "TGCGTGCA" (middle position only) and "TGCATGAA" (end position only)
    assert "TGCgTGCA" in template_df.at[1, "wide_mutated_target_intermediate"]
    assert "TGCATGaA" in template_df.at[1, "wide_mutated_target_intermediate"]

    # Position lists should also be populated
    assert len(template_df.at[0, "deepeditposition_intermediate"]) > 0
    assert len(template_df.at[1, "deepeditposition_intermediate"]) > 0

    # Check values for non-MultibpReplacement entries
    assert template_df.at[2, "wide_mutated_target_intermediate"] is None
    assert template_df.at[2, "deepeditposition_intermediate"] is None


def test_initialize_template_counts():
    """Test the initialize_template_counts function."""
    # Create a sample template dataframe
    template_df = pd.DataFrame({
        "wide_initial_target": ["ACGT", "TGCA"],
        "wide_mutated_target": ["ACGa", "TGCt"],
    })

    # Initialize counts
    result_df = initialize_template_counts(template_df)

    # Check that new columns were added with correct initial values
    assert "editedcount" in result_df.columns
    assert "uneditedcount" in result_df.columns
    assert "indelcount" in result_df.columns
    assert "nickindelcount" in result_df.columns
    assert "intermediatecount" in result_df.columns
    assert "intermediatepositions" in result_df.columns
    assert "beforeflapindelcount" in result_df.columns
    assert "totalreads" in result_df.columns

    # Check initial values
    assert result_df["editedcount"].tolist() == [0, 0]
    assert result_df["uneditedcount"].tolist() == [0, 0]
    assert result_df["indelcount"].tolist() == [0, 0]
    assert result_df["totalreads"].tolist() == [0, 0]

    # Check that intermediatepositions is initialized as empty lists
    assert isinstance(result_df["intermediatepositions"][0], list)
    assert len(result_df["intermediatepositions"][0]) == 0


def test_preprocess_template():
    """Test the preprocess_template function with a minimal template."""
    # Create a minimal CSV file in memory
    csv_data = io.StringIO(
        """uniquename,barcode,RTrevcomp,PBSrevcomp,Mutation_Type,Correction_Type,"""
        """Correction_Length,RTlength,wide_initial_target,wide_mutated_target
    sample1,ACGTAC,gTATGC,TTCATG,SinglebpReplacement,Replacement,1,15,ACGTACGTACGTACGT,ACGTACGTACGcACGT
    sample2,TGCATG,cCGTAC,GCATGC,MultibpReplacement,Replacement,2,20,TGCATGCATGCATGCA,TGCATGCgTaCATGCA
    """
    )

    # Mock the pd.read_csv to return our data
    with patch('pandas.read_csv', return_value=pd.read_csv(csv_data)):
        template_df = preprocess_template("dummy_path.csv")

    # Check that new columns were created
    assert "barcode_rev" in template_df.columns
    assert "RTTrevcomp_original" in template_df.columns
    assert "RTTrevcomp" in template_df.columns
    assert "extension" in template_df.columns
    assert "extension_15bp" in template_df.columns
    assert "Name" in template_df.columns

    # Check barcode reverse complement
    assert template_df.at[0, "barcode_rev"] != template_df.at[0, "barcode"]

    # Check RTT processing
    assert template_df.at[0, "RTTrevcomp_original"] == "gTATGC"
    assert template_df.at[0, "RTTrevcomp"] == "GTATGC"  # Should be uppercase

    # Check extension is concatenation of RTT and PBS
    assert template_df.at[0, "extension"] == "GTATGCTTCATG"

    # Check that intermediate columns were added (verifying process_multibp_replacements ran)
    assert "wide_mutated_target_intermediate" in template_df.columns
    assert "deepeditposition_intermediate" in template_df.columns

    # Check specific intermediate values for MultibpReplacement rows
    multibp_rows = template_df[template_df["Mutation_Type"] == "MultibpReplacement"]
    assert len(multibp_rows) > 0, "Expected at least one MultibpReplacement row"
    for _, row in multibp_rows.iterrows():
        assert row["wide_mutated_target_intermediate"] is not None


def test_extract_sequences():
    """Test the extract_sequences function."""
    # Define test inputs
    targetend = "ACGTACGTACGTACGTACGTCTCGATAG"
    variantindex = 0

    # Define sequence parameters
    rtt_length = 5
    correction_length = 2
    basebeforenick = 2
    baseafterflap = 5
    target_start_to_nick = 10
    nick_pos_in_template = 10

    # Calculate expected slice positions based on parameters
    nick_start = target_start_to_nick

    # Set up column indices for test
    col_wide_initial = 0
    col_wide_mutated = 1

    # Set up a mock numpy array for templatenumpy - just use 2 columns for simplicity
    templatenumpy = np.zeros((1, 2), dtype=object)
    templatenumpy[0, col_wide_initial] = "ACGTACGTACGTACGTACGTCTCGATAG"  # wide_initial_target
    templatenumpy[0, col_wide_mutated] = "ACGTACGTACtgACGTACGTCTCGATAG"  # wide_mutated_target (change at pos 10-11)

    # Test case 1: Replacement
    sequences = extract_sequences(
        targetend=targetend,
        variantindex=variantindex,
        templatenumpy=templatenumpy,
        correction_type="Replacement",
        RTTlength=rtt_length,
        correction_length=correction_length,
        target_start_to_nick_distance=target_start_to_nick,
        nick_pos_in_wide_initial_target=nick_pos_in_template,
        basesbeforenick=basebeforenick,
        baseafterflap=baseafterflap,
        wide_initial_idx=col_wide_initial,
        wide_mutated_idx=col_wide_mutated
    )

    # Check that all expected keys are present
    expected_keys = [
        "sequenceWT", "sequenceWT_beforeflap", "sequenceWT_nick",
        "sequenceMUT", "sequenceMUT_beforeflap", "sequenceMUT_nick",
        "controlWT", "controlWT_beforeflap", "controlWT_nick",
        "controlMUT", "controlMUT_beforeflap", "controlMUT_nick"
    ]
    for key in expected_keys:
        assert key in sequences

    # Check some specific values
    assert sequences["sequenceWT_nick"] == targetend[nick_start-basebeforenick:nick_start + 2]

    # For replacement:
    expected_length = rtt_length + baseafterflap
    assert sequences["sequenceWT"] == targetend[nick_start-basebeforenick:nick_start+expected_length]

    # controlWT should come from template at correct position
    expected_template_wt = templatenumpy[0, col_wide_initial][
          nick_pos_in_template-basebeforenick:nick_pos_in_template+expected_length]
    assert sequences["controlWT"] == expected_template_wt.upper()

    assert sequences["controlWT"] == sequences["sequenceWT"]

    wt_expected_length = basebeforenick + rtt_length + baseafterflap
    mut_expected_length = basebeforenick + rtt_length + baseafterflap
    assert len(sequences["controlWT"]) == wt_expected_length

    # Test case 2: Deletion

    # Set up a mock numpy array for templatenumpy - just use 2 columns for simplicity
    templatenumpy = np.zeros((1, 2), dtype=object)
    templatenumpy[0, col_wide_initial] = "ACGTACGTACGTACGTACGTCTCGATAG"  # wide_initial_target
    templatenumpy[0, col_wide_mutated] = "ACGTACGTACTACGTACGTCTCGATAG"  # wide_mutated_target (delete at position 10)

    rtt_length = 5
    correction_length = 1

    sequences = extract_sequences(
        targetend=targetend,
        variantindex=variantindex,
        templatenumpy=templatenumpy,
        correction_type="Deletion",
        RTTlength=rtt_length,
        correction_length=correction_length,
        target_start_to_nick_distance=target_start_to_nick,
        nick_pos_in_wide_initial_target=nick_pos_in_template,
        basesbeforenick=basebeforenick,
        baseafterflap=baseafterflap,
        wide_initial_idx=col_wide_initial,
        wide_mutated_idx=col_wide_mutated
    )

    # Check some specific values
    assert sequences["sequenceWT_nick"] == targetend[nick_start-basebeforenick:nick_start + 2]

    expected_length = rtt_length + correction_length + baseafterflap
    assert sequences["sequenceWT"] == targetend[nick_start-basebeforenick:nick_start+expected_length]

    expected_template_wt = templatenumpy[0, col_wide_initial][
          nick_pos_in_template-basebeforenick:nick_pos_in_template+expected_length]
    assert sequences["controlWT"] == expected_template_wt.upper()
    assert sequences["controlWT"] == sequences["sequenceWT"]

    wt_expected_length = basebeforenick + rtt_length + correction_length + baseafterflap
    mut_expected_length = basebeforenick + rtt_length + baseafterflap
    assert len(sequences["controlWT"]) == wt_expected_length
    assert len(sequences["controlMUT"]) == mut_expected_length
    assert sequences["controlWT"] == sequences["sequenceWT"]

    # Test case 3: Insertion
    templatenumpy = np.zeros((1, 2), dtype=object)
    templatenumpy[0, col_wide_initial] = "ACGTACGTACGTACGTACGTCTCGATAG"  # wide_initial_target
    templatenumpy[0, col_wide_mutated] = "ACGTACGTACAGTACGTACGTCTCGATAG"  # wide_mutated_target (delete at position 10)

    rtt_length = 5
    correction_length = 1

    sequences = extract_sequences(
        targetend=targetend,
        variantindex=variantindex,
        templatenumpy=templatenumpy,
        correction_type="Insertion",
        RTTlength=rtt_length,
        correction_length=correction_length,
        target_start_to_nick_distance=target_start_to_nick,
        nick_pos_in_wide_initial_target=nick_pos_in_template,
        basesbeforenick=basebeforenick,
        baseafterflap=baseafterflap,
        wide_initial_idx=col_wide_initial,
        wide_mutated_idx=col_wide_mutated
    )

    # Check some specific values
    assert sequences["sequenceWT_nick"] == targetend[nick_start-basebeforenick:nick_start + 2]

    expected_length = rtt_length - correction_length + baseafterflap
    assert sequences["sequenceWT"] == targetend[nick_start-basebeforenick:nick_start+expected_length]

    expected_template_wt = templatenumpy[0, col_wide_initial][
          nick_pos_in_template-basebeforenick:nick_pos_in_template+expected_length]
    assert sequences["controlWT"] == expected_template_wt.upper()
    assert sequences["controlWT"] == sequences["sequenceWT"]

    wt_expected_length = basebeforenick + rtt_length - correction_length + baseafterflap
    mut_expected_length = basebeforenick + rtt_length + baseafterflap
    assert len(sequences["controlWT"]) == wt_expected_length
    assert len(sequences["controlMUT"]) == mut_expected_length
    assert sequences["controlWT"] == sequences["sequenceWT"]

    # Test error cases
    # Test case: Unknown correction type
    with pytest.raises(ValueError, match="Unknown correction type"):
        extract_sequences(
            targetend=targetend,
            variantindex=variantindex,
            templatenumpy=templatenumpy,
            correction_type="InvalidType",  # Invalid correction type
            RTTlength=rtt_length,
            correction_length=correction_length,
            target_start_to_nick_distance=target_start_to_nick,
            nick_pos_in_wide_initial_target=nick_pos_in_template,
            basesbeforenick=basebeforenick,
            baseafterflap=baseafterflap,
            wide_initial_idx=col_wide_initial,
            wide_mutated_idx=col_wide_mutated
        )

    # Test case: Target sequence too short
    short_target = "ACGT"  # Very short target sequence
    with pytest.raises(ValueError, match="Target sequence is too short"):
        extract_sequences(
            targetend=short_target,
            variantindex=variantindex,
            templatenumpy=templatenumpy,
            correction_type="Replacement",
            RTTlength=rtt_length,
            correction_length=correction_length,
            target_start_to_nick_distance=target_start_to_nick,
            nick_pos_in_wide_initial_target=nick_pos_in_template,
            basesbeforenick=basebeforenick,
            baseafterflap=baseafterflap,
            wide_initial_idx=col_wide_initial,
            wide_mutated_idx=col_wide_mutated
        )

    # Test case: WT template sequence too short
    short_template = np.zeros((1, 2), dtype=object)
    short_template[0, col_wide_initial] = "ACGT"  # Short WT template
    short_template[0, col_wide_mutated] = "ACGTACGTACGTACGTACGTCTCGATAG"  # Normal MUT template

    with pytest.raises(ValueError, match="WT template sequence too short"):
        extract_sequences(
            targetend=targetend,
            variantindex=variantindex,
            templatenumpy=short_template,
            correction_type="Replacement",
            RTTlength=rtt_length,
            correction_length=correction_length,
            target_start_to_nick_distance=target_start_to_nick,
            nick_pos_in_wide_initial_target=nick_pos_in_template,
            basesbeforenick=basebeforenick,
            baseafterflap=baseafterflap,
            wide_initial_idx=col_wide_initial,
            wide_mutated_idx=col_wide_mutated
        )

    # Test case: MUT template sequence too short
    short_template = np.zeros((1, 2), dtype=object)
    short_template[0, col_wide_initial] = "ACGTACGTACGTACGTACGTCTCGATAG"  # Normal WT template
    short_template[0, col_wide_mutated] = "ACGT"  # Short MUT template

    with pytest.raises(ValueError, match="MUT template sequence too short"):
        extract_sequences(
            targetend=targetend,
            variantindex=variantindex,
            templatenumpy=short_template,
            correction_type="Replacement",
            RTTlength=rtt_length,
            correction_length=correction_length,
            target_start_to_nick_distance=target_start_to_nick,
            nick_pos_in_wide_initial_target=nick_pos_in_template,
            basesbeforenick=basebeforenick,
            baseafterflap=baseafterflap,
            wide_initial_idx=col_wide_initial,
            wide_mutated_idx=col_wide_mutated
        )


def test_process_matches():
    """Test the process_matches function."""
    # Create a mock diseasedict with various match scenarios
    diseasedict = {
        "read1": {
            "spacermatch": [0, 1],
            "barcodematch": [0, 2],
            "exttemplatematch": [0, 3],
            "tevopreq_in_extension": True
        },
        "read2": {
            "spacermatch": [1],
            "barcodematch": [1],
            "exttemplatematch": [1],
            "tevopreq_in_extension": True
        },
        "read3": {
            "spacermatch": [2],
            "barcodematch": [2],
            "exttemplatematch": [3],  # No common match with spacer and barcode for extension
            "tevopreq_in_extension": True
        },
        "read4": {
            "spacermatch": None,  # Missing spacermatch
            "barcodematch": [0],
            "exttemplatematch": [0],
            "tevopreq_in_extension": True
        },
        "read5": {
            "spacermatch": [1, 2],
            "barcodematch": [2, 3],
            "exttemplatematch": [2, 4],
            "tevopreq_in_extension": True
        },
        "read6": {
            "spacermatch": [5],
            "barcodematch": [6],  # No intersection with spacermatch
            "exttemplatematch": [7],
            "tevopreq_in_extension": True
        },
        "read7": {
            "spacermatch": [0, 1],
            "barcodematch": [0, 2],
            "exttemplatematch": [0, 3],
            "tevopreq_in_extension": False
        }
    }

    # Process the matches
    diseasedf, diseasedf_tevo, final_diseasedf = process_matches(diseasedict)

    # Check filtering of rows with missing values
    assert "read4" not in diseasedf.index, "Rows with None values should be filtered out"

    # Check that read7 is in the initial diseasedf but not in diseasedf_tevo
    assert "read7" in diseasedf.index, "read7 should be in diseasedf"
    assert "read7" not in diseasedf_tevo.index, "read7 should not be in diseasedf_tevo"
    assert "read7" not in final_diseasedf.index, "read7 should not be in final_diseasedf"

    # Check that match column is created correctly (intersection of spacermatch and barcodematch)
    assert diseasedf.loc["read1", "match"] == [0], "Intersection of spacermatch and barcodematch should be [0]"
    assert diseasedf.loc["read2", "match"] == [1], "Intersection should be [1]"
    assert diseasedf.loc["read3", "match"] == [2], "Intersection should be [2]"
    assert diseasedf.loc["read5", "match"] == [2], "Intersection should be [2]"
    assert len(diseasedf.loc["read6", "match"]) == 0, "Intersection should be empty"
    assert diseasedf.loc["read7", "match"] == [0], "Intersection for read7 should be [0]"

    # Check that matchwithextension column is created correctly (intersection of all three)
    assert diseasedf.loc["read1", "matchwithextension"] == [0], "Three-way intersection should be [0]"
    assert diseasedf.loc["read2", "matchwithextension"] == [1], "Three-way intersection should be [1]"
    assert not diseasedf.loc["read3", "matchwithextension"], "Three-way intersection should be empty"
    assert diseasedf.loc["read5", "matchwithextension"] == [2], "Three-way intersection should be [2]"
    assert not diseasedf.loc["read6", "matchwithextension"], "Three-way intersection should be empty"
    assert diseasedf.loc["read7", "matchwithextension"] == [0], "Three-way intersection for read7 should be [0]"

    # Check recombination flag
    assert not diseasedf.loc["read1", "recombination"], "read1 should not be a recombination"
    assert not diseasedf.loc["read2", "recombination"], "read2 should not be a recombination"
    assert diseasedf.loc["read3", "recombination"], "read3 should be a recombination"
    assert not diseasedf.loc["read5", "recombination"], "read5 should not be a recombination"
    assert diseasedf.loc["read6", "recombination"], "read6 should be a recombination"
    assert not diseasedf.loc["read7", "recombination"], "read7 should not be a recombination"

    # Check final_diseasedf filtering
    assert "read4" not in final_diseasedf.index, "Reads with missing values should be excluded"
    assert "read6" not in final_diseasedf.index, "Reads with empty match list should be excluded"
    assert "read7" not in final_diseasedf.index, "read7 should not be in final_diseasedf"

    # Verify reads with non-empty match lists are included
    assert "read1" in final_diseasedf.index, "read1 should be in final_diseasedf"
    assert "read2" in final_diseasedf.index, "read2 should be in final_diseasedf"
    assert "read3" in final_diseasedf.index, "read3 should be in final_diseasedf"
    assert "read5" in final_diseasedf.index, "read5 should be in final_diseasedf"

    # Check the contents of final_diseasedf
    assert "match" in final_diseasedf.columns, "final_diseasedf should have match column"
    assert "matchwithextension" in final_diseasedf.columns, "final_diseasedf should have matchwithextension column"
    assert list(final_diseasedf.columns) == [
        "match", "matchwithextension"], "final_diseasedf should only have these two columns"


def test_process_matches_with_read_length():
    """Test process_matches function with different read length settings."""
    # Create mock diseasedict for testing
    diseasedict = {
        "read1": {
            "spacermatch": [0],
            "barcodematch": [0],
            "exttemplatematch": [0],
            "tevopreq_in_extension": True
        },
        "read2": {
            "spacermatch": [1],
            "barcodematch": [1],
            "exttemplatematch": [1],
            "tevopreq_in_extension": False
        },
        "read3": {
            "spacermatch": [2],
            "barcodematch": [2],
            "exttemplatematch": [2],
            "tevopreq_in_extension": True
        }
    }

    # Test with short read length (150bp) - should keep all reads
    diseasedf, diseasedf_tevo, final_diseasedf = process_matches(diseasedict, read_length=150)

    # Check that diseasedf_tevo contains all reads regardless of tevopreq
    assert len(diseasedf) == 3, "Should have 3 total reads"
    assert len(diseasedf_tevo) == 3, "With short reads, should keep all reads"
    assert "read2" in diseasedf_tevo.index, "With short reads, should keep read2 despite missing tevopreq"

    # Test with long read length (250bp) - should filter reads without tevopreq
    diseasedf, diseasedf_tevo, final_diseasedf = process_matches(diseasedict, read_length=250)

    # Check that diseasedf_tevo filters reads without tevopreq
    assert len(diseasedf) == 3, "Should have 3 total reads"
    assert len(diseasedf_tevo) == 2, "With long reads, should filter out read2"
    assert "read2" not in diseasedf_tevo.index, "With long reads, should filter out read2"


def test_sample_read_lengths():
    """Test the read length sampling function."""
    # Create a temporary FASTQ file with known read lengths
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a simple FASTQ file
        test_file = os.path.join(tmpdir, "test.fastq.gz")

        with gzip.open(test_file, 'wt') as f:
            # Write 10 reads with varying lengths (100-109bp)
            for i in range(10):
                seq = "A" * (100 + i)
                f.write(f"@read{i}\n{seq}\n+\n{'F' * len(seq)}\n")

        # Test with default max_reads
        lengths = sample_read_lengths(tmpdir, "test.fastq.gz", max_reads=5)
        assert len(lengths) == 5, "Should get 5 read lengths"
        assert all(100 <= length <= 104 for length in lengths), "Lengths should be 100-104"

        # Test with larger max_reads
        lengths = sample_read_lengths(tmpdir, "test.fastq.gz", max_reads=20)
        assert len(lengths) == 10, "Should get all 10 read lengths"
        assert lengths == list(range(100, 110)), "Lengths should be 100-109"

        # Test with non-existent file
        lengths = sample_read_lengths(tmpdir, "nonexistent.fastq.gz")
        assert lengths == [], "Should return empty list for missing file"


def test_print_summary_stats_short_reads(capsys):
    """Test print_summary_stats with short read length."""
    # Create mock DataFrames
    final_diseasedf = pd.DataFrame({'match': [0, 1, 2]})
    diseasedf = pd.DataFrame(index=range(100))  # 100 total reads
    diseasedf_tevo = pd.DataFrame(index=range(90))  # 90 reads kept

    # Call with short read length
    print_summary_stats(final_diseasedf, diseasedf, diseasedf_tevo, read_length=150)

    # Capture printed output
    captured = capsys.readouterr()
    output = captured.out

    # Check for expected messages
    assert "NOTE: Extension recombination statistics not available (short read mode)" in output
    assert "Recombined reads at extension:" not in output
    assert "Percentage extension recombined:" not in output
    assert "Recombined reads at scaffold: 87" in output  # 90 - 3


def test_print_summary_stats_long_reads(capsys):
    """Test print_summary_stats with long read length."""
    # Create mock DataFrames
    final_diseasedf = pd.DataFrame({'match': [0, 1, 2]})
    diseasedf = pd.DataFrame(index=range(100))  # 100 total reads
    diseasedf_tevo = pd.DataFrame(index=range(90))  # 90 reads kept

    # Call with long read length
    print_summary_stats(final_diseasedf, diseasedf, diseasedf_tevo, read_length=250)

    # Capture printed output
    captured = capsys.readouterr()
    output = captured.out

    # Check for expected messages
    assert "NOTE: Extension recombination statistics not available" not in output
    assert "Recombined reads at extension: 10" in output  # 100 - 90
    assert "Percentage extension recombined: 10.0%" in output
    assert "Recombined reads at scaffold: 87" in output  # 90 - 3


def test_importfiles_shortreads(test_data_path):
    """Simplified test for importfiles function using real test files."""
    # Path to test FASTQ files
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")

    read_length = 150

    if read_length == 150:
        # Sample name matching your test files
        sample_name = "test_MC4R-lib-control1"
    elif read_length == 300:
        sample_name = "test_long_MC4R-lib-control1"

    # Create small lookup dictionaries for testing
    # These are example sequences - replace with actual data from your files
    spacer_sequences = ["GTCATTCAAACCATAGTTT", "GAAGGCCAGCAGATGGAGG", "AGAGGATTGCTGTCCTCCC"]
    rtt_sequences = ["GGAGGGAAACAAACCTAAA", "TCTCCGAGACCTCCT", "GATGGTCCCTGTACTGGGA"]
    barcode_sequences = ["GGACTT", "CGTCCA", "AGTTGA"]

    # Generate lookup dictionaries
    spacerlookup = {}
    for i, seq in enumerate(spacer_sequences):
        spacerlookup["G" + seq] = [i]

    rttlookup = {}
    for i, seq in enumerate(rtt_sequences):
        rttlookup[seq[:6]] = [i]

    barcodelookup = {}
    for i, seq in enumerate(barcode_sequences):
        barcodelookup[seq] = [i]

    # Call importfiles with the real files and the test lookup dictionaries
    # This will also expose the match functions
    result, read_length = importfiles(
        shortname=sample_name,
        path=fastq_path,
        spacerlookup=spacerlookup,
        rttlookup=rttlookup,
        barcodelookup=barcodelookup,
        barcode_start=0,
        barcode_length=6,
        rtt_lookup_length=6
    )

    # Basic verification that we got a dictionary back
    assert isinstance(result, dict), "Result should be a dictionary"
    print(f"importfiles returned a dictionary with {len(result)} entries")

    # Print a few sample entries to help you understand what's in the result
    samples = list(result.keys())[:3] if result else []  # Get up to 3 sample IDs
    for sample_id in samples:
        print(f"Sample entry for {sample_id}: {result[sample_id]}")

    # Verify the structure - each entry should have a diseasedict with lookup results
    for read_id, read_data in list(result.items())[:5]:  # Check first 5 entries
        assert isinstance(read_data, dict), f"Entry for {read_id} should be a dictionary"

        # Check for expected keys, but don't fail if some are missing (real data might not match all lookups)
        for key in ["spacermatch", "barcodematch", "exttemplatematch"]:
            if key in read_data:
                # Just verify the type, but allow None values for non-matches
                assert read_data[key] is None or isinstance(read_data[key], list), f"{key} should be None or a list"

    # Test the exposed match functions
    # Test match_spacer function
    diseasedict = {}
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[1]), id="test_read1"), diseasedict)
    assert "test_read1" in diseasedict
    assert "spacermatch" in diseasedict["test_read1"]
    assert diseasedict["test_read1"]["spacermatch"] == [1]

    # Test with a sequence that won't match
    importfiles.match_spacer(SeqRecord(Seq("NONMATCHING"), id="test_read2"), diseasedict)
    assert "test_read2" in diseasedict
    assert "spacermatch" in diseasedict["test_read2"]
    assert diseasedict["test_read2"]["spacermatch"] is None

    # Test updating an existing spacermatch
    diseasedict = {
        "test_read1": {"spacermatch": [99]}  # Existing entry with spacermatch
    }
    # This should overwrite the existing spacermatch value
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[1]), id="test_read1"), diseasedict)
    assert diseasedict["test_read1"]["spacermatch"] == [1], "match_spacer should update existing spacermatch"

    # Test match_barcode function
    # Test with a sequence that should match (only first 6 bases are checked)
    diseasedict = {"test_read1": {}}
    importfiles.match_barcode(SeqRecord(Seq(barcode_sequences[0] + "EXTRA"), id="test_read1"), diseasedict)
    assert "barcodematch" in diseasedict["test_read1"]
    assert diseasedict["test_read1"]["barcodematch"] == [0]

    # Test with a sequence that won't match
    diseasedict = {"test_read2": {}}
    importfiles.match_barcode(SeqRecord(Seq("NONMATCHING"), id="test_read2"), diseasedict)
    assert "barcodematch" in diseasedict["test_read2"]
    assert diseasedict["test_read2"]["barcodematch"] is None

    # Test match_extension function
    # Test with a sequence that should match
    diseasedict = {"test_read1": {}}
    importfiles.match_extension(SeqRecord(Seq(rtt_sequences[0] + "EXTRA"), id="test_read1"), diseasedict)
    assert "exttemplatematch" in diseasedict["test_read1"]
    assert diseasedict["test_read1"]["exttemplatematch"] == [0]

    # Test with a sequence that won't match
    diseasedict = {"test_read2": {}}
    importfiles.match_extension(SeqRecord(Seq("NONMATCHING"), id="test_read2"), diseasedict)
    assert "exttemplatematch" in diseasedict["test_read2"]
    assert diseasedict["test_read2"]["exttemplatematch"] is None

    # Test match_barcode with a non-existent identifier (should do nothing)
    diseasedict = {}  # Empty dictionary
    importfiles.match_barcode(SeqRecord(Seq(barcode_sequences[0]), id="new_read"), diseasedict)
    assert len(diseasedict) == 0, "match_barcode should not add new entries to diseasedict"

    # Test tevopreq_in_extension with sequence containing CGCGGTTCTATCTA
    diseasedict = {}  # Start with empty dictionary
    sequence_with_tevopreq = rtt_sequences[0] + "CGCGGTTCTATCTA" + "EXTRA"
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[0]), id="test_read_tevo"), diseasedict)
    importfiles.match_extension(SeqRecord(Seq(sequence_with_tevopreq), id="test_read_tevo"), diseasedict)
    assert "tevopreq_in_extension" in diseasedict["test_read_tevo"]
    assert diseasedict["test_read_tevo"]["tevopreq_in_extension"] is True, "Should detect tevopreq sequence"

    # Test tevopreq_in_extension with sequence NOT containing CGCGGTTCTATCTA
    diseasedict = {}  # Start with empty dictionary
    sequence_without_tevopreq = rtt_sequences[0] + "ACGTACGTACGTACGT"
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[0]), id="test_read_no_tevo"), diseasedict)
    importfiles.match_extension(SeqRecord(Seq(sequence_without_tevopreq), id="test_read_no_tevo"), diseasedict)
    assert "tevopreq_in_extension" in diseasedict["test_read_no_tevo"]
    assert diseasedict["test_read_no_tevo"]["tevopreq_in_extension"] is False, "Should not detect tevopreq sequence"

    # Test match_extension with a non-existent identifier (should do nothing)
    diseasedict = {}  # Empty dictionary
    importfiles.match_extension(SeqRecord(Seq(rtt_sequences[0]), id="new_read"), diseasedict)
    assert len(diseasedict) == 0, "match_extension should not add new entries to diseasedict"


def test_importfiles_longreads(test_data_path):
    """Simplified test for importfiles function using real test files."""
    # Path to test FASTQ files
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")

    read_length = 300

    if read_length == 150:
        # Sample name matching your test files
        shortname = "test_MC4R-lib-control1"
    elif read_length == 300:
        shortname = "test_long_MC4R-lib-control1"

    # Create small lookup dictionaries for testing
    # These are example sequences - replace with actual data from your files
    spacer_sequences = ["GTCATTCAAACCATAGTTT", "GAAGGCCAGCAGATGGAGG", "AGAGGATTGCTGTCCTCCC"]
    rtt_sequences = ["GGAGGGAAACAAACCTAAA", "TCTCCGAGACCTCCT", "GATGGTCCCTGTACTGGGA"]
    barcode_sequences = ["GGACTT", "CGTCCA", "AGTTGA"]

    # Generate lookup dictionaries
    spacerlookup = {}
    for i, seq in enumerate(spacer_sequences):
        spacerlookup["G" + seq] = [i]

    rttlookup = {}
    for i, seq in enumerate(rtt_sequences):
        rttlookup[seq[:6]] = [i]

    barcodelookup = {}
    for i, seq in enumerate(barcode_sequences):
        barcodelookup[seq] = [i]

    # Call importfiles with the real files and the test lookup dictionaries
    # This will also expose the match functions
    result, read_length = importfiles(
        shortname=shortname,
        path=fastq_path,
        spacerlookup=spacerlookup,
        rttlookup=rttlookup,
        barcodelookup=barcodelookup,
        barcode_start=0,
        barcode_length=6,
        rtt_lookup_length=6
    )

    # Basic verification that we got a dictionary back
    assert isinstance(result, dict), "Result should be a dictionary"
    print(f"importfiles returned a dictionary with {len(result)} entries")

    # Print a few sample entries to help you understand what's in the result
    samples = list(result.keys())[:3] if result else []  # Get up to 3 sample IDs
    for sample_id in samples:
        print(f"Sample entry for {sample_id}: {result[sample_id]}")

    # Verify the structure - each entry should have a diseasedict with lookup results
    for read_id, read_data in list(result.items())[:5]:  # Check first 5 entries
        assert isinstance(read_data, dict), f"Entry for {read_id} should be a dictionary"

        # Check for expected keys, but don't fail if some are missing (real data might not match all lookups)
        for key in ["spacermatch", "barcodematch", "exttemplatematch"]:
            if key in read_data:
                # Just verify the type, but allow None values for non-matches
                assert read_data[key] is None or isinstance(read_data[key], list), f"{key} should be None or a list"

    # Test the exposed match functions
    # Test match_spacer function
    diseasedict = {}
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[1]), id="test_read1"), diseasedict)
    assert "test_read1" in diseasedict
    assert "spacermatch" in diseasedict["test_read1"]
    assert diseasedict["test_read1"]["spacermatch"] == [1]

    # Test with a sequence that won't match
    importfiles.match_spacer(SeqRecord(Seq("NONMATCHING"), id="test_read2"), diseasedict)
    assert "test_read2" in diseasedict
    assert "spacermatch" in diseasedict["test_read2"]
    assert diseasedict["test_read2"]["spacermatch"] is None

    # Test updating an existing spacermatch
    diseasedict = {
        "test_read1": {"spacermatch": [99]}  # Existing entry with spacermatch
    }
    # This should overwrite the existing spacermatch value
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[1]), id="test_read1"), diseasedict)
    assert diseasedict["test_read1"]["spacermatch"] == [1], "match_spacer should update existing spacermatch"

    # Test match_barcode function
    # Test with a sequence that should match (only first 6 bases are checked)
    diseasedict = {"test_read1": {}}
    importfiles.match_barcode(SeqRecord(Seq(barcode_sequences[0] + "EXTRA"), id="test_read1"), diseasedict)
    assert "barcodematch" in diseasedict["test_read1"]
    assert diseasedict["test_read1"]["barcodematch"] == [0]

    # Test with a sequence that won't match
    diseasedict = {"test_read2": {}}
    importfiles.match_barcode(SeqRecord(Seq("NONMATCHING"), id="test_read2"), diseasedict)
    assert "barcodematch" in diseasedict["test_read2"]
    assert diseasedict["test_read2"]["barcodematch"] is None

    # Test match_extension function
    # Test with a sequence that should match
    diseasedict = {"test_read1": {}}
    importfiles.match_extension(SeqRecord(Seq(rtt_sequences[0] + "EXTRA"), id="test_read1"), diseasedict)
    assert "exttemplatematch" in diseasedict["test_read1"]
    assert diseasedict["test_read1"]["exttemplatematch"] == [0]

    # Test with a sequence that won't match
    diseasedict = {"test_read2": {}}
    importfiles.match_extension(SeqRecord(Seq("NONMATCHING"), id="test_read2"), diseasedict)
    assert "exttemplatematch" in diseasedict["test_read2"]
    assert diseasedict["test_read2"]["exttemplatematch"] is None

    # Test match_barcode with a non-existent identifier (should do nothing)
    diseasedict = {}  # Empty dictionary
    importfiles.match_barcode(SeqRecord(Seq(barcode_sequences[0]), id="new_read"), diseasedict)
    assert len(diseasedict) == 0, "match_barcode should not add new entries to diseasedict"

    # Test tevopreq_in_extension with sequence containing CGCGGTTCTATCTA
    diseasedict = {}  # Start with empty dictionary
    sequence_with_tevopreq = rtt_sequences[0] + "CGCGGTTCTATCTA" + "EXTRA"
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[0]), id="test_read_tevo"), diseasedict)
    importfiles.match_extension(SeqRecord(Seq(sequence_with_tevopreq), id="test_read_tevo"), diseasedict)
    assert "tevopreq_in_extension" in diseasedict["test_read_tevo"]
    assert diseasedict["test_read_tevo"]["tevopreq_in_extension"] is True, "Should detect tevopreq sequence"

    # Test tevopreq_in_extension with sequence NOT containing CGCGGTTCTATCTA
    diseasedict = {}  # Start with empty dictionary
    sequence_without_tevopreq = rtt_sequences[0] + "ACGTACGTACGTACGT"
    importfiles.match_spacer(SeqRecord(Seq(spacer_sequences[0]), id="test_read_no_tevo"), diseasedict)
    importfiles.match_extension(SeqRecord(Seq(sequence_without_tevopreq), id="test_read_no_tevo"), diseasedict)
    assert "tevopreq_in_extension" in diseasedict["test_read_no_tevo"]
    assert diseasedict["test_read_no_tevo"]["tevopreq_in_extension"] is False, "Should not detect tevopreq sequence"

    # Test match_extension with a non-existent identifier (should do nothing)
    diseasedict = {}  # Empty dictionary
    importfiles.match_extension(SeqRecord(Seq(rtt_sequences[0]), id="new_read"), diseasedict)
    assert len(diseasedict) == 0, "match_extension should not add new entries to diseasedict"


def test_analyze_sequence():
    """Test the analyze_sequence function."""
    # Create a mock templatenumpy array
    templatenumpy = np.zeros((1, 7), dtype=object)
    templatenumpy[0, 0] = 0  # uneditedcount
    templatenumpy[0, 1] = 0  # editedcount
    templatenumpy[0, 2] = 0  # indelcount
    templatenumpy[0, 3] = 0  # nickindelcount
    templatenumpy[0, 4] = 0  # intermediatecount
    templatenumpy[0, 5] = []  # intermediatepositions
    templatenumpy[0, 6] = 0  # beforeflapindel

    # Set up column indices
    col_indices = {
        'uneditedcount': 0,
        'editedcount': 1,
        'indelcount': 2,
        'nickindelcount': 3,
        'intermediatecount': 4,
        'intermediatepos': 5,
        'beforeflapindel': 6
    }

    # Test Case 1: WT sequence matches control WT
    sequences = {
        'sequenceWT': 'ACGT',
        'controlWT': 'ACGT',
        'sequenceMUT': 'ACGT',
        'controlMUT': 'ACGA',
        'sequenceWT_nick': 'AC',
        'controlWT_nick': 'AC',
        'sequenceMUT_nick': 'AC',
        'controlMUT_nick': 'AC',
        'sequenceWT_beforeflap': 'ACG',
        'controlWT_beforeflap': 'ACG',
        'sequenceMUT_beforeflap': 'ACG',
        'controlMUT_beforeflap': 'ACG'
    }

    analyze_sequence(sequences, templatenumpy, 0, col_indices, 'SinglebpReplacement')
    assert templatenumpy[0, 0] == 1, "uneditedcount should be incremented"
    assert templatenumpy[0, 1] == 0, "editedcount should not change"
    assert templatenumpy[0, 2] == 0, "indelcount should not change"

    # Reset templatenumpy for next test
    templatenumpy[0, 0] = 0

    # Test Case 2: MUT sequence matches control MUT
    sequences = {
        'sequenceWT': 'ACGA',  # Not matching controlWT
        'controlWT': 'ACGT',
        'sequenceMUT': 'ACGA',  # Matching controlMUT
        'controlMUT': 'ACGA',
        'sequenceWT_nick': 'AC',
        'controlWT_nick': 'AC',
        'sequenceMUT_nick': 'AC',
        'controlMUT_nick': 'AC',
        'sequenceWT_beforeflap': 'ACG',
        'controlWT_beforeflap': 'ACG',
        'sequenceMUT_beforeflap': 'ACG',
        'controlMUT_beforeflap': 'ACG'
    }

    analyze_sequence(sequences, templatenumpy, 0, col_indices, 'SinglebpReplacement')
    assert templatenumpy[0, 0] == 0, "uneditedcount should not change"
    assert templatenumpy[0, 1] == 1, "editedcount should be incremented"
    assert templatenumpy[0, 2] == 0, "indelcount should not change"

    # Reset templatenumpy for next test
    templatenumpy[0, 1] = 0

    # Test Case 3: No exact matches (indel case)
    sequences = {
        'sequenceWT': 'ACGC',  # Not matching controlWT
        'controlWT': 'ACGT',
        'sequenceMUT': 'ACGC',  # Not matching controlMUT
        'controlMUT': 'ACGA',
        'sequenceWT_nick': 'AC',
        'controlWT_nick': 'AC',
        'sequenceMUT_nick': 'AC',
        'controlMUT_nick': 'AC',
        'sequenceWT_beforeflap': 'ACG',
        'controlWT_beforeflap': 'ACG',
        'sequenceMUT_beforeflap': 'ACG',
        'controlMUT_beforeflap': 'ACG'
    }

    analyze_sequence(sequences, templatenumpy, 0, col_indices, 'SinglebpReplacement')
    assert templatenumpy[0, 0] == 0, "uneditedcount should not change"
    assert templatenumpy[0, 1] == 0, "editedcount should not change"
    assert templatenumpy[0, 2] == 1, "indelcount should be incremented"

    # Reset templatenumpy for next test
    templatenumpy[0, 2] = 0

    # Test Case 4: MultibpReplacement with intermediate matching
    sequences = {
        'sequenceWT': 'ACGA',  # Should match an intermediate state
        'controlWT': 'ACGT',
        'sequenceMUT': 'ACGA',  # Should match an intermediate state
        'controlMUT': 'ACTA',
        'sequenceWT_nick': 'AC',
        'controlWT_nick': 'AC',
        'sequenceMUT_nick': 'AC',
        'controlMUT_nick': 'AC',
        'sequenceWT_beforeflap': 'ACG',
        'controlWT_beforeflap': 'ACG',
        'sequenceMUT_beforeflap': 'ACG',
        'controlMUT_beforeflap': 'ACT'
    }

    # Create intermediate data with a match
    intermediate_data = {
        'editlist': ['ACGA', 'ACTT'],  # First element matches the sequence
        'positions': [[3], [2, 3]]  # Position info for each intermediate
    }

    analyze_sequence(sequences, templatenumpy, 0, col_indices, 'MultibpReplacement', intermediate_data)
    assert templatenumpy[0, 2] == 1, "indelcount should be incremented"
    assert templatenumpy[0, 4] == 1, "intermediatecount should be incremented"
    assert templatenumpy[0, 5] == [[3]], "intermediatepos should contain the position"

    # Reset templatenumpy for next test
    templatenumpy[0, 2] = 0
    templatenumpy[0, 4] = 0
    templatenumpy[0, 5] = []

    # Test Case 5: Indel around nick site
    sequences = {
        'sequenceWT': 'ACGTATGTACGTACGTACGT',  # Full sequence
        'controlWT': 'ACGTACGTACGTACGTACGT',
        'sequenceMUT': 'ACGTATGTACGTACGTACGT',
        'controlMUT': 'ACGTACGTACGTACGTACGT',
        'sequenceWT_nick': 'ACGTAT',  # 4bp region around nick site (modified)
        'controlWT_nick': 'ACGTAC',   # Original control sequence
        'sequenceMUT_nick': 'ACGTAT',  # Same modification in MUT
        'controlMUT_nick': 'ACGTAC',
        'sequenceWT_beforeflap': 'ACGTACGTACGT',
        'controlWT_beforeflap': 'ACGTACGTACGT',
        'sequenceMUT_beforeflap': 'ACGTACGTACGT',
        'controlMUT_beforeflap': 'ACGTACGTACGT'
    }

    analyze_sequence(sequences, templatenumpy, 0, col_indices, 'SinglebpReplacement')
    assert templatenumpy[0, 2] == 1, "indelcount should be incremented"
    assert templatenumpy[0, 3] == 1, "nickindelcount should be incremented"
    assert templatenumpy[0, 6] == 0, "beforeflapindel should not change"

    # Reset templatenumpy for next test
    templatenumpy[0, 2] = 0
    templatenumpy[0, 3] = 0

    # Test Case 6: Indel before flap
    sequences = {
        'sequenceWT': 'ACGTACGTACCTACGTACGT',  # Changed G to C in the middle
        'controlWT': 'ACGTACGTACGTACGTACGT',
        'sequenceMUT': 'ACGTACGTACCTACGTACGT',
        'controlMUT': 'ACGTACGTACGTACGTACGT',
        'sequenceWT_nick': 'ACGTAC',  # Matching controlWT_nick
        'controlWT_nick': 'ACGTAC',
        'sequenceMUT_nick': 'ACGTAC',  # Matching controlMUT_nick
        'controlMUT_nick': 'ACGTAC',
        'sequenceWT_beforeflap': 'ACGTACGTACCT',  # Different from controlWT_beforeflap (mutation here)
        'controlWT_beforeflap': 'ACGTACGTACGT',
        'sequenceMUT_beforeflap': 'ACGTACGTACCT',  # Different from controlMUT_beforeflap
        'controlMUT_beforeflap': 'ACGTACGTACGT'
    }

    analyze_sequence(sequences, templatenumpy, 0, col_indices, 'SinglebpReplacement')
    assert templatenumpy[0, 2] == 1, "indelcount should be incremented"
    assert templatenumpy[0, 3] == 0, "nickindelcount should not change"
    assert templatenumpy[0, 6] == 1, "beforeflapindel should be incremented"


def test_calculate_statistics():
    """Test the calculate_statistics function."""
    # Create a sample template dataframe with some test data
    data = {
        "uneditedcount": [10, 5, 0, 20],
        "editedcount": [5, 10, 0, 20],
        "indelcount": [5, 5, 0, 10],
        "intermediatepositions": [
            [1, 2, 1, 3],  # List with duplicates to test sorted_counter
            [4, 5, 4],
            [],
            [10, 11, 12]
        ]
    }
    templatedf = pd.DataFrame(data)

    # Call the function
    result_df = calculate_statistics(templatedf)

    # Test total reads calculation
    assert 'totalreads' in result_df.columns, "totalreads column should be added"
    assert result_df.loc[0, 'totalreads'] == 20, "Row 0 should have 20 total reads"
    assert result_df.loc[1, 'totalreads'] == 20, "Row 1 should have 20 total reads"
    assert pd.isna(result_df.loc[2, 'totalreads']), "Row 2 should have NaN for totalreads"
    assert result_df.loc[3, 'totalreads'] == 50, "Row 3 should have 50 total reads"

    # Test percentage calculations
    assert 'percentageedited' in result_df.columns, "percentageedited column should be added"
    assert 'percentageunedited' in result_df.columns, "percentageunedited column should be added"
    assert 'percentageindel' in result_df.columns, "percentageindel column should be added"

    # Check specific percentage values
    assert result_df.loc[0, 'percentageedited'] == 25, "Row 0 should have 25% edited"
    assert result_df.loc[0, 'percentageunedited'] == 50, "Row 0 should have 50% unedited"
    assert result_df.loc[0, 'percentageindel'] == 25, "Row 0 should have 25% indel"

    assert result_df.loc[1, 'percentageedited'] == 50, "Row 1 should have 50% edited"
    assert result_df.loc[1, 'percentageunedited'] == 25, "Row 1 should have 25% unedited"
    assert result_df.loc[1, 'percentageindel'] == 25, "Row 1 should have 25% indel"

    assert pd.isna(result_df.loc[2, 'percentageedited']), "Row 2 should have NaN percentages"

    assert result_df.loc[3, 'percentageedited'] == 40, "Row 3 should have 40% edited"
    assert result_df.loc[3, 'percentageunedited'] == 40, "Row 3 should have 40% unedited"
    assert result_df.loc[3, 'percentageindel'] == 20, "Row 3 should have 20% indel"

    # Test intermediate positions processing
    # Should be converted to sorted counters
    assert isinstance(result_df.loc[0, 'intermediatepositions'], dict), "Should be converted to a dictionary"

    # Check that count is correct
    assert result_df.loc[0, 'intermediatepositions'] == {1: 2, 3: 1, 2: 1}, "Counts should be correct"
    assert result_df.loc[1, 'intermediatepositions'] == {4: 2, 5: 1}, "Counts should be correct"
    assert result_df.loc[2, 'intermediatepositions'] == {}, "Empty list should give empty dict"
    assert result_df.loc[3, 'intermediatepositions'] == {10: 1, 11: 1, 12: 1}, "Counts should be correct"

    # Verify that the original dataframe has been modified (not a copy)
    assert id(templatedf) == id(result_df), "Function should modify and return the input dataframe"


def test_fix_column_names():
    """Test the fix_column_names function."""
    # Create a sample template dataframe with columns to be renamed
    data = {
        'RToverhanglength': [10, 20],
        'RTlength': [30, 40],
        'RTseqoverhangrevcomp': ['ACGT', 'TGCA'],
        'RTTrevcomp_original': ['gTATGC', 'aCGTAC'],
        'RTrevcomp': ['GTATGC', 'ACGTAC'],
        'RTmt': [0.5, 0.7],
        'RToverhangmt': [0.8, 0.9],
        'RToverhangmatches': [3, 4],
        'Spacer-Sequence': ['GTATCGA', 'GCTAGCT'],
        'protospacermt': [0.2, 0.3],
        'other_column': ['value1', 'value2']  # Column that shouldn't be renamed
    }
    templatedf = pd.DataFrame(data)

    # Call the function
    result_df = fix_column_names(templatedf)

    # Test that all new columns were added with correct values
    assert 'RTToverhanglength' in result_df.columns
    assert 'RTTlength' in result_df.columns
    assert 'RTTseqoverhangrevcomp' in result_df.columns
    assert 'RTTrevcomp_lowercase' in result_df.columns
    assert 'RTTrevcomp' in result_df.columns
    assert 'RTTmt' in result_df.columns
    assert 'RTToverhangmt' in result_df.columns
    assert 'RTToverhangmatches' in result_df.columns
    assert 'Spacer' in result_df.columns
    assert 'spacermt' in result_df.columns

    # Test that original columns still exist
    assert 'RToverhanglength' in result_df.columns
    assert 'RTlength' in result_df.columns
    assert 'RTseqoverhangrevcomp' in result_df.columns
    assert 'RTTrevcomp_original' in result_df.columns
    assert 'RTrevcomp' in result_df.columns
    assert 'RTmt' in result_df.columns
    assert 'RToverhangmt' in result_df.columns
    assert 'RToverhangmatches' in result_df.columns
    assert 'Spacer-Sequence' in result_df.columns
    assert 'protospacermt' in result_df.columns
    assert 'other_column' in result_df.columns

    # Test specific values to ensure the mapping was correct
    assert result_df.loc[0, 'RTToverhanglength'] == 10
    assert result_df.loc[1, 'RTTlength'] == 40
    assert result_df.loc[0, 'RTTseqoverhangrevcomp'] == 'ACGT'
    assert result_df.loc[1, 'RTTrevcomp_lowercase'] == 'aCGTAC'
    assert result_df.loc[0, 'RTTrevcomp'] == 'GTATGC'
    assert result_df.loc[1, 'RTTmt'] == 0.7
    assert result_df.loc[0, 'RTToverhangmt'] == 0.8
    assert result_df.loc[1, 'RTToverhangmatches'] == 4
    assert result_df.loc[0, 'Spacer'] == 'GTATCGA'
    assert result_df.loc[1, 'spacermt'] == 0.3

    # Test that the original dataframe is modified and returned
    assert id(templatedf) == id(result_df)


def test_print_summary_stats(capsys):
    """Test print_summary_stats prints the correct information."""
    # Create mock DataFrames
    final_diseasedf = pd.DataFrame({
        'match': [0, 1, 2, 0]  # 3 unique matches
    })
    diseasedf = pd.DataFrame(index=range(100))  # 100 total reads
    diseasedf_tevo = pd.DataFrame(index=range(80))  # 80 reads after extension filtering

    # Call the function
    print_summary_stats(final_diseasedf, diseasedf, diseasedf_tevo)

    # Capture printed output
    captured = capsys.readouterr()
    output = captured.out

    # Check that the output contains expected information
    assert "Unique matched variants: 3" in output
    assert "Total reads: 100" in output
    assert "Successfully assigned reads: 4" in output
    assert "Recombined reads at extension: 20" in output  # 100 - 80
    assert "Percentage extension recombined: 20.0%" in output
    assert "Recombined reads at scaffold: 76" in output  # 80 - 4
    assert "Percentage scaffold recombined: 95.0%" in output  # 76/80 * 100


def test_print_summary_stats_no_reads():
    """Test that print_summary_stats raises an exception when there are no reads."""
    # Create empty DataFrames
    final_diseasedf = pd.DataFrame({'match': []})
    diseasedf = pd.DataFrame()
    diseasedf_tevo = pd.DataFrame()

    # Call the function with empty DataFrames
    with pytest.raises(ValueError, match="No reads found for analysis"):
        print_summary_stats(final_diseasedf, diseasedf, diseasedf_tevo)


def test_print_summary_stats_no_tevo_reads():
    """Test that print_summary_stats raises an exception when there are no tevo reads."""
    # Create DataFrames with diseasedf but empty diseasedf_tevo
    final_diseasedf = pd.DataFrame({'match': []})
    diseasedf = pd.DataFrame(index=range(10))  # 10 reads in diseasedf
    diseasedf_tevo = pd.DataFrame()  # No reads in diseasedf_tevo

    # Call the function
    with pytest.raises(ValueError, match="No reads found after extension filtering"):
        print_summary_stats(final_diseasedf, diseasedf, diseasedf_tevo)


def test_save_summary_stats():
    """Test save_summary_stats creates a CSV file with correct statistics."""
    # Create mock DataFrames
    final_diseasedf = pd.DataFrame({
        'match': [0, 1, 2, 0]  # 3 unique matches
    })
    diseasedf = pd.DataFrame(index=range(100))  # 100 total reads
    diseasedf_tevo = pd.DataFrame(index=range(80))  # 80 reads after extension filtering

    # Create a temporary directory for output
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = tmpdir + '/'
        date_prefix = "20250101"
        shortname = "test_sample"

        # Call the function
        stats = save_summary_stats(
            final_diseasedf,
            diseasedf,
            diseasedf_tevo,
            output_dir,
            date_prefix,
            shortname
        )

        # Check that the file was created
        stats_file = f"{output_dir}{date_prefix}_{shortname}_summary_stats.csv"
        assert os.path.exists(stats_file), f"Stats file not created: {stats_file}"

        # Load the file and check its contents
        stats_df = pd.read_csv(stats_file)
        assert len(stats_df) == 1, "Expected a single row in the stats file"

        # Check specific values
        assert stats_df.iloc[0]['sample'] == "test_sample"
        assert stats_df.iloc[0]['unique_matched_variants'] == 3
        assert stats_df.iloc[0]['total_reads'] == 100
        assert stats_df.iloc[0]['successfully_assigned_reads'] == 4
        assert stats_df.iloc[0]['recombined_reads_extension'] == 20
        assert stats_df.iloc[0]['percentage_extension_recombined'] == 20.0
        assert stats_df.iloc[0]['recombined_reads_scaffold'] == 76
        assert stats_df.iloc[0]['percentage_scaffold_recombined'] == 95.0

        # Check that the returned stats match the file contents
        assert stats['unique_matched_variants'] == 3
        assert stats['total_reads'] == 100
        assert stats['successfully_assigned_reads'] == 4


def test_save_summary_stats_no_reads():
    """Test save_summary_stats raises an exception when there are no reads."""
    # Create empty DataFrames
    final_diseasedf = pd.DataFrame({'match': []})
    diseasedf = pd.DataFrame()
    diseasedf_tevo = pd.DataFrame()

    # Call the function with empty DataFrames
    with pytest.raises(ValueError, match="No reads found for analysis"):
        save_summary_stats(
            final_diseasedf,
            diseasedf,
            diseasedf_tevo,
            "output/",
            "20250101",
            "test_sample"
        )


def test_save_summary_stats_no_tevo_reads():
    """Test save_summary_stats raises an exception when there are no tevo reads."""
    # Create DataFrames with reads in diseasedf but not in diseasedf_tevo
    final_diseasedf = pd.DataFrame(columns=['match'])
    diseasedf = pd.DataFrame(index=range(10))  # 10 reads
    diseasedf_tevo = pd.DataFrame()  # 0 reads after extension filtering

    # Call the function
    with pytest.raises(ValueError, match="No reads found after extension filtering"):
        save_summary_stats(
            final_diseasedf,
            diseasedf,
            diseasedf_tevo,
            "output/",
            "20250101",
            "test_sample"
        )


def test_process_target_file_increments_skippedcount_realistic(tmp_path):
    """Integration-style test: a too-short read triggers skippedcount increment."""
    template_df = pd.DataFrame({
        "uneditedcount": [0],
        "editedcount": [0],
        "indelcount": [0],
        "nickindelcount": [0],
        "intermediatecount": [0],
        "intermediatepositions": [[]],
        "beforeflapindelcount": [0],
        "skippedcount": [0],
        "wide_initial_target": ["A" * 50],  # template expects long sequence
        "wide_mutated_target": ["A" * 50],
        "wide_mutated_target_intermediate": [[]],
        "deepeditposition_intermediate": [[]],
        "RTlength": [20],
        "Correction_Type": ["Replacement"],
        "Mutation_Type": ["Replacement"],
        "Correction_Length": [10],
    })
    templatenumpy = template_df.to_numpy()
    col_indices = get_column_indices(template_df)

    diseasedf = pd.DataFrame({"match": [(0,)]}, index=["read1"])
    final_diseasedf = diseasedf.copy()

    fastq_data = "@read1\nACG\n+\nIII\n"  # only 3 nt long
    target_file = tmp_path / "shorttest_target_5trim.fastq.gz"
    with gzip.open(target_file, "wt") as f:
        f.write(fastq_data)

    process_target_file(
        path=tmp_path,
        shortname="shorttest",
        diseasedf=diseasedf,
        final_diseasedf=final_diseasedf,
        templatenumpy=templatenumpy,
        col_indices=col_indices,
        filtered="",  # still needed
        barcode_start=0,
        barcode_length=6,
        min_quality=0,
        basesbeforenick=2,
        baseafterflap=5,
        target_start_to_nick_distance=27,
        nick_pos_in_wide_initial_target=26,
    )
    # representing "totalreads":
    assert (templatenumpy[0, col_indices["uneditedcount"]] + templatenumpy[0, col_indices["editedcount"]] +
            templatenumpy[0, col_indices["indelcount"]] == 0)
    assert templatenumpy[0, col_indices["skippedcount"]] == 1


def test_process_target_file_integrated_short(test_data_path, template_file, sample_config):
    """Test process_target_file in an integrated manner with real test files."""
    # Setup paths
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")

    read_length = 150

    if read_length == 150:
        # Sample name matching your test files
        shortname = "test_MC4R-lib-control1"
    elif read_length == 300:
        shortname = "test_long_MC4R-lib-control1"

    # Run the steps leading up to process_target_file
    templatedf = preprocess_template(template_file)

    # Extract template sequences for lookup
    spacertemplate = templatedf['Spacer-Sequence'].values.flatten()
    exttemplate = templatedf['extension_15bp'].values.flatten()
    barcodetemplate = templatedf['barcode_rev'].values.flatten()

    # Build lookup dictionaries
    spacerlookup, rttlookup, barcodelookup = lookup(
        spacertemplate, exttemplate, barcodetemplate,
        sample_config["rtt_lookup_length"]
    )

    # Import and match files
    diseasedict, read_lengths = importfiles(
        shortname, fastq_path, spacerlookup, rttlookup, barcodelookup,
        filtered='', barcode_start=sample_config["barcode_start"],
        barcode_length=sample_config["barcode_length"],
        rtt_lookup_length=sample_config["rtt_lookup_length"]
    )

    # Process matches
    diseasedf, diseasedf_tevo, final_diseasedf = process_matches(diseasedict, read_length)

    # Initialize template for analysis
    templatedf = initialize_template_counts(templatedf)

    # Get column indices
    col_indices = get_column_indices(templatedf)

    # Convert to numpy
    templatenumpy = templatedf.to_numpy()

    # Now call process_target_file
    read_count = process_target_file(
        fastq_path, shortname, diseasedf, final_diseasedf,
        templatenumpy, col_indices, filtered='',
        barcode_start=sample_config["barcode_start"],
        barcode_length=sample_config["barcode_length"],
        min_quality=sample_config["min_quality"],
        basesbeforenick=sample_config["basesbeforenick"],
        baseafterflap=sample_config["baseafterflap"],
        target_start_to_nick_distance=sample_config["target_start_to_nick_distance"],
        nick_pos_in_wide_initial_target=sample_config["nick_pos_in_wide_initial_target"]
    )

    # Check that reads were processed
    assert read_count > 0, "Expected to process some reads"

    # Convert back to dataframe for easier checking
    templatedf_result = pd.DataFrame(
        data=templatenumpy,
        index=templatedf.index.tolist(),
        columns=templatedf.columns.tolist()
    )

    # Check that some counts were updated - these depend on the specific test data
    total_counts = (templatedf_result["uneditedcount"].sum() + templatedf_result["editedcount"].sum()
                    + templatedf_result["indelcount"].sum())
    assert total_counts > 0, "Expected some reads to be counted"

    # Check that at least one variant has counts
    variants_with_counts = (templatedf_result["uneditedcount"] > 0).sum()
    + (templatedf_result["editedcount"] > 0).sum() + (templatedf_result["indelcount"] > 0).sum()
    assert variants_with_counts > 0, "Expected at least one variant to have counts"


def test_process_target_file_integrated_long(test_data_path, template_file, sample_config):
    """Test process_target_file in an integrated manner with real test files."""
    # Setup paths
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")

    read_length = 300

    if read_length == 150:
        # Sample name matching your test files
        shortname = "test_MC4R-lib-control1"
    elif read_length == 300:
        shortname = "test_long_MC4R-lib-control1"

    # Run the steps leading up to process_target_file
    templatedf = preprocess_template(template_file)

    # Extract template sequences for lookup
    spacertemplate = templatedf['Spacer-Sequence'].values.flatten()
    exttemplate = templatedf['extension_15bp'].values.flatten()
    barcodetemplate = templatedf['barcode_rev'].values.flatten()

    # Build lookup dictionaries
    spacerlookup, rttlookup, barcodelookup = lookup(
        spacertemplate, exttemplate, barcodetemplate,
        sample_config["rtt_lookup_length"]
    )

    # Import and match files
    diseasedict, read_lengths = importfiles(
        shortname, fastq_path, spacerlookup, rttlookup, barcodelookup,
        filtered='', barcode_start=sample_config["barcode_start"],
        barcode_length=sample_config["barcode_length"],
        rtt_lookup_length=sample_config["rtt_lookup_length"]
    )

    # Process matches
    diseasedf, diseasedf_tevo, final_diseasedf = process_matches(diseasedict, read_length)

    # Initialize template for analysis
    templatedf = initialize_template_counts(templatedf)

    # Get column indices
    col_indices = get_column_indices(templatedf)

    # Convert to numpy
    templatenumpy = templatedf.to_numpy()

    # Now call process_target_file
    read_count = process_target_file(
        fastq_path, shortname, diseasedf, final_diseasedf,
        templatenumpy, col_indices, filtered='',
        barcode_start=sample_config["barcode_start"],
        barcode_length=sample_config["barcode_length"],
        min_quality=sample_config["min_quality"],
        basesbeforenick=sample_config["basesbeforenick"],
        baseafterflap=sample_config["baseafterflap"],
        target_start_to_nick_distance=sample_config["target_start_to_nick_distance"],
        nick_pos_in_wide_initial_target=sample_config["nick_pos_in_wide_initial_target"]
    )

    # Check that reads were processed
    assert read_count > 0, "Expected to process some reads"

    # Convert back to dataframe for easier checking
    templatedf_result = pd.DataFrame(
        data=templatenumpy,
        index=templatedf.index.tolist(),
        columns=templatedf.columns.tolist()
    )

    # Check that some counts were updated - these depend on the specific test data
    total_counts = (templatedf_result["uneditedcount"].sum() + templatedf_result["editedcount"].sum()
                    + templatedf_result["indelcount"].sum())
    assert total_counts > 0, "Expected some reads to be counted"

    # Check that at least one variant has counts
    variants_with_counts = (templatedf_result["uneditedcount"] > 0).sum()
    + (templatedf_result["editedcount"] > 0).sum() + (templatedf_result["indelcount"] > 0).sum()
    assert variants_with_counts > 0, "Expected at least one variant to have counts"


def test_process_target_file_with_mocks():
    """Test process_target_file with mocked inputs to verify unedited sequence logic."""
    # Create mocked dataframes
    diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    final_diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]],
        'matchwithextension': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    # Create a simple template array for testing
    templatenumpy = np.zeros((3, 15), dtype=object)

    # Set up template values for the first variant
    templatenumpy[0, 11] = 10  # RTTlength
    templatenumpy[0, 12] = "Replacement"  # correction_type
    templatenumpy[0, 13] = "SinglebpReplacement"  # mutation_type
    templatenumpy[0, 14] = 1  # correction_length
    templatenumpy[0, 8] = "ACGTACGTACGTACGTACGTACGTACGT"  # wide_initial_target
    templatenumpy[0, 9] = "ACGTACGTACGTACGTACGTACaTACGT"  # wide_mutated_target

    # Column indices for lookups
    col_indices = {
        'uneditedcount': 0,
        'editedcount': 1,
        'indelcount': 2,
        'nickindelcount': 3,
        'intermediatecount': 4,
        'intermediatepos': 5,
        'beforeflapindel': 6,
        'wide_initial': 8,
        'wide_mutated': 9,
        'wide_intermediate': 10,
        'deepeditpos_inter': 7,
        'RTTlength': 11,
        'correction_type': 12,
        'mutation_type': 13,
        'correction_length': 14
    }

    # Create a mock sequence record that should match variant 0 (unedited)
    seq_record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGT"), id="id1")
    seq_record.letter_annotations["phred_quality"] = [40] * 28  # High quality

    # Mock the file opening and parsing
    with patch('gzip.open'), \
         patch('Bio.SeqIO.parse', return_value=[seq_record]):

        # Call the function with our mocked data
        read_count = process_target_file(
            "dummy_path", "test_sample", diseasedf, final_diseasedf,
            templatenumpy, col_indices, filtered='', barcode_start=0,
            barcode_length=6, min_quality=20, basesbeforenick=2,
            baseafterflap=5, target_start_to_nick_distance=10,
            nick_pos_in_wide_initial_target=10
        )

    # Check that one read was processed
    assert read_count == 1, "Expected to process 1 read"

    # Check that the uneditedcount was incremented for variant 0
    # Since our mock sequence matches the wide_initial_target exactly
    assert templatenumpy[0, 0] == 1, "Expected uneditedcount to be incremented"
    assert templatenumpy[0, 1] == 0, "Expected editedcount to remain 0"
    assert templatenumpy[0, 2] == 0, "Expected indelcount to remain 0"


def test_process_target_file_mocked_with_edited():
    """Test process_target_file with a mock that matches an edited sequence."""
    # Create mocked dataframes
    diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    final_diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]],
        'matchwithextension': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    # Create a simple template array for testing
    templatenumpy = np.zeros((3, 15), dtype=object)

    # Set up template values for the first variant
    templatenumpy[0, 11] = 10  # RTTlength
    templatenumpy[0, 12] = "Replacement"  # correction_type
    templatenumpy[0, 13] = "SinglebpReplacement"  # mutation_type
    templatenumpy[0, 14] = 1  # correction_length
    templatenumpy[0, 8] = "ACGTACGTACGTACGTACGTACGTACGT"  # wide_initial_target
    templatenumpy[0, 9] = "ACGTACGTACGTACGTACGTACaTACGT"  # wide_mutated_target

    # Column indices for lookups
    col_indices = {
        'uneditedcount': 0,
        'editedcount': 1,
        'indelcount': 2,
        'nickindelcount': 3,
        'intermediatecount': 4,
        'intermediatepos': 5,
        'beforeflapindel': 6,
        'wide_initial': 8,
        'wide_mutated': 9,
        'wide_intermediate': 10,
        'deepeditpos_inter': 7,
        'RTTlength': 11,
        'correction_type': 12,
        'mutation_type': 13,
        'correction_length': 14
    }

    # Create a mock sequence record that should match the mutated variant
    seq_record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACATACGT"), id="id1")  # This matches the mutated sequence
    seq_record.letter_annotations["phred_quality"] = [40] * 28  # High quality

    # Mock the file opening and parsing
    with patch('gzip.open'), \
         patch('Bio.SeqIO.parse', return_value=[seq_record]):

        # Call the function with our mocked data
        read_count = process_target_file(
            "dummy_path", "test_sample", diseasedf, final_diseasedf,
            templatenumpy, col_indices, filtered='', barcode_start=0,
            barcode_length=6, min_quality=20, basesbeforenick=2,
            baseafterflap=5, target_start_to_nick_distance=10,
            nick_pos_in_wide_initial_target=10
        )

    # Check that one read was processed
    assert read_count == 1, "Expected to process 1 read"

    # Check that the editedcount was incremented for variant 0
    # Since our mock sequence matches the wide_mutated_target exactly
    assert templatenumpy[0, 0] == 0, "Expected uneditedcount to remain 0"
    assert templatenumpy[0, 1] == 1, "Expected editedcount to be incremented"
    assert templatenumpy[0, 2] == 0, "Expected indelcount to remain 0"


def test_process_target_file_mocked_with_indel():
    """Test process_target_file with a mock that contains an indel."""
    # Create mocked dataframes
    diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    final_diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]],
        'matchwithextension': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    # Create a simple template array for testing
    templatenumpy = np.zeros((3, 15), dtype=object)

    # Set up template values for the first variant
    templatenumpy[0, 11] = 10  # RTTlength
    templatenumpy[0, 12] = "Replacement"  # correction_type
    templatenumpy[0, 13] = "SinglebpReplacement"  # mutation_type
    templatenumpy[0, 14] = 1  # correction_length
    templatenumpy[0, 8] = "ACGTACGTACGTACGTACGTACGTACGT"  # wide_initial_target
    templatenumpy[0, 9] = "ACGTACGTACGTACGTACGTACaTACGT"  # wide_mutated_target

    # Column indices for lookups
    col_indices = {
        'uneditedcount': 0,
        'editedcount': 1,
        'indelcount': 2,
        'nickindelcount': 3,
        'intermediatecount': 4,
        'intermediatepos': 5,
        'beforeflapindel': 6,
        'wide_initial': 8,
        'wide_mutated': 9,
        'wide_intermediate': 10,
        'deepeditpos_inter': 7,
        'RTTlength': 11,
        'correction_type': 12,
        'mutation_type': 13,
        'correction_length': 14
    }

    # Create a mock sequence record with an indel (doesn't match either template)
    seq_record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACCTACGT"), id="id1")  # Different mutation
    seq_record.letter_annotations["phred_quality"] = [40] * 28  # High quality

    # Mock the file opening and parsing
    with patch('gzip.open'), \
         patch('Bio.SeqIO.parse', return_value=[seq_record]):

        # Call the function with our mocked data
        read_count = process_target_file(
            "dummy_path", "test_sample", diseasedf, final_diseasedf,
            templatenumpy, col_indices, filtered='', barcode_start=0,
            barcode_length=6, min_quality=20, basesbeforenick=2,
            baseafterflap=5, target_start_to_nick_distance=10,
            nick_pos_in_wide_initial_target=10
        )

    # Check that one read was processed
    assert read_count == 1, "Expected to process 1 read"

    # Check that the indelcount was incremented for variant 0
    # Since our mock sequence has an indel compared to both templates
    assert templatenumpy[0, 0] == 0, "Expected uneditedcount to remain 0"
    assert templatenumpy[0, 1] == 0, "Expected editedcount to remain 0"
    assert templatenumpy[0, 2] == 1, "Expected indelcount to be incremented"


def test_process_target_file_multibp_replacement_intermediate():
    """Test process_target_file with a MultibpReplacement that matches an intermediate state."""
    # Create mocked dataframes
    diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    final_diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]],
        'matchwithextension': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    # Create a simple template array for testing
    templatenumpy = np.zeros((3, 15), dtype=object)

    # Set up template values for MultibpReplacement
    templatenumpy[0, 11] = 10  # RTTlength
    templatenumpy[0, 12] = "Replacement"  # correction_type
    templatenumpy[0, 13] = "MultibpReplacement"  # mutation_type
    templatenumpy[0, 14] = 3  # correction_length spanning 3 bases
    templatenumpy[0, 8] = "ACGTACGTACGTACGTACGTACGTACGT"  # wide_initial_target
    templatenumpy[0, 9] = "ACGTACGTACGTACGTACGTtCaTACGT"  # wide_mutated_target (2 changes)

    # Set up intermediate states for MultibpReplacement
    intermediate1 = "ACGTACGTACGTACGTACGTtCGTACGT"  # First change only
    intermediate2 = "ACGTACGTACGTACGTACGTACaTACGT"  # Second change only
    templatenumpy[0, 10] = [intermediate1, intermediate2]  # wide_intermediate
    templatenumpy[0, 7] = [[20], [22]]  # deepeditpos_inter - positions of changes

    # Column indices for lookups
    col_indices = {
        'uneditedcount': 0,
        'editedcount': 1,
        'indelcount': 2,
        'nickindelcount': 3,
        'intermediatecount': 4,
        'intermediatepos': 5,
        'beforeflapindel': 6,
        'wide_initial': 8,
        'wide_mutated': 9,
        'wide_intermediate': 10,
        'deepeditpos_inter': 7,
        'RTTlength': 11,
        'correction_type': 12,
        'mutation_type': 13,
        'correction_length': 14
    }

    # Initialize intermediate positions list
    templatenumpy[0, 5] = []

    # Create a mock sequence record that matches an intermediate state
    seq_record = SeqRecord(Seq(intermediate1), id="id1")  # Matches first intermediate
    seq_record.letter_annotations["phred_quality"] = [40] * 28  # High quality

    # Mock the file opening and parsing
    with patch('gzip.open'), \
         patch('Bio.SeqIO.parse', return_value=[seq_record]):

        # Call the function with our mocked data
        read_count = process_target_file(
            "dummy_path", "test_sample", diseasedf, final_diseasedf,
            templatenumpy, col_indices, filtered='', barcode_start=0,
            barcode_length=6, min_quality=20, basesbeforenick=2,
            baseafterflap=5, target_start_to_nick_distance=10,
            nick_pos_in_wide_initial_target=10
        )

    # Check that one read was processed
    assert read_count == 1, "Expected to process 1 read"

    # For intermedia sequences, they are counted as indels but also increment intermediatecount
    assert templatenumpy[0, 0] == 0, "Expected uneditedcount to remain 0"
    assert templatenumpy[0, 1] == 0, "Expected editedcount to remain 0"
    assert templatenumpy[0, 2] == 1, "Expected indelcount to be incremented for intermediate"
    assert templatenumpy[0, 4] == 1, "Expected intermediatecount to be incremented"
    assert templatenumpy[0, 5] == [[20]], "Expected intermediatepos to record the position of change"


def test_process_target_file_low_quality():
    """Test process_target_file skips low quality reads."""
    # Create mocked dataframes
    diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    final_diseasedf = pd.DataFrame({
        'match': [[0], [1], [2]],
        'matchwithextension': [[0], [1], [2]]
    }, index=['id1', 'id2', 'id3'])

    # Create a simple template array for testing
    templatenumpy = np.zeros((3, 15), dtype=object)

    # Column indices for lookups
    col_indices = {
        'uneditedcount': 0,
        'editedcount': 1,
        'indelcount': 2,
        'nickindelcount': 3,
        'intermediatecount': 4,
        'intermediatepos': 5,
        'beforeflapindel': 6,
        'wide_initial': 8,
        'wide_mutated': 9,
        'wide_intermediate': 10,
        'deepeditpos_inter': 7,
        'RTTlength': 11,
        'correction_type': 12,
        'mutation_type': 13,
        'correction_length': 14
    }

    # Create a mock sequence record with low quality scores
    seq_record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGT"), id="id1")
    seq_record.letter_annotations["phred_quality"] = [10] * 28  # Low quality scores

    # Mock the file opening and parsing
    with patch('gzip.open'), \
         patch('Bio.SeqIO.parse', return_value=[seq_record]):

        # Call the function with high min_quality threshold
        read_count = process_target_file(
            "dummy_path", "test_sample", diseasedf, final_diseasedf,
            templatenumpy, col_indices, filtered='', barcode_start=0,
            barcode_length=6, min_quality=30,  # High threshold to skip the read
            basesbeforenick=2, baseafterflap=5,
            target_start_to_nick_distance=10,
            nick_pos_in_wide_initial_target=10
        )

    # Check that read was counted but not analyzed
    assert read_count == 0, "Expected to count 0 read"
    assert templatenumpy[0, 0] == 0, "Low quality read should be skipped (no uneditedcount increment)"
    assert templatenumpy[0, 1] == 0, "Low quality read should be skipped (no editedcount increment)"
    assert templatenumpy[0, 2] == 0, "Low quality read should be skipped (no indelcount increment)"


def test_main_with_real_data_shortreads(monkeypatch, test_data_path, tmpdir):
    """Test the main function with real test data files instead of extensive mocking."""
    # Set up paths and arguments
    template_file = os.path.join(test_data_path, "20240624_mc4r_library_designed_final.csv")
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")
    output_dir = str(tmpdir)
    date_prefix = "20250101"
    library_name = "MC4R"
    shortname = "test_MC4R-lib-control1"

    # Mock sys.argv
    monkeypatch.setattr('sys.argv', [
        'self_targeting_analysis.py',
        '--shortname', shortname,
        '--path', fastq_path,
        '--template_file', template_file,
        '--output_dir', output_dir,
        '--date_prefix', date_prefix,
        '--library_name', library_name,
        '--barcode_start', '0',
        '--barcode_length', '6',
        '--min_quality', '25',
        '--basesbeforenick', '2',
        '--baseafterflap', '5',
        '--target_start_to_nick_distance', '27',
        '--nick_pos_in_wide_initial_target', '26',
        '--rtt_lookup_length', '6',
        '--read_length', '150'
    ])

    # Only mock a few functions that require real resources or that we want to verify
    with patch('scripts.self_targeting_analysis.print_summary_stats') as mock_print, \
         patch('scripts.self_targeting_analysis.save_summary_stats', return_value={}) as mock_save:

        # Call the main function
        main()

        # Verify output file was created
        output_file = f"{output_dir}/{date_prefix}_{library_name}_{shortname}_analysisdf_focused.csv"
        assert os.path.exists(output_file), f"Output file not created: {output_file}"

        # Verify summary stats functions were called
        assert mock_print.called, "print_summary_stats should be called"
        assert mock_save.called, "save_summary_stats should be called"

        # Optionally load and check content of output file
        df = pd.read_csv(output_file)
        assert not df.empty, "Output DataFrame should not be empty"


def test_main_with_real_data_longreads(monkeypatch, test_data_path, tmpdir):
    """Test the main function with real test data files instead of extensive mocking."""
    # Set up paths and arguments
    template_file = os.path.join(test_data_path, "20240624_mc4r_library_designed_final.csv")
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")
    output_dir = str(tmpdir)
    date_prefix = "20250101"
    library_name = "MC4R"
    shortname = "test_long_MC4R-lib-control1"

    # Mock sys.argv
    monkeypatch.setattr('sys.argv', [
        'self_targeting_analysis.py',
        '--shortname', shortname,
        '--path', fastq_path,
        '--template_file', template_file,
        '--output_dir', output_dir,
        '--date_prefix', date_prefix,
        '--library_name', library_name,
        '--barcode_start', '0',
        '--barcode_length', '6',
        '--min_quality', '25',
        '--basesbeforenick', '2',
        '--baseafterflap', '5',
        '--target_start_to_nick_distance', '27',
        '--nick_pos_in_wide_initial_target', '26',
        '--rtt_lookup_length', '6',
        '--read_length', '300'
    ])

    # Only mock a few functions that require real resources or that we want to verify
    with patch('scripts.self_targeting_analysis.print_summary_stats') as mock_print, \
         patch('scripts.self_targeting_analysis.save_summary_stats', return_value={}) as mock_save:

        # Call the main function
        main()

        # Verify output file was created
        output_file = f"{output_dir}/{date_prefix}_{library_name}_{shortname}_analysisdf_focused.csv"
        assert os.path.exists(output_file), f"Output file not created: {output_file}"

        # Verify summary stats functions were called
        assert mock_print.called, "print_summary_stats should be called"
        assert mock_save.called, "save_summary_stats should be called"

        # Optionally load and check content of output file
        df = pd.read_csv(output_file)
        assert not df.empty, "Output DataFrame should not be empty"


def test_main_with_real_data_longreads_shortconfig(monkeypatch, test_data_path, tmpdir, capsys):
    """Test the main function with real test data files instead of extensive mocking."""
    # Using long reads but short read config for testing ###
    # Should print a warning but run through without error ###
    # Set up paths and arguments
    template_file = os.path.join(test_data_path, "20240624_mc4r_library_designed_final.csv")
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")
    output_dir = str(tmpdir)
    date_prefix = "20250101"
    library_name = "MC4R"
    shortname = "test_long_MC4R-lib-control1"

    # Mock sys.argv
    monkeypatch.setattr('sys.argv', [
        'self_targeting_analysis.py',
        '--shortname', shortname,
        '--path', fastq_path,
        '--template_file', template_file,
        '--output_dir', output_dir,
        '--date_prefix', date_prefix,
        '--library_name', library_name,
        '--barcode_start', '0',
        '--barcode_length', '6',
        '--min_quality', '25',
        '--basesbeforenick', '2',
        '--baseafterflap', '5',
        '--target_start_to_nick_distance', '27',
        '--nick_pos_in_wide_initial_target', '26',
        '--rtt_lookup_length', '6',
        '--read_length', '150'
    ])

    # Only mock a few functions that require real resources or that we want to verify
    with patch('scripts.self_targeting_analysis.print_summary_stats') as mock_print, \
         patch('scripts.self_targeting_analysis.save_summary_stats', return_value={}) as mock_save:

        # Call the main function
        main()

        # Capture stdout
        captured = capsys.readouterr()

        # Check for the warning message
        assert "WARNING: Read length configuration mismatch!" in captured.out
        assert "Config specifies 150bp, but median read length is" in captured.out
        assert "Tevopreq filtering is disabled, but could be enabled" in captured.out

        # Verify output file was created
        output_file = f"{output_dir}/{date_prefix}_{library_name}_{shortname}_analysisdf_focused.csv"
        assert os.path.exists(output_file), f"Output file not created: {output_file}"

        # Verify summary stats functions were called
        assert mock_print.called, "print_summary_stats should be called"
        assert mock_save.called, "save_summary_stats should be called"

        # Optionally load and check content of output file
        df = pd.read_csv(output_file)
        assert not df.empty, "Output DataFrame should not be empty"


def test_main_with_real_data_shortreads_longconfig(monkeypatch, test_data_path, tmpdir, capsys):
    """Test the main function with real test data files instead of extensive mocking."""
    # Using short reads but long read config for testing ###
    # Should raise an error ###
    # Set up paths and arguments
    template_file = os.path.join(test_data_path, "20240624_mc4r_library_designed_final.csv")
    fastq_path = os.path.join(test_data_path, "fastq", "trimmed")
    output_dir = str(tmpdir)
    date_prefix = "20250101"
    library_name = "MC4R"
    shortname = "test_MC4R-lib-control1"

    # Mock sys.argv
    monkeypatch.setattr('sys.argv', [
        'self_targeting_analysis.py',
        '--shortname', shortname,
        '--path', fastq_path,
        '--template_file', template_file,
        '--output_dir', output_dir,
        '--date_prefix', date_prefix,
        '--library_name', library_name,
        '--barcode_start', '0',
        '--barcode_length', '6',
        '--min_quality', '25',
        '--basesbeforenick', '2',
        '--baseafterflap', '5',
        '--target_start_to_nick_distance', '27',
        '--nick_pos_in_wide_initial_target', '26',
        '--rtt_lookup_length', '6',
        '--read_length', '250'
    ])

    # Only mock a few functions that require real resources or that we want to verify
    with patch('scripts.self_targeting_analysis.print_summary_stats') as mock_print, \
         patch('scripts.self_targeting_analysis.save_summary_stats', return_value={}) as mock_save:

        # Call the main function
        with pytest.raises(ValueError, match="Read length configuration mismatch."):
            main()

            # Capture stdout
            captured = capsys.readouterr()

            # Check for the warning message
            assert "WARNING: Read length configuration mismatch!" in captured.out
            assert "Config specifies 250bp, but median read length is" in captured.out

            # Verify output file was created
            output_file = f"{output_dir}/{date_prefix}_{library_name}_{shortname}_analysisdf_focused.csv"
            assert not os.path.exists(output_file), f"Output file wrongfully created: {output_file}"

            # Verify summary stats functions were called
            assert not mock_print.called, "print_summary_stats should be called"
            assert not mock_save.called, "save_summary_stats should be called"
