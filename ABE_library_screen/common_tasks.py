import pandas as pd
import numpy as np 
from prefect import task
from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path
import gzip


def open_fastq(fastq_path, mode='r'):
    if isinstance(fastq_path, Path): 
        fastq_path = str(fastq_path)
    if fastq_path.endswith('.gz'):
        return gzip.open(fastq_path, mode+'t')
    else:
        return open(fastq_path, mode)


def decompose_read_into_parts(read_id, r, found_anchor_pos, anchor_idx, part_lens, part_names):
    """Splits read into parts accoring to the layout."""

    layout_anchor_start = sum(part_lens[:anchor_idx])
    part_seqs = []
    read_len = len(r)
    for i in range(len(part_lens)):
        layout_part_start = sum(part_lens[:i])
        offset = layout_part_start - layout_anchor_start
        pos_in_read = found_anchor_pos + offset
        part_len = part_lens[i]
        if pos_in_read < 0:
            print(f'Warning: read {read_id} part {part_names[i]} start is outside of the read, make sure the part_offset is correct')
            return None
        elif pos_in_read + part_len > read_len:
            print(f'Warning: read {read_id} part {part_names[i]} end is outside of the read, make sure the part_len is correct')
            return None
        part_seq = r[pos_in_read: pos_in_read + part_len]
        part_seqs.append(part_seq)
    return part_seqs


def anchor_reads(sample, layout, anchor_seq, fastq: Path, out_csv: Path, reverse_complement=False, part_anchor='anchor'):
    """
    Anchors reads by finding the anchor sequence and extracting sequence parts.

    params: 
        parts: tuples (part_len, part_name)
        anchor_seq: anchor sequence to find in the read.
        reverse_complement: if True, reverse complement the read before anchoring.

        part_anchor: name of the part to use an an anchor.

    """
    num_anchor_error = 0 
    num_reads = 0 
    part_names = [part_name for _, part_name in layout]
    part_lens = [part_len for part_len, _ in layout]
    try:
        anchor_idx = part_names.index(part_anchor)
    except ValueError:
        raise ValueError(f"Anchor part '{part_anchor}' not found in layout.")
    
    with open(out_csv, 'w') as out_f, open_fastq(fastq, 'r') as in_f:
         out_f.write(','.join(['read_id'] + part_names) + '\n')
         for read in SeqIO.parse(in_f, 'fastq'): 
            num_reads += 1
            r = str(read.seq) if not reverse_complement else str(read.seq.reverse_complement())
            found_anchor_pos = r.find(anchor_seq)

            if found_anchor_pos == -1 or r.count(anchor_seq) > 1:
                num_anchor_error += 1
                continue
            part_seqs = decompose_read_into_parts(read.id, r, found_anchor_pos, anchor_idx, part_lens, part_names)
            if part_seqs is None: 
                num_anchor_error += 1
                continue
            out_f.write(','.join([read.id] + part_seqs) + '\n')
           

    print(f"Processed {sample} (reverse={reverse_complement}): reads {num_reads}, anchor errors {num_anchor_error}")

    frac = num_anchor_error / num_reads if num_reads > 0 else np.nan
    stats = pd.Series({
        'sample': sample, 
        'num_reads': num_reads,
        'frac_anchor_error': frac,
    })
    return stats



@task(log_prints=True)
def process_sample(sample, r1_fastq, config, r2_fastq=None):
    print(f'Processing {sample}')
    """Anchor reads and extract read parts."""
    r1_fastq = config.fastq_r1_file(r1_fastq)
    assert r1_fastq.exists(), f'R1 fastq file {r1_fastq} does not exist'
    if not r2_fastq: 
        r2_fastq = config.fastq_r2_file(r1_fastq)
    else: 
        r2_fastq = config.fastq_r2_file(r2_fastq)
    assert r2_fastq.exists(), f'R2 fastq file {r2_fastq} does not exist'
    
    r1_stats = anchor_reads(sample, config.LIB_LAYOUT['r1_layout'], config.LIB_LAYOUT['r1_anchor'], 
        r1_fastq, config.parts_r1(sample), part_anchor='r1_anchor'
    )
    r2_stats = anchor_reads(sample, config.LIB_LAYOUT['r2_layout'], config.LIB_LAYOUT['r2_anchor'],  r2_fastq, config.parts_r2(sample), 
        reverse_complement=True, part_anchor='r2_anchor')
    
    
    d = pd.merge(pd.read_csv(config.parts_r1(sample), index_col=0), pd.read_csv(config.parts_r2(sample), index_col=0), left_index=True, right_index=True)
    
    # the whole target site is added to the amplicon as a reverse complement, so reverse it back.
    d.loc[:, 'target_site'] = d.target_site_rc.map(lambda x: str(Seq(x).reverse_complement()))
    
    d.to_csv(config.parts(sample))
    return pd.Series({
        'sample_name': sample, 'num_reads': r1_stats['num_reads'], 'r1_frac_anchor_error': r1_stats['frac_anchor_error'], 'r2_frac_anchor_error': r2_stats['frac_anchor_error'],
    })


@task(log_prints=True)
def filter_reads(sample: str, lib: pd.DataFrame, config) -> pd.Series: 
    """ Filter out reads that don't match library members."""
    d = pd.read_csv(config.parts(sample), index_col=0)
    num_reads = len(d.index)
    frac_sublib_barcode_match = d.sublibrary_barcode.isin(lib['sublibrary_barcode'].unique()).sum() / num_reads
    frac_member_barcode_match = d.member_barcode.isin(lib['member_barcode'].unique()).sum() / num_reads
    
    # retain reads with known sub-library and member barcodes
    d = pd.merge(d, lib[['guide', 'sublibrary_barcode', 'member_barcode', 'site_id', 'sublibrary', 'target_site_rc']], 
                 left_on=['sublibrary_barcode', 'member_barcode'], 
                 right_on=['sublibrary_barcode', 'member_barcode'], 
                 how='inner', suffixes=('', '_lib'))
    
    d.loc[:, 'target_site_lib'] = d.target_site_rc_lib.map(lambda x: str(Seq(x).reverse_complement()))

    d['guide_match'] = d.guide == d.guide_lib
    d['guide_recombination'] = (d.guide != d.guide_lib) & (d.guide.isin(lib['guide'].unique()))
    num_reads_barcode_match = len(d.index)
    frac_guide_match = d.guide_match.sum() / num_reads_barcode_match
    frac_guid_recombination = d.guide_recombination.sum() / num_reads_barcode_match

    # retain reads with matching guide 
    d = d[d.guide_match]
    d.drop(columns=['guide_match', 'guide_recombination', 'guide_lib', 'guide'], inplace=True)
    num_reads_guide_match = len(d.index)    

    d.to_csv(config.filtered(sample))

    return pd.Series({
        'sample_name': sample, 'num_reads': num_reads, 'frac_sublib_barcode_match': frac_sublib_barcode_match, 
        'frac_member_barcode_match': frac_member_barcode_match, 'frac_guide_match': frac_guide_match, 
        'frac_guid_recombination': frac_guid_recombination, 'num_reads_barcode_match': num_reads_barcode_match, 
        'num_reads_guide_match': num_reads_guide_match
    })



    