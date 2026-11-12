from prefect import task, unmapped

import pandas as pd 
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt 
from itertools import product

import screen_utils
import common_tasks


import pandas as pd

def get_library(config):  # config unused but kept for compatibility
    df = pd.read_excel(config.lib_filename, index_col=0)  # assuming pd.read_excel as the load function
    return df.rename(columns={
        'concatsequence-250bp': 'amplicon',
        'part3-guidespacer-20bp': 'guide',
        'part5-targetinreverse-64bp-setup20plus4plus20plus20': 'target_site_rc',
        'part6-BC1-7bp': 'member_barcode',
        'part7-BC2-15bp': 'sublibrary_barcode',
    })


def get_sublibraries(config):
    return ['MC4R-ABE']

def load_sublibrary(config):  # config unused but kept for compatibility
    return pd.read_excel(config.sublib_filename)

@task(log_prints=True)
def split_sublibraries(sample: str, config):
    print(f'Splitting sublibraries for {sample}')
    print(f'Parsing target sites for {sample}')
    d = pd.read_csv(config.filtered(sample), index_col=0)

    def extract_spacer_pam_for_sublibrary(d_sublib, sublibrary):
        """Extract spacer and PAM for a given sublibrary DataFrame."""
        if sublibrary == 'MC4R-ABE':
            # flatnk -- spacer (20bp) -- pam (4bp) -- flank (20bp)
            d_sublib['spacer'] = d_sublib['target_site'].str.slice(20, 40)
            d_sublib['PAM'] = d_sublib['target_site'].str.slice(40, 44)
            d_sublib['lib_spacer'] = d_sublib['target_site_lib'].str.slice(20, 40)
            d_sublib['lib_PAM'] = d_sublib['target_site_lib'].str.slice(40, 44)
        else:
            raise ValueError(f'Unknown sublibrary: {sublibrary}')
        return d_sublib
    
    # Split reads per sublibrary and apply the extraction function to the whole sub-DataFrame
    for sublibrary in d['sublibrary'].unique():
        d_sublib = d[d['sublibrary'] == sublibrary].copy()
        d_sublib = extract_spacer_pam_for_sublibrary(d_sublib, sublibrary)
        d_sublib.to_csv(config.sublib_reads(sample, sublibrary))


@task(log_prints=True) 
def calculate_library_distribution(sample, lib, config): 
    d = pd.read_csv(config.filtered(sample), index_col=0)
    member_counts = d.groupby(['sublibrary_barcode', 'member_barcode']).size().reset_index()
    member_counts = pd.merge(member_counts, lib[['sublibrary_barcode', 'member_barcode', 'sublibrary', 'site_id']], 
                             on=['sublibrary_barcode', 'member_barcode'], how='left').drop(columns=['sublibrary_barcode', 'member_barcode'])
    member_counts = member_counts.rename(columns={0: 'reads'})
    member_counts['sample_name'] = sample
    return member_counts


@task
def process_and_filter_reads(config):
    samples = config.get_samples()
    res  = common_tasks.process_sample.map(samples['sample_name'], samples['file_name_R1'], unmapped(config), samples['file_name_R2'])
    res = [res.result() for res in res]
    res = pd.DataFrame(res)
    res.to_excel(config.WORKDIR / 'results/read_stats.xlsx')
    
    lib = get_library(config)
    res = common_tasks.filter_reads.map(samples['sample_name'], unmapped(lib), unmapped(config))
    res = [res.result() for res in res]
    res = pd.DataFrame(res)
    res.to_excel(config.WORKDIR / f'results/library_stats.xlsx')

    res = calculate_library_distribution.map(samples['sample_name'], unmapped(lib), unmapped(config))
    res = [res.result() for res in res]
    res = pd.concat(res).reset_index(drop=True)
    res.to_excel(config.WORKDIR / f'results/library_member_stats.xlsx')


@task(log_prints=True)
def filter_target_sites(sample, sublib: str, config):
    d = pd.read_csv(config.sublib_reads(sample, sublib), index_col=0)
    num_reads = len(d.index)

    d = d[d['PAM'] == d['lib_PAM']].drop(columns=['lib_PAM'])
    frac_pam_match = len(d.index) / num_reads
    d.to_csv(config.filtered_sublib_reads(sample, sublib))
    return pd.Series({'sample_name': sample, 'lib': sublib, 'num_reads': num_reads, 'frac_pam_match': frac_pam_match})


@task(log_prints=True)
def count_outcomes(sample: str, sublib: str, config):
    """Count editing outcomes for a given sample."""
    print(f'Counting outcomes for {sample} {sublib}')
    d_sublib = pd.read_csv(config.filtered_sublib_reads(sample, sublib), index_col=0)
    # TODO: do this in a library function, always take the first elem from the rest.
    outcomes = d_sublib.groupby(['site_id', 'sublibrary', 'PAM', 'spacer', 'lib_spacer']).size().reset_index(name='num_reads')
    outcomes.to_csv(config.outcomes(sample, sublib))


@task(log_prints=True)
def process_target_sites(config):
    samples = config.get_samples()
    
    res = split_sublibraries.map(samples['sample_name'], unmapped(config))
    _ = [res.result() for res in res]
    
    sample_sublib_combinations = product(samples['sample_name'], get_sublibraries(config))
    
    stats = filter_target_sites.map(*zip(*sample_sublib_combinations), unmapped(config))
    stats = [res.result() for res in stats]
    stats = pd.DataFrame(stats)
    stats.to_excel(config.WORKDIR / 'results/sublib_stats.xlsx')
    
    sample_sublib_combinations = product(samples['sample_name'], get_sublibraries(config))
    res = count_outcomes.map(*zip(*sample_sublib_combinations), unmapped(config))
    _ = [res.result() for res in res]
    
    
@task(log_prints=True)
def calculate_precise_abe_edit_rate(sample, sublib, config):
    """"precise means only on-target bases are different"""
    print(f'Analyzing bystander rates for {sample} {sublib}')
    o = pd.read_csv(config.outcomes(sample, sublib), index_col=0)
    o = screen_utils.extract_abe_edit_positions(o, col_ref_spacer='lib_spacer')
    o = o[o['edit_position'].notna()]
    o['edit_position'] = o['edit_position'].astype(int)
    o['is_target_base_edited'] = o.apply(lambda row: screen_utils.is_precise_ontarget_a_to_g_conversion(row.spacer, row.lib_spacer, row.edit_position), axis=1)
    e = screen_utils.calculate_edit_rate_per_group(o, col_is_outcome_edited='is_target_base_edited', cols_groupby=['site_id', 'edit_position'])
    e.to_csv(config.edit_rate(f'{sample}-precise-unfiltered', sublib))
    
@task(log_prints=True)
def calculate_abe_edit_rate(sample, sublib, config):
    print(f'Analyzing bystander rates for {sample} {sublib}')
    o = pd.read_csv(config.outcomes(sample, sublib), index_col=0)
    o = screen_utils.extract_abe_edit_positions(o, col_ref_spacer='lib_spacer')
    o = o[o['edit_position'].notna()]
    o['edit_position'] = o['edit_position'].astype(int)
    o['is_target_base_edited'] = o.apply(lambda row: screen_utils.is_ontarget_a_to_g_conversion(row.spacer, row.lib_spacer, row.edit_position), axis=1)
    e = screen_utils.calculate_edit_rate_per_group(o, col_is_outcome_edited='is_target_base_edited', cols_groupby=['site_id', 'edit_position'])
    e.to_csv(config.edit_rate(f'{sample}-unfiltered', sublib))


@task(log_prints=True)
def filter_edit_rates(sample, lib_name, config,in_suffix='unfiltered',out_suffix='.edit_rate'):
    print(f'Filtering edit rates for {sample} {lib_name}')
    d = pd.read_csv(config.edit_rate(f'{sample}-{in_suffix}', lib_name), index_col=0)
    # filter out low reads.
    d.loc[d.num_reads < config.LOW_READ_TH, 'edit_rate'] = pd.NA
    # fill out missing site ids
    sites = load_sublibrary(config)
    missing_sites = list(set(sites['site_id']) - set(d['site_id']))
    if missing_sites:
        missing_vals = sites.set_index('site_id').loc[missing_sites, ['edit_position', 'PAM']].reset_index()
        missing_vals['num_reads'] = 0
        missing_vals['edit_rate'] = pd.NA
        missing_vals['sublibrary'] = lib_name
        d = pd.concat([d, missing_vals], ignore_index=True)
  
    frac_low_reads = d['edit_rate'].isna().sum() / len(d.index)
    d.to_csv(config.edit_rate(sample, lib_name,suffix=out_suffix))
    return pd.Series({'sample_name': sample, 'lib': lib_name, 'frac_low_reads': frac_low_reads})
    
    

@task(log_prints=True)
def collect_outcomes(config): 
    print('Collecting all outcomes')
    samples = config.get_samples()
    sublibs = get_sublibraries(config)
    res = pd.concat([pd.read_csv(config.outcomes(sample, sublib)).assign(sample_name=sample) for sample in samples['sample_name'] for sublib in sublibs])
    res.to_csv(config.outcomes('all-samples', 'all-libs'))


@task(log_prints=True)
def collect_abe_edit_rates(config,in_suffix='.edit_rate',out_suffix='.edit_rate'):
    """ Collect edit rates for all samples and sublibraries. """

    samples = config.get_samples()
    sublibs = get_sublibraries(config)
    
    all_sample_rates = pd.concat(
        [
            pd.read_csv(config.edit_rate(sample.sample_name, sublib,in_suffix), index_col=0).assign(
                sublibrary=sublib, 
                editor=sample.editor, 
                sample_name=sample.sample_name, 
                bioreplicate=sample.bioreplicate,
                is_ctrl=sample.is_ctrl, 
                replicate=sample.replicate,
            ) 
            for sample in samples.itertuples() 
            for sublib in sublibs
        ]
    ).reset_index(drop=True)
    all_sample_rates.loc[:, 'is_ctrl'] = all_sample_rates.is_ctrl == 'yes'
    
    
    # annotate target position.
    sublibs = get_sublibraries(config)
    sites = pd.concat(
        [
            load_sublibrary(config)[['site_id', 'edit_position']]
            .assign(sublibrary=sublib)
            .rename(columns={'edit_position': 'target_position'})
            for sublib in sublibs
        ],
        ignore_index=True
    ).reset_index(drop=True)
    all_sample_rates = pd.merge(all_sample_rates, sites, on=['site_id', 'sublibrary'], how='left')
    
    all_sites = all_sample_rates[['site_id', 'editor', 'sublibrary', 'target_position']].drop_duplicates()
    missing_positions = all_sites[all_sites.target_position.isna()].drop_duplicates()
    if not missing_positions.empty:
        raise ValueError(f"Missing library entries for the following entries:\n{missing_positions}")
    
    all_sample_rates.to_csv(config.edit_rate('all-abe-samples-not-normalized', 'all-libs',out_suffix))


@task(log_prints=True)
def normalize_abe_editing_rate(config,in_suffix='.edit_rate'):
    d = pd.read_csv(config.edit_rate('all-abe-samples-not-normalized', 'all-libs',in_suffix), index_col=0)
    d = screen_utils.normalize_editing_rate(d, 
        col_normalize_over=['site_id', 'sublibrary', 'edit_position'],
        col_control='is_ctrl',
        th_ctrl_missing=0.7
    ).drop(columns=['is_ctrl'])
    d.to_csv(config.edit_rate('all-abe-samples-normalized', 'all-libs',in_suffix))


@task(log_prints=True)
def aggregate_abe_replicates(config,in_suffix='.edit_rate'):
    d = pd.read_csv(config.edit_rate('all-abe-samples-normalized', 'all-libs',in_suffix), index_col=0)
    d = screen_utils.aggregate_replicates(d, 
        cols_group=['editor', 'site_id', 'sublibrary', 'edit_position'],
        col_rate='edit_rate'
    ).drop(columns=['sample_name', 'replicate'])
    d.to_csv(config.edit_rate('all-abe-editors', 'all-libs',in_suffix))
    

@task(log_prints=True)
def extract_ontarget_editing_rates(config,in_suffix='.edit_rate'): 
    """ 
    Extract on-target editing rates for all editors and sublibraries.
    
    We calculate editing rates for all As in the spacer. Here we extract on the on-target As.
    """
    d = pd.read_csv(config.edit_rate('all-abe-editors', 'all-libs',in_suffix), index_col=0)
    d_sites = d[['site_id', 'editor', 'sublibrary']].drop_duplicates()
    
    d = d[d.target_position == d.edit_position]
    
    d_sites_found_target_pos = d[['site_id', 'editor', 'sublibrary']].drop_duplicates()
    missing_sites = d_sites[~d_sites.site_id.isin(d_sites_found_target_pos.site_id)]
    if not missing_sites.empty:
        raise ValueError(f"Missing target positions for the following entries:\n{missing_sites}")
    
    d.to_csv(config.edit_rate('all-abe-editors-ontarget', 'all-libs',in_suffix))



@task(log_prints=True)
def calculate_abe_edit_rates(config):
    abe_samples = config.get_samples()
    sublibs = get_sublibraries(config)
    abe_sample_sublib_combinations = list(product(abe_samples['sample_name'], sublibs))
    
    res = calculate_abe_edit_rate.map(*zip(*abe_sample_sublib_combinations), unmapped(config))
    res = [r.result() for r in res]
    
    res = filter_edit_rates.map(*zip(*abe_sample_sublib_combinations), unmapped(config))
    res = [r.result() for r in res]
    
    
    
@task(log_prints=True)
def calculate_precise_abe_edit_rates(config):
    abe_samples = config.get_samples()
    sublibs = get_sublibraries(config)
    abe_sample_sublib_combinations = list(product(abe_samples['sample_name'], sublibs))
    
    res = calculate_precise_abe_edit_rate.map(*zip(*abe_sample_sublib_combinations), unmapped(config))
    res = [r.result() for r in res]
    
    res = filter_edit_rates.map(*zip(*abe_sample_sublib_combinations), unmapped(config),unmapped('precise-unfiltered'),unmapped('.precise-edit_rate'))
    res = [r.result() for r in res]
    



