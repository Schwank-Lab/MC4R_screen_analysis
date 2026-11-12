import pandas as pd
import numpy as np
from Bio.Seq import Seq


def extract_abe_edit_positions(df, col_ref_spacer='ref_spacer', col_edit_position='edit_position'):
    """ Extract all possible A->G editing positions from spacers in the df. """ 
    if not df[col_ref_spacer].str.contains('A').all():
        print("Warning: Some reference spacers do not contain 'A':", 
              df[~df[col_ref_spacer].str.contains('A')][col_ref_spacer].tolist())
        
    df[col_edit_position] = df[col_ref_spacer].apply(lambda s: [i + 1 for i, b in enumerate(s) if b == "A"])
    return df.explode(col_edit_position)


def is_ontarget_a_to_g_conversion(spacer, ref_spacer, edit_position):
    """Check if edit is an A to G conversion at target position.
    
    Args:
        spacer (str): spacer sequence
        ref_spacer (str): reference spacer sequence
        edit_position (int): position of edit in spacer
    Returns:
        bool: True if edit is an A to G conversion at target position
    """
    if edit_position is None:
        return False
    if len(spacer) != len(ref_spacer):
        return False
    if edit_position < 1 or edit_position > len(spacer):
        return False
    return ref_spacer[edit_position - 1] == 'A' and spacer[edit_position - 1] == 'G'

def is_precise_ontarget_a_to_g_conversion(spacer, ref_spacer, edit_position):
    """Check if edit is an A to G conversion at target position.
    
    Args:
        spacer (str): spacer sequence
        ref_spacer (str): reference spacer sequence
        edit_position (int): position of edit in spacer
    Returns:
        bool: True if edit is an A to G conversion at target position
    """
    if edit_position is None:
        return False
    if len(spacer) != len(ref_spacer):
        return False
    if edit_position < 1 or edit_position > len(spacer):
        return False
    # Check if reference base at edit position is A and edited base is G
    if ref_spacer[edit_position - 1] != 'A' or spacer[edit_position - 1] != 'G':
        return False
    
    # Check that all other positions are the same (no bystander edits)
    for i in range(len(spacer)):
        if i != edit_position - 1 and spacer[i] != ref_spacer[i]:
            return False
    
    return True

def calculate_edit_rate_per_group(outcomes: pd.DataFrame, col_is_outcome_edited: str, cols_groupby: list[str], 
                                    col_num_reads: str = 'num_reads', col_edit_rate: str = 'edit_rate') -> pd.DataFrame:
    """ Aggregate outcomes into edit rates.
    
    Args:
        col_is_outcome_edited (str): column name for outcome flag. If true, outcomes is counted as edited.
        cols_groupby (List[str]): columns by which to group outcomes together.
        col_num_reads (str): column name for number of reads
        col_edit_rate (str): column name for editing rate
    Returns:
        pd.DataFrame: df with columns
            - col_num_reads: number of total reads per group
            - col_edit_rate: editing rate per group 
            - rest of the columns in ooutocmes
    """

    outcomes['num_outcome_is_edited'] = outcomes[col_num_reads].where(outcomes[col_is_outcome_edited], 0)

    # calculate editing rate per group
    d_grouped = outcomes.groupby(cols_groupby).agg({
        col_num_reads: 'sum',
        'num_outcome_is_edited': 'sum',
        **{col: 'first' for col in outcomes.columns if col not in cols_groupby + [col_num_reads, 'num_outcome_is_edited']}
    }).reset_index()
    d_grouped[col_edit_rate] = d_grouped.num_outcome_is_edited / d_grouped[col_num_reads]
    d_grouped = d_grouped.drop(columns=['num_outcome_is_edited', col_is_outcome_edited])
    return d_grouped

def normalize_editing_rate(data, col_normalize_over=['site_id', 'sublibrary'], col_rate='edit_rate', col_control='is_ctrl', th_ctrl_missing=0.05) -> pd.DataFrame:
    """ Normalize editing rate by substracting background editing rate 

    Args:
        data (pd.DataFrame): data to normalize
        data (pd.DataFrame): data to normalize
        col_normalize_over (List[str]): columns to normalize over
        col_rate (str): column name for editing rate
        col_control (str): column name for control flag.
        th_ctrl_missing (float): this fraction of sites is allowed to be missing control editing rates. Otherwise the script will complain.
        col_rate (str): column name for editing rate
        col_control (str): column name for control flag.
        th_ctrl_missing (float): this fraction of sites is allowed to be missing control editing rates. Otherwise the script will complain.
    Returns:
        pd.DataFrame: normalized data
        pd.DataFrame: normalized data
    """
    d_ctrl = data[data[col_control]]
    d_ctrl = d_ctrl[col_normalize_over + [col_rate]].rename(columns={col_rate: 'edit_rate_ctrl'})
    d_ctrl['edit_rate_ctrl'] = d_ctrl['edit_rate_ctrl'].fillna(0.0)
    # there might be multiple control points per group, aggregate
    d_ctrl = d_ctrl.groupby(col_normalize_over).edit_rate_ctrl.mean().reset_index()
    
    d_no_ctrl = data[~data[col_control]].copy()

    d_norm = pd.merge(d_no_ctrl, d_ctrl, on=col_normalize_over, how='left')
    ctrl_missing = d_norm[d_norm['edit_rate_ctrl'].isna()]
    ctrl_missing = ctrl_missing.groupby(col_normalize_over).first()
    n_missing = len(ctrl_missing) 
    n_total = len(d_ctrl) + n_missing
    if n_missing / n_total <= th_ctrl_missing:
        print(f'{n_missing} out of {n_total} sites are missing control editing rate. Using 0.0 for {n_missing} points.')
        d_norm['edit_rate_ctrl'] = d_norm['edit_rate_ctrl'].fillna(0.0)
    else:
        raise ValueError(f'{n_missing} out of {n_total} sites are missing control editing rate.')

    d_norm[col_rate] = (d_norm[col_rate] - d_norm['edit_rate_ctrl'])/(1- d_norm['edit_rate_ctrl'])
    d_norm[col_rate] = d_norm[col_rate].clip(lower=0)
    d_norm = d_norm.drop(columns=['edit_rate_ctrl'])
    return d_norm


def aggregate_replicates(data: pd.DataFrame, cols_group = ['sublibrary','site_id', 'PAM', 'bioreplicate'], col_rate='edit_rate') -> pd.DataFrame:
    """ Aggregate replicates by taking the mean editing rate. 

    Args:
        data (pd.DataFrame): data to aggregate
        data (pd.DataFrame): data to aggregate
        cols_group (List[str]): columns to group by
        col_rate (str): column name for editing rate
        col_rate (str): column name for editing rate
    Returns:
        pd.DataFrame: aggregated data
    """    
    agg_funcs = {col: 'first' for col in data.columns if col not in cols_group + [col_rate]}
    agg_funcs[col_rate] = 'mean'
    return data.groupby(cols_group).agg(agg_funcs).reset_index()