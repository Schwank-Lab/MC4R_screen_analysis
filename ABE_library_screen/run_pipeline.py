from prefect import flow, task

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import seaborn as sns

import pipeline
import config


@flow
def run_pipeline():
    pipeline.process_and_filter_reads(config)
    pipeline.process_target_sites(config)
    pipeline.collect_outcomes(config)
    
    
    pipeline.calculate_abe_edit_rates(config)
    pipeline.collect_abe_edit_rates(config)
    pipeline.normalize_abe_editing_rate(config)
    pipeline.aggregate_abe_replicates(config)
    pipeline.extract_ontarget_editing_rates(config)
    
    pipeline.calculate_precise_abe_edit_rates(config)
    pipeline.collect_abe_edit_rates(config, in_suffix='.precise-edit_rate', out_suffix='.edit_rate_precise')
    pipeline.normalize_abe_editing_rate(config,in_suffix='.edit_rate_precise')
    pipeline.aggregate_abe_replicates(config,in_suffix='.edit_rate_precise')
    pipeline.extract_ontarget_editing_rates(config,in_suffix='.edit_rate_precise')

    
    

if __name__ == '__main__': 
    run_pipeline()
   