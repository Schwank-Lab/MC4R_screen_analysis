import pandas as pd 
from pathlib import Path 
import os
import library_layouts



LIB_LAYOUT = library_layouts.LAYOUTS['libnb-010-MC4R']
lib_filename = 'LibNB-010_MC4R.xlsx'
sublib_filename = 'SubLib_MC4R-ABE.xlsx'


LOW_READ_TH = 100 # minimum number of reads to be considered a valid point in any analysis


WORKDIR = Path("/home/tao/t4-volume/tao/screens/MC4R")

assert WORKDIR.exists(), f'WORKDIR {WORKDIR} does not exist'
FASTQ_DIR = WORKDIR / 'data/fastqs'

def get_samples(): 
    return pd.read_excel(WORKDIR / 'meta/sample_sheet.xlsx', index_col=False)

def fastq_r1_file(fastq_file): 
    return FASTQ_DIR / fastq_file

def fastq_r2_file(fastq_file):
    return FASTQ_DIR / fastq_file





parts_dir = WORKDIR  / 'out/01_processed'
parts_dir.mkdir(exist_ok=True, parents=True)
def parts(sample):
    return parts_dir / (sample + '.parts.csv')
def parts_r1(sample): 
    return parts_dir / (sample + '.parts_R1.csv')
def parts_r2(sample): 
    return parts_dir / (sample + '.parts_R2.csv')

filtered_dir = WORKDIR / 'out/02_filtered'
filtered_dir.mkdir(exist_ok=True, parents=True)
def filtered(sample): 
    return filtered_dir / (sample + '.filtered.csv')

sublibs_dir = WORKDIR / 'out/03_sublibs'
sublibs_dir.mkdir(exist_ok=True, parents=True)
def sublib_reads(sample, sublib):
    return sublibs_dir / f'{sample}_{sublib}.filtered.csv'

filtered_sublibs_dir = WORKDIR / 'out/04_filtered_sublibs'
filtered_sublibs_dir.mkdir(exist_ok=True, parents=True)
def filtered_sublib_reads(sample, sublib):
    return filtered_sublibs_dir / f'{sample}_{sublib}.filtered.csv'

outcome_dir = WORKDIR / 'out/05_outcomes'
outcome_dir.mkdir(exist_ok=True, parents=True)
def outcomes(sample, lib): 
    return outcome_dir / f'{sample}_{lib}.outcomes.csv'

edit_rate_dir = WORKDIR / 'out/06_edit_rate'
edit_rate_dir.mkdir(exist_ok=True, parents=True)
def edit_rate(sample, lib,suffix='.edit_rate'):
    return edit_rate_dir / f'{sample}_{lib}{suffix}.csv'


FASTQC_DIR = WORKDIR / 'out/fastqc'
FASTQC_DIR.mkdir(exist_ok=True, parents=True)
MULTIQC_DIR = WORKDIR / 'out/multiqc'
MULTIQC_DIR.mkdir(exist_ok=True, parents=True) 