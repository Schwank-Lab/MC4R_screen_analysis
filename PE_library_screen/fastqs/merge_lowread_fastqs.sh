#!/bin/bash
zcat MC4R_lib_PE7_sample1_S169_L001_R1_001.fastq.gz MC4R_lib_PE71_S7_L001_R1_001.fastq.gz | gzip > MC4R_lib_PE71_merge_L001_R1_001.fastq.gz
zcat MC4R_lib_PE7_sample1_S169_L001_R2_001.fastq.gz MC4R_lib_PE71_S7_L001_R2_001.fastq.gz | gzip > MC4R_lib_PE71_merge_L001_R2_001.fastq.gz
