#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import pandas as pd
import numpy as np
import os
import glob
from joblib import Parallel, delayed


# %%


#Read in hmmer results
hmmerHits = pd.read_csv('/home/gridsan/lbrenes/hhpred/hhblts/hhmscan_multis_nr/multi_extended_gene_systems_NR_hits_filtered_dom_nodupes.csv',usecols=[2,3,5],names=['Domain','DomainAcc','Query'],header=0)
hmmerHits



# %%


mEp506_small_sys_DUF2594s = set(hmmerHits[(hmmerHits['Domain'] == 'mEp506_small_sys_DUF2594.extendedcov.blts')]['Query'].values.tolist())

# %%


#Assembly taxid and species taxid information for every genbank and refseq genome and store in dictionary
assemblySums = pd.read_csv('/home/gridsan/lbrenes/hhpred/assembly_summary.txt',sep = '\t',skiprows = 1,usecols = ['#assembly_accession','species_taxid'])
Acc2Taxid = {row[0]:row[1] for index, row in assemblySums.iterrows()}


# %%


def checkFT(FT_file):
    mEp506_small_sys_DUF2594_system = []
    
    FT = pd.read_csv(FT_file, sep='\t')
    FT = FT[FT['# feature'] == 'CDS'].reset_index()
    
    Acc = FT.at[0, 'assembly']
    ID = Acc2Taxid.get(Acc, 'Unknown')
    FTloc = os.path.basename(FT_file)
        
    all_prots = set(FT['non-redundant_refseq'])
    mEp506_small_sys_DUF2594_IDs = all_prots.intersection(mEp506_small_sys_DUF2594s)
    
    if mEp506_small_sys_DUF2594_IDs:
        mEp506_small_sys_DUF2594_indices = FT[FT['non-redundant_refseq'].isin(mEp506_small_sys_DUF2594_IDs)].index
        for mEp506_small_sys_DUF2594_index in mEp506_small_sys_DUF2594_indices:
                mEp506_small_sys_DUF2594_system.append((FTloc, Acc, ID, FT.at[mEp506_small_sys_DUF2594_index, 'non-redundant_refseq']))
    
    return mEp506_small_sys_DUF2594_system

    

# Collect all feature table files
all_FT_files = glob.glob('/home/gridsan/lbrenes/LaubLab_shared/all_refseq_feature_tables/*.txt')

# Run checkFT function on the files in parallel
results = Parallel(n_jobs=47, verbose=10)(delayed(checkFT)(x) for x in all_FT_files)



# Flatten results into a list of tuples
mEp506_small_sys_DUF2594_out = [item for sublist in results for item in sublist]

# Filter out items that do not have exactly 4 elements and print warnings
filtered_results = []
for item in mEp506_small_sys_DUF2594_out:
    if len(item) != 4:
        print(f"Warning: Found item with length {len(item)} - {item}")
    else:
        filtered_results.append(item)

mEp506_small_sys_DUF2594s_out = pd.DataFrame(filtered_results, columns=['Featuretable', 'Genome', 'TaxID', 'Prot_ID'])
mEp506_small_sys_DUF2594s_out.to_csv('TaxIDswmEp506_small_sys_DUF2594s_nr_dom_filtered_STRINGENT.csv')


# %%




