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
hmmerHits = pd.read_csv('/home/gridsan/lbrenes/hhpred/hhblts/hhmscan_singles_nr/kilAN_extendedcov_NR_hits_filtered_stringent_dom_nodupes.csv',usecols=[2,3,5],names=['Domain','DomainAcc','Query'],header=0)
hmmerHits



# %%


mEpX1_KilANs = set(hmmerHits[(hmmerHits['Domain'] == 'mEpX1_KilAN.extendedcov.blts')]['Query'].values.tolist())
PhiR49_1_KilANs = set(hmmerHits[(hmmerHits['Domain'] == 'PhiR49.1_KilAN.extendedcov.blts')]['Query'].values.tolist())

# %%


#Assembly taxid and species taxid information for every genbank and refseq genome and store in dictionary
assemblySums = pd.read_csv('/home/gridsan/lbrenes/hhpred/assembly_summary.txt',sep = '\t',skiprows = 1,usecols = ['#assembly_accession','species_taxid'])
Acc2Taxid = {row[0]:row[1] for index, row in assemblySums.iterrows()}


# %%


def checkFT(FT_file):
    mEp332_cIoperonic_homolog_system = []
    mEpX1_Gp15 = []
    mEpX1_KilAN = []
    mEpX1_Nun = []
    PhiR12_1_cIoperonic_system = []
    PhiR26_1_SNIPE_system = []
    PhiR46_1_spacklelike = []
    PhiR49_1_KilAN = []
    
    FT = pd.read_csv(FT_file,sep = '\t')
    FT = FT[FT['# feature'] == 'CDS'].reset_index()
    
    Acc = FT.at[0,'assembly']
    ID = Acc2Taxid[Acc]
    FTloc = os.path.basename(FT_file)
    
    all_prots = set(FT['non-redundant_refseq'])
    mEpX1_KilAN_IDs = Prot =  all_prots.intersection(mEpX1_KilANs)
    
    if mEpX1_KilAN_IDs:
        mEpX1_KilAN = [(FTloc,Acc,ID,Prot)]

    PhiR49_1_KilAN_IDs = Prot = all_prots.intersection(PhiR49_1_KilANs)
    
    if PhiR49_1_KilAN_IDs:
        PhiR49_1_KilAN = [(FTloc,Acc,ID,Prot)]

    return mEpX1_KilAN, PhiR49_1_KilAN 


# %%


all_FT_files = glob.glob('/home/gridsan/lbrenes/LaubLab_shared/all_refseq_feature_tables/*.txt')


# %%


results=Parallel(n_jobs=47)(delayed(checkFT)(x) for x in all_FT_files)


# %%


results = pd.DataFrame(results)
mEpX1_KilAN_out = []
PhiR49_1_KilAN_out = []
for index, row in results.iterrows():
    if row[0]:
        mEpX1_KilAN_out = mEpX1_KilAN_out + row[0]
    if row[1]:
        PhiR49_1_KilAN_out = PhiR49_1_KilAN_out + row[1]

mEpX1_KilAN_out = pd.DataFrame(mEpX1_KilAN_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
mEpX1_KilAN_out.to_csv('TaxIDswmEpX1_KilAN_nr_dom_filtered_STRINGENT_prot.csv')

PhiR49_1_KilAN_out = pd.DataFrame(PhiR49_1_KilAN_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
PhiR49_1_KilAN_out.to_csv('TaxIDswPhiR49_1_KilAN_nr_dom_filtered_STRINGENT_prot.csv')



# %%




