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
hmmerHits = pd.read_csv('/home/gridsan/lbrenes/hhpred/hhblts/hhmscan_singles_nr/single_gene_systems_NR_hits_filtered_dom_nodupes.csv',usecols=[2,3,5],names=['Domain','DomainAcc','Query'],header=0)
hmmerHits



# %%


mEp332_cIoperonic_homolog_systems = set(hmmerHits[(hmmerHits['Domain'] == 'mEp332_cIoperonic_homolog_system.blts')]['Query'].values.tolist())
mEpX1_Gp15s = set(hmmerHits[(hmmerHits['Domain'] == 'mEpX1_Gp15.blts')]['Query'].values.tolist())
mEpX1_KilANs = set(hmmerHits[(hmmerHits['Domain'] == 'mEpX1_KilAN.blts')]['Query'].values.tolist())
mEpX1_Nuns = set(hmmerHits[(hmmerHits['Domain'] == 'mEpX1_Nun.blts')]['Query'].values.tolist())
PhiR12_1_cIoperonic_systems = set(hmmerHits[(hmmerHits['Domain'] == 'PhiR12.1_cIoperonic_system.blts')]['Query'].values.tolist())
PhiR26_1_SNIPE_systems = set(hmmerHits[(hmmerHits['Domain'] == 'PhiR26.1_SNIPE_system.blts')]['Query'].values.tolist())
PhiR46_1_spacklelikes = set(hmmerHits[(hmmerHits['Domain'] == 'PhiR46.1_spacklelike.blts')]['Query'].values.tolist())
PhiR49_1_KilANs = set(hmmerHits[(hmmerHits['Domain'] == 'PhiR49.1_KilAN.blts')]['Query'].values.tolist())

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
    mEp332_cIoperonic_homolog_system_IDs = Prot = all_prots.intersection(mEp332_cIoperonic_homolog_systems)
    
    if mEp332_cIoperonic_homolog_system_IDs:
        mEp332_cIoperonic_homolog_system = [(FTloc,Acc,ID,Prot)]

    mEpX1_Gp15_IDs = Prot = all_prots.intersection(mEpX1_Gp15s)
    
    if mEpX1_Gp15_IDs:
        mEpX1_Gp15 = [(FTloc,Acc,ID,Prot)]
        

    mEpX1_KilAN_IDs = Prot = all_prots.intersection(mEpX1_KilANs)
    
    if mEpX1_KilAN_IDs:
        mEpX1_KilAN = [(FTloc,Acc,ID,Prot)]

    mEpX1_Nun_IDs = Prot = all_prots.intersection(mEpX1_Nuns)
    
    if mEpX1_Nun_IDs:
        mEpX1_Nun = [(FTloc,Acc,ID,Prot)]

    PhiR12_1_cIoperonic_system_IDs = Prot = all_prots.intersection(PhiR12_1_cIoperonic_systems)
    
    if PhiR12_1_cIoperonic_system_IDs:
        PhiR12_1_cIoperonic_system = [(FTloc,Acc,ID,Prot)]

    PhiR26_1_SNIPE_system_IDs = Prot = all_prots.intersection(PhiR26_1_SNIPE_systems)
    
    if PhiR26_1_SNIPE_system_IDs:
        PhiR26_1_SNIPE_system = [(FTloc,Acc,ID,Prot)]

    PhiR46_1_spacklelike_IDs = Prot = all_prots.intersection(PhiR46_1_spacklelikes)
    
    if PhiR46_1_spacklelike_IDs:
        PhiR46_1_spacklelike = [(FTloc,Acc,ID,Prot)]

    PhiR49_1_KilAN_IDs = Prot = all_prots.intersection(PhiR49_1_KilANs)
    
    if PhiR49_1_KilAN_IDs:
        PhiR49_1_KilAN = [(FTloc,Acc,ID,Prot)]

    return mEp332_cIoperonic_homolog_system, mEpX1_Gp15, mEpX1_KilAN, mEpX1_Nun, PhiR12_1_cIoperonic_system, PhiR26_1_SNIPE_system, PhiR46_1_spacklelike, PhiR49_1_KilAN 


# %%


all_FT_files = glob.glob('/home/gridsan/lbrenes/LaubLab_shared/all_refseq_feature_tables/*.txt')


# %%


results=Parallel(n_jobs=47)(delayed(checkFT)(x) for x in all_FT_files)


# %%


results = pd.DataFrame(results)
mEp332_cIoperonic_homolog_system_out = []
mEpX1_Gp15_out = []
mEpX1_KilAN_out = []
mEpX1_Nun_out = []
PhiR12_1_cIoperonic_system_out = []
PhiR26_1_SNIPE_system_out = []
PhiR46_1_spacklelike_out = []
PhiR49_1_KilAN_out = []
for index, row in results.iterrows():
    if row[0]:
        mEp332_cIoperonic_homolog_system_out = mEp332_cIoperonic_homolog_system_out + row[0]
    if row[1]:
        mEpX1_Gp15_out = mEpX1_Gp15_out + row[1]
    if row[2]:
        mEpX1_KilAN_out = mEpX1_KilAN_out + row[2]
    if row[3]:
        mEpX1_Nun_out = mEpX1_Nun_out + row[3]
    if row[4]:
        PhiR12_1_cIoperonic_system_out = PhiR12_1_cIoperonic_system_out + row[4]
    if row[5]:
        PhiR26_1_SNIPE_system_out = PhiR26_1_SNIPE_system_out + row[5]
    if row[6]:
        PhiR46_1_spacklelike_out = PhiR46_1_spacklelike_out + row[6]
    if row[7]:
        PhiR49_1_KilAN_out = PhiR49_1_KilAN_out + row[7]

mEp332_cIoperonic_homolog_system_out = pd.DataFrame(mEp332_cIoperonic_homolog_system_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
mEp332_cIoperonic_homolog_system_out.to_csv('TaxIDswmEp332_cIoperonic_homolog_system_nr_dom_filtered_prot.csv')

mEpX1_Gp15_out = pd.DataFrame(mEpX1_Gp15_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
mEpX1_Gp15_out.to_csv('TaxIDswmEpX1_Gp15_nr_dom_filtered_prot.csv')

mEpX1_KilAN_out = pd.DataFrame(mEpX1_KilAN_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
mEpX1_KilAN_out.to_csv('TaxIDswmEpX1_KilAN_LAX_nr_dom_filtered_prot.csv')

mEpX1_Nun_out = pd.DataFrame(mEpX1_Nun_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
mEpX1_Nun_out.to_csv('TaxIDswmEpX1_Nun_nr_dom_filtered_prot.csv')

PhiR12_1_cIoperonic_system_out = pd.DataFrame(PhiR12_1_cIoperonic_system_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
PhiR12_1_cIoperonic_system_out.to_csv('TaxIDswPhiR12_1_cIoperonic_system_nr_dom_filtered_prot.csv')

PhiR26_1_SNIPE_system_out = pd.DataFrame(PhiR26_1_SNIPE_system_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
PhiR26_1_SNIPE_system_out.to_csv('TaxIDswPhiR26_1_SNIPE_system_nr_dom_filtered_prot.csv')

PhiR46_1_spacklelike_out = pd.DataFrame(PhiR46_1_spacklelike_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
PhiR46_1_spacklelike_out.to_csv('TaxIDswPhiR46_1_spacklelike_nr_dom_filtered_prot.csv')

PhiR49_1_KilAN_out = pd.DataFrame(PhiR49_1_KilAN_out,columns = ['Featuretable','Genome','TaxID','Prot_ID'])
PhiR49_1_KilAN_out.to_csv('TaxIDswPhiR49_1_KilAN_LAX_nr_dom_filtered_prot.csv')



# %%




