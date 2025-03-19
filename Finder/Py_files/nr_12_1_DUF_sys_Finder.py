import pandas as pd
import numpy as np
import os
import glob
from joblib import Parallel, delayed

# Read in hmmer results
hmmerHits = pd.read_csv('/home/gridsan/lbrenes/hhpred/hhblts/hhmscan_multis_nr/multi_gene_systems_NR_hits.txt', sep=' ', comment='#', header=None, skipinitialspace=True, usecols=[0, 1, 2, 3, 4], names=['Domain', 'DomainAcc', 'Query', 'Empty', 'EVal'])

# Define the sets of protein IDs
PhiR12_1_DUF_system_DUF1327s = set(hmmerHits[hmmerHits['Domain'] == 'PhiR12.1_DUF_system_DUF1327.blts']['Query'].values.tolist())
PhiR12_1_DUF_system_DUF1823s = set(hmmerHits[hmmerHits['Domain'] == 'PhiR12.1_DUF_system_DUF1823.blts']['Query'].values.tolist())

# Assembly taxid and species taxid information for every genbank and refseq genome
assemblySums = pd.read_csv('/home/gridsan/lbrenes/hhpred/assembly_summary.txt', sep='\t', skiprows=1, usecols=['#assembly_accession', 'species_taxid'])
Acc2Taxid = {row[0]: row[1] for index, row in assemblySums.iterrows()}

# Function to check feature tables and find matches
def checkFT(FT_file):
    PhiR12_1_DUF_system = []
    
    FT = pd.read_csv(FT_file, sep='\t')
    FT = FT[FT['# feature'] == 'CDS'].reset_index()
    
    Acc = FT.at[0, 'assembly']
    ID = Acc2Taxid.get(Acc, 'Unknown')
    FTloc = os.path.basename(FT_file)
    
    all_prots = set(FT['non-redundant_refseq'])                        

    PhiR12_1_DUF_system_DUF1327_IDs = all_prots.intersection(PhiR12_1_DUF_system_DUF1327s)
    PhiR12_1_DUF_system_DUF1823_IDs = all_prots.intersection(PhiR12_1_DUF_system_DUF1823s)
    
    if PhiR12_1_DUF_system_DUF1327_IDs and PhiR12_1_DUF_system_DUF1823_IDs:
        PhiR12_1_DUF_system_DUF1327_indices = FT[FT['non-redundant_refseq'].isin(PhiR12_1_DUF_system_DUF1327_IDs)].index
        for PhiR12_1_DUF_system_DUF1327_index in PhiR12_1_DUF_system_DUF1327_indices:
            if any(FT.iloc[PhiR12_1_DUF_system_DUF1327_index-2:PhiR12_1_DUF_system_DUF1327_index+3]['non-redundant_refseq'].isin(PhiR12_1_DUF_system_DUF1823_IDs)):
                PhiR12_1_DUF_system.append((FTloc, Acc, ID, FT.at[PhiR12_1_DUF_system_DUF1327_index, 'non-redundant_refseq']))
            break
    
    return PhiR12_1_DUF_system


# Collect all feature table files
all_FT_files = glob.glob('/home/gridsan/lbrenes/LaubLab_shared/all_refseq_feature_tables/*.txt')

# Run checkFT function on the files in parallel
results = Parallel(n_jobs=47, verbose=10)(delayed(checkFT)(x) for x in all_FT_files)

# Flatten results into a list of tuples
PhiR12_1_DUF_system_out = [item for sublist in results for item in sublist]

# Filter out items that do not have exactly 4 elements and print warnings
filtered_results = []
for item in PhiR12_1_DUF_system_out:
    if len(item) != 4:
        print(f"Warning: Found item with length {len(item)} - {item}")
    else:
        filtered_results.append(item)

# Convert to DataFrame
PhiR12_1_DUF_system_df = pd.DataFrame(filtered_results, columns=['Featuretable', 'Genome', 'TaxID', 'Prot_ID'])
PhiR12_1_DUF_system_df.to_csv('TaxIDs_nr_PhiR12_1_DUF_system.csv', index=False)
