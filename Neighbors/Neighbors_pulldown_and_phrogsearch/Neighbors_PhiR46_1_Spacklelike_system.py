import pandas as pd
import numpy as np
import os
import glob
from joblib import Parallel, delayed



#Read in hmmer results 
hmmerHits = pd.read_csv('/home/gridsan/lbrenes/hhpred/hhblts/hhmscan_singles_nr/single_gene_systems_NR_hits_filtered_dom_nodupes.csv',usecols=[2,3,5],names=['Domain','DomainAcc','Query'],header=0)
hmmerHits

PhiR46_1_spacklelikes = set(hmmerHits[(hmmerHits['Domain'] == 'PhiR46.1_spacklelike.blts')]['Query'].values.tolist())

# Step 1: Read the DataFrame
PhiR46_1_spacklelike_hits=pd.read_csv('/home/gridsan/lbrenes/hhpred/taxids/single_gene_systems_nr/TaxIDswPhiR46_1_spacklelike_nr_dom_filtered_prot.csv')

# Step 2: Make a dictionary connection between the FT and Prot ID column values
FT_Prot_ID_match = dict(zip(PhiR46_1_spacklelike_hits['Featuretable'], PhiR46_1_spacklelike_hits['Prot_ID']))

# Step 3: Define the directory path
all_FT_files_dir = '/home/gridsan/lbrenes/LaubLab_shared/all_refseq_feature_tables/'

# Step 4: Initialize a dictionary to hold file-protein ID matches
PhiR46_1_spacklelike_matches = {}

# Step 5: Iterate over each file in the directory and match with DataFrame
for Featuretable in os.listdir(all_FT_files_dir):
    if Featuretable in FT_Prot_ID_match:
        Prot_ID = FT_Prot_ID_match[Featuretable]
        PhiR46_1_spacklelike_matches[Featuretable] = Prot_ID

def neighbors(directory, matches):
    # Initialize a dictionary to store neighboring rows for each matched protein ID
    neighbors_dict = {}

    # Process each entry in the matches dictionary
    for Featuretable, Prot_ID in PhiR46_1_spacklelike_matches.items():
        # Construct the full path to the feature table file
        FT_file_path = os.path.join(directory, Featuretable)
        
        # Check if the file exists
        if os.path.isfile(FT_file_path):
            # Read the feature table file into a DataFrame
            FT = pd.read_csv(FT_file_path, sep='\t')
            FT = FT[FT['# feature'] == 'CDS'].reset_index()
            FT = FT[FT['class'] == 'with_protein'].reset_index()

            all_prots = set(FT['non-redundant_refseq'])                        

            PhiR46_1_spacklelike_system_IDs = all_prots.intersection(PhiR46_1_spacklelikes)

            # Ensure 'non-redundant_refseq' column exists
            if 'non-redundant_refseq' not in FT.columns:
                raise ValueError(f"The 'non-redundant_refseq' column is missing from the file: {Featuretable}")
            
            # Find rows where the Prot_ID matches in the 'non-redundant_refseq' column
            matched_rows = FT[FT['non-redundant_refseq'] == Prot_ID]
            
            if not matched_rows.empty:
                # Extract the genomic_accession value for matched rows
                genomic_accession_value = matched_rows['genomic_accession'].iloc[0]

                # Get the integer-based index positions of the matched rows
                matched_indices = matched_rows.index

                # Retrieve neighboring rows for each occurrence
                for idx in matched_indices:
                    # Define the range for neighbors (20 rows above and below)
                    start_idx = max(0, idx - 20)
                    end_idx = min(len(FT) - 1, idx + 20)
            
                    # Extract the rows in the range
                    neighbor_rows = FT.iloc[start_idx:end_idx + 1]

                    # Filter neighbor_rows to include only those with the same genomic_accession
                    neighbor_rows = neighbor_rows[neighbor_rows['genomic_accession'] == genomic_accession_value]

                    neighbor_rows = neighbor_rows[~neighbor_rows['non-redundant_refseq'].isin(PhiR46_1_spacklelike_system_IDs)]
                    # Only keep neighbor rows if there are at least 10
                    if len(neighbor_rows) >= 10:
                        # Use a composite key based on both protein_id and file_name
                        key = (Prot_ID, Featuretable)

                        # Combine results into the dictionary
                        if key in neighbors_dict:
                            neighbors_dict[key] = pd.concat([neighbors_dict[key], neighbor_rows]).drop_duplicates()
                        else:
                            neighbors_dict[key] = neighbor_rows

    return neighbors_dict

def save_results(results):
    # Ensure the output directory exists
    output_dir = '/home/gridsan/lbrenes/hhpred/neighbors/PhiR46_1_spacklelike_system'
    os.makedirs(output_dir, exist_ok=True)

    # Initialize a list to collect all Prot_IDs for the summary CSV
    refseq_summary = []
    all_rows = []

    # Save each set of neighbor rows to a separate CSV file
    for Prot_ID, neighbor_rows in results.items():
        output_file = os.path.join(output_dir, f'{Prot_ID}_neighbors.csv')
        neighbor_rows.to_csv(output_file, index=False)
        
        # Collect non-redundant_refseq for each row
        refseq_summary.extend(neighbor_rows['non-redundant_refseq'].tolist())
        # Collect all rows for detailed summary
        all_rows.append(neighbor_rows)

    # Create a DataFrame for refseq summary
    refseq_df = pd.DataFrame({'non_redundant_refseq': refseq_summary})

    # Save the refseq summary to a CSV file
    refseq_summary_file = os.path.join(output_dir, 'PhiR46_1_spacklelike_system_neighbors_prot_refseq_summary.csv')
    refseq_df.to_csv(refseq_summary_file, index=False)
  
    # Combine all rows into a single DataFrame
    all_rows_df = pd.concat(all_rows, ignore_index=True)

    # Save the detailed summary to a CSV file
    detailed_summary_file = os.path.join(output_dir, 'PhiR46_1_spacklelike_system_neighbors_detailed_summary.csv')
    all_rows_df.to_csv(detailed_summary_file, index=False)

    

# Define the directory containing the feature table files
directory = '/home/gridsan/lbrenes/LaubLab_shared/all_refseq_feature_tables/'

matches = PhiR46_1_spacklelike_matches

results = neighbors(directory, matches)

# Save results to files
save_results(results)
