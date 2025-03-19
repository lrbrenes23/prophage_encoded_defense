import pandas as pd
import numpy as np
import os
import glob
from joblib import Parallel, delayed

def gethhsrch(faa):
    base_name = os.path.basename(faa)
    gene_oid = os.path.splitext(base_name)[0]
    os.system(f"hhsearch -i {faa} -d ~/hhpred/phrogs_curated/phrogs -o ~/hhpred/neighbors/mEp506_small_sys_DUF2594_system/hhsrch/{gene_oid}.hhsrch.hhr -p 20 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000")


# Get all .faa files recursively
all_faa_files = glob.glob('/home/gridsan/lbrenes/hhpred/neighbors/mEp506_small_sys_DUF2594_system/fastas/*')


# Run hhsrch function on the files in parallel
results = Parallel(n_jobs=47, verbose=10)(delayed(gethhsrch)(x) for x in all_faa_files)
