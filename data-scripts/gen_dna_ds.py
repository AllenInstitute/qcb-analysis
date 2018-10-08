import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import io as skio
from aicsimage import io, processing
from aicsfeature.extractor import dna

#
# Loading metadata
#

with open('tmp_meta_drug_ds.pkl', 'rb') as fp:
    df_meta = pickle.load(fp)

#
# Looping over all cells
#

rootf = '/allen/aics/assay-dev/Analysis/QCB_DRUG_database/cell_info'

tablef = pd.DataFrame([])
for cid in tqdm(range(df_meta.shape[0])):
    seg_path = os.path.join(rootf,df_meta.cell_id[cid],'seg.ome.tif')
    raw_path = os.path.join(rootf,df_meta.cell_id[cid],'raw.ome.tif')
    SEG = skio.imread(seg_path)
    RAW = skio.imread(raw_path)
    tablef = tablef.append(dna.GetFeatures(None,seg=RAW[1,:,:,:]*SEG[1,:,:,:]), ignore_index=True)

#
# Save as pickle
#

with open('tmp_dna_drug_ds.pkl', 'wb') as fp:
    pickle.dump(tablef,fp)

