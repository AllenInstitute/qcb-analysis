import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import io as skio
import datasetdatabase as dsdb
from aicsimage import io, processing
from aicsfeature.extractor import structure

rootf = '/allen/aics/assay-dev/Analysis/QCB_DRUG_database/cell_info'

#
# Loading metadata
#

with open('tmp_meta_drug_ds.pkl', 'rb') as fp:
    df_meta = pickle.load(fp)

#
# Structures to be analyzed
#

Struct_To_Process = [
        {'df_name': 'tubulin', 'ex_name': 'er', 'save_as': 'tmp_str_tubulin_drug_ds.pkl'},
        {'df_name': 'sec61b', 'ex_name': 'er', 'save_as': 'tmp_str_er_drug_ds.pkl'},
        {'df_name': 'golgi', 'ex_name': 'golgi', 'save_as': 'tmp_str_golgi_drug_ds.pkl'}]

#
# Looping over all cells
#

df_meta = df_meta.set_index('cell_id')

for struct in Struct_To_Process:

    print('structure:', struct['df_name'])

    df_meta_struct = df_meta.loc[df_meta.structure_name==struct['df_name']]

    tablef = pd.DataFrame([])
    for row in tqdm(range(df_meta_struct.shape[0])):
        cid = df_meta_struct.index[row]
        seg_path = os.path.join(rootf,cid,'seg.ome.tif')
        raw_path = os.path.join(rootf,cid,'raw.ome.tif')
        SEG = skio.imread(seg_path)
        RAW = skio.imread(raw_path)
        tablef = tablef.append(structure.GetFeatures(None,seg=RAW[2,:,:,:]*SEG[2,:,:,:], structure_name=struct['ex_name']), ignore_index=True)

    tablef.index = df_meta_struct.index

    #
    # Save as pickle
    #

    with open(struct['save_as'], 'wb') as fp:
        pickle.dump(tablef,fp)

