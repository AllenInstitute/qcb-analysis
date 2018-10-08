import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import io as skio
from aicsimage import io, processing
from aicsfeature.extractor import fov

#
# Loading metadata
#

with open('tmp_meta_drug_ds.pkl', 'rb') as fp:
    df_meta = pickle.load(fp)

#
# Creating a stack data frame
#

df_meta['czi_folder'] = ''
czi_formats = {'2':12, '3': 10}
for fs, fd in czi_formats.items():
    valids = df_meta['czi_filename'].str.startswith(fs)
    df_meta.loc[valids,'czi_folder'] = df_meta.loc[valids,'czi_filename'].str[0:fd]

valid_suffix = ['_seg_coarse.ome.tif','_seg_default.ome.tif','.czi_seg_default.ome.tif']
base_folder = '/allen/aics/assay-dev/Analysis/QCB_DRUG_database/structure_segmentation/'
df_stack = df_meta[['czi_filename', 'czi_folder']].copy()
df_stack = df_stack.drop_duplicates().reset_index(drop=True)
df_stack['suffix'] = None

for idx, row in df_stack.iterrows():
    for suffix in valid_suffix:
        full_path = os.path.join(base_folder,row.czi_folder,row.czi_filename.replace('.czi',suffix))
        if os.path.isfile(full_path):
            df_stack.loc[idx,'suffix'] = suffix
    if df_stack.suffix[idx] == None:
        print('[fail]',full_path)
        break

#
# Looping over all stacks
#

tablef = pd.DataFrame([])
for sid in tqdm(range(df_stack.shape[0])):
    seg_path = os.path.join(base_folder,df_stack.czi_folder[sid],df_stack.czi_filename[sid].replace('.czi',df_stack.suffix[sid]))
    SEG = skio.imread(seg_path)
    tablef = tablef.append(fov.GetFeatures(None,seg=SEG), ignore_index=True)

tablef = tablef.join(df_stack['czi_filename'])

#
# Save as pickle
#

with open('tmp_fov_drug_ds.pkl', 'wb') as fp:
    pickle.dump(tablef,fp)

