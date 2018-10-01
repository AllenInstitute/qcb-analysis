import pickle
import numpy as np
import pandas as pd
import datasetdatabase as dsdb

#
# Load meta
#

with open('../data-raw/tmp_ds_meta.pkl', 'rb') as fp:
	df_meta = pickle.load(fp)

print(df_meta.head())

df_meta.to_csv('../engine/data-processed/tmp-meta.csv')

#
# Load mem and dna
#

with open('../data-raw/tmp_ds_mem.pkl', 'rb') as fp:
	df_mem = pickle.load(fp)

with open('../data-raw/tmp_ds_dna.pkl', 'rb') as fp:
	df_dna = pickle.load(fp)

df = pd.concat([df_mem.reset_index(drop=True), df_dna], axis=1)

print(df.head())

df.to_csv('../engine/data-processed/tmp-feature.csv')


