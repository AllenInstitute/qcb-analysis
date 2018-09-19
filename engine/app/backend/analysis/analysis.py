import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from uuid import getnode
from keras.models import Model, load_model
from sklearn.manifold import TSNE

sys.path.insert(0,'../aux/')
from auxFunctions import loadData, myGenerator, backupFiles

pd.set_option('display.max_columns', None)
np.set_printoptions(threshold=np.nan)

#
# Arguments
#

parser = argparse.ArgumentParser(description='Exploration of QCB data')
parser.add_argument('-p','--path', help='Path to tables', required=True)
args = vars(parser.parse_args())

#
# Feature extraction
#

METADATA = pd.read_csv(os.path.join(args['path'],'meta.csv'))
FEATURES = pd.read_csv(os.path.join(args['path'],'features.csv'))

print('METADATA',METADATA.shape)
print('FEATURES',FEATURES.shape)

TABLE = METADATA.merge(FEATURES, on='cell_id', how='inner')

print('TABLE',TABLE.shape)

for col_id, col in enumerate(TABLE.columns.tolist()):
	print(col_id, col)

TABLE = TABLE[TABLE['mitosis']==0]

#
# Save
#

DATA = []

factor_name = 'structure_name'

factor = np.unique(TABLE[factor_name].values)

print(factor)

for fac in factor:

	print(fac)

	ids = TABLE[factor_name] == fac

	ncells = ids.sum()

	subDATA = dict()
	subDATA['x'] = TABLE['cell_seg_vol'][ids].values.tolist()
	subDATA['y'] = TABLE['dna_seg_vol'][ids].values.tolist()
	subDATA['label'] = [fac] * ncells
	subDATA['name']  = fac
	subDATA['id']  = TABLE['cell_id'][ids].values.tolist()

	DATA.append(subDATA)

#
# Save JSON
#

with open('../../static/results/features.json','w') as fj:
	json.dump(DATA, fj, indent=4, sort_keys=True)

#
# Backup
#

#backupFiles('./featurex.py', 'featurex.py', '../../static/results/config-'+suffix+'.json')
