import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from uuid import getnode
from keras.models import Model, load_model
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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
METADATA = METADATA.drop(columns=['row_id'])
FEATURES = FEATURES.drop(columns=['row_id'])

print('METADATA',METADATA.shape)
print('FEATURES',FEATURES.shape)

#
# PCA
#

TABLE = FEATURES

TABLE = StandardScaler().fit(FEATURES).transform(TABLE)
TABLE = PCA(n_components=2, whiten=True).fit_transform(TABLE)
TABLE = pd.DataFrame(data=TABLE, columns=['PC1','PC2'])

print('TABLE',TABLE.shape)

for col_id, col in enumerate(TABLE.columns.tolist()):
	print(col_id, col)

#
# LABEL
#

TABLE = pd.concat([TABLE, METADATA], axis=1)

print(TABLE.head())

#
# Save
#

DATA = []
xvar = 'PC1'
yvar = 'PC2'
factor_name = 'mitosis'

factor = np.unique(TABLE[factor_name].values)

print(factor)

for fac in factor:

	print(fac)

	ids = TABLE[factor_name] == fac

	ncells = ids.sum()

	subDATA = dict()
	subDATA['x'] = TABLE[xvar][ids].values.tolist()
	subDATA['y'] = TABLE[yvar][ids].values.tolist()
	subDATA['label'] = [np.str(fac)] * ncells
	subDATA['name']  = np.str(fac)
	subDATA['id']  = 'id'#TABLE['cell_id'][ids].values.tolist()

	DATA.append(subDATA)

#
# Save JSON
#

with open('../../static/results/features.json','w') as fj:
	json.dump({'DATA': DATA, 'LABS': [xvar,yvar]}, fj, indent=4, sort_keys=True)

os.system('rm ../../static/results/analysis.jpg')
print('Rscript cond.R fix/free',xvar,yvar)

# #
# # Backup
# #

# #backupFiles('./featurex.py', 'featurex.py', '../../static/results/config-'+suffix+'.json')
