import os
import sys
import json
import keras
import argparse
import numpy as np
from uuid import getnode
from keras.models import Model, load_model
from sklearn.manifold import TSNE

sys.path.insert(0,'../aux/')
from auxFunctions import loadData, myGenerator, backupFiles

#
# Arguments
#

parser = argparse.ArgumentParser(description='Multimodal Epilepsy Classifier')
parser.add_argument('-s','--suffix', help='Model suffix', required=True)
args = vars(parser.parse_args())

#
# Load model
#

MODEL = load_model('../../static/results/model-'+args['suffix']+'.h5')

MODEL = Model(inputs=MODEL.input, outputs=MODEL.get_layer('dense_1').output)

MODEL.summary()

#
# Feature extraction
#

with open('../../static/resources/deepheader.json','r') as fj:
	HEADER = json.load(fj)

DATA = []
nsamples = 0

for group in HEADER['body']:

	for mode in ['train','valid']:

		for sample in group[mode]:

			X, Y = loadData(Config=HEADER['head'], Instances=[sample])

			Z = MODEL.predict(X[1])

			DATA.append(Z)

			nsamples = nsamples + 1

DATA = np.array(DATA).reshape(nsamples,-1)

#
# TSNE
#

tsne = TSNE(n_components=2)
tsne_result = tsne.fit_transform(DATA)

# print(tsne_result.shape)

#
# Save
#

# tsne_result = np.random.rand(735,2)

PROJ = []

nsamples = 0

for mode in ['train','valid']:

	for group in HEADER['body']:

		files = []
		x_coo = []
		y_coo = []
		label = []

		for sample in group[mode]:

			files.append(os.path.basename(sample['file'][0]).replace('_HR','').replace('.tif',''))
			x_coo.append(tsne_result[nsamples,0].tolist())
			y_coo.append(tsne_result[nsamples,1].tolist())
			label.append(group['class']+'-'+mode)
		
			nsamples = nsamples + 1

		DATA = dict()
		DATA['files'] = files
		DATA['x_coo'] = x_coo
		DATA['y_coo'] = y_coo
		DATA['label'] = label
		DATA['name']  = group['class']+'-'+mode

		PROJ.append(DATA)

#
# Save JSON
#

with open('../../static/results/features-'+args['suffix']+'.json','w') as fj:
	json.dump(PROJ, fj, indent=4, sort_keys=True)

#
# Backup
#

backupFiles('./featurex.py', 'featurex.py', '../../static/results/config-'+suffix+'.json')
