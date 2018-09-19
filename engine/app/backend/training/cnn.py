import os
import sys
import json
import keras
import numpy as np
from skimage import io
from uuid import getnode
from keras.models import Model
from keras.utils import plot_model
from keras.callbacks import CSVLogger

sys.path.insert(0,'../aux/')
from auxFunctions import loadData, myGenerator, getModel

#
# Folders setup
#

mac = getnode()
print('Mac Address:',mac)
if mac == 279278385726721 or mac == 154505274806339:
	Folder = '../../static/'
else:
	Folder = "./"

for dirname in ['imgs/','results/']:
	if not os.path.exists(Folder+dirname):
	    os.makedirs(Folder+dirname)

print('Folder:', Folder)


#
# Load JSON
#

with open('../../static/resources/deepheader.json','r') as fj:
	HEADER = json.load(fj)

Parameters = {'batch_size': 32, 'header': HEADER}

W = []
for group_id, group in enumerate(HEADER):
	W.append(group['weight'])
W = 1.0 - np.array(W)
W /= np.min(1-W)

WEIGHTS = dict()
for group_id, group in enumerate(HEADER):
	WEIGHTS[group_id] = W[group_id]

X, Y = next(myGenerator(Parameters,'train'))

X = X[:,:,:,0]

n, h, w = X.shape

X = np.transpose(X,(0,2,1))

io.imsave(Folder+'imgs/batch_example.jpg', np.transpose(X.reshape(-1,h)))

#
# External aurguments
#

import argparse
parser = argparse.ArgumentParser(description='Multimodal Epilepsy Classifier')
parser.add_argument('-m','--mode',   help='Available modes: train, test and check', required=True)
parser.add_argument('-s','--suffix', help='Model suffix', required=True)
parser.add_argument('-e','--epochs', help='Number of epochs', required=True)
args = vars(parser.parse_args())

nepochs = np.int(args['epochs'])

#
# CSV Logger
#

csv_logger = CSVLogger(Folder+'results/loss-'+args['suffix']+'.csv', separator=',', append=False)

#
# Training
#

MODEL = getModel()

MODEL.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

MODEL.summary()

try:
	plot_model(MODEL, to_file=Folder+'imgs/model-'+args['suffix']+'.png', show_shapes=True, show_layer_names=True)
except: pass

MODEL.fit_generator(myGenerator(Parameters, mode='train'), steps_per_epoch=16, epochs=nepochs, validation_data=myGenerator(Parameters, mode='valid'), validation_steps=4, verbose=1, callbacks=[csv_logger])

#
# Save trained model
#

MODEL.save(Folder+'results/model-'+args['suffix']+'.h5')

