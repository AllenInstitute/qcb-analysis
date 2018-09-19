import os
import sys
import json
import keras
import numpy as np
from skimage import io
from keras.models import Model
from keras.utils import plot_model
from keras.callbacks import CSVLogger, ModelCheckpoint

sys.path.insert(0,'../aux/')
from auxFunctions import loadData, myGenerator, getModel_3D_ext, backupFiles

#
# Folders setup
#

Folder = '../../static/'

for dirname in ['imgs/','results/']:
	if not os.path.exists(Folder+dirname):
	    os.makedirs(Folder+dirname)

#
# Load JSON
#

with open('../../static/resources/deepheader.json','r') as fj:
	HEADER = json.load(fj)

Parameters = {'batch_size': 32, 'header': HEADER}

#
# Testing the image generator
#

χ, γ = next(myGenerator(Parameters,'train'))

for x_id, X in enumerate(χ):

	b, d, h, w, c = X.shape

	x = X[:,np.int(0.5*d),:,:,0]

	x = np.transpose(x,(1,0,2)).reshape(h,-1)

	io.imsave(Folder+'imgs/batch_example_'+str(x_id)+'.jpg', (0.5*(1+x)*255).astype(np.uint8))

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
# Callbacks
#

csv_logger = CSVLogger(Folder+'results/loss-'+args['suffix']+'.csv', separator=',', append=False)

chk_point = ModelCheckpoint(Folder+'results/model-'+args['suffix']+'.h5', monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=False, mode='min', period=1)

#
# Training
#

MODEL = getModel_3D_ext(Config=HEADER['head'])

MODEL.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

MODEL.summary()

plot_model(MODEL, to_file=Folder+'imgs/model-'+args['suffix']+'.png', show_shapes=True, show_layer_names=True)

MODEL.fit_generator(myGenerator(Parameters, mode='train'), steps_per_epoch=2, epochs=nepochs, validation_data=myGenerator(Parameters, mode='valid'), validation_steps=4, verbose=1, callbacks=[csv_logger,chk_point])

#
# Backup
#

backupFiles('./cnn3D.py', 'cnn3D.py', '../../static/results/config-'+suffix+'.json')
backupFiles('../aux/auxFunctions.py', 'auxFunctions.py', '../../static/results/config-'+suffix+'.json')
