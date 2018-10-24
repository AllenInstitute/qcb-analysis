import os
import json
import keras
import numpy as np
from skimage import io
from keras.models import Model, load_model
from keras.layers import Conv2D, MaxPooling2D, GlobalMaxPooling2D
from keras.layers import Conv3D, MaxPooling3D, GlobalMaxPooling3D
from keras.layers import Input, Dense, Dropout, Concatenate, Lambda, Multiply
from keras.preprocessing.image import ImageDataGenerator
from keras.utils import plot_model
from keras import backend as K
from keras.callbacks import CSVLogger
from keras.optimizers import SGD
import tensorflow as tf

#
# Data Augmentation
#

AugmentationEngine = ImageDataGenerator(featurewise_center=False,
	samplewise_center=False,
	featurewise_std_normalization=False,
	samplewise_std_normalization=False,
	zca_whitening=False,
	rotation_range=5.0,
	width_shift_range=0.03,
	height_shift_range=0.03,
	shear_range=0.0,
	zoom_range=0.1,
	channel_shift_range=0.0,
	fill_mode='constant',
	cval=0.0,
	horizontal_flip=False,
	vertical_flip=False)

#
# Loading function
#

def loadData(Config, Instances):

	#
	# Template for input data
	#

	γ = []
	χ = [[]]

	for i in range(1,Config['ntypes']):

		χ.append([])


	#
	# Load input according to its type
	#

	for type_id in range(0,Config['ntypes']):

		batch = []

		for sample_id, sample in enumerate(Instances):

			fname = sample['file'][type_id]

			if sample['type'][type_id] == 'image2d':

				ζ = io.imread(fname)

			elif sample['type'][type_id] == 'image3d':

				ζ = io.imread(fname)
			
			else:
			
				print('Reader for this data type has not been implemented.')

			batch.append(ζ)

		batch = np.array(batch, dtype=np.float32).reshape(-1,*Config['data']['dim'][type_id],1)

		χ[type_id] = batch

	for sample_id, sample in enumerate(Instances):

		γ.append(sample['target'])
	
	γ = np.array(γ).reshape(-1,1)
	
	return χ, γ

#
# Custom Generator
#

def myGenerator(Parameters, mode):

	n_samples_per_group = np.int(Parameters['batch_size'] / len(Parameters['header']['body']))

	while True:

		Instances = []

		for group in Parameters['header']['body']:

			samples = np.random.choice(a=len(group[mode]), size=n_samples_per_group, replace=False)

			for s in samples:

				Instances.append(group[mode][s])

		np.random.shuffle(Instances)

		χ, γ = loadData(Parameters['header']['head'], Instances)

		#
		# Normalization and noise
		#

		for x_id, X in enumerate(χ):

			# X *= 0.5
			# X += 0.5

			# W = np.random.choice([0, 1], size=X.shape, p=[0.1,0.9])

			# X *= W

			#
			# Data augmentation goes here
			#

			# AugmentationEngine.fit(X)

			#Xt, _ = AugmentationEngine.flow(X, γ, seed=171, shuffle=False).next()

			χ[x_id] = np.copy(X)

		yield χ, γ

#
# Backup files
#

def backupFiles(fname, flabel, fbackup):

	with open(fname, 'r') as fc:
		content = fc.read()

	try:
		with open(fbackup,'r') as fb:
			backup = json.load(fb)
		backup[flabel] = content
		with open(fbackup,'w') as fb:
			json.dump(backup, fb, indent=4, sort_keys=True)
	except:
		with open(fbackup,'w') as fb:
			json.dump({flabel: content}, fb, indent=4, sort_keys=True)

#
# Models
#

def getModel_2D():

	χ1 = Input((240, 180, 4))

	φ1 = Conv2D(filters=16, kernel_size=(3,3), strides=(2,2), activation='relu')(χ1)
	φ1 = Conv2D(filters=32, kernel_size=(3,3), strides=(2,2), activation='relu')(φ1)
	φ1 = MaxPooling2D(pool_size=(2,2))(φ1)
	φ1 = Dropout(rate=0.3)(φ1)

	φ1 = Conv2D(filters=32, kernel_size=(3,3), strides=(2,2), activation='relu')(φ1)
	φ1 = Conv2D(filters=32, kernel_size=(3,3), strides=(2,2), activation='relu')(φ1)
	φ1 = MaxPooling2D(pool_size=(2,2))(φ1)
	φ1 = Dropout(rate=0.3)(φ1)

	φ1 = GlobalMaxPooling2D()(φ1)

	γ2 = Dense(units=64,  activation='relu')(φ1)

	γ1 = Dense(units=1,  activation='sigmoid')(γ1)

	MODEL = Model(inputs=[χ1],outputs=[γ1])

	return MODEL

def getModel_3D(Config):

	x1_dim = Config['data']['dim'][0]
	x2_dim = Config['data']['dim'][1]

	χ1 = Input((*x1_dim,1))
	χ2 = Input((*x2_dim,1))

	φ2 = Conv3D(filters=16, kernel_size=(3,3,3), strides=(2,2,2), activation='relu')(χ2)
	φ2 = Conv3D(filters=32, kernel_size=(3,3,3), strides=(1,1,1), activation='relu')(φ2)
	φ2 = MaxPooling3D(pool_size=(2,2,2))(φ2)
	φ2 = Dropout(rate=0.3)(φ2)

	φ2 = Conv3D(filters=32, kernel_size=(3,3,3), strides=(2,2,2), activation='relu')(φ2)
	φ2 = Conv3D(filters=32, kernel_size=(3,3,3), strides=(1,1,1), activation='relu')(φ2)
	φ2 = MaxPooling3D(pool_size=(2,2,2))(φ2)
	φ2 = Dropout(rate=0.3)(φ2)

	φ2 = GlobalMaxPooling3D()(φ2)

	γ2 = Dense(units=64,  activation='relu')(φ2)

	γ2 = Dense(units=1,  activation='sigmoid')(γ2)

	MODEL = Model(inputs=[χ1,χ2],outputs=[γ2])

	return MODEL


def getModel_3D_ext(Config):

	x1_dim = Config['data']['dim'][0]
	x2_dim = Config['data']['dim'][1]

	χ1 = Input((*x1_dim,1))
	χ2 = Input((*x2_dim,1))

	φ2 = Conv3D(filters=16, kernel_size=(3,3,3), strides=(2,2,2), activation='relu')(χ2)
	φ2 = Conv3D(filters=16, kernel_size=(3,3,3), strides=(1,1,1), activation='relu')(φ2)
	φ2 = MaxPooling3D(pool_size=(2,2,2))(φ2)
	φ2 = Dropout(rate=0.3)(φ2)

	φ2 = Conv3D(filters=32, kernel_size=(3,3,3), strides=(2,2,2), activation='relu')(φ2)
	φ2 = Conv3D(filters=32, kernel_size=(3,3,3), strides=(1,1,1), activation='relu')(φ2)
	φ2 = MaxPooling3D(pool_size=(2,2,2))(φ2)
	φ2 = Dropout(rate=0.3)(φ2)

	φ2 = GlobalMaxPooling3D()(φ2)

	γ2 = Dense(units=x1_dim[1],  activation='softmax')(φ2)

	m = Lambda(lambda x: tf.transpose(x,[0,4,3,2,1]))(χ1)
	m = Multiply()([m,γ2])
	m = Lambda(lambda x: tf.transpose(x,[0,4,3,2,1]))(m)
	χ3 = Lambda(lambda x: tf.reduce_sum(x,1))(m)

	φ3 = Conv2D(filters=16, kernel_size=(3,3), strides=(2,2), activation='relu')(χ3)
	φ3 = Conv2D(filters=16, kernel_size=(3,3), strides=(2,2), activation='relu')(φ3)
	φ3 = MaxPooling2D(pool_size=(2,2))(φ3)
	φ3 = Dropout(rate=0.3)(φ3)

	φ3 = Conv2D(filters=32, kernel_size=(3,3), strides=(2,2), activation='relu')(φ3)
	φ3 = Conv2D(filters=32, kernel_size=(3,3), strides=(2,2), activation='relu')(φ3)
	φ3 = MaxPooling2D(pool_size=(2,2))(φ3)
	φ3 = Dropout(rate=0.3)(φ3)

	φ3 = GlobalMaxPooling2D()(φ3)

	γ3 = Dense(units=256,  activation='relu')(φ3)
	γ3 = Dense(units=128,  activation='relu')(γ3)

	γ3 = Dense(units=1,  activation='sigmoid')(γ3)

	MODEL = Model(inputs=[χ1,χ2],outputs=[γ3])

	return MODEL


