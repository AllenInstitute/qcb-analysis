import os
import sys
import json
import glob
import getopt
import numpy as np

sys.path.insert(0,'../aux/')
from auxFunctions import backupFiles


opts, _ = getopt.getopt(sys.argv[1:],'p:s:S:',['path=','split=','suffix='])

for opt, arg in opts:
	if opt in ('-p','--path'):
		rootfolder = arg
	if opt in ('-S','--suffix'):
		suffix = arg
	if opt in ('-s','--split'):
		split = np.int(arg)

#
# List all available directories
#

Directories = [ x for x in os.listdir(rootfolder) if os.path.isdir(os.path.join(rootfolder,x)) ]

#
# Number of samples inside each directory
#

Classes = []

total_samples = 0

for d in Directories:

	Filenames = os.listdir(os.path.join(rootfolder,d))

	# We might have to test here if the image is corrupted

	samples = len(Filenames)

	total_samples += samples

	Classes = np.append(Classes, d)

#
# Config
#

with open("config.json","r") as fj:
	CONFIG = json.load(fj)

#
# For each file
#

BODY = []

for d in Directories:

	Data = []

	target = np.int(np.where(Classes==d)[0][0])

	Filenames = glob.glob(os.path.join(rootfolder,d,CONFIG['data']['prefix'][0])+'*')

	np.random.shuffle(Filenames)

	for f in Filenames:

		fnames = []
		ftypes = []
		fdimen = []

		for data_id in range(CONFIG['ntypes']):

			fname = f.replace(CONFIG['data']['prefix'][0],CONFIG['data']['prefix'][data_id])
			fnames.append(os.path.join(rootfolder,d,fname))
			ftypes.append(CONFIG['data']['type'][data_id])
			fdimen.append(CONFIG['data']['dim'][data_id])

		Data = np.append(Data,{'file': fnames, 'type': ftypes, 'size': fdimen, 'target': target})

	valid_samples = np.random.choice(a=len(Data), size=split, replace=False).tolist()

	train_samples = list(set(np.array(range(len(Data))))-set(valid_samples))

	DataTrain = Data[train_samples].tolist()
	
	DataValid = Data[valid_samples].tolist()

	BODY = np.append(BODY,{'class': d, 'weight': 1-len(Filenames)/total_samples, 'train': DataTrain, 'valid': DataValid})

BODY = BODY.tolist()

HEADER = {'head': CONFIG, 'body': BODY}

#
# Save JSON
#

with open('../../static/resources/deepheader.json','w') as fj:
	json.dump(HEADER, fj, indent=4, sort_keys=True)

#
# Backup
#

backupFiles('./config.json', 'config.json', '../../static/results/config-'+suffix+'.json')
backupFiles('../../static/resources/deepheader.json', 'deepheader.json', '../../static/results/config-'+suffix+'.json')
backupFiles('./'+os.path.basename(__file__), os.path.basename(__file__), '../../static/results/config-'+suffix+'.json')