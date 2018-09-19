#
# Load data using the header created by deepheader
#

import numpy as np
from skimage import io, transform

def loaddata(Instances,Parameters):

	Target = []
	BatchIM2D_1 = []
	BatchIM2D_2 = []
	BatchAR1D_1 = []

	#
	# Number of classes from the JSON
	#

	nclasses = len(Parameters['header'])

	#
	# Load data according to its type
	#

	for sample in Instances:

		img = io.imread(sample['file'][0])

		BatchIM2D_1 = np.append(BatchIM2D_1,img)

		img = io.imread(sample['file'][1])

		BatchIM2D_2 = np.append(BatchIM2D_2,img)

		arr = np.loadtxt(sample['file'][2])

		BatchAR1D_1 = np.append(BatchAR1D_1,arr)

		target = sample['target']

		Target = np.append(Target,target)

	BatchIM2D_1 = BatchIM2D_1.reshape(-1,240,180,4)

	BatchIM2D_2 = BatchIM2D_2.reshape(-1,510,510,3)

	BatchAR1D_1 = BatchAR1D_1.reshape(-1,128,1)

	Target = Target.reshape(-1,1)

	return BatchIM2D_1, BatchIM2D_2, BatchAR1D_1, Target