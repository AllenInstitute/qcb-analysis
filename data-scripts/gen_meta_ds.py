import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import io as skio
import datasetdatabase as dsdb
from aicsimage import io, processing

#
# Multiprocessing
#

os.environ["DSDB_PROCESS_LIMIT"] = "16"

#
# Loading metadata
#

prod = dsdb.DatasetDatabase(config='/allen/aics/assay-dev/Analysis/QCB_database/prod_config.json')
ds_meta = prod.get_dataset(name='QCB_drug_cell_meta')

with open('tmp_meta_drug_ds.pkl', 'wb') as fp:
    pickle.dump(ds_meta.ds,fp)

