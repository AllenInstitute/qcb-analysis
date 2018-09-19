import numpy as np
import pandas as pd
import datasetdatabase as dsdb

# [ 2] 18339 - MitoEval20180807
# [ 4]  2436 - mito predictions                     (Nuclear lamin B1)
# [ 5]  3839 - 3451d1a1-0a99-4265-b170-106d3de0a46  (Cell and DNA features)
# [12]     7 - MitoEval20180821                     (Golgi)
# [15]  4481 - 51e405e2-5617-4b70-8f8e-7b999316f677 (Meta)
# [16]  4000 - a5bd9d48-27bb-4769-aa03-46f549a6db20
# [17]   230 - 6c0d3484-f0ba-4b86-aa35-77de7c1c16d2 (Random)
# [18]  4000 - 2701363a-3783-4fc3-9bc9-d2cfa737c3bf
# [19]  2000 - 63ca65dd-ac7e-4641-a460-0db36344b5c9
# [24]  1000 - 8f7a51bf-9cd4-4148-a838-d5ecf5af605c
# [25]  2000 - 285f60b8-2476-4309-b889-e06741f8b7c6
# [26]  2000 - 110e8878-a2cc-4278-8112-0d0b2cc91419
#
# From [15] (Meta): FBL,ZO1,ACTB,ACTN1,MYOSIN,DSP,TUBA,LMNB,TOM20,SEC61


prod = dsdb.DatasetDatabase(config='//allen/aics/assay-dev/Analysis/QCB_database/prod_config.json')

ds = prod.get_dataset(12).ds

print(ds.describe())

ds.to_csv('features.csv')


