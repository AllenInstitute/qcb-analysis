import os
import sys
import json
import pickle
import argparse
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

sys.path.insert(0,"../aux/")
from auxFunctions import backupFiles

pd.set_option("display.max_columns", None)
np.set_printoptions(threshold=np.nan)

#
# Arguments
#

parser = argparse.ArgumentParser(description="Exploration of QCB data")
parser.add_argument("-m","--meta", help="Path to metadata df", required=True)
parser.add_argument("-1","--df1", help="Path to df 1", required=True)
parser.add_argument("-2","--df2", help="Path to df 2", required=True)
parser.add_argument("-x","--x", help="variable for x axis", required=True)
parser.add_argument("-y","--y", help="variable for y axis", required=True)
parser.add_argument("-c","--color", help="variable for color", required=True)
args = vars(parser.parse_args())

#
# Feature extraction
#

with open(args["meta"], "rb") as fp:
	df_meta = pickle.load(fp)
with open(args["df1"], "rb") as fp:
	df1 = pickle.load(fp)
with open(args["df2"], "rb") as fp:
	df2 = pickle.load(fp)

df = df_meta.merge(df1.merge(df2, left_index=True, right_index=True), left_index=True, right_index=True)

print("Shapes:", df1.shape, df2.shape, df.shape)

for col_id, col in enumerate(df.columns.tolist()):
	print(col_id, col)

print(df.mem_volume[0],df.dna_volume[0])

#
# Save
#

DATA = []

xvar = args["x"]
yvar = args["y"]
factor_name = args["color"]

factor = np.unique(df[factor_name].values)

print(factor)

for fac in factor:

	print(fac)

	ids = df[factor_name] == fac

	ncells = ids.sum()

	subDATA = dict()
	subDATA["x"] = df[xvar][ids].values.tolist()
	subDATA["y"] = df[yvar][ids].values.tolist()
	subDATA["label"] = [fac] * ncells
	subDATA["name"]  = fac
	subDATA["id"]  = df["cell_id"][ids].values.tolist()

	DATA.append(subDATA)

#
# Save JSON
#

with open("../../static/results/features.json","w") as fj:
	json.dump({"DATA": DATA, "LABS": [xvar,yvar]}, fj, indent=4, sort_keys=True)

os.system("rm ../../static/results/analysis.jpg")
print("Rscript cond.R fix/free",xvar,yvar)

#
# Backup
#

#backupFiles("./featurex.py", "featurex.py", "../../static/results/config-"+suffix+".json")
