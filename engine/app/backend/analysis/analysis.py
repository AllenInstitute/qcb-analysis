import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

#sys.path.insert(0,"../aux/")
# from auxFunctions import backupFiles

pd.set_option("display.max_columns", None)
np.set_printoptions(threshold=np.nan)

#
# Arguments
#

parser = argparse.ArgumentParser(description="Exploration of QCB data")
parser.add_argument("-df","--df", help="Path to data frame", required=True)
parser.add_argument("-x","--x", help="variable for x axis", required=True)
parser.add_argument("-y","--y", help="variable for y axis", required=True)
parser.add_argument("-c","--color", help="variable for color", required=True)
args = vars(parser.parse_args())

#
# Feature extraction
#

df = pd.read_csv(os.path.join("../../../../data-raw",args["df"]+".csv"))

for col_id, col in enumerate(df.columns.tolist()):
	print(col_id, col)

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
