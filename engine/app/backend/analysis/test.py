import numpy as np
import pandas as pd

U1 = pd.Series([  1,  2,  3,  4,  5,  6])
X = pd.Series([0.1,0.4,0.2,0.5,0.8,0.3])
Y = pd.Series(['A','A','A','B','B','B'])

df1 = pd.DataFrame({'id':U1, 'x':X, 'y':Y})

U2 = pd.Series([ 6, 4, 2])
Z = pd.Series([10,53,62])
W = pd.Series([16,23,62])

df2 = pd.DataFrame({'id':U2, 'z':Z, 'w':W})

print(df1)

print(df2)

df = df1.merge(df2, on='id', how='inner')

print(df)

ids = (df['y']=='B')

print(ids)

print(df[ids])

print(np.unique(df1['y'].values))