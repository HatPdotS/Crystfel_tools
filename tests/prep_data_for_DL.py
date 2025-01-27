import crystfel_tools.handling.scaling as scaling
import crystfel_tools.handling.fast_math as fast_math
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
import numpy as np


def addpadding(array,size=1000):
    if array.ndim == 1:
        new = np.zeros(size)
        new[:array.shape[0]] = array[:size]
        return new
    if array.ndim == 2:
        new = np.zeros((size,array.shape[1]))
        new[:array.shape[0]] = array[:size]
        return new
    print('array has wrong dimensions')
vectors = []
partialities = []
df = pd.read_feather('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/test_df.feather')
cell = [61.51,91.01,151.11,90,90,90]
scaler = scaling.scale_it(df,14000,0.01,cell,pg='mmm')
scaler.set_reference('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/partialator_3.hkl')
rec = fast_math.calc_rec_space_vector(cell)
for i in tqdm(range(10000)):
    scaler.load_crystal(i)
    sorting = np.argsort(scaler.current_crystal.I.values)[::-1]
    latticevectors = fast_math.convert_to_rec_space(scaler.current_crystal[['h','k','l']].values,rec)[sorting]
    latticevectors = addpadding(latticevectors).reshape(1,-1,3)
    par = (scaler.current_crystal.I.values / scaler.current_ref.I.values)[sorting]
    par = addpadding(par).reshape(1,-1)
    par = par/par.max()
    vectors.append(latticevectors)
    partialities.append(par)
vectors = np.concatenate(vectors,axis=0)
partialities = np.concatenate(partialities,axis=0)
np.save('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/vectors.npy',vectors)
np.save('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/partialities.npy',partialities)