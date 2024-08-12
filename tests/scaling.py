import sparse
from handling import fast_math 
import numpy as np
from time import time
from handling.scaling import scale_it
import pandas as pd

x = sparse.load_npz('/das/work/p17/p17490/Peter/Library/crystfel_tools/tests/data.npz')
weigths = sparse.load_npz('/das/work/p17/p17490/Peter/Library/crystfel_tools/tests/weights.npz')
df = pd.read_feather('/das/work/p17/p17490/Peter/Library/crystfel_tools/tests/data.feather')

scale = scale_it(df)

# half1, half2 = scale.write_out_half_datasets()
# print(fast_math.calc_cc(half1,half2))

# for i in range(10):
scale.split_crystal()