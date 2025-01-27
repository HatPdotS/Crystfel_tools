import crystfel_tools.handling.scaling as scaling
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
df = pd.read_feather('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/test_df.feather')
cell = [61.51,91.01,151.11,90,90,90]
scaler = scaling.scale_it(df,12000,0.01,cell,pg='mmm')
res = defaultdict(list)
scaler.shuffle_crystals()
scaler.set_reference('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/partialator_3.hkl')
for i in tqdm(range(5000)):
    scaler.load_crystal(i)

    scaler.scale_current_crystal_simple()
    
    scaler.apply_scaling_simple()

    cc_old,cc_new,rsplit_unweighted_old,rsplit_unweighted_new,rsplit_weighted_old,rsplit_weighted_new = scaler.check_scaling(print_results=False)
    res['cc_old'].append(cc_old)
    res['cc_new'].append(cc_new)
    res['rsplit_unweighted_old'].append(rsplit_unweighted_old)
    res['rsplit_unweighted_new'].append(rsplit_unweighted_new)
    res['rsplit_weighted_old'].append(rsplit_weighted_old)
    res['rsplit_weighted_new'].append(rsplit_weighted_new)
df = pd.DataFrame(res)
print(df.mean())

