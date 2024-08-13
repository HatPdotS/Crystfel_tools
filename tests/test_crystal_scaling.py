import crystfel_tools.handling.scaling as scaling
import pandas as pd
import crystfel_tools.handling.data_conversion as data_conversion
from crystfel_tools.handling import fast_math
import numpy as np
from scipy.optimize import minimize
from time import time

df = pd.read_feather('/home/esrf/hans1507/library/crystfel_tools/test_data/test_df.feather')


ref = data_conversion.load_reference('/home/esrf/hans1507/library/crystfel_tools/test_data/partialator_11.hkl',df)
ref = ref.loc[ref['I'] > 0]

ref.sort_values('idx_hkl',inplace=True)
ref['idx_hkl_new'] = range(len(ref))
ref.set_index('idx_hkl',inplace=True)

df = df.loc[df['idx'].isin(ref.index)]
df['idx_hkl'] = df['idx_hkl'].map(ref['idx_hkl_new'])
df.dropna(inplace=True)
ref.reset_index(inplace=True)
ref.idx_hkl = ref.idx_hkl_new
ref.drop(columns='idx_hkl_new',inplace=True)
ref.set_index('idx_hkl',inplace=True)
S = scaling.scale_it(df)
t = time()
for i in range(1000):
    crystal = S.df_split[i]
    crystal = S.sanitize_crystal(crystal)
    I_ref = ref.iloc[crystal['idx_hkl'].values]

    before = fast_math.calc_cc(crystal['I'].values,I_ref['I'].values)

    def calc_residual(factors):
        residual = np.sum((np.log(I) - np.log(I_ref) - np.log(factors[0]) - factors[1] * hkl[:,0]**2 - factors[2] * hkl[:,1]**2 - factors[3] * hkl[:,2]**2)**2)
        return residual

    I = crystal['I'].values
    hkl = crystal[['h','k','l']].values 
    I_ref = I_ref['I'].values
    s = 1/crystal.resolution.values


    res = minimize(calc_residual,np.array((1,0,0,0)),method='L-BFGS-B',bounds=[(1e-20,None),(None,None),(None,None),(None,None)])
    I_scaled = I * res.x[0] * np.exp(-(res.x[1] * hkl[:,0]**2 + res.x[2] * hkl[:,1]**2 + res.x[3] * hkl[:,2]**2))

    after = fast_math.calc_cc(I_scaled,I_ref)

print(t -time())