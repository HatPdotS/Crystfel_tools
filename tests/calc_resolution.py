from crystfel_tools import crystfel_tools as cft 
import numpy as np
import reciprocalspaceship as rs
from crystfel_tools.handling.fast_math import get_resolution


cell = [14.98,18.87,18.95,88.36,84.87,67.86] 
hkl = '/das/work/p17/p17489/Peter/2025-01-31_SwissFEL_SFX/Indexing_all/indexing_per_run_dd_2/work/semifinal_merge_all_separated/dark_0/dark_0.hkl'


hkl = cft.read_partialator_hkl(hkl)


hkl_dataset = rs.read_mtz('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/indexing_all_0_0.mtz')
hkl_dataset.reset_index(inplace=True)
cell = hkl_dataset.cell.parameters

hkl_dataset['res'] = get_resolution(hkl_dataset[['H','K','L']].values,cell)


hkl['res'] = get_resolution(hkl[['h','k','l']].values,cell)

print(hkl.res.min(),hkl.res.max())
print(hkl_dataset.res.min(),hkl_dataset.res.max())
