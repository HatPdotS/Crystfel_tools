import crystfel_tools.crystfel_tools as cft


cell = [14.98,18.87,18.95,88.36,84.87,67.86]

hkl = cft.read_partialator_hkl('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/indexing_all_0_0.hkl')

hkl = cft.bin_hkl(hkl,cell,bins=50)
print(hkl.groupby('bin').res.mean())