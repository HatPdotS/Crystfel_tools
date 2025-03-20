import crystfel_tools.crystfel_tools as cft


hklin ='/das/work/p17/p17489/Peter/2025-01-31_SwissFEL_SFX/Indexing_all/indexing_per_run_dd_2/work/semifinal_merge_all/all_dark_0/all_dark_0.hkl'
outfile = '/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/indexing_all_0_0.mtz'
space_group = 'P1'
unit_cell = [14.98,18.87,18.95,88.36,84.87,67.86]
cft.convert_crystfel_to_mtz_new(hklin, outfile, unit_cell,space_group,max_res=0.9)