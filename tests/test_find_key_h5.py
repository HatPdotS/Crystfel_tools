import crystfel_tools.crystfel_tools as cft

file_esrf = '/das/work/p17/p17488/seidel_h/2024-06-22_ESRF_SFX/S2_new_embedded_0_days/run_01_ssx_injector_collection/S2_cage-S2_cage_dense_2_00009.h5'



d = cft.find_key_h5_files(file_esrf)

print(d)

