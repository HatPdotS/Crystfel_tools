import crystfel_tools.crystfel_tools as cft

cf = '/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/config.json'




Experiment = cft.Experiment(cf,load_experiment=False,experiment_name='test',Cellfile='/das/work/p17/p17488/seidel_h/2024-06-22_ESRF_SFX/S2.cell',Geomin='/das/work/p17/p17488/seidel_h/2024-06-22_ESRF_SFX/geom/refined3.geom',
                            data_dir='/das/work/p17/p17488/seidel_h/2024-06-22_ESRF_SFX/S2_new_embedded_0_days/run_01_ssx_injector_collection',partition='day')

Experiment.optimize_geometry()

Experiment.run_align_detector()



