import crystfel_tools.crystfel_tools as cft

configpath = ''           # set the path to the config file

list_in = None            # list of images to optimize geometry for, if None, h5 files from given data directory are used


Experiment = cft.Experiment(configpath=configpath)

Experiment.optimize_geometry(list_in=list_in)          # runs indexing with mille flag, uses default indexing parameters, takes indexing dict like the one in run_bulk_indexing.py

Experiment.run_align_detector()         # runs align detector with created bin files creates a new geom and sets that as the new geometry file in the config file



