import crystfel_tools.crystfel_tools as cft

configpath = ''           # set the path to the config file

list_in = None            # list of images to optimize geometry for, if None, h5 files from given data directory are used


indexamajiq_config = None

Experiment = cft.Experiment(configpath=configpath)

Experiment.run_idexamajig_for_mille(list_in=list_in,indexamajig_config=indexamajiq_config)          # runs indexing with mille flag, uses default indexing parameters, takes indexing dict like the one in run_bulk_indexing.py
Experiment.wait_until_done()                                  # waits until indexing is done
Experiment.run_align_detector()                               # runs align detector with created bin files creates a new geom and sets that as the new geometry file in the config file



