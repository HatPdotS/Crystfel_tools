import crystfel_tools.crystfel_tools as cft




configpath = ''            # Path to the config file, has to be set, work directory is created in the same directory as the config file 




experiment_name = None      # Name of the experiment, if None, the user is prompted to enter it later

Cellfile = None             # Path to the cell file, if None, the user is prompted to enter it later   

regex_data = None             # Regex matching h5 files, if None, the user is prompted to enter it later

geometry_file = None        # Path to the geometry file, if None, the user is prompted to enter it later



Experiment = cft.Experiment(configpath,load_experiment=False,experiment_name=experiment_name,Cellfile=Cellfile,Geomin=geometry_file,
                            regex_data=regex_data)