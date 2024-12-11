import crystfel_tools.crystfel_tools as cft

configpath = ''             # set the path to the config file

str_to_match = ''           # string to match for h5 files 


# Creates a list in a standardized location for this experiment based on the provided regular expression

# lists are created automatically based on the provided data directory 
# if no list is present in stadardized location

# This script is for the scenario when the directory structure is more
# complex and the user wants to create a list of all h5 files in multiple directories


Experiment = cft.Experiment(configpath=configpath)


Experiment.setup_list(str_to_match)