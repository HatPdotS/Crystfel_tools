import crystfel_tools.crystfel_tools as cft

configpath = ''    # set the path to the config file



Experiment = cft.Experiment(configpath=configpath)

pg = '' # set the pointgroup of your crystals

partialator_dict = dict()
partialator_dict['--model'] = 'xsphere' # 'unity' or 'xsphere' default is 'xsphere'
# partialator_dict['--model'] = 'unity' # comment this in if you want to use the unity model


sbatch_config = dict()
sbatch_config['-c'] = '56' # number of cores to use
sbatch_config['--mem'] = '450G' # amount of memory to use

partialator_id = Experiment.setup_partialator('Run_1',pg,partialator_config=partialator_dict,sbatch_config=sbatch_config)
Experiment.execute_partialator(partialator_id)