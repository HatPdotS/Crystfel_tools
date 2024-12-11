import crystfel_tools as cft


configpath = '' # set the path to the config file


nchunks = 20   # number of chunks to split the job into got number is 20 for ra and 100 for esrf

param = dict()

# Example of a parameter dictionary for indexing

param['--indexing'] = 'xgandalf'
param['--min-snr'] = '4'
param['--min-pix-count'] = '2'
param['--threshold'] = '800'
param['--max-res'] = '2000'
param['--min-res'] = '20'
param['--int-radius'] = '2,4,6'




Experiment = cft.Experiment(configpath=configpath)


runid = Experiment.setup_run(indexamajig_config=param)


Experiment.execute_job_split(runid,nchunks=nchunks) # execute the job split in chunks be careful, if nchunk = 1000, 1000 jobs will be submitted