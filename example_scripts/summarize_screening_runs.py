import crystfel_tools.crystfel_tools as cft


configpath = ''   # set the path to the config file

settings_to_consider = ['--min-pix-count','--threshold','n_indexed','--min-snr'] # list of settings to consider for listing in df





Experiment = cft.Experiment(configpath=configpath)
df = Experiment.summarize_runs()
df = df[settings_to_consider]
df.sort_values(by='n_indexed',inplace=True,ascending=False)

print(df.head(10))