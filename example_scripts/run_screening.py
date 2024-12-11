import crystfel_tools as cft

configpath = ''           # set the path to the config file

list_in = None            # list of images to optimize geometry for, if None, h5 files from given data directory are used


Experiment = cft.Experiment(configpath=configpath)

# example screening dictionary

screen = dict()
screen['--indexing'] = ['xgandalf']
screen['--min-snr'] = ['3.5','4','4.5','5']
screen['--min-pix-count'] = ['1','2']
screen['--threshold'] = ['700','800','900','1000','1100']
screen['--max-res'] = ['1200']
screen['--min-res'] = ['20']
# screen['--xgandalf-sampling-pitch'] = ['6','7']
# screen['--xgandalf-grad-desc-iterations'] = ['6','7']
# screen['--xgandalf-min-lattice-vector-length'] = ['5','10','20','40','80']
# screen['--xgandalf-max-lattice-vector-length'] = ['30','50','100','200','400']


Experiment.screen_indexing_parameters(screen,list_in=list_in) # list files optional 