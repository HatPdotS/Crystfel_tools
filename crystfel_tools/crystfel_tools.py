import os
from itertools import product
import json
import slurm_tools
import pandas as pd
from glob import glob
import h5py
import shutil
import numpy as np
from collections import defaultdict
import datetime
import copy
import reciprocalspaceship
import gemmi
from crystfel_tools.handling import fast_math

def edit_geom(template: str, to_change: dict,outfile):
    with open(template) as f:
        with open(outfile,'w') as w:
            for line in f:
                if line[0] == ';':
                    continue
                linenew = line
                for key in to_change.keys():
                    if key == line[:len(key)]:
                        linenew = key + ' = ' + to_change[key] + '\n'

                        break
                w.write(linenew)

def read_runtime_from_log(logfile):
    with open(logfile) as f:
        lines = f.readlines()
        firstline = lines[0]
        lastline = lines[-1]
        start = datetime.datetime.strptime(firstline.strip(),"%a %b %d %H:%M:%S %Z %Y")
        try:
            end = datetime.datetime.strptime(lastline.strip(),"%a %b %d %H:%M:%S %Z %Y")
        except ValueError:
            print('Returning NAN, could not read end time, process likely still running or quit unexpectedly for log:',logfile)
            return np.nan
    return end - start

def get_len_text_file(file):
    with open(file) as f:
        return len(f.readlines())

def bin_hkl(hkl,cell,bins=20):
    hkl['res'] = fast_math.get_resolution(hkl[['h','k','l']].values,cell)
    sorted_res = np.sort(hkl.res)[::-1]
    sorted_res[0] = sorted_res[0] + 10
    sorted_res[-1] = 0
    reflection_to_target = sorted_res.shape[0] // bins
    bins = sorted_res[::reflection_to_target]
    bins[-1] = 0
    hkl['bin'] = np.digitize(hkl.res,bins)
    return hkl

def convert_crystfel_to_mtz(file,outfile,cell,symm):
    if isinstance(cell,list):
        cell = ' '.join([str(p) for p in cell])
    os.system(f"sed -n '/End of reflections/q;p' {file} > create-mtz.temp.hkl")
    cmd = f"""module load ccp4; f2mtz HKLIN create-mtz.temp.hkl HKLOUT {outfile} > out.html << EOF
TITLE Reflections from CrystFEL
NAME PROJECT wibble CRYSTAL wibble DATASET wibble
CELL {cell}
SYMM {symm}
SKIP 3
LABOUT H K L IMEAN SIGIMEAN
CTYPE  H H H J     Q
FORMAT '(3(F4.0,1X),F10.2,10X,F10.2)'
EOF"""
    os.system(cmd)
    os.system('rm create-mtz.temp.hkl')

def convert_crystfel_to_mtz_new(file,outfile,cell,symm,max_res=None):
    h = []
    k = []
    l = []
    I = []
    SIGMA = []
    nmeas = []
    with open(file) as f:
        next(f)
        next(f)
        next(f)
        for line in f:
            if line.strip('\n') == 'End of reflections':
                break
            _h,_k,_l,_I,_,_sigma,_nmeas = line.strip('\n').split()
            h.append(np.int32(_h))
            k.append(np.int32(_k))
            l.append(np.int32(_l))
            I.append(float(_I))
            SIGMA.append(float(_sigma))
            nmeas.append(int(_nmeas))
    dataset = reciprocalspaceship.DataSet({'H':h,'K':k,'L':l,'I':I,'sigma(I)':SIGMA,'NMEAS':nmeas})
    dataset.set_index(['H','K','L'],inplace=True)
    dataset.cell = cell
    dataset.spacegroup = symm
    if max_res != None:
        resolution = fast_math.get_resolution(np.array([h,k,l]).T,cell)
        dataset = dataset.loc[resolution > max_res]
    dataset.infer_mtz_dtypes(inplace=True)
    dataset.write_mtz(outfile)

def parse_sinfo_to_dataframe() -> pd.DataFrame:
    import subprocess
    data = str(subprocess.check_output(["sinfo"]))
    data = data.replace("\\n", "\n")
    lines = [line.strip("'").strip() for line in data.split("\n") if line.strip()][1:]  
    lines = [line for line in lines if line]
    parsed_data = []
    for line in lines:
        parsed_data.append(line.split())
    df = pd.DataFrame(parsed_data, columns=["Partition", "Avail", "TimeLimit", "Nodes", "State", "Nodelist"])
    timlimits = []
    for i in df.TimeLimit:
        if i == "infinite":
            timlimits.append(pd.to_timedelta("1000 days"))
            continue
        timlimits.append(pd.to_timedelta(i.replace("-", " days ")))
    df["TimeLimit"] = timlimits
    df["Nodes"] = df["Nodes"].astype(int)
    default_partition = [i for i in df.Partition if '*' in i][1]
    df.Partition = df.Partition.str.replace("*", "")
    df = df.loc[df.Partition != "admin"]
    print('I found the following partitions:')
    df_dense = df.groupby("Partition").first()
    df_dense.Nodes = df.groupby("Partition").Nodes.sum()
    print(df_dense)
    print('The default partition is:' , default_partition.strip('*'))   
    df_dense = df_dense.loc[df_dense.Nodes > 5]
    df_dense = df_dense.sort_values(by="TimeLimit", ascending=False)
    return list(df_dense.index)

def sorted_dict(d):
    return dict(sorted(d.items(), key=lambda item: item[1], reverse=True))

def recursive_inventory_search(open_h5_file,paths):
    try: 
        keys = open_h5_file.keys()
        for key in keys:
            recursive_inventory_search(open_h5_file[key],paths)
    except AttributeError:
        paths[open_h5_file.name] = open_h5_file.size

def find_key_h5_files(h5file_path):
    with h5py.File(h5file_path,'r') as f:
        paths = dict()
        recursive_inventory_search(f,paths)
    for key in list(paths.keys()):
        if 'Simulator' in key:
            del paths[key]
    print(sorted_dict(paths))
    likely_path = list(sorted_dict(paths).keys())[0]
    print('found the most likely path to data:',likely_path)
    return likely_path

def get_all_events_smart(pathin,h5py_path):
        pathout = pathin.replace('.lst','_all_events.lst')
        with open(pathin,'r') as f:
            with open(pathout,'w') as fnew:
                for line in f:
                    line = line.strip()
                    id = line.split('/')[-1].split('.')[-2].split('_')[0]
                    try: 
                        with h5py.File(line,'r') as f:
                            n  = f[h5py_path].shape[0]
                    except: continue
                    for i in range(n):

                        to_write = line + ' //' + str(i) + '\n'
                        fnew.write(to_write)
        return pathout

def make_difference_map(mtz1,mtz2,outdir,dark_model,restraints,resmax=100,resmin=0):
    from functions.io_functions import make_diff_map_Xtrapol8
    os.makedirs(outdir,exist_ok=True)
    make_diff_map_Xtrapol8(mtz1, mtz2, dark_model,outdir, restraints ,res_max=resmax,res_min=resmin)

def make_indexamajiq_cmd(parameters: dict,IO_config: dict):
    base = 'indexamajig '
    base += make_cmd_from_dict(parameters)
    base += make_cmd_from_dict(IO_config)
    return base

def make_partialator_cmd(parameters: dict):
    base = 'partialator '
    base += make_cmd_from_dict(parameters) 
    return base

def make_cmd_from_dict(dic):
    base = ''
    for key in dic.keys():
        if dic[key] == None:
            continue
        if key.startswith('--'):
            if dic[key] == '':
                base += str(key) +' '
                continue
            base += str(key) + '='+ str(dic[key]) + ' '
        elif key.startswith('-'):
            base += str(key) + ' ' + str(dic[key]) + ' '
    return base

def make_sbatch_cmd(parameters: dict):
    base = 'sbatch '
    base += make_cmd_from_dict(parameters)
    return base

def config_equal(config1,config2):
    for key in config1.keys():
        if config1[key] == None:
            continue
        if not key in config2.keys():
            return False
        if config1[key] != config2[key]:
            return False
    for key in config2.keys():
        if config2[key] == None:
            continue
        if not key in config1.keys():
            return False
        if config2[key] != config1[key]:
            return False
    return True

def run_single_indexamajiq(batch_dict, indexamajiq_dict,IO_config):
    sbatch = make_sbatch_cmd(batch_dict)
    indexer = make_indexamajiq_cmd(indexamajiq_dict,IO_config)
    cmd_template = '''{sbatch} --wrap="
module purge
module load crystfel/0.11.1
date
echo {indexer}
{indexer}
date"'''
    cmd = cmd_template.format(sbatch = sbatch,indexer = indexer)
    return slurm_tools.submit_string(cmd)

def run_partialator(batch_dict, partialator_config):
    sbatch = make_sbatch_cmd(batch_dict)
    partialator = make_partialator_cmd(partialator_config)
    cmd_template = '''{sbatch} --wrap="
module purge
module load crystfel/0.11.1
echo {partialator}
{partialator}"'''
    cmd = cmd_template.format(sbatch = sbatch,partialator = partialator)
    return slurm_tools.submit_string(cmd)

def shuffle_lines_list(filein,n=None):
    import random
    if n == None:
        fileout = filein.replace('.lst','.shuffled.lst')
    else:
        fileout = filein.replace('.lst',f'.shuffled_sample_{n}.lst')
    with open(filein,'r') as f:
        d = f.readlines()
    if n == None or n >= len(d):
        random.shuffle(d)
    else:
        d = random.sample(d,n)
    with open(fileout,'w') as g:
        g.writelines(d)
    return fileout

def cut_stream_after_ncrystals(streamfile,outfile,n):
    counter = 0
    with open(streamfile,'r') as f:
        with open(outfile,'w') as g:
            for line in f:
                if '--- Begin crystal' == line.strip('\n'):
                    counter += 1
                if line.strip('\n') == '----- Begin chunk -----':
                    if counter >= n:
                        break
                g.write(line)

def produce_1000_events_per_file(pathin):
    pathout = pathin.replace('.lst','_all_events.lst')
    with open(pathin,'r') as f:
        with open(pathout,'w') as n:
            for line in f:
                line = line.strip()

                for i in range(1000):
                    to_write = line + ' //' + str(i) + '\n'
                    n.write(to_write)

def split_list_file(pathin,nchunks,list_sub_dir = None):
    if list_sub_dir == None:
        name_list = [pathin.replace('.lst','_' + str(i)+'.lst') for i in range(nchunks)]
    else:
        path_list = '/'.join(pathin.split('/')[:-1])
        os.makedirs(path_list + '/' + list_sub_dir,exist_ok=True)
        name_list = [path_list + '/' + list_sub_dir + '/' + pathin.split('/')[-1].replace('.lst','_' + str(i)+'.lst') for i in range(nchunks)]
    with open(pathin,'r') as f:
        chunksize = len(f.readlines())//nchunks
        print(chunksize)
    with open(pathin,'r') as f:
        for file in name_list:
            with open(file,'w') as w:
                for _ in range(chunksize):
                    l = next(f)
                    w.write(l)
    return name_list

def produce_5000_events_per_file(pathin):
    pathout = pathin.replace('.lst','_all_events.lst')
    with open(pathin,'r') as f:
        with open(pathout,'w') as n:
            for line in f:
                line = line.strip()
                for i in range(5000):

                    to_write = line + ' //' + str(i) + '\n'
                    n.write(to_write)
    return pathout

def get_n_indexed_from_stream(streamfile):
    c = 0
    try:
        with open(streamfile,'r') as f:
            for line in f:
                if 'Cell' == line[:4]:
                    c += 1
    except FileNotFoundError:
        print('could not find file',streamfile)
    return c

def setup_standard_config():
    sbatch_parameters = dict()
    sbatch_parameters['-p'] = 'low'
    sbatch_parameters['-t'] = '1-00:00:00'
    sbatch_parameters['-c'] = '32'
    parameters = dict()
    parameters['-i'] = None
    parameters['-g'] = None
    parameters['-j'] = '32'
    parameters['-p'] = None
    parameters['--peaks'] = 'peakfinder8'
    parameters['--min-snr'] = '5'
    parameters['--min-peaks'] = '5'
    parameters['--int-radius'] = '5,7,10'
    parameters['--indexing'] = 'xgandalf'
    parameters['--min-pix-count'] = '2'
    parameters['--multi'] = ''
    parameters['--threshold'] = '700'
    parameters['--no-non-hits-in-stream'] = ''
    sbatch_parameters['--output'] = None
    sbatch_parameters['--error'] = None
    parameters['-o'] = None
    return parameters,sbatch_parameters

def mask_maker(h5_in,n=1000,id=None,boxsize=7,mask_inner_edges=True,mask_outer_edges=True):
    import h5py
    from numpy.lib.stride_tricks import sliding_window_view
    with h5py.File(h5_in,'r') as f:
        if id == None:
            id = h5_in.split('/')[-1].split('.')[-2]
        
        data = f[f'data/{id}/data'][:n].mean(axis=0)
        data_mean = np.pad(sliding_window_view(data,(boxsize,boxsize),axis=(0,1)).mean(axis=(2,3)),boxsize//2)
        data_std= np.pad(sliding_window_view(data,(boxsize,boxsize),axis=(0,1)).mean(axis=(2,3)),boxsize//2)
        mask = data < 2*data_std + data_mean
    if mask_inner_edges:
        edges = data == 0
        edges_all = np.logical_or(edges,np.roll(edges,1,axis=1))
        edges_all = np.logical_or(edges_all,np.roll(edges,1,axis=0))
        edges_all = np.logical_or(edges_all,np.roll(edges,-1,axis=1))
        edges_all = np.logical_or(edges_all,np.roll(edges,-1,axis=0))
        mask = np.logical_and(mask,~edges_all)
    if mask_outer_edges:
        mask[0] = False
        mask[-1] = False
        mask[:,0] = False
        mask[:,-1] = False
        mask[0] = False
        mask[513] = False
        mask[514] = False
        mask[1027] = False
        mask[1028] = False
        mask[1541] = False
        mask[1542] = False
        mask[2055] = False
        mask[2056] = False
        mask[2569] = False
        mask[2570] = False
        mask[3083] = False
        mask[3084] = False
        mask[3597] = False
        mask[3598] = False
        mask[4111] = False
    return mask,data

def cat_files(list_files,outpath,remove_group_statements = True):
    with open(outpath,'w') as f:
        for file in list_files:
            try:
                with open(file,'r') as g:
                    for line in g:
                        if remove_group_statements:
                            if line.strip('\n') == 'bandwidth = 1.000000e-08':
                                continue
                            if line[:6] == 'group_':
                                continue
                        f.write(line)
            except: print('could not find file',file)

def make_list(str_to_match,list_out):
    files = glob(str_to_match)
    with open(list_out,'w') as f:
        for file in files:
            f.write(file + '\n')
    return files

def save_mask(mask,outpath):
    import h5py
    with h5py.File(outpath,'w') as f:
        f.create_dataset('pixel_mask',data=mask)
    
def plot_mask(mask,mean_data,outpath):
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(1,2)
    ax[0].imshow(mask)
    ax[1].imshow(mean_data)
    plt.savefig(outpath)

def read_hits(file_open):
    Event = 0
    filename = ''
    fs = []
    ss = []
    res = []
    I = []
    Panel = []
    while True: 
        line = next(file_open)
        if line[:14] == 'Image filename':
            filename = line.split()[-1].strip('\n')
        if line[:5] == 'Event':
            Event = line.split()[-1].strip('//').strip('\n')
        if line[:5] == 'Peaks':
            line = next(file_open)

            while True:
                line = next(file_open)
                if line[:16] == 'End of peak list':
                    break
                _fs, _ss, _res , _I, _Panel = line.strip('\n').split()
                fs.append(float(_fs))
                ss.append(float(_ss))
                res.append(float(_res))
                I.append(float(_I))
                Panel.append(_Panel)
            break
    df = pd.DataFrame({'fs': fs,'ss': ss, 'res': res, 'I': I, 'Panel': Panel, 'filename': filename, 'Event': Event, })
    return df, filename , Event

def remove_geometry_file(filein):
    with open(filein) as f:
        with open(filein.replace('.stream','_no_geom.stream'),'w') as g:
            geom = False
            for line in f:
                if line.strip('\n') == '----- Begin geometry file -----':
                    geom = True
                if line.strip('\n') == '----- End geometry file -----':
                    geom = False
                if not geom:
                    g.write(line)

def read_crystal(file_open):
    lines = ['--- Begin crystal\n']
    cell = ''
    while True: 
        line = next(file_open)
        lines.append(line)
        if line[:15] == 'Cell parameters':
            cell = line[15:].strip('\n').replace('nm,','').replace('deg','')
            continue
        if line.strip('\n') == 'Reflections measured after indexing':
            lines.append(line)
            line = next(file_open)
            h = []
            k = []
            l = []
            I = []
            sigmaI = []
            peak = []
            background = []
            fs = []
            ss = []
            panel = []
            while True: 
                line = next(file_open)
                if line.strip('\n') == 'End of reflections':
                    break
                _h,_k,_l,_I,_sigmaI,_peak,_background,_fs,_ss,_panel = line.strip('\n').split()
                h.append(int(_h))
                k.append(int(_k))
                l.append(int(_l))
                I.append(float(_I))
                sigmaI.append(float(_sigmaI))
                peak.append(float(_peak))
                background.append(float(_background))
                fs.append(float(_fs))
                ss.append(float(_ss))
                panel.append(_panel)
            break
    df = pd.DataFrame({'h': h, 'k': k, 'l': l, 'I': I, 'sigmaI': sigmaI, 'peak': peak, 'background': background, 'fs': fs, 'ss': ss, 'panel': panel, 'cell': cell})
    return df, lines

def read_only_indexed(filein, n_hits=100000):
    with open(filein) as f:
        hits = []
        crystals = []
        for line in f: 
            if line.strip('\n') == '----- Begin chunk -----':
                df_hit, _ , _ = read_hits(f)
                if next(f).strip('\n') == '--- Begin crystal':
                    df_crystal, _ = read_crystal(f)
                    hits.append(df_hit)
                    crystals.append(df_crystal)
                    if len(hits) == n_hits:
                        break
    return hits, crystals

def get_crystals(filein):
    with open(filein) as f:
        crystals = []
        for line in f: 
            if line.strip('\n') == '----- Begin chunk -----':
                _, filename , event = read_hits(f)
                if next(f).strip('\n') == '--- Begin crystal':
                    df_crystal, _ = read_crystal(f)
                    df_crystal['filename'] = filename
                    df_crystal['event'] = event
                    crystals.append(df_crystal)
    crystals = pd.concat(crystals)
    return crystals

def get_hits(filein,only_indexed = False):
    with open(filein) as f:
        hits = []
        for line in f: 
            if line.strip('\n') == '----- Begin chunk -----':
                df_hit, filename , event = read_hits(f)
                acq = int(filename.split('acq')[1].split('.')[0])
                runnr = int(filename.split('run')[1].split('-')[0])
                df_hit['acq'] = acq
                df_hit['runnr'] = runnr
                df_hit['event'] = event
                if next(f).strip('\n') == '--- Begin crystal' or not only_indexed:
                    hits.append(df_hit)
    return hits

def read_intensities_from_stream(stream):
    with open(stream) as f:
        intensities = []
        for line in f:
            if line == 'Reflections measured after indexing\n':
                next(f)
                for line in f:
                    if line == 'End of reflections\n':
                        break
                    intensities.append(float(line.split()[3]))
    return np.array(intensities)
    
def get_cells(filein):
    with open(filein) as f:
        cells = []
        for line in f: 
            if line.startswith('Cell parameters'):
                _,_,a,b,c,_,alpha,beta,gamma,_ =line.split()
                cells.append((float(a)*10,float(b)*10,float(c)*10,float(alpha),float(beta),float(gamma)))
    cells = np.array(cells)
    return cells

def get_max_peak_resolution(filein):
    with open(filein) as f:
        res = []
        for line in f: 
            if line.startswith('peak_resolution'):
                res.append(float(line.split()[-2]))
    res = np.array(res)
    return res

def get_profile_radii_from_stream(streamin):
    radii = []
    with open(streamin) as f:
        for line in f:
            if line.startswith('profile_radius'):
                radii.append(float(line.split()[2]))
    return radii

def align_spots(df_hits,df_crystal):
    panels = df_hits['Panel'].unique()
    df_hits['deviation'] = np.nan
    for i in panels:
        hits_panel = df_hits[df_hits['Panel'] == i]
        crystal_panel = df_crystal[df_crystal['panel'] == i]
        for f in hits_panel.index:
            df_hits.loc[f,'deviation'] = np.min(np.sqrt((df_hits.loc[f,'fs'] - crystal_panel.fs)**2 + (df_hits.loc[f,'ss'] - crystal_panel.ss)**2))
    return df_hits

def calc_deviations(streamin):
    path_out  = streamin.replace('.stream','_deviations.csv')
    if os.path.exists(path_out):
        return pd.read_csv(path_out)
    hits,crystals = read_only_indexed(streamin)
    [align_spots(a,b) for a,b in zip(hits,crystals)]
    hits = pd.concat(hits)
    hits.to_csv(path_out)
    return hits

def get_indexamajiq_parameters():
    parameters = dict()
    parameters['-j'] = '32'
    parameters['--peaks'] = 'peakfinder8'
    parameters['--min-snr'] = '5'
    parameters['--min-peaks'] = '5'
    parameters['--int-radius'] = '2,4,6'
    parameters['--indexing'] = 'xgandalf'
    parameters['--min-pix-count'] = '2'
    parameters['--multi'] = ''
    parameters['--threshold'] = '700'
    parameters['--no-non-hits-in-stream'] = ''
    parameters['--max-res'] = '3000'
    return parameters

def get_sbatch_standard():
    sbatch_parameters = dict()
    sbatch_parameters['-p'] = 'day'
    sbatch_parameters['-t'] = '1-00:00:00'
    sbatch_parameters['-c'] = '32'
    return sbatch_parameters

def create_config(configpath,cell_file, geompath, outdir=None):
    params, sbatch_params = setup_standard_config()
    with open(configpath,'r') as f:
        config = json.load(f)
    params['-p'] = cell_file
    params['-g'] = geompath
    if outdir == None:
        outdir = os.path.join(os.path.dirname(configpath),'work')
    params['-o'] = outdir
    return config

def make_standard_screen():
    screen = dict()
    screen['--min-snr'] = ['4','4.5','5','5.5','6']
    screen['--threshold'] = ['5000']
    screen['--indexing'] = ['xgandalf','pinkIndexer']
    screen['--int-radius'] = ['5,7,10']
    screen['--min-pix-count'] = ['2','3']
    return screen

def get_statistics(hklin,cell,resmax=None,resmin=None,nshells=20,
                   fom_list = ['rsplit','cc','"cc*"']):
    import subprocess
    dfs = []
    check_hkl_out = hklin.replace('.hkl','_stats.dat')
    cmd = '''module use MX; module load crystfel/0.11.1; check_hkl {hklin} -p {cell} --nshells={nshells} --shell-file={check_hkl_out}'''
    cmd = cmd.format(hklin=hklin,cell=cell,nshells=nshells,check_hkl_out=check_hkl_out)
    if resmax != None:
        cmd += ' --highres={rescut}'.format(rescut=resmax)
    if resmin != None:
        cmd += ' --lowres={rescut}'.format(rescut=resmin)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print(f"Error: {result.stderr}")
    df = pd.read_csv(check_hkl_out,sep=r'\s+',skiprows=1,names=['1/d centre','nrefs_possible','nrefs_observed','Compl','Meas','Red','SNR','Mean I','d/A','Min 1/nm','Max 1/nm'])
    df.set_index('d/A',inplace=True)
    df = df[['nrefs_observed','Compl','Meas','Red','SNR','Mean I']]
    dfs.append(df)
    for fom in fom_list:
        hkl1 = hklin.replace('.hkl',f'.hkl1')
        hkl2 = hklin.replace('.hkl',f'.hkl2')
        path_out = hklin.replace('.hkl',f'_{fom}.dat').replace('"cc*"','ccstar')
        cmd = '''module use MX; module load crystfel/0.11.1; compare_hkl {hkl1} {hkl2} -p {cell} --fom={fom} --nshells={nshells} --shell-file={path_out}'''
        cmd = cmd.format(hkl1=hkl1,hkl2=hkl2,cell=cell,fom=fom,nshells=nshells,path_out=path_out)
        if resmax != None:
            cmd += ' --highres={rescut}'.format(rescut=resmax)
        if resmin != None:
            cmd += ' --lowres={rescut}'.format(rescut=resmin)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(f"Error: {result.stderr}")
        if fom == '"cc*"':
            fom = 'ccstar'
        df = pd.read_csv(path_out,sep=r'\s+',skiprows=1,names=['1/d centre',fom,'nref','d/A','Min 1/nm','Max 1/nm'])
        df.set_index('d/A',inplace=True)
        dfs.append(df[fom])
    df = pd.concat(dfs,axis=1)
    print(df.to_markdown())
    with open(hklin.replace('.hkl','.statistics.dat'),'w') as f:
        f.write(df.to_markdown())
    df.to_csv(hklin.replace('.hkl','.statistics.csv'))
    return df

def convert_crystfel_to_mtz(file,outfile,cell,symm):
    os.makedirs(os.path.dirname(outfile),exist_ok=True)
    if isinstance(cell,list) or isinstance(cell,tuple):
        cell = ' '.join([str(p) for p in cell])
    cmd = r"""module load ccp4; sed -n '/End\ of\ reflections/q;p' {file} > create-mtz.temp.hkl;
    f2mtz HKLIN create-mtz.temp.hkl HKLOUT {outfile} > out.html << EOF
TITLE Reflections from CrystFEL
NAME PROJECT wibble CRYSTAL wibble DATASET wibble
CELL {cell}
SYMM {symm}
SKIP 3
LABOUT H K L IMEAN SIGIMEAN
CTYPE  H H H J     Q
FORMAT '(3(F4.0,1X),F10.2,10X,F10.2)'
EOF"""
    cmd = cmd.format(file=file,outfile=outfile,cell=cell,symm=symm)
    print(cmd)
    os.system(cmd)
    os.system('rm create-mtz.temp.hkl')

def read_partialator_hkl(file):
    res = defaultdict(list)
    with open(file) as f:
        next(f)
        next(f)
        next(f)
        for line in f:
            if line.strip('\n') == 'End of reflections':
                break
            h,k,l,I,phase,sigma,nmeas = line.strip('\n').split()
            res['h'].append(int(h))
            res['k'].append(int(k))
            res['l'].append(int(l))
            res['I'].append(float(I))
            if phase == '-':
                phase = np.nan
            res['sigma'].append(float(sigma))
            res['phase'].append(float(phase))
            res['nmeas'].append(int(nmeas))
    df = pd.DataFrame(res)
    df.h = df.h.astype(int)
    df.k = df.k.astype(int)
    df.l = df.l.astype(int)
    df.I = df.I.astype(float)
    df.sigma = df.sigma.astype(float)
    df.nmeas = df.nmeas.astype(int)
    # print(df)
    return df

def write_df_as_hklf4(df,outfile):
    df.sigma = df.sigma / df.I.mean() * 100
    df.I = df.I / df.I.mean() * 100
    with open(outfile,'w') as f:
        # f.write('H K L I SIGMA\n')
        for i in range(df.shape[0]):
            f.write(f'{int(df.iloc[i].h):>4}{int(df.iloc[i].k):>4}{int(df.iloc[i].l):>4}{df.iloc[i].I:>8.2f}{df.iloc[i].sigma:>8.2f}\n')

def confert_hkl_to_hklf4(filein,fileout):
    df = read_partialator_hkl(filein)
    write_df_as_hklf4(df,fileout)
    
def read_crystfel_cell_file(file):
    with open(file) as f:
        for line in f:
            if line[:3] == 'a =':
                a = float(line.split()[2])
            if line[:3] == 'b =':
                b = float(line.split()[2])
            if line[:3] == 'c =':
                c = float(line.split()[2])
            if line[:4] == 'al =':
                alpha = float(line.split()[2])
            if line[:4] == 'be =':
                beta = float(line.split()[2])
            if line[:4] == 'ga =':
                gamma = float(line.split()[2])
    return a,b,c,alpha,beta,gamma

def setup_experiment(configpath,experiment_name):
    config = dict()
    config['path'] = configpath
    config['experiment'] = experiment_name
    return config

def save_config(config: dict, path_out):
    with open(path_out,'w') as f:
        json.dump(config,f, indent=4)

def fill_name_config(name, indexamajiq_config, sbatch_parameters):
    indexamajiq_config['-i'] = name 
    indexamajiq_config['-o'] = name + '.stream'
    sbatch_parameters['--output'] = name + '.log'
    sbatch_parameters['--error'] = name + '.log'
    return indexamajiq_config, sbatch_parameters

def ask_for_input_until_file_exists(string_in):
    while True:
        path = input(string_in)
        if os.path.exists(path):
            break
        else:
            print('File does not exist')
    return path

def run_align_detector(geom,geom_out,mille_file_match,sbatch_parameters=None,level=1):
    if sbatch_parameters == None:
        sbatch_parameters = get_sbatch_standard()
    sbatch_parameters['-c'] = '1'
    sbatch_parameters['--chdir'] = os.path.dirname(geom_out)
    sbatch = make_sbatch_cmd(sbatch_parameters)
    cmd_template = '''{sbatch} --wrap="
module purge
module load crystfel/0.11.1
echo {angle_analyser}
{angle_analyser}"'''
    angle_analyser = f'align_detector -i {geom} -o {geom_out} -l {level} {mille_file_match}'
    cmd = cmd_template.format(sbatch = sbatch,angle_analyser = angle_analyser)
    return slurm_tools.submit_string(cmd)

def sanitize_stream(stream_in,Stream_out,write_hits=True):
    with open(stream_in,'r') as f:
        with open(Stream_out,'w') as g:
            for line in f:
                g.write(line)
                if line.startswith('----- Begin chunk -----'):
                    chunk = [line]
                    crystals = []
                    crystals_ok_beamradius = []
                    has_crystals = False
                    for line in f:
                        chunk.append(line)
                        if line.startswith('----- End chunk -----'):
                            break
                        if line.startswith('--- Begin crystal'):
                            has_crystals = True
                            crystal = []
                            crystal_beamradius_ok = True
                            for line in f:
                                crystal.append(line)
                                if line.startswith('--- End crystal'):
                                    break
                                if line.startswith('profile_radius'):
                                    if float(line.split('=')[-1].split()[0].strip()) == 0:
                                        print('found a crystal with beamradius 0')
                                        crystal_beamradius_ok = False
                            crystals_ok_beamradius.append(crystal_beamradius_ok)
                            crystals.append(crystal)
                    if write_hits or has_crystals:
                        g.writelines(chunk[:-1])
                        if has_crystals:
                            for i,crystal in enumerate(crystals):
                                if crystals_ok_beamradius[i]:
                                    g.writelines(crystal)
                        g.write(chunk[-1])
    return Stream_out

class Experiment:
    def __init__(self,configpath, load_experiment=True,h5py_path = None,workdir = 'work' ,experiment_name=None, Cellfile=None, Geomin = None, regex_data = None,partition=None,ncores=32) -> None:
        self.jobs = []   
        if load_experiment:
            with open(configpath,'r') as f:
                self.config = json.load(f)
        else:
            if experiment_name == None:
                experiment_name = input('Enter experiment name: ')
            self.config = setup_experiment(configpath,experiment_name)
            if Cellfile == None:
                Cellfile = ask_for_input_until_file_exists('enter cell file: ')
            if Geomin == None:
                Geomin = ask_for_input_until_file_exists('enter geom file: ')
            if regex_data == None:
                regex_data = input('enter regex matching data: ')
            if partition == None:
                valid_partitions = parse_sinfo_to_dataframe()
                print('Partition suggestions are:',valid_partitions)
                f = False
                while not f:
                    partition = input('Enter partition: ')
                    if partition in valid_partitions:
                        f = True
                    else:
                        print('Partition not in suggestions, are you sure? [y/n]')
                        if input() == 'y':
                            f = True
            self.config['sbatch_default'] = get_sbatch_standard()
            self.config['sbatch_default']['-p'] = partition
            self.config['sbatch_default']['-c'] = str(ncores)
            self.config['cell'] = Cellfile
            self.config['regex_data'] = regex_data
            self.config['geom'] = Geomin
            self.config['configpath'] = configpath
            self.config['workpath'] = os.path.join(os.path.dirname(configpath),workdir)
            self.config['h5py_path'] = h5py_path    
            self.config['runs'] = []
            os.makedirs(self.config['workpath'],exist_ok=True)
            save_config(self.config,configpath)
    
    def cat_runs(self,regex,prefix='concatenated_stream',outpath=None):
        import re
        if outpath == None:
            
            outpath = os.path.join(self.config['workpath'],'concatenated_streams', prefix + '_{i}.stream')
            i = 0
            while os.path.exists(outpath.format(i=i)):
                i += 1
            outpath = outpath.format(i=i)
            run_id = prefix + f'_{i}.stream'
            outdir = os.path.dirname(outpath)
        os.makedirs(outdir,exist_ok=True)
        run_names = self.config['runs']
        if isinstance(regex,str):
            regex = [regex]
        assert isinstance(regex,list)
        all_runs = []
        for reg in regex:
            reg = re.compile(reg)
            runs = [run for run in run_names if reg.match(run)]
            all_runs += runs
        all_runs = list(set(all_runs))
        for run in all_runs:
            if self.config[run]['Ran_distributed'] and not os.path.exists(self.config[run]['Stream_name']):
                self.cat_files(run)
        stream_names = [self.config[run]['Stream_name'] for run in all_runs]
        cat_files(stream_names,outpath)
        self.config['runs'].append(run_id)
        self.config[run_id] = {'stream_in': stream_names, 'Stream_name': outpath}
        save_config(self.config,self.config['configpath'])
        print('concatenated streams:',stream_names)
        print('output:',outpath)

    def cut_stream(self,streamfile,n,runnumber = None):
        if runnumber == None:
            i = 0 
            while True:
                outdir = os.path.join(self.config['workpath'],f'stream_cut_{i}')
                if not os.path.exists(outdir):
                    break
                i += 1
            runnumber = i
        outdir = os.path.join(self.config['workpath'],f'stream_cut_{runnumber}')
        runid = f'stream_cut_{runnumber}'
        os.makedirs(outdir,exist_ok=True)
        streamout = os.path.join(outdir,f'stream_cut_{runnumber}.stream')
        self.config['runs'].append(runid)
        self.config[runid] = {'stream_in': streamfile, 'Stream_name': streamout, 'n_crystals': n}
        cut_stream_after_ncrystals(streamfile,streamout,n)
        save_config(self.config,self.config['configpath'])
        return runid

    def setup_run(self,list_in=None,runnumber=None,cell=None,geom=None,indexamajig_config=None,sbatch_parameters=None,prefix='run',copy_geom=True,sub_folder=None):
        if cell == None:
            cell = self.config['cell']
        if geom == None:
            geom = self.config['geom']
        if runnumber == None:
            i = 0 
            while True:
                if sub_folder != None:
                    outdir = os.path.join(self.config['workpath'],sub_folder,f'{prefix}_{i}')
                else:
                    outdir = os.path.join(self.config['workpath'],f'{prefix}_{i}')
                if not os.path.exists(outdir):
                    break
                i += 1
            runnumber = i
        run_id = f'{prefix}_{runnumber}'
        if list_in == None:
            list_in = self.check_if_standard_list_exists()
        
        os.makedirs(outdir,exist_ok=True)
        list_local = os.path.join(outdir,f'{run_id}.lst')
        shutil.copy(list_in,list_local)
        streamout = os.path.join(outdir,f'{run_id}.stream')
        self.indexamajig_config = get_indexamajiq_parameters()
        if indexamajig_config != None:
            self.indexamajig_config.update(indexamajig_config)
        self.sbatch_parameters = self.config['sbatch_default']
        self.sbatch_parameters['--chdir'] = outdir
        if sbatch_parameters != None:
            self.sbatch_parameters.update(sbatch_parameters)
        if copy_geom:
            geom_internal = os.path.join(outdir,geom.split('/')[-1])
            shutil.copy(geom,geom_internal)
            geom = geom_internal
        if 'detector_distance' in self.indexamajig_config.keys():
            new_detector_distance = self.indexamajig_config['detector_distance']
            new_geom = geom.replace('.geom',f'_distance_{new_detector_distance}.geom')
            to_edit =  {'clen': '{i:.4f}'.format(i=new_detector_distance/1000)}
            edit_geom(geom, to_edit, new_geom)
            geom = new_geom
        IO_config = dict()
        IO_config['-i'] = list_local
        IO_config['-o'] = streamout
        IO_config['-p'] = cell
        IO_config['-g'] = geom

        log = os.path.join(outdir,f'{run_id}.log')
        self.sbatch_parameters['--output'] = log
        self.sbatch_parameters['--error'] = log
        run = {'Run_started': False, 'Run_finished': False, 'Stream_name': streamout, 'Run_path': outdir, 
        'Run_config': self.indexamajig_config, 'Run_sbatch': self.sbatch_parameters,'IO_config': IO_config,
        'Ran_distributed': False}
        run['n_images'] = get_len_text_file(list_local)
        self.config[run_id] = copy.deepcopy(run)
        self.config['runs'].append(run_id)
        save_config(self.config,self.config['configpath'])
        return run_id
    
    def execute_run(self,run_id):
        run = self.config[run_id]
        self.jobs.append(run_single_indexamajiq(run['Run_sbatch'],run['Run_config'],run['IO_config']))
        run['Run_started'] = True
    
    def execute_run_split(self,run_id,nchunks=None,images_per_chunk=None,wait_until_done=True,chunk_subfolder='chunks'):
        if nchunks == None and images_per_chunk == None:
            raise ValueError('Either nchunks or images_per_chunk must be given')
        run = self.config[run_id]
        if nchunks == None:
            nchunks = int(np.ceil(run['n_images']/images_per_chunk))
            print(f'run {run_id} will be split into {nchunks} chunks, it has {run["n_images"]} images')
        run = self.config[run_id]
        run_path = run['Run_path']
        lists = split_list_file(run['IO_config']['-i'],nchunks,list_sub_dir='lists')
        run['Run_started'] = True
        run['Ran_distributed'] = True
        run['nchunks'] = nchunks
        stream_parts = []
        for i,list_ in enumerate(lists):
            outdir = os.path.join(run_path,chunk_subfolder,f'chunk_{i}')
            os.makedirs(outdir,exist_ok=True)
            streamout = os.path.join(outdir,f'chunk_{i}.stream')
            log = os.path.join(outdir,f'chunk_{i}.log')
            stream_parts.append(streamout)
            if '--mille-dir' in run['Run_config']:
                run['Run_config']['--mille-dir'] = outdir + '/mille'
            run['IO_config']['-i'] = list_
            run['IO_config']['-o'] = streamout
            run['Run_sbatch']['--output'] = log
            run['Run_sbatch']['--error'] = log
            self.jobs.append(run_single_indexamajiq(run['Run_sbatch'],run['Run_config'],run['IO_config']))
        run['stream_parts'] = stream_parts
        self.config[run_id] = copy.deepcopy(run)
        save_config(self.config,self.config['configpath'])
        if wait_until_done:
            self.wait_until_done()
            self.cat_files(run_id)
    
    def cat_files(self,run_id):
        run = self.config[run_id]
        cat_files(run['stream_parts'],run['Stream_name'])
    
    def edit_sbatch_default(self,sbatch_parameters: dict):
        self.config['sbatch_default'].update(sbatch_parameters)
        save_config(self.config,self.config['configpath'])
    
    def add_sbatch_default(self):
        self.config['sbatch_default'] = get_sbatch_standard()
        save_config(self.config,self.config['configpath'])

    def get_last_run(self):
        return self.config['runs'][-1]
    
    def setup_list(self,string_to_match=None):
        if string_to_match == None:
            string_to_match = self.config['regex_data']
        list_out = os.path.join(self.config['workpath'],'list_files.lst')
        make_list(string_to_match,list_out)
        if self.config['h5py_path'] == None:
            h5py_file = glob(string_to_match)[0]
            self.config['h5py_path'] = find_key_h5_files(h5py_file)
            save_config(self.config,self.config['configpath'])
        list_out_all = get_all_events_smart(list_out,self.config['h5py_path'])
        return list_out_all

    def setup_partialator(self,runid = None, stream_in=None, sbatch_config = None, partialator_config = None,save_config_=True,prefix='partialator',outpath=None,subfolder=None):  
        if runid == None and stream_in == None:
            stream_in = input('No runid or stream_in given, please enter stream_path:')
        elif stream_in == None:
            stream_in = self.config[runid]['Stream_name']
        workpath = self.config['workpath']
        if outpath != None:
            workpath = outpath
        i = 0 
        if subfolder != None:
            workpath = os.path.join(workpath,subfolder)
        while os.path.exists(workpath + f'/{prefix}_{i}'):
            i += 1
        if not '-y' in partialator_config.keys():
            try: 
                pg = self.config['pg']
                partialator_config['-y'] = pg
            except:
                pg = input('Enter partialator pg: ')
                partialator_config['-y'] = pg
                self.config['pg'] = pg
        outpath = workpath + f'/{prefix}_{i}'
        outhkl = outpath + f'/{prefix}_{i}.hkl'
        self.sbatch_config = self.config['sbatch_default']
        self.sbatch_config['--output'] = outpath + f'/{prefix}_{i}.log'
        self.sbatch_config['--error'] = outpath + f'/{prefix}_{i}.log'
        if sbatch_config != None:
            self.sbatch_config.update(sbatch_config)
        os.makedirs(outpath,exist_ok=True)
        self.partialator_config = dict()
        self.partialator_config['-i'] = stream_in
        self.partialator_config['-o'] = outhkl
        self.partialator_config['-j'] = self.sbatch_config['-c']
        if partialator_config != None:
            self.partialator_config.update(partialator_config)
        run = {'Run_started': False, 'Run_finished': False, 'Partialator_out': outhkl, 
               'Run_path': outhkl, 'Partialator_config': self.partialator_config, 'Partialator_sbatch': self.sbatch_config}
        run_id = f'{prefix}_{i}'
        self.config[run_id] = run.copy()
        self.config['runs'].append(run_id)
        if save_config_:
            save_config(self.config,self.config['configpath'])
        return run_id

    def execute_partialator(self,run_id):
        run = self.config[run_id]
        self.jobs.append(run_partialator(run['Partialator_sbatch'],run['Partialator_config']))
        run['Run_started'] = True
        
    def check_if_standard_list_exists(self):
        list_out = os.path.join(self.config['workpath'],'list_files_all_events.lst')
        if not os.path.exists(list_out):
            list_out = self.setup_list()
        return list_out
    
    def screen_indexing_parameters(self,screen,sbatch_parameters=None,split_job=False,list_in=None,sample_size=5000,subfolder=None,wait_until_done=True,ignore_repeats=True):
        if list_in == None:
            list_in = self.check_if_standard_list_exists()
        if os.path.exists(list_in.replace('.lst',f'.shuffled_sample_{sample_size}.lst')):
            list_in = list_in.replace('.lst',f'.shuffled_sample_{sample_size}.lst')
        else:
            list_in = shuffle_lines_list(list_in,sample_size)
        product_values = product(*[v if isinstance(v, (list, tuple)) else [v] for v in screen.values()])
        params = [dict(zip(screen.keys(), values)) for values in product_values]
        for param in params:
            base_config = get_indexamajiq_parameters()
            base_config.update(param)
            indexamajig_config = base_config
            if ignore_repeats:
                if any([config_equal(base_config,self.config[run]['Run_config']) for run in self.config['runs']]):
                    continue
            run_id = self.setup_run(list_in,indexamajig_config=indexamajig_config,sbatch_parameters=sbatch_parameters,prefix='screening',sub_folder=subfolder)
            if split_job:
                self.execute_job_split(run_id)
            else:
                self.execute_run(run_id)
            # sleep(1)
        if wait_until_done: 
            self.wait_until_done()

    def run_idexamajig_for_mille(self,indexamajig_config=None,sbatch_parameters=None,nframes=5000,list_in=None,split_chunk=None,sub_folder=None):   
        if indexamajig_config == None:
            indexamajig_config = dict()
        indexamajig_config['--mille'] = ''
        if list_in == None:
            list_in = self.check_if_standard_list_exists()
        list_in = shuffle_lines_list(list_in,nframes)
        if sub_folder == None:
            mille_out_base = os.path.join(self.config['workpath'],'mille_{i}')
        else:
            mille_out_base = os.path.join(self.config['workpath'],sub_folder,'mille_{i}')
        i = 0
        while os.path.exists(mille_out_base.format(i=i)):
            i += 1
        indexamajig_config['--mille-dir'] = mille_out_base.format(i=i)
        runnr = self.setup_run(list_in=list_in,indexamajig_config=indexamajig_config,sbatch_parameters=sbatch_parameters,sub_folder=sub_folder)
        if split_chunk == None:
            self.execute_run(runnr)
        else:
            self.execute_run_split(runnr,split_chunk)


    def run_align_detector(self,mille_files=None,geom_in=None):
        i = 0 
        while os.path.exists(self.config['workpath'] + f'/align_detector_{i}'):
            i += 1
        outdir = os.path.join(self.config['workpath'],f'align_detector_{i}')
        os.makedirs(outdir,exist_ok=True)
        geom_out = os.path.join(outdir,f'align_detector_{i}.geom')
        if mille_files == None:
            i = 0
            while os.path.exists(self.config['workpath'] + f'/mille_{i+1}'):
                i += 1
            mille_files = self.config['workpath'] + f'/mille_{i}/*.bin'
            if len(glob(mille_files)) == 0:
                print('No mille files found')
                return
        log_file = os.path.join(outdir,f'align_detector_{i}.log')
        if geom_in == None:
            geom_in = self.config['geom']
        shutil.copy(geom_in,outdir)
        sbatch_params = self.config['sbatch_default'].copy()
        sbatch_params['--output'] = log_file
        run_align_detector(geom_in,geom_out,mille_files,sbatch_parameters=sbatch_params)
        if not 'old_geoms'in self.config:
            self.config['old_geoms'] = []
        self.config['old_geoms'].append(geom_in)
        self.config['geom'] = geom_out
        save_config(self.config,self.config['configpath'])

    def wait_until_done(self):
        slurm_tools.wait_until_queue_empty(self.jobs)
    
    def summarize_runs(self,key_to_summarize='screening'):
        run_ids = [key for key in self.config.keys() if key_to_summarize in key]
        res = []
        for run_id in run_ids:
            run = self.config[run_id]
            n_indexed = get_n_indexed_from_stream(run['Stream_name'])
            config = run['Run_config']
            log = run['Run_sbatch']['--output']
            time = read_runtime_from_log(log)
            config['run_time'] = time
            config['run_id'] = run_id
            config['n_indexed'] = [n_indexed]
            df = pd.DataFrame(config)
            res.append(df)
        df = pd.concat(res)
        df = df.set_index('run_id')
        return df

class scan_clen:
    def __init__(self, run_list: str,init_geom: str,range_to_scan_mm: tuple,cell_file: str,
                 indexamajiq_config= None, sbatch_parameters=None, optimize_single_runs = False,
                min_sample=1000,path_out= './') -> None:
        self.run_list = run_list
        self.all_images = get_all_events_smart(self.run_list)
        self.geom = init_geom
        self.path_out = path_out
        self.range_to_scan = range_to_scan_mm
        self.stream_dir = os.path.join(self.path_out,'streams')
        ref_configs = setup_standard_config()
        self.jobs = []
        if isinstance(indexamajiq_config,dict):
            pass
        elif isinstance(indexamajiq_config,str):
            with open(indexamajiq_config) as f:
                indexamajiq_config = json.load(f)
        else:
            indexamajiq_config = ref_configs[0]
            save_config(indexamajiq_config,os.path.join(self.path_out,'indexamajiq_config.json')) 
        if isinstance(sbatch_parameters,dict):
            pass
        elif isinstance(sbatch_parameters,str):
            with open(sbatch_parameters) as f:
                sbatch_parameters = json.load(f)
        else:
            sbatch_parameters = ref_configs[1]
            save_config(sbatch_parameters,os.path.join(self.path_out,'sbatch_parameters.json'))
        self.indexamajiq_config = indexamajiq_config
        self.sbatch_parameters = sbatch_parameters
        if not optimize_single_runs:
            self.sample = shuffle_lines_list(self.all_images,min_sample)
        self.indexamajiq_config['-p'] = cell_file
        
    def setup_geom(self):
        self.dir_geom_out = os.path.join(self.path_out,'geoms')
        os.makedirs(self.dir_geom_out,exist_ok=True)
        for i in range(*self.range_to_scan):
            outfile = os.path.join(self.dir_geom_out,f'clen_{i}.geom')
            to_change = {'clen': '{i:.4f}'.format(i=i/1000)}
            edit_geom(self.geom,to_change,outfile)
    
    def run_indexing(self):
        self.stream_files = []
        os.makedirs(self.stream_dir,exist_ok=True)
        for geom in glob(self.dir_geom_out + '/*.geom'):
            name = os.path.join(self.stream_dir,geom.split('/')[-1].replace('.geom',''))
            self.indexamajiq_config,self.sbatch_parameters = fill_name_config(name,self.indexamajiq_config,self.sbatch_parameters)
            self.indexamajiq_config['-g'] = os.path.join(self.dir_geom_out,geom)
            self.indexamajiq_config['-i'] = self.sample
            self.stream_files.append(self.indexamajiq_config['-o'])
            print(self.indexamajiq_config)
            self.jobs.append(run_single_indexamajiq(self.sbatch_parameters,self.indexamajiq_config))
    
    def wait_until_done(self):
        slurm_tools.wait_until_queue_empty(self.jobs)
    
    def get_results(self,outpath=None):
        import matplotlib.pyplot as plt
        clen = []
        indexed = []
        for file in self.stream_files:
            cell = get_n_indexed_from_stream(file)
            print(file,cell)
            clen.append(float(file.split('/')[-1].split('_')[-1].split('.')[0]))
            indexed.append(cell)
        df = pd.DataFrame({'clen': clen, 'indexed': indexed})
        df.to_csv(os.path.join(self.path_out,'clen_scan_results.csv'))
        plt.plot(clen,indexed)
        if outpath == None:
            outpath = os.path.join(self.path_out,'clen_scan_results.png')
        plt.savefig(outpath)

