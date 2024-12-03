from collections import defaultdict
import pandas as pd
from multiprocessing import Pool


def read_single_hitlist(openfile):
    """
    Read a single hit from a file moves iterator to the endpostion of the hit
    """
    hits = defaultdict(list)
    for line in openfile:
        if line.strip('\n') == 'End of peak list':
            break
        fs, ps, res, I, panel = line.split()
        hits['fs'].append(float(fs))
        hits['ps'].append(float(ps))
        hits['res'].append(float(res))
        hits['I'].append(float(I))
        hits['panel'].append(panel)
    return hits

def read_crystal_reflections(open_file):
    res =  defaultdict(list)
    next(open_file)
    for line in open_file:
        if line.strip('\n') == 'End of reflections':
            break
        h,k,l,I,sigma,peak,background,fs,ps,panel = line.split()
        res['h'].append(int(h))
        res['k'].append(int(k))
        res['l'].append(int(l))
        res['I'].append(float(I))
        res['sigma'].append(float(sigma))
        res['peak'].append(float(peak))
        res['background'].append(float(background))
        res['fs'].append(float(fs))
        res['ps'].append(float(ps))
        res['panel'].append(panel)

    return res



def write_hkl(file_out,df,Symmetry):
    triplet = df[['h','k','l']].values
    mean = df['I'].values
    std = df['sigma'].values
    nmeas = df['nmeas'].values
    with open(file_out,'w') as f:
        # header = '!ITEM_H=1\n!ITEM_K=2\n!ITEM_L=3\n!ITEM_IOBS=4\n!ITEM_SIGMA(IOBS)=5\n!END_OF_HEADER\n'
        # f.write(header)
        f.write('CrystFEL reflection list version 2.0\n')
        f.write(f'Symmetry: {Symmetry}\n')
        f.write('   h    k    l          I    phase   sigma(I)   nmeas\n')
        for i in range(len(triplet)):
            if str(mean[i]).lower() == 'nan':
                continue
            s = f'{triplet[i][0]:>4}{triplet[i][1]:>5}{triplet[i][2]:>5}{mean[i]:>11.2f}        -{std[i]:>11.2f}{nmeas[i]:>8}\n'
            f.write(s)
        f.write('End of reflections')

def read_crystal(open_file):
    crystal = {}
    cell_line = next(open_file).split()
    crystal['cell'] = [cell_line[idx] for idx in [2,3,4,6,7,8]]
    for line in open_file:
        if line.strip('\n') == '--- End crystal':
            break
        if line.strip('\n') == 'Reflections measured after indexing':
            crystal['reflections'] = read_crystal_reflections(open_file)
            continue
        if line[:24] == 'predict_refine/det_shift':
            split = line.split()
            crystal['det_shift'] = (split[3],split[6])
            continue
        try:
            key, value = line.split('=')
            crystal[key.strip()] = value.strip()
        except ValueError:
            print(line)
    return crystal

def read_chunk(open_file):
    chunk = {}
    ncryst = 0
    chunk['crystals'] = []
    for _ in range(3):
        key, value = next(open_file).split(':')
        chunk[key.strip()] = value.strip()
    for line in open_file:
        if line.strip('\n') == '----- End chunk -----':
            break
        if line.strip('\n') == 'Peaks from peak search':
            next(open_file)
            chunk['hits'] = read_single_hitlist(open_file)
            continue
        if '--- Begin crystal' == line.strip('\n'):
            ncryst += 1
            chunk['crystals'].append(read_crystal(open_file))
            continue
        try:
            key, value = line.split('=')
            chunk[key.strip()] = value.strip()
        except ValueError:
            print(line)
    return chunk

def prep_lines_crystal(crystal):
    lines = []
    lines.append('--- Begin crystal\n')
    cell = crystal['cell']
    cell_line = f'Cell parameters {cell[0]} {cell[1]} {cell[2]} nm, {cell[3]} {cell[4]} {cell[5]} deg'  + '\n'
    lines.append(cell_line)
    lines.append(' = '.join(['astar',crystal['astar']]) + '\n')
    lines.append(' = '.join(['bstar',crystal['bstar']]) + '\n')
    lines.append(' = '.join(['cstar',crystal['cstar']]) + '\n')
    lines.append(' = '.join(['lattice_type',crystal['lattice_type']]) + '\n')
    lines.append(' = '.join(['centering',crystal['centering']]) + '\n')
    lines.append(' = '.join(['unique_axis',crystal['unique_axis']]) + '\n')
    lines.append(' = '.join(['profile_radius',crystal['profile_radius']]) + '\n')
    lines.append(' = '.join(['predict_refine/final_residual',crystal['predict_refine/final_residual']]) + '\n')
    lines.append(' = '.join(['predict_refine/det_shift x',crystal['det_shift'][0]]) + ' = '.join([' y',crystal['det_shift'][1]]) + ' mm' + '\n')
    lines.append(' = '.join(['diffraction_resolution_limit',crystal['diffraction_resolution_limit']]) + '\n')
    lines.append(' = '.join(['num_reflections',crystal['num_reflections']]) + '\n')
    lines.append(' = '.join(['num_saturated_reflections',crystal['num_saturated_reflections']]) + '\n')
    lines.append(' = '.join(['num_implausible_reflections',crystal['num_implausible_reflections']]) + '\n')
    lines.append('Reflections measured after indexing\n')
    lines.append('   h    k    l          I   sigma(I)       peak background  fs/px  ss/px panel\n')
    for i in range(len(crystal['reflections']['h'])):
        h = crystal['reflections']['h'][i]
        k = crystal['reflections']['k'][i]
        l = crystal['reflections']['l'][i]
        I = crystal['reflections']['I'][i]
        sigma = crystal['reflections']['sigma'][i]
        peak = crystal['reflections']['peak'][i]
        background = crystal['reflections']['background'][i]
        fs = crystal['reflections']['fs'][i]
        ps = crystal['reflections']['ps'][i]
        panel = crystal['reflections']['panel'][i]
        lines.append(f'{h:>4}{k:>5}{l:>5}{I:>11.2f}{sigma:>11.2f}{peak:>11.2f}{background:>11.2f}{fs:>7.1f}{ps:>7.1f}{panel:>3}\n')
    lines.append('End of reflections\n')
    lines.append('--- End crystal\n')
    return lines

def load_preamble_stream(file):
    lines = []
    with open(file) as f:
        for line in f:
            if line.strip('\n') == '----- Begin chunk -----':
                break
            lines.append(line)
    return lines

def prep_lines_chunk(chunk,write_all_crystals=True,crystals_to_write=[]):
    lines = []
    lines.append('----- Begin chunk -----\n')
    lines.append('Image filename: ' + chunk['Image filename'] + '\n')
    lines.append('Event: ' + chunk['Event'] + '\n')
    lines.append('Image serial number: ' + chunk['Image serial number'] + '\n')
    lines.append('hit = ' + chunk['hit'] + '\n')
    lines.append('indexed_by = ' + chunk['indexed_by'] + '\n')
    if 'n_indexing_tries' in chunk:
        lines.append('n_indexing_tries = ' + chunk['n_indexing_tries'] + '\n')
    lines.append('photon_energy_eV = ' + chunk['photon_energy_eV'] + '\n')
    lines.append('beam_divergence = ' + chunk['beam_divergence'] + '\n')
    lines.append('beam_bandwidth = ' + chunk['beam_bandwidth'] + '\n')
    lines.append('average_camera_length = ' + chunk['average_camera_length'] + '\n')
    lines.append('num_peaks = ' + chunk['num_peaks'] + '\n')
    lines.append('peak_resolution = ' + chunk['peak_resolution'] + '\n')
    lines.append('Peaks from peak search\n')
    lines.append('  fs/px   ss/px (1/d)/nm^-1   Intensity  Panel\n')
    for i in range(len(chunk['hits']['fs'])):
        fs = chunk['hits']['fs'][i]
        ps = chunk['hits']['ps'][i]
        res = chunk['hits']['res'][i]
        I = chunk['hits']['I'][i]
        panel = chunk['hits']['panel'][i]
        lines.append(f'{fs:>7.2f}{ps:>8.2f}{res:>11.2f}{I:>12.2f}{panel:>5}\n') 
    lines.append('End of peak list\n')
    if write_all_crystals:
        for crystal in chunk['crystals']:
            lines.extend(prep_lines_crystal(crystal))
    else:
        for i in crystals_to_write:
            if i >= len(chunk['crystals']):
                print('Crystal index out of range, please check which crystals you want to write')
                continue
            lines.extend(prep_lines_crystal(chunk['crystals'][int(i)]))
    lines.append('----- End chunk -----\n')
    return lines

def read_cell(file,readable=True):
    """
    Read the cell from a stream.
    """

    cell = {}
    if readable:
        with open(file) as f:
            for line in f:
                if line.strip('\n') == '----- Begin unit cell -----' or line.strip('\n') == 'CrystFEL unit cell file version 1.0':
                    for line in f:
                        if line.strip('\n') == '----- End unit cell -----':
                            break
                        try: 
                            key, value = line.split('=')
                            cell[key.strip()] = float(value.strip().split()[0])
                        except: pass
                    break
        return list(cell.values())
    else:
        with open(file) as f:
            for line in f:
                if line.strip('\n') == '----- Begin unit cell -----':
                    for line in f:
                        if line.strip('\n') == '----- End unit cell -----':
                            break
                        try: 
                            key, value = line.split('=')
                            cell[key.strip()] = value.strip()
                        except: pass
                    break
        return cell

def prep_lines_cell(cell):
    lines = []
    lines.append('----- Begin unit cell -----\n')
    lines.append('CrystFEL unit cell file version 1.0\n')
    for key in cell.keys():
        lines.append(f'{key} = {cell[key]} \n')
    lines.append('----- End unit cell -----\n')
    return lines
            
def read_geometry_stream(file):
    """
    Read the geometry from a stream.
    """
    with open(file) as f:
        geom = {}
        for line in f:
            if line.strip('\n') == '----- Begin geometry file -----':
                for line in f: 
                    if line.strip('\n') == '----- End geometry file -----':
                        break
                    if line[0] == ';' or not line.strip():
                        continue
                    key, value = line.split('=')
                    geom[key.strip()] = value.strip()

                break        
    return geom

def prep_lines_geometry(geom):
    lines = []
    lines.append('----- Begin geometry file ----- \n')
    for key in geom.keys():
        lines.append(f'{key:<10} = {geom[key]:>10} \n')
    lines.append('----- End geometry file ----- \n')
    return lines

def read_unitcell_stream(file):
    """
    Read the unitcell from a stream.
    """
    with open(file) as f:
        for line in f:
            if line.strip('\n') == '----- Begin unit cell -----':
                return read_cell(f)


def read_entire_stream(file):
    """
    Read the entire stream.
    """
    chunks = []
    with open(file) as f:
        for line in f:
            if line.strip('\n') == '----- Begin chunk -----':
                chunks.append(read_chunk(f))
    return chunks

def get_number_of_crystals(chunks):
    return sum([len(chunk['crystals']) for chunk in chunks])

def read_only_crystals(file):
    """read only crystals from stream"""
    crystals = []
    with open(file) as f:
        for line in f:
            if line.strip('\n') == '----- Begin crystal':
                crystal = read_crystal(f)
                crystals.append(crystal)


def filter_for_crystals(stream):
    """
    Filter the stream for crystals.
    """
    return [chunk for chunk in stream if len(chunk['crystals']) > 0]

def convert_chunk_to_df(args):
    chunk, crystal, idx = args
    crystal_data = pd.DataFrame(crystal['reflections'])
    crystal_data['idx'] = idx
    data = pd.DataFrame({ 'filename': chunk['Image filename'],'event': chunk['Event'],'idx':[idx]})
    return crystal_data, data




def build_dataframe_generator(chunks):
    counter = 0
    for chunk in chunks:
        for crystal in chunk['crystals']:
            yield chunk, crystal, counter
            counter += 1

def build_dataframe(stream):
    """
    Convert a stream to a DataFrame.
    """
    with Pool(10) as p:
        data, df = zip(*p.map(convert_chunk_to_df, build_dataframe_generator(stream)))
    data = pd.concat(data)
    df = pd.concat(df)
    return data, df

def read_crystfel_hkl(path):
    ref = pd.read_csv(path, sep=r'\s+', names=['h','k','l','I','phase','sigma','nmeas'],skiprows=3,on_bad_lines='skip').dropna()
    ref.h = ref.h.astype(int)
    ref.k = ref.k.astype(int)
    ref.l = ref.l.astype(int)
    ref.I = ref.I.astype(float)
    ref.sigma = ref.sigma.astype(float) 
    ref.nmeas = ref.nmeas.astype(int)
    return ref



                    
