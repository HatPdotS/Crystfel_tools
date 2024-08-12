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
        h,k,l,I,sigma,peak,background,_,_,_ = line.split()
        res['h'].append(int(h))
        res['k'].append(int(k))
        res['l'].append(int(l))
        res['I'].append(float(I))
        res['sigma'].append(float(sigma))
        res['peak'].append(float(peak))
        res['background'].append(float(background))
    return res

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




                    
