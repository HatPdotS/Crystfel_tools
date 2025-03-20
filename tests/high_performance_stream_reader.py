import pandas as pd
from time import time 
from multiprocessing import pool as Pool
from multiprocessing.pool import ThreadPool
import atexit
from collections import defaultdict

from multiprocessing import shared_memory, Lock ,Value

import pyarrow as pa
import orjson
import numpy as np


stream_file = '/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/chunk_1.stream'

def get_starting_positions(stream_file, nthreads):
    with open(stream_file, 'r') as f:
        f.seek(0, 2)  # Move to end of file
        file_size = f.tell()  # Get current position (which is the file size)
        starting_points = range(0,file_size,file_size//nthreads)
        starting_points_clean = []
        for starting_point in starting_points:
            f.seek(starting_point)
            while True:
                line = f.readline()
                if not line:
                    break
                if line.startswith('----- Begin chunk -----'):
                    starting_points_clean.append(f.tell())
                    break
    return sorted(list(set(starting_points_clean)))

def get_filesize(stream_file):
    with open(stream_file, 'r') as f:
        f.seek(0, 2)  # Move to end of file
        file_size = f.tell()  # Get current position (which is the file size)
    return file_size

def read_file_segment(args):
    file_path, start, end = args
    with open(file_path, "rb") as f:  # Open in binary mode to handle all file types
        f.seek(start)  # Move to the start position
        data = f.read(end - start)  # Read the specified range
    return data

def read_file_in_segments(stream_file,nthreads):
    starting_pos = get_starting_positions(stream_file, nthreads)
    file_size = get_filesize(stream_file)
    ends = starting_pos[1:] + [file_size]
    args = [(stream_file, start, end) for start, end in zip(starting_pos, ends)]
    with ThreadPool(nthreads) as p:
        data = p.map(read_file_segment, args)
    return data

def is_valid_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def read_chunks(chunk_of_data):
    chunks = []
    lines = chunk_of_data.decode().splitlines()
    iter_lines = iter(lines)
    while True:
        chunk = dict()
        try:
            line = next(iter_lines)
        except StopIteration:
            break
        chunk['filename'] = line.split()[-1]
        chunk['Event'] = int(next(iter_lines).split()[-1].strip('//'))
        chunk['Image serial number'] = int(next(iter_lines).split()[-1])
        for line in iter_lines:
            if line.startswith('Peaks from peak search'):
                break
            key, value = line.split('=')
            value = value.strip().split()[0]    
            try: 
                value = int(value)
            except:
                value = str(value)
            chunk[key.strip()] = value
        hits = defaultdict(list)
        line = next(iter_lines)
        for line in iter_lines:
            if line.startswith('End of peak list'):
                break
            line_split = line.split()
            hits['fs/px'].append(float(line_split[0]))
            hits['ss/px'].append(float(line_split[1]))
            hits['(1/d)/nm^-1'].append(float(line_split[2]))
            hits['Intensity'].append(float(line_split[3]))
            hits['Panel'].append(line_split[4])
        hits = pa.serialize_pandas(pd.DataFrame(hits))
        chunk['hits'] = hits
        chunk['crystals'] = []
        for line in iter_lines:
            if line.startswith('----- End chunk -----'):
                break
            if line.startswith('--- Begin crystal'):
                crystal = dict()
                line = next(iter_lines)
                cell = [float(i) for i in line.split() if is_valid_float(i)]
                crystal['cell'] = cell
                line = next(iter_lines)
                crystal['astar'] = [float(i) for i in line.split() if is_valid_float(i)]
                line = next(iter_lines)
                crystal['bstar'] = [float(i) for i in line.split() if is_valid_float(i)]
                line = next(iter_lines)
                crystal['cstar'] = [float(i) for i in line.split() if is_valid_float(i)]
                line = next(iter_lines)
                for line in iter_lines:
                    if line.startswith('Reflections measured after indexing'):
                        break
                    try:
                        key, value = line.split('=')
                        value = value.strip().split()[0]
                        try:
                            value = int(value)
                        except:
                            value = str(value)
                        crystal[key.strip()] = value
                    except:
                        key, valuex, valuey = line.split('=')
                        valuex = valuex.strip().split()[0]
                        valuey = valuey.strip().split()[0]
                        crystal[key.strip()] = [valuex, valuey]
                reflections = defaultdict(list)
                line = next(iter_lines)
                for line in iter_lines:
                    if line.startswith('End of reflections'):
                        break
                    line_split = line.split()
                    reflections['h'].append(int(line_split[0]))
                    reflections['k'].append(int(line_split[1]))
                    reflections['l'].append(int(line_split[2]))
                    reflections['I'].append(float(line_split[3]))
                    reflections['sigma'].append(float(line_split[4]))
                    reflections['peak'].append(float(line_split[5]))
                    reflections['background'].append(float(line_split[6]))
                    reflections['fs/px'].append(float(line_split[7]))
                    reflections['ss/px'].append(float(line_split[8]))
                    reflections['panel'].append(line_split[9])
                reflections = pd.DataFrame(reflections)
                crystal['reflections'] = pa.serialize_pandas(reflections)
                chunk['crystals'].append(crystal)
        chunk['ncrystals'] = len(chunk['crystals'])
        chunks.append(chunk)
        for line in iter_lines:
            if line.startswith('----- Begin chunk -----'):
                break
    return chunks

def read_chunks_shared_mem(args):
    chunk_of_data = args
    chunks = []
    lines = chunk_of_data.decode().splitlines()
    iter_lines = iter(lines)
    while True:
        chunk = dict()
        try:
            line = next(iter_lines)
        except StopIteration:
            break
        chunk['filename'] = line.split()[-1]
        chunk['Event'] = int(next(iter_lines).split()[-1].strip('//'))
        chunk['Image serial number'] = int(next(iter_lines).split()[-1])
        for line in iter_lines:
            if line.startswith('Peaks from peak search'):
                break
            key, value = line.split('=')
            value = value.strip().split()[0]    
            try: 
                value = int(value)
            except:
                value = str(value)
            chunk[key.strip()] = value
        hits = defaultdict(list)
        line = next(iter_lines)
        for line in iter_lines:
            if line.startswith('End of peak list'):
                break
            line_split = line.split()
            hits['fs/px'].append(float(line_split[0]))
            hits['ss/px'].append(float(line_split[1]))
            hits['(1/d)/nm^-1'].append(float(line_split[2]))
            hits['Intensity'].append(float(line_split[3]))
            hits['Panel'].append(float(line_split[4].strip('m')))
        hits_array = np.array((hits['fs/px'],hits['ss/px'],hits['(1/d)/nm^-1'],hits['Intensity'],hits['Panel']),dtype = np.float64)
        len_hits = len(hits['fs/px'])
        with lock:
            start = current_idx_hits.value
            current_idx_hits.value += len_hits
            end = current_idx_hits.value
        hitarray[start:end] = hits_array.T
        chunk['hits'] = slice(start,end)
        chunk['crystals'] = []
        for line in iter_lines:
            if line.startswith('----- End chunk -----'):
                break
            if line.startswith('--- Begin crystal'):
                crystal = dict()
                line = next(iter_lines)
                cell = [float(i) for i in line.split() if is_valid_float(i)]
                crystal['cell'] = cell
                line = next(iter_lines)
                crystal['astar'] = [float(i) for i in line.split() if is_valid_float(i)]
                line = next(iter_lines)
                crystal['bstar'] = [float(i) for i in line.split() if is_valid_float(i)]
                line = next(iter_lines)
                crystal['cstar'] = [float(i) for i in line.split() if is_valid_float(i)]
                line = next(iter_lines)
                for line in iter_lines:
                    if line.startswith('Reflections measured after indexing'):
                        break
                    try:
                        key, value = line.split('=')
                        value = value.strip().split()[0]
                        try:
                            value = int(value)
                        except:
                            value = str(value)
                        crystal[key.strip()] = value
                    except:
                        key, valuex, valuey = line.split('=')
                        valuex = valuex.strip().split()[0]
                        valuey = valuey.strip().split()[0]
                        crystal[key.strip()] = [valuex, valuey]
                reflections = defaultdict(list)
                line = next(iter_lines)
                for line in iter_lines:
                    if line.startswith('End of reflections'):
                        break
                    line_split = line.split()
                    reflections['h'].append(int(line_split[0]))
                    reflections['k'].append(int(line_split[1]))
                    reflections['l'].append(int(line_split[2]))
                    reflections['I'].append(float(line_split[3]))
                    reflections['sigma'].append(float(line_split[4]))
                    reflections['peak'].append(float(line_split[5]))
                    reflections['background'].append(float(line_split[6]))
                    reflections['fs/px'].append(float(line_split[7]))
                    reflections['ss/px'].append(float(line_split[8]))
                    reflections['panel'].append(int(line_split[9].strip('m')))
                reflections_array = np.array((reflections['h'],reflections['k'],reflections['l'],reflections['I'],reflections['sigma'],reflections['peak'],reflections['background'],reflections['fs/px'],reflections['ss/px'],reflections['panel']),dtype = np.float64)
                len_reflections = len(reflections['h'])
                with lock:
                    start = current_idx_indexed.value
                    current_idx_indexed.value += len_reflections
                    end = current_idx_indexed.value
                    indexed_array[start:end] = reflections_array.T
                crystal['reflections'] = slice(start,end)
                chunk['crystals'].append(crystal)
        chunk['ncrystals'] = len(chunk['crystals'])
        chunks.append(chunk)
        for line in iter_lines:
            if line.startswith('----- Begin chunk -----'):
                break
    return chunks

def read_stream(stream_file,nthreads = 56):
    data = read_file_in_segments(stream_file,nthreads)
    with Pool.Pool(processes = nthreads) as p:
        chunks = p.map(read_chunks, data)
    chunks = [item for sublist in chunks for item in sublist]
    return chunks

def read_stream_shared_array(stream_file,nthreads = 56):
    data = read_file_in_segments(stream_file,nthreads)
    file_size = get_filesize(stream_file)
    max_array_size = file_size // 5
    shared_array_shape = (max_array_size//5,5)
    dtype = np.dtype('float64')
    global hitarray,indexed_array,current_idx_hits,current_idx_indexed,lock,hitarray_shm,indexed_array_shm
    def cleanup_shared_memory():
        """Cleanup shared memory on exit."""
        global hitarray_shm, indexed_array_shm
        if hitarray_shm:
            hitarray_shm.close()
            hitarray_shm.unlink()
        if indexed_array_shm:
            indexed_array_shm.close()
            indexed_array_shm.unlink()
    atexit.register(cleanup_shared_memory)
    hitarray_shm = shared_memory.SharedMemory(create=True, size=np.prod(shared_array_shape) * np.dtype(dtype).itemsize)
    hitarray = np.ndarray(shared_array_shape, dtype=dtype, buffer=hitarray_shm.buf)
    hitarray.fill(0)
    shared_array_shape = (max_array_size//10,10)
    indexed_array_shm = shared_memory.SharedMemory(create=True, size=np.prod(shared_array_shape) * np.dtype(dtype).itemsize)
    indexed_array = np.ndarray(shared_array_shape, dtype=dtype, buffer=indexed_array_shm.buf)
    indexed_array.fill(0)
    current_idx_hits = Value('i', 0)
    current_idx_indexed = Value('i', 0)
    lock = Lock()
    with Pool.Pool(processes = nthreads) as p:
        chunks = p.map(read_chunks_shared_mem, data)
    chunks = [item for sublist in chunks for item in sublist]
    hitarray = hitarray[:current_idx_hits.value]
    indexed_array = indexed_array[:current_idx_indexed.value]

    return chunks, hitarray, indexed_array

def unpack_pds(chunks):
    for chunk in chunks:
        for crystal in chunk['crystals']:
            crystal['reflections'] = pa.deserialize_pandas(crystal['reflections'])
    pass

class pystream:
    def __init__(self,stream_file=None,nthreads = 56):
        self.nthreads = nthreads

        if isinstance(stream_file, str):
            self.stream_file = stream_file
            self.read_stream(stream_file)
        
    def read_stream(self,streamfile):
        self.stream, self.hit_array, self.indexed_array = read_stream_shared_array(streamfile,self.nthreads)
        self.hits = pd.DataFrame(self.hit_array,columns = ['fs/px','ss/px','(1/d)/nm^-1','Intensity','Panel'])
        self.indexed = pd.DataFrame(self.indexed_array,columns = ['h','k','l','I','sigma','peak','background','fs/px','ss/px','panel'])
        print('stream read')
        print('Chunks length',len(self.stream))
        print('Crystals length',sum([len(chunk['crystals']) for chunk in self.stream]))

    def universal_assigner(self, df, column, values, slices):
        df[column] = -1
        column_index = df.columns.get_loc(column)
        global assigner
        def assigner(args):
            column, value, slice = args
            df.iloc[slice,column] = value

        args = [(column_index, value, slice) for slice,value in zip(slices,values)]
        for arg in args:    
            assigner(arg)

    def universal_selector(self, df, slices):
        mega_slice = np.r_[tuple(range(s.start, s.stop) for s in slices)]
        return df.iloc[mega_slice]

    def move_to_main_process_memory(self):
        self.hits = self.hits.copy()
        self.indexed = self.indexed.copy()

    def label_reflections_runnr(self):
        runnrs = []
        s = []
        for chunk in self.stream:
            for crystal in chunk['crystals']:
                try:
                    runnrs.append(int(chunk['filename'].split('run')[-1].lstrip('_').split('.')[0].split('_')[0].split('-')[0]))
                    s.append(crystal['reflections'])
                except:
                    print(chunk['filename'])
        self.universal_assigner(self.indexed,'runnr',runnrs,s)
    
    def label_reflections_event(self):
        events = []
        s = []
        for chunk in self.stream:
            for crystal in chunk['crystals']:
                events.append(chunk['Event'])
                s.append(crystal['reflections'])
        self.universal_assigner(self.indexed,'Event',events,s)
                    
    def label_reflections_acqnr(self):
        acqnrs = []
        s = []
        for chunk in self.stream:
            for crystal in chunk['crystals']:
                acqnrs.append(int(chunk['filename'].split('acq')[-1].split('.')[0].split('_')[0].split('-')[0]))
                s.append(crystal['reflections'])
        self.universal_assigner(self.indexed,'acqnr',acqnrs,s)

    def cat_all_reflections(self):
        reflections = []
        for chunk in self.stream:
            for crystal in chunk['crystals']:
                reflections.append(crystal['reflections'])
        return pd.concat(reflections)
    
    def select_based_on_chunk_level_key(self,key,lam):
        slices = []
        for chunk in self.stream:
            if lam(chunk[key]):
                for crystal in chunk['crystals']:
                    slices.append(crystal['reflections'])
        return self.universal_selector(self.indexed,slices)
