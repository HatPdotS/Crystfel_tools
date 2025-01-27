from crystfel_tools.handling import io_functions
from collections import defaultdict 
import numpy as np
import pandas as pd
import multiprocessing as mp

class pystream:

    def __init__(self,stream: str) -> None:
        self.chunks = io_functions.read_entire_stream(stream)
        self.cell = io_functions.read_cell(stream,readable=False)
        self.geometry = io_functions.read_geometry_stream(stream)
        self.cell_readable = io_functions.read_cell(stream)
        self.lookup_table = self.build_lookup_table()
        self.preamble = io_functions.load_preamble_stream(stream)

    def build_lookup_table(self):
        lookup_table = defaultdict(list)
        for i,chunk in enumerate(self.chunks):
            id = i
            file = chunk['Image filename']
            event = chunk['Event']
            ncrystals = len(chunk['crystals'])
            crystal_id = np.nan
            if ncrystals > 0:
                for crystal_id,crystal in enumerate(chunk['crystals']):
                    lookup_table['id'].append(id)
                    lookup_table['file'].append(file)
                    lookup_table['event'].append(event)
                    lookup_table['ncrystals'].append(ncrystals)
                    lookup_table['crystal_id'].append(crystal_id)
                    lookup_table['a'].append(float(crystal['cell'][0]))
                    lookup_table['b'].append(float(crystal['cell'][1]))
                    lookup_table['c'].append(float(crystal['cell'][2]))
                    lookup_table['alpha'].append(float(crystal['cell'][3]))
                    lookup_table['beta'].append(float(crystal['cell'][4]))
                    lookup_table['gamma'].append(float(crystal['cell'][5]))
                    lookup_table['lattice_type'].append(crystal['lattice_type'])  
                    lookup_table['indexed_by'].append(chunk['indexed_by'])
                    lookup_table['n_indexing_tries'].append(chunk['n_indexing_tries'])
                    lookup_table['det_shift_x'].append(float(crystal['det_shift'][0]))
                    lookup_table['det_shift_y'].append(float(crystal['det_shift'][1]))
                    lookup_table['diffraction_resolution_limit'].append(float(crystal['diffraction_resolution_limit'].split()[-2]))
            else:
                lookup_table['id'].append(id)
                lookup_table['file'].append(file)
                lookup_table['event'].append(event)
                lookup_table['ncrystals'].append(ncrystals)
                lookup_table['crystal_id'].append(np.nan)
                lookup_table['a'].append(np.nan)
                lookup_table['b'].append(np.nan)  
                lookup_table['c'].append(np.nan)
                lookup_table['alpha'].append(np.nan)
                lookup_table['beta'].append(np.nan)
                lookup_table['gamma'].append(np.nan)
                lookup_table['lattice_type'].append('none')
                lookup_table['indexed_by'].append(chunk['indexed_by'])
                lookup_table['n_indexing_tries'].append('none')
                lookup_table['det_shift_x'].append(np.nan)
                lookup_table['det_shift_y'].append(np.nan)
                lookup_table['diffraction_resolution_limit'].append(np.nan)
        df = pd.DataFrame(lookup_table)
        self.lookup_table = df
        return df

    def add_stream(self,stream):
        self.chunks.extend(io_functions.read_entire_stream(stream))

    def add_streams_mp(self,streams):
        pool = mp.Pool(32)
        chunks = pool.map(io_functions.read_entire_stream,streams)
        self.chunks.extend(sum(chunks,[]))

    def clean_lookup_table(self):
        self.lookup_table = self.lookup_table.dropna()

    def write_stream(self,out: str):
        lines = self.preamble
        for chunk in self.lookup_table.id.unique():
            crystals_to_write = self.lookup_table.loc[self.lookup_table.id == chunk,'crystal_id'].tolist()
            lines.extend(io_functions.prep_lines_chunk(self.chunks[chunk],crystals_to_write=crystals_to_write,write_all_crystals=False))
        with open(out,'w') as f:
            f.writelines(lines)







        


        

