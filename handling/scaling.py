import sparse
import pandas as pd
import numpy as np
from handling import data_conversion
from handling import fast_math
import multiprocessing as mp

def select_crystal(i):
    return sp[:,i]

class scale_it:
    def __init__(self,df):
        self.sparse, self.weigths = data_conversion.make_sparse(df)
        self.res = df.resolution
    
    def split_crystal(self):
        sp = self.sparse
        with mp.Pool(56) as p:
            self.crystals = p.map(select_crystal,range(sp.shape[1]))
    
    def sanitize_crystal(self,crystal):
        crystal = crystal.copy()
        idx = np.where(crystal.data > 0)
        crystal.coords = crystal.coords[:,idx]
        crystal.data = crystal.data[idx]
        return crystal

        
    def get_mean(self):
        return fast_math.calc_mean(self.data * self.brightness,self.weigths/self.brightness)
    


    def scale_frames_brightness(self):
        mean = self.get_mean()
        normed = self.data/mean.reshape(-1,1)
        idx = np.where(normed.data < 0)
        normed.data[idx] = 0
        data_there = self.data_there.copy()
        data_there.data[idx] = 0
        scale_factor = np.sum(normed,axis=0).todense() / np.sum(data_there,axis=0).todense()
        scale_factor /= scale_factor.mean()
        self.brightness *= scale_factor
    
    def write_out_half_datasets(self):
        return fast_math.write_out_half_datasets(self.data * self.brightness,self.weigths / self.brightness)






    
    
    


