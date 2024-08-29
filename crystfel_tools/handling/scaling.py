import pandas as pd
import numpy as np
from crystfel_tools.handling import data_conversion
from crystfel_tools.handling import fast_math
from scipy.optimize import minimize



class scale_it:
    def __init__(self,df,photon_energy,bandwidth,cell,pg='-1',startalpha=0,startbeta=0):
        self.df_raw = df
        self.df_split = data_conversion.split_df(df)
        c = 3e18 # speed of light in Angstrom/s
        h = 6.62607015e-34 # Planck constant in Js
        self.cell = cell
        e = photon_energy * 1.60217662e-19 # photon energy in Joules
        self.ech = e / (h * c) 
        self.deltaech = bandwidth * self.ech
        self.startalpha = startalpha
        self.startbeta = startbeta
        self.pg = pg

    
    def load_crystal(self,nr):
        self.current_crystal = self.sanitize_crystal(self.df_split[nr])
        self.current_crystal = self.current_crystal.loc[self.current_crystal.idx_hkl.isin(self.reference.idx_hkl)]
        self.current_ref = self.reference.loc[self.reference.idx_hkl.isin(self.current_crystal.idx_hkl)]
        self.current_ref = self.current_ref.sort_values('idx_hkl')
        self.current_crystal = self.current_crystal.sort_values('idx_hkl')

    def sanitize_crystal(self,crystal):
        crystal = crystal.copy()
        crystal = crystal.loc[crystal['I']>0]
        crystal = crystal.groupby('idx_hkl').mean().reset_index()
        return crystal

    def get_mean(self):
        return fast_math.calc_mean(self.data * self.brightness,self.weigths/self.brightness)

    def set_reference(self,reference_hkl):
        self.reference = data_conversion.load_reference(reference_hkl,self.df_raw,pg=self.pg)
        self.reference = self.reference.loc[self.reference['I']>0]
        

    def scale_current_crystal(self):
        factor0 = np.array([1,self.startalpha,self.startbeta,1])
        bounds = [(1e-32,None),(-1*np.pi,3*np.pi),(-1*np.pi,3*np.pi),(1e-5,1e5)]
        self.res = minimize(fast_math.calc_residual_new,factor0,args=(self.current_crystal[['h','k','l']].values,self.current_crystal['I'].values,self.current_ref['I'].values,self.cell,self.ech,self.deltaech)
                            ,bounds=bounds,method='L-BFGS-B')
        return self.res
    
    def scale_current_crystal_xsphere(self):
        factor0 = np.array([1,self.startalpha,self.startbeta,1])
        bounds = [(1e-32,None),(-1*np.pi,2*np.pi),(-1*np.pi,2*np.pi),(1e-5,1e5)]
        self.res = minimize(fast_math.residual_xpshere,factor0,args=(self.current_crystal[['h','k','l']].values,self.current_crystal['I'].values,self.current_ref['I'].values,self.cell,self.ech,self.deltaech)
                            ,bounds=bounds,method='L-BFGS-B')
        return self.res
    
    def apply_scaling_xsphere(self):
        self.I_scaled, self.weights, self.mask = fast_math.apply_scaling_xsphere(self.res.x,self.current_crystal[['h','k','l']].values,self.current_crystal['I'].values,self.cell,self.ech,self.deltaech)
    
    def apply_scaling(self):
        self.I_scaled, self.weights, self.mask = fast_math.evaluate_crystal_new(self.res.x,self.current_crystal[['h','k','l']].values,self.current_crystal['I'].values,self.cell,self.ech,self.deltaech)

    def screen_angles(self,size=100):
        angles = np.linspace(0,2*np.pi,size)
        res = np.zeros((size,size))
        for i, alpha in enumerate(angles):
            for j, beta in enumerate(angles):
                angles_current = np.array([alpha,beta]) 
                res[i,j] = fast_math.residual_simple(angles_current,self.current_crystal[['h','k','l']].values,self.current_crystal['I'].values,self.current_ref['I'].values,self.cell,self.ech,self.deltaech)
        min_idx = np.unravel_index(np.argmin(res, axis=None), res.shape)
        self.startalpha = angles[min_idx[0]]
        self.startbeta = angles[min_idx[1]]
        return self.startalpha, self.startbeta

    def check_scaling(self,print_results=True):
        cc_old = fast_math.calc_cc(self.current_crystal['I'].values[self.mask],self.current_ref['I'].values[self.mask])
        cc_new = fast_math.calc_cc(self.I_scaled,self.current_ref['I'].values[self.mask])
        rsplit_unweighted_old = fast_math.calc_rsplit_weighted(self.current_crystal['I'].values[self.mask],self.current_ref['I'].values[self.mask])
        rsplit_unweighted_new = fast_math.calc_rsplit_weighted(self.I_scaled,self.current_ref['I'].values[self.mask])
        rsplit_weighted_old = fast_math.calc_rsplit_weighted(self.current_crystal['I'].values[self.mask],self.current_ref['I'].values[self.mask],self.weights)
        rsplit_weighted_new = fast_math.calc_rsplit_weighted(self.I_scaled,self.current_ref['I'].values[self.mask],self.weights)
        if print_results:
            print('CC old:',cc_old)
            print('CC new:',cc_new)
            print('Rsplit unweighted old:',rsplit_unweighted_old)
            print('Rsplit unweighted new:',rsplit_unweighted_new)
            print('Rsplit weighted old:',rsplit_weighted_old)
            print('Rsplit weighted new:',rsplit_weighted_new)
        return cc_old,cc_new,rsplit_unweighted_old,rsplit_unweighted_new,rsplit_weighted_old,rsplit_weighted_new

    def write_out_half_datasets(self):
        return fast_math.write_out_half_datasets(self.data * self.brightness,self.weigths / self.brightness)






    
    
    


