import pandas as pd
import sparse
import numpy as np
from crystfel_tools.handling import io_functions

def load_reference(hkl,ref_for_relabeling,pg='-1'):
    """
    Load the reference hkl file
    """
    ref_for_relabeling = ref_for_relabeling.groupby(['h','k','l']).first().reset_index()
    ref_for_relabeling['triplet'] = ref_for_relabeling['h'].astype(str) + '_' + ref_for_relabeling['k'].astype(str) + '_' + ref_for_relabeling['l'].astype(str)
    
    ref = io_functions.read_crystfel_hkl(hkl)
    ref = extend_reference_to_full_ewaldsphere(ref,pg=pg)
    ref['triplet'] = ref['h'].astype(str) + '_' + ref['k'].astype(str) + '_' + ref['l'].astype(str)
    ref['idx_hkl'] = np.nan
    ref['idx_hkl'] = ref['triplet'].map(ref_for_relabeling.set_index('triplet')['idx_hkl']).astype('Float64')
    ref = ref.loc[~ref['idx_hkl'].isna()]

    ref.idx_hkl = ref.idx_hkl.astype(int)
    return ref

def relabel_hkl(df):
    """
    Relabel the hkl columns in the reflections dataframe to match the
    hkl columns in the crystal dataframe
    """
    df.sort_values(by=['h','k','l'], inplace=True)
    df['triplet'] = df['h'].astype(str) + '_' + df['k'].astype(str) + '_' + df['l'].astype(str)
    hkl = df.groupby('triplet').first()
    hkl['id'] = range(len(hkl))
    df['idx_hkl'] = df['triplet'].map(hkl['id'])
    df.drop(columns='triplet', inplace=True)
    return df

def apply_pg(df, pg,merge_friedel=True):
    """
    Apply the point group symmetry to the reflections dataframe
    """
    if merge_friedel:
        print(df.loc[df.h < 0,['h','k','l']])
        df.loc[df.h < 0,['h','k','l']] *= -1
    if pg == '-1':
        df.loc[df.h < 0,['h','k','l']] *= -1 
        return df
    if pg == '222':
        df.k = np.abs(df.k)
        return df
    if pg == 'mmm':
        df.k = np.abs(df.k)
        df.l = np.abs(df.l)
        return df
    print('Symmetry not implemented')

def extend_reference_to_full_ewaldsphere(ref,pg):
    if pg == '-1':
        refnew = ref.copy()
        refnew['h'] *= -1 
        refnew['k'] *= -1 
        refnew['l'] *= -1 
        ref = pd.concat([ref,refnew])
        return ref

def get_resolution(df,cell):
    """
    Get the resolution of the reflections dataframe
    """
    df['resolution'] = np.sqrt(1 / ((df['h'].astype(int)/cell[0])**2 + (df['k'].astype(int)/cell[1])**2 + (df['l'].astype(int)/cell[2])**2))
    return df

def make_sparse(df):
    """
    Make a sparse matrix from the reflections dataframe
    """
    df = df[['idx_hkl','idx', 'I', 'sigma']]
    data = sparse.COO(coords=(df['idx_hkl'], df['idx']), data=df['I'], shape=(df['idx_hkl'].max()+1, df['idx'].max()+1))
    weights = sparse.COO(coords=(df['idx_hkl'], df['idx']), data=np.ones_like(df['I']), shape=(df['idx_hkl'].max()+1, df['idx'].max()+1))
    return data, weights



def split_df(df):
    """
    Split the reflections dataframe into a list of dataframes
    """
    dfs = [x for _, x in df.groupby('idx')]
    return dfs