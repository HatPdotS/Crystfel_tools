import pandas as pd
import numpy as np

def relabel_hkl(df):
    """
    Relabel the hkl columns in the reflections dataframe to match the
    hkl columns in the crystal dataframe
    """
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
    print(df)
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