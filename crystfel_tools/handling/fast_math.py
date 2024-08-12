import numba
import numpy as np



def calc_cc(set1,set2):
    """
    Calculate the correlation coefficient between two sets of data
    """
    mean1 = set1.mean()
    mean2 = set2.mean()
    normed1 = set1 - mean1
    normed2 = set2 - mean2
    return np.sum(normed1 * normed2) / np.sqrt(np.sum(normed1**2) * np.sum(normed2**2))

def calc_mean(data,weights):
    sum_weigths = np.sum(weights,axis=1).todense()
    sum_weigths[sum_weigths == 0] = 1
    return np.sum(data,axis=1).todense() / sum_weigths

def write_out_half_datasets(data, weigths):
    half1 = data[:,:data.shape[1]//2]

    half2 = data[:,data.shape[1]//2:]

    weigths1 = weigths[:,:weigths.shape[1]//2]
    weigths2 = weigths[:,weigths.shape[1]//2:]
    half1 = calc_mean(half1,weigths1)
    half2 = calc_mean(half2,weigths2)
    return half1,half2