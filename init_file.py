###3rd August 2018
## Code author SP

###########
### Python libraries
############


import numpy as np
np.set_printoptions(threshold=np.nan)

import os

from numpy.lib.stride_tricks import as_strided

#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.close('all')


############
#### define variables
############

#make folders for outputs
os.mkdir(path_python+'/superpixel_MIBI_python/count_data')
os.mkdir(path_python+'/figures_python')

#clusterChannels = ['CD20','Pan-Keratin','CD45','HLA-DR','CD3','CD68']
clusterChannels = ['CD20','Pan-Keratin','CD45','HLA-DR','CD3','CD68','FoxP3','CD4','CD8','CD16','MPO','CD56','CD209','CD11c','CD11b','Keratin6','Keratin17','p53','Beta catenin','EGFR','Vimentin','SMA','CD31']

neighbourhood_conn = 2 # 4 is also an option
scale = 10 # to see the markers overlaid with the mask
threshold_pixels = 8 # to be considered on the nucleus, was 8

#se = np.one(10) # for finding neighbours via 8-connectivity
divisor_patch = 6

points = 1 # number of patients

window_neigh = 4
dp = 2 #patch for finding neighbours

max_number_MCMC = 3

# The nuclear markers have worked well in this cohort
# Plugin what are the markers that have worked well
certain_markers = ['dsDNA', 'H3K9ac', 'H3K27me3']

step_count = 20 # for number of grid squares width wise = image_width/step_count

# also place the name of the segmentation mask here

def intersection(lst1, lst2):
    # Use of hybrid method
    #temp = set(lst2)
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

# [list(filter(lambda x: x in sure_coords, sublist)) for sublist in regions_revised[0]['coords']]


def intersection1(array1, array2):
    intersection = np.empty([array1.shape[0], array2.shape[0]], dtype=np.int8)
    array2_sets = map(set, array2)
    for i,row1 in enumerate(array1):
        set1 = set(row1)
        for j, set2 in enumerate(array2_sets):
            intersection[i,j] = len(set1.intersection(set2))
    return intersection

def intersect2d(lst1,lst2):
    aset = set([tuple(x) for x in lst1])
    bset = set([tuple(x) for x in lst2])
    l = np.array([x for x in aset & bset])
    return l

def setdiff2d(lst1,lst2):
    aset = set([tuple(x) for x in lst1])
    bset = set([tuple(x) for x in lst2])
    l = np.array([x for x in aset - bset])
    return l



def sliding_window(arr, window_size):
    """ Construct a sliding window view of the array"""
    arr = np.asarray(arr)
    window_size = int(window_size)
    if arr.ndim != 2:
        raise ValueError("need 2-D input")
    if not (window_size > 0):
        raise ValueError("need a positive window size")
    shape = (arr.shape[0] - window_size + 1,
             arr.shape[1] - window_size + 1,
             window_size, window_size)
    if shape[0] <= 0:
        shape = (1, shape[1], arr.shape[0], shape[3])
    if shape[1] <= 0:
        shape = (shape[0], 1, shape[2], arr.shape[1])
    strides = (arr.shape[1]*arr.itemsize, arr.itemsize,
               arr.shape[1]*arr.itemsize, arr.itemsize)
    return as_strided(arr, shape=shape, strides=strides)

def cell_neighbors(arr, i, j, d):
    """Return d-th neighbors of cell (i, j)"""
    w = sliding_window(arr, 2*d+1)

    ix = np.clip(i - d, 0, w.shape[0]-1)
    jx = np.clip(j - d, 0, w.shape[1]-1)

    i0 = max(0, i - d - ix)
    j0 = max(0, j - d - jx)
    i1 = w.shape[2] - max(0, d - i + ix)
    j1 = w.shape[3] - max(0, d - j + jx)

    return w[ix, jx][i0:i1,j0:j1].ravel()


execfile(os.path.join(path_python+"/MIBI_readin_2.py"))
execfile(os.path.join(path_python+"/call_R.py"))
execfile(os.path.join(path_python+"/MIBI_model_postprocess.py"))
