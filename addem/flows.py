"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Mon Mar 14, 2016  05:15PM

"""

import numpy as np
from progressbar import ProgressBar

def hello():
    """says hello"""
    print "Hello! I am the 'addem.flows' module!"
    return None

def sinks(arr, verbose=False, pbar=False):
    """
    Identify sinks in given DEM array.
    """
    import sys
    nx, ny = arr.shape
    if verbose: print("testing for loop...")
    if pbar: prog = ProgressBar(maxval=nx).start()
    count = 0
    for i in range(1, nx - 1):
        sinks_curr = np.ones(len(nx - 2)).astype("bool")
#         dif_nn = arr[i, :] - arr[i - 1, :]
#         dif_ss = arr[i + 1, :] - arr[i, :]
#         dif_ww = np.r_[arr[i, 0], arr[i, :-1] - arr[i, 1:]]
#         dif_ee = np.r_[arr[i, 1:] - arr[i, :-1], arr[i, -1]]
#         dif_nw = np.r_[arr[i, 0], arr[i - 1, :-1] - arr[i, 1:]]
#         dif_ne = np.r_[arr[i, 1:] - arr[i, :-1], arr[i, -1]]
#         dif_sw = np.r_[arr[i, 0], arr[i + 1, :-1] - arr[i, 1:]]
#         dif_se = np.r_[arr[i + 1, 1:] - arr[i, :-1], arr[i, -1]]
#         tmp = np.c_[
#                     dif_nn,
#                     dif_ss,
#                     dif_ww,
#                     dif_ee,
#                     dif_nw,
#                     dif_ne,
#                     dif_sw,
#                     dif_se,
#                     ]
#         tmp[tmp>0] = 0
#         tmp[tmp<0] = 1
#         tmp = tmp.astype("bool")
#         idx = tmp.all(axis=1)
#         count += sum(idx)
        if pbar: prog.update(i)
    if pbar: prog.finish()    
    return count
