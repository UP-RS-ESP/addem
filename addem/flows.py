"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Tue Mar 15, 2016  01:35PM

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
    nx, ny = arr.shape
    if verbose: print("testing for loop...")
    if pbar: prog = ProgressBar(maxval=nx-2).start()
    count = 0
    span = range(1, nx - 1)
    for i in span:
        # indices for different positions / axes
        c = span                    # centre
        n = i - 1                   # north
        s = i + 1                   # south
        w = range(2, nx)            # west
        e = range(0, nx - 2)        # east
        # arrays for different positions / axes
        cc = arr[i, span]
        nn = arr[n, span]
        ss = arr[s, span]
        ww = arr[c, w]
        ee = arr[c, e]
        nw = arr[n, w]
        ne = arr[n, e]
        sw = arr[s, w]
        se = arr[s, e]
        # differences between arrays at different locations
        # the array 'diff' has shape (8, nx - 2)
        locs = [nn, ss, ww, ee, nw, ne, sw, se]
        diff = np.array([loc - cc for loc in locs])
        # sinks := those columns in diff with only positive entries
        count += (diff > 0).all(axis=0).sum()
        if pbar: prog.update(i)
    if pbar: prog.finish()    
    return count
