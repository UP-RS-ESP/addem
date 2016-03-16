"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Wed Mar 16, 2016  05:34PM

"""

import numpy as np
from progressbar import ProgressBar
from addem.helpers import sinks

def hello():
    """says hello"""
    print "Hello! I am the 'addem.flows' module!"
    return None

def locate_sinks(arr, verbose=False, pbar=False):
    """
    Identify sinks in given DEM array.
    """
    ns, row, col = sinks.locate(arr.astype("float64"))
    row = row[:ns]
    col = col[:ns]
    return ns, row, col

def locate_sinks_purepython(arr, verbose=False, pbar=False):
    """
    Identify sinks in given DEM array.
    """
    nx, ny = arr.shape
    if verbose: print("testing for loop...")
    if pbar: prog = ProgressBar(maxval=nx-2).start()
    span = range(1, nx - 1)
    # only 1 out of 9 locations can be a sink
    # => max. possible no. of sinks = (nx * ny) / 9
    max_ns = (nx * ny) / 9
    row = np.zeros(max_ns, dtype="int")
    col = np.zeros(max_ns, dtype="int")
    # core loop for detecting sinks
    ns = 0                          # no. of sinks detected
    # indices for west & east positions
    w = range(2, nx)                # west
    e = range(0, nx - 2)            # east
    for i in span:
        # indices for center, northm & south positions
        c = i                       # centre
        n = i - 1                   # north
        s = i + 1                   # south
        # arrays for different positions / axes
        cc = arr[c, span]
        nn = arr[n, span]
        ss = arr[s, span]
        ww = arr[c, w]
        ee = arr[c, e]
        nw = arr[n, w]
        ne = arr[n, e]
        sw = arr[s, w]
        se = arr[s, e]
        # construct array of all 8 neighbours
        # the array 'arr_nbrs' has shape (8, nx - 2)
        arr_nbrs = np.array([nn, ss, ww, ee, nw, ne, sw, se])
        min_nbrs = arr_nbrs.min(axis=0)
        # sinks := where value is less than minimum of nbrs
        sinks = cc < min_nbrs
        num = sinks.sum()
        if num > 0:
            sink_col_idx = np.where(sinks)[0] + 1
            row[ns:ns+num] = i
            col[ns:ns+num] = sink_col_idx
            ns += num
        if pbar: prog.update(i)
    if pbar: prog.finish()    
    row = row[:ns]
    col = col[:ns]
    return ns, row, col
