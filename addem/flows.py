"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Tue Mar 29, 2016  04:22PM

"""

import numpy as np
from addem.helpers import sinks

def hello():
    """says hello"""
    print "Hello! I am the 'addem.flows' module!"
    return None

def sink_filling(arr, verbose=False):
    """
    Finds and fills sinks in given DEM array.
    """
    if type(arr) is np.ma.core.MaskedArray:
        fv = arr.fill_value
        arr = arr.data
    else:
        fv = -9999.
    arr = arr.astype("float64")
    if verbose: print("find single pixel sinks...")
    ns, row, col, fval = sinks.find_single_pixel(arr, fv)
    print "%d single-pixel sinks detected."%ns
    #import sys
    #import matplotlib.pyplot as pl
    #arr[arr==fv] = np.nan
    #pl.imshow(arr)
    #pl.plot(col, row, "ks", ms=1, mec="none", alpha=0.5)
    #pl.show()
    #sys.exit()
    if verbose: print("fill single pixel sinks...")
    filled_arr = sinks.fill_single_pixel(arr, ns, row, col, fval)
    if verbose: print("find multi pixel sinks...")
    s_list, c_list, fval, outlet = sinks.find_multi_pixel(filled_arr, fv)
    if verbose: print("fill multi pixel sinks...")
    filled_arr_final = sinks.fill_multi_pixel(filled_arr, s_list, fval)
    return filled_arr_final, ns, row, col, fval

