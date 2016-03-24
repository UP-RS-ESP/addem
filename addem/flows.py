"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Thu Mar 24, 2016  05:13PM

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
        arr = arr.data
    arr = arr.astype("float64")
    if verbose: print("find single pixel sinks...")
    ns, row, col, fval = sinks.find_single_pixel(arr)
    if verbose: print("fill single pixel sinks...")
    filled_arr = sinks.fill_single_pixel(arr, ns, row, col, fval)
    if verbose: print("find multi pixel sinks...")
    sink_list, c_list, fval, outlet = sinks.find_multi_pixel(filled_arr)
    import sys
    #sys.exit()
    test = np.zeros(arr.shape, dtype="float64")
    for loc in c_list[0]:
        test[loc] = arr[loc]#np.random.randint(1000)
    test[test==0] = np.nan
    import matplotlib.pyplot as pl
    pl.imshow(test, cmap="jet", interpolation="none")
    pl.colorbar()
    pl.show()
    sys.exit()
    if verbose: print("fill multi pixel sinks...")
    filled_arr_final = sinks.fill_multi_pixel(filled_arr, sink_list, fval)
    if verbose: print("done.")
    return filled_arr_final, ns, row, col, fval

