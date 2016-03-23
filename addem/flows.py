"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Wed Mar 23, 2016  04:37PM

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
    ns, row, col, fval = sinks.find_single_pixel(arr)
    filled_arr = sinks.fill_single_pixel(arr, ns, row, col, fval)
    sink_list, c_list, fval, outlet = sinks.find_multi_pixel(filled_arr)
    filled_arr_final = sinks.fill_multi_pixel(filled_arr, sink_list, fval)
    return filled_arr_final, ns, row, col, fval

