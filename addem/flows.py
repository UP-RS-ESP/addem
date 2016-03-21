"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Mon Mar 21, 2016  03:00PM

"""

import numpy as np
from addem.helpers import sinks

def hello():
    """says hello"""
    print "Hello! I am the 'addem.flows' module!"
    return None

def locate_sinks(arr, verbose=False):
    """
    Identify sinks in given DEM array.
    """
    if type(arr) is np.ma.core.MaskedArray:
        arr = arr.data
    arr = arr.astype("float64")
    ns, row, col = sinks.locate(arr)
    return ns, row, col

def fill_sinks(arr, ns=None, row=None, col=None, verbose=False):
    """
    Fill sinks in given DEM array.
    """
    if type(arr) is np.ma.core.MaskedArray:
        arr = arr.data
    arr = arr.astype("float64")
    if ns is None:
        if row is None or col is None:
            raise "Row / column indices are needed!"
        ns, row, col = sinks.locate(arr)
    sink_filled_arr = sinks.fill(arr, ns, row, col)
    return sink_filled_arr

