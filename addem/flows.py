"""
FLOWS.PY

Created: Thu Mar 10, 2016  01:55PM
Last modified: Fri Apr 15, 2016  04:27PM

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
    si, np_si, nsi = sinks.find(arr, fv)
    print "%d sinks detected."%nsi
    arr = sinks.fill(arr, fv)
    si, np_si, nsi = sinks.find(arr, fv)
    print "%d sinks remained after sink filling routine."%nsi
    return arr

