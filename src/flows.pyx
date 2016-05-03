"""
FLOWS.PYX

Created: Tue May 03, 2016  03:25PM
Last modified: Tue May 03, 2016  03:39PM

"""

import numpy as np
cimport numpy as np
cimport cython

INT = np.int
ctypedef np.int_t INT_t
FLOAT = np.float64
ctypedef np.float64_t FLOAT_t


def hello():
    """says hello"""
    print "Hello! I am the 'addem.flows' module!"
    return None


@cython.boundscheck(False) # turn off bounds-checking for this function
def direction(np.ndarray[FLOAT_t, ndim=2] arr not None,
              float missing_value, int verbose=1):
    """
    Estimates the flow direction for each pixel in given DEM raster.
    """
    if verbose:
        print("Coming soon...")
    return None
