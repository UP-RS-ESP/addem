#! /usr/bin/env python
"""
SINKS.PYX

Created: Wed Mar 16, 2016  02:41PM
Last modified: Wed Mar 16, 2016  06:02PM

"""

import numpy as np
cimport numpy as np


INT = np.int
ctypedef np.int_t INT_t
FLOAT = np.float64
ctypedef np.float64_t FLOAT_t

def hello():
    """says hello"""
    print "Hello! I am the 'addem.helpers.sinks' module!"
    return None

def locate(np.ndarray[FLOAT_t, ndim=2] arr not None):
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef int max_ns = (nx * ny) // 9
    cdef int ns = 0
    cdef Py_ssize_t i, j, k
    cdef float m
    cdef np.ndarray[INT_t, ndim=1] row = np.zeros([max_ns], dtype=INT)
    cdef np.ndarray[INT_t, ndim=1] col = np.zeros([max_ns], dtype=INT)
    k = 0
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            m = min(# minimum of neighbours
                    arr[i-1, j], arr[i+1, j],       # north & south
                    arr[i, j-1], arr[i, j+1],       # east & west
                    arr[i-1, j-1], arr[i-1, j+1],   # nw & ne
                    arr[i+1, j-1], arr[i+1, j+1]    # sw & se
                    )
            if arr[i, j] < m:
                row[k] = i
                col[k] = j
                k += 1
                ns += 1
    row = row[:k]
    col = col[:k]
    return ns, row, col

