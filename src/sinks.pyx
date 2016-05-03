"""
SINKS.PYX

Created: Wed Mar 16, 2016  02:41PM
Last modified: Fri Apr 29, 2016  04:27PM

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
    print "Hello! I am the 'addem.sinks' module!"
    return None

@cython.boundscheck(False) # turn off bounds-checking for this function
def find(np.ndarray[FLOAT_t, ndim=2] arr not None,
         float missing_value, int verbose=1):
    """
    Finds multi-pixel depressions in given DEM array.
    """
    cdef np.ndarray[FLOAT_t, ndim=2] arr_out = arr
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef Py_ssize_t i, j, k
    cdef FLOAT_t x, m
    # get data pixels
    cdef int max_ndp = nx * ny
    cdef np.ndarray[INT_t, ndim=2] dp_ind = np.zeros([2, max_ndp], dtype=INT)
    cdef int dp_num
    dp_ind, dp_num = data_pixels(arr_out[1:nx-1, 1:ny-1],
                                 nx-2, ny-2, missing_value)
    dp_ind = dp_ind + 1
    cdef np.ndarray[INT_t, ndim=2] visited = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=1] row = dp_ind[0]
    cdef np.ndarray[INT_t, ndim=1] col = dp_ind[1]
    cdef int ns = 0         # no. of single-pixel sinks
    cdef int nm = 0         # no. of multi-pixel sinks
    cdef np.ndarray[INT_t, ndim=2] ind
    cdef bint is_sink
    for k in range(dp_num):
        i = <Py_ssize_t>row[k]
        j = <Py_ssize_t>col[k]
        x = <FLOAT_t>arr_out[i, j]
        m = <FLOAT_t>min(
                         arr_out[i - 1, j],
                         arr_out[i + 1, j],
                         arr_out[i, j - 1],
                         arr_out[i, j + 1],
                         arr_out[i - 1, j - 1],
                         arr_out[i - 1, j + 1],
                         arr_out[i + 1, j - 1],
                         arr_out[i + 1, j + 1],
                        )
        if m != missing_value:
            if visited[i, j] == 0:
                if x < m:
                    ns += 1
                    visited[i, j] = 1
                elif x == m:
                    is_sink, visited = check_if_sink(arr_out, visited,
                                                     i, j,
                                                     nx, ny,
                                                     dp_num,
                                                     missing_value)
                    if is_sink:
                        nm += 1
    return ns, nm

@cython.boundscheck(False) # turn off bounds-checking for this function
def fill(np.ndarray[FLOAT_t, ndim=2] arr not None,
         float missing_value, int verbose=1):
    """
    Finds and fills depressions in given DEM array.
    """
    cdef np.ndarray[FLOAT_t, ndim=2] arr_out = arr
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef Py_ssize_t i, j, k
    cdef FLOAT_t x
    cdef FLOAT_t m
    # get data pixels
    cdef int max_ndp = nx * ny
    cdef np.ndarray[INT_t, ndim=2] dp_ind = np.ndarray([2, max_ndp], dtype=INT)
    cdef int dp_num
    dp_ind, dp_num = data_pixels(arr_out[1:nx-1, 1:ny-1],
                                 nx-2, ny-2, missing_value)
    dp_ind = dp_ind + 1
    # fill single pixel sinks...
    cdef np.ndarray[INT_t, ndim=1] row = dp_ind[0]
    cdef np.ndarray[INT_t, ndim=1] col = dp_ind[1]
    for k in range(dp_num):
        i = <Py_ssize_t>row[k]
        j = <Py_ssize_t>col[k]
        x = <FLOAT_t>arr_out[i, j]
        m = <FLOAT_t>min(
                         arr_out[i - 1, j],
                         arr_out[i + 1, j],
                         arr_out[i, j - 1],
                         arr_out[i, j + 1],
                         arr_out[i - 1, j - 1],
                         arr_out[i - 1, j + 1],
                         arr_out[i + 1, j - 1],
                         arr_out[i + 1, j + 1],
                        )
        if m != missing_value:
            if x < m:
                arr_out[i, j] = m
    # fill multi-pixel sinks iteratively
    cdef np.ndarray[INT_t, ndim=2] visited = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] vi_tmp = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ind
    cdef bint is_sink
    for k in range(dp_num):
        i = <Py_ssize_t>row[k]
        j = <Py_ssize_t>col[k]
        x = <FLOAT_t>arr_out[i, j]
        m = <FLOAT_t>min(
                         arr_out[i - 1, j],
                         arr_out[i + 1, j],
                         arr_out[i, j - 1],
                         arr_out[i, j + 1],
                         arr_out[i - 1, j - 1],
                         arr_out[i - 1, j + 1],
                         arr_out[i + 1, j - 1],
                         arr_out[i + 1, j + 1],
                        )
        if m != missing_value:
            if visited[i, j] == 0:
                if x == m:
                    is_sink, vi_tmp = check_if_sink(arr_out,
                                                    vi_tmp,
                                                    i, j,
                                                    nx, ny,
                                                    dp_num,
                                                    missing_value)
                    if is_sink:
                        arr_out, visited = floodfill(arr_out,
                                                     visited,
                                                     i, j,
                                                     nx, ny,
                                                     dp_num)
    return arr_out

@cython.boundscheck(False) # turn off bounds-checking for this function
def data_pixels(np.ndarray[FLOAT_t, ndim=2] arr not None,
                int nx, int ny, float missing_value):
    """
    Returns indices of data pixels in DEM array.
    """
    cdef Py_ssize_t i, j
    cdef int max_ndp = nx * ny
    cdef np.ndarray[INT_t, ndim=2] dp_ind = np.zeros([2, max_ndp], dtype=INT)
    cdef int dp_num = 0
    for i in range(nx):
        for j in range(ny):
            if arr[i, j] is not missing_value:
                dp_ind[0, dp_num] = i
                dp_ind[1, dp_num] = j
                dp_num += 1
    dp_ind = dp_ind[:, :dp_num]
    return dp_ind, dp_num

@cython.boundscheck(False) # turn off bounds-checking for this function
cdef check_if_sink(np.ndarray[FLOAT_t, ndim=2] arr,
                   np.ndarray[INT_t, ndim=2] va_global,
                   Py_ssize_t i, Py_ssize_t j,
                   int nx, int ny, int dp_num,
                   float missing_value):
    """
    Returns flatland at (i, j); min. of bordering pixels; and if way down.
    """
    cdef FLOAT_t y
    cdef FLOAT_t x
    cdef int nbr_i_min = 0
    cdef int nbr_i_max = 0
    cdef int nbr_j_min = 0
    cdef int nbr_j_max = 0
    cdef int max_nbrs = 8
    cdef Py_ssize_t idx = 0
    cdef bint way_down = 0
    cdef bint at_border = 0
    cdef bint trapped = 0
    cdef Py_ssize_t ii, jj, kk
    cdef Py_ssize_t k
    cdef Py_ssize_t num = 1
    cdef np.ndarray[INT_t, ndim=2] va_local = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ind = np.zeros([2, dp_num], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] nbr = np.zeros([2, max_nbrs], dtype=INT)
    ####
    ind[0, idx] = <INT_t>i
    ind[1, idx] = <INT_t>j
    x = arr[i, j]
    va_local[i, j] = <INT_t>1
    while not trapped and not way_down:
        i = ind[0, idx]
        j = ind[1, idx]
        # get neighbour list
        nbr_i_min = i - 1
        nbr_i_max = i + 2
        if i == 0:
            nbr_i_min = 0
        elif i == nx - 1:
            nbr_i_max = nx
        nbr_j_min = j - 1
        nbr_j_max = j + 2
        if j == 0:
            nbr_j_min = 0
        elif j == ny - 1:
            nbr_j_max = ny
        ### go through the range of i, j values and create a nbr list.
        k = 0
        for ii in range(nbr_i_min, nbr_i_max):
            for jj in range(nbr_j_min, nbr_j_max):
                if not (ii == i and jj == j):
                    nbr[0, k] = <INT_t>ii
                    nbr[1, k] = <INT_t>jj
                    k += 1
        # go through the nbr list and add to flatland if necessary
        for kk in range(k):
            ii = nbr[0, kk]
            jj = nbr[1, kk]
            y = arr[ii, jj]
            if va_local[ii, jj] == 0:   # i.e. pixel is not yet visited
                if y == x:
                    ind[0, num] = <INT_t>ii
                    ind[1, num] = <INT_t>jj
                    num += 1
                    va_local[ii, jj] = <INT_t>1
                elif y < x:
                    way_down = <bint>1
                if y == missing_value:
                    at_border = <bint>1
        idx += 1
        if idx == num:
            trapped = <bint>1
    # decide if floodplain is sink or not
    cdef bint is_sink = 0
    if not at_border:
        if not way_down:
            is_sink = <bint>1
            for k in range(num):
                ii = ind[0, k]
                jj = ind[1, k]
                va_global[ii, jj] = <INT_t>1
    return is_sink, va_global

@cython.boundscheck(False) # turn off bounds-checking for this function
def floodfill(np.ndarray[FLOAT_t, ndim=2] arr, # not None,
              np.ndarray[INT_t, ndim=2] va_global, # not None,
              Py_ssize_t i, Py_ssize_t j, int nx, int ny, int dp_num):
    """
    blabla
    """
    cdef FLOAT_t y
    cdef FLOAT_t x
    cdef FLOAT_t m
    cdef int max_nbrs = 8
    cdef int nbr_i_min = 0
    cdef int nbr_i_max = 0
    cdef int nbr_j_min = 0
    cdef int nbr_j_max = 0
    cdef Py_ssize_t idx = 0
    cdef Py_ssize_t num = 1
    cdef Py_ssize_t ii, jj, kk
    cdef Py_ssize_t k
    cdef np.ndarray[INT_t, ndim=2] va_local = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ind = np.zeros([2, dp_num], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] nbrs = np.zeros([2, max_nbrs], dtype=INT)
    cdef bint way_down = 0
    cdef bint trapped = 0
    ind[0, idx] = <INT_t>i
    ind[1, idx] = <INT_t>j
    va_local[i, j] = <INT_t>1
    while not way_down:
        trapped = <bint>0
        m = <FLOAT_t>100000
        x = arr[i, j]
        idx = 0
        while not trapped:
            i = ind[0, idx]
            j = ind[1, idx]
            # get neighbour list
            nbr_i_min = i - 1
            nbr_i_max = i + 2
            if i == 0:
                nbr_i_min = 0
            elif i == nx - 1:
                nbr_i_max = nx
            nbr_j_min = j - 1
            nbr_j_max = j + 2
            if j == 0:
                nbr_j_min = 0
            elif j == ny - 1:
                nbr_j_max = ny
            ### go through the range of i, j values and create a nbr list.
            k = 0
            for ii in range(nbr_i_min, nbr_i_max):
                for jj in range(nbr_j_min, nbr_j_max):
                    if not (ii == i and jj == j):
                        nbrs[0, k] = <INT_t>ii
                        nbrs[1, k] = <INT_t>jj
                        k += 1
            # go through the nbr list and add to flatland if necessary
            for kk in range(k):
                ii = nbrs[0, kk]
                jj = nbrs[1, kk]
                y = arr[ii, jj]
                if va_local[ii, jj] == 0:   # pixel not yet visited locally
                    if y == x:
                        ind[0, num] = <INT_t>ii
                        ind[1, num] = <INT_t>jj
                        num += 1
                        va_local[ii, jj] = <INT_t>1
                    elif y < x:
                        way_down = <bint>1
                        trapped = <bint>1
                    elif y > x:
                        if y < m:
                            m = y
            idx += 1
            if idx == num:
                trapped = <bint>1
        if not way_down:
            for k in range(num):
                ii = ind[0, k]
                jj = ind[1, k]
                arr[ii, jj] = m
    for k in range(num):
        ii = ind[0, k]
        jj = ind[1, k]
        va_global[ii, jj] = <INT_t>1
    return arr, va_global
