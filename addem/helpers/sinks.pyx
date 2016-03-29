"""
SINKS.PYX

Created: Wed Mar 16, 2016  02:41PM
Last modified: Tue Mar 29, 2016  06:13PM

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

def find_single_pixel(np.ndarray[FLOAT_t, ndim=2] arr not None,
                      float fill_val):
    """
    Finds single-pixel depressions in given DEM array.
    """
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef int max_ns = (nx * ny) // 9
    cdef int ns = 0
    cdef Py_ssize_t i, j, k
    cdef float m
    cdef np.ndarray[INT_t, ndim=1] row = np.zeros([max_ns], dtype=INT)
    cdef np.ndarray[INT_t, ndim=1] col = np.zeros([max_ns], dtype=INT)
    cdef np.ndarray[FLOAT_t, ndim=1] fval = np.zeros([max_ns], dtype=FLOAT)
    cdef undefined = False
    k = 0
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            if arr[i, j] is not fill_val:
                m = min(# minimum of neighbours
                        arr[i-1, j], arr[i+1, j],       # north & south
                        arr[i, j-1], arr[i, j+1],       # east & west
                        arr[i-1, j-1], arr[i-1, j+1],   # nw & ne
                        arr[i+1, j-1], arr[i+1, j+1]    # sw & se
                        )
                if m is not fill_val: # i.e. (i, j) not at border
                    # sinks_1px := array entries less than min. of nbrs
                    if arr[i, j] < m:
                        row[k] = i
                        col[k] = j
                        fval[k] = m
                        k += 1
                        ns += 1
    row = row[:k]
    col = col[:k]
    fval = fval[:k]
    return ns, row, col, fval

def fill_single_pixel(np.ndarray[FLOAT_t, ndim=2] arr not None,
                      int ns,
                      np.ndarray[INT_t, ndim=1] row not None,
                      np.ndarray[INT_t, ndim=1] col not None,
                      np.ndarray[FLOAT_t, ndim=1] fval not None,
                      ):
    """
    Fills single-pixel depressions in given DEM array.
    """
    cdef Py_ssize_t k
    cdef float m
    cdef np.ndarray[FLOAT_t, ndim=1] nbrs = np.zeros([8], dtype=FLOAT)
    for k in range(ns):
        arr[row[k], col[k]] = fval[k]
    return arr

def find_multi_pixel(np.ndarray[FLOAT_t, ndim=2] arr not None,
                     float fill_val):
    """
    Finds multi-pixel depressions in given DEM array.
    """
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef int max_ca = arr.size
    cdef int max_ns = max_ca // 9
    cdef int ns = 0
    cdef Py_ssize_t i, j, k
    cdef float m
    cdef list sinks = [0] * max_ns
    cdef list catchments = [0] * max_ns
    cdef np.ndarray[FLOAT_t, ndim=1] fval = np.zeros([max_ns], dtype=FLOAT)
    cdef np.ndarray[INT_t, ndim=1] ca_row = np.zeros([max_ca], dtype=INT)
    cdef np.ndarray[INT_t, ndim=1] ca_rol = np.zeros([max_ca], dtype=INT)
    cdef list visited = []
    cdef list sink = []
    cdef float f
    k = 0
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            if arr[i, j] is not fill_val:
                if (i, j) not in visited:
                    m = min(# minimum of neighbours
                            arr[i-1, j], arr[i+1, j],       # north & south
                            arr[i, j-1], arr[i, j+1],       # east & west
                            arr[i-1, j-1], arr[i-1, j+1],   # nw & ne
                            arr[i+1, j-1], arr[i+1, j+1]    # sw & se
                            )
                    if m is not fill_val:
                        if arr[i, j] == m:
                            sink = [(i, j)]
                            visited = [(i, j)]
                            catch = []
                            print "plateau at (%d, %d) ..."%(i, j)
                            s, c, v = plateau_at(arr, i, j, m,
                                                 sink, catch, visited)
                            print "...identified with %d pixels"%len(s)
                            print "catchment at (%d, %d)..."%(i, j)
                            catchment, visited = sink_catchment(arr, m, c,
                                                                visited,
                                                                fill_val)
                            f, outlet = fill_value_multi_pixel(arr,
                                                               catchment,
                                                               sink)
                            sinks[k] = sink
                            catchments[k] = catchment
                            fval[k] = f
                            k += 1
    sinks = sinks[:k]
    catchments = catchments[:k]
    fval = fval[:k]
    return sinks, catchments, fval, outlet

def plateau_at(np.ndarray[FLOAT_t, ndim=2] arr not None, 
               int i, int j, float m,
               list sink, list catchment, list visited,
               ):
    """
    Returns the plateau region in given array starting at position i, j.
    """
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef list i_range, j_range
    cdef int u, v
    cdef list locs
    cdef Py_ssize_t idx_locs, k
    trapped = False
    while not trapped:
        k = 0
        locs = [0] * 9
        i_range, j_range = neighbour_indices(i, j, nx, ny)
        for u in i_range:
            for v in j_range:
                locs[k] = (u, v)
                k += 1
        locs = locs[:k]
        idx_locs = 0
        found_new_loc = False
        while not found_new_loc:
            if idx_locs < k:
                loc = locs[idx_locs]
                if loc not in visited:
                    if arr[loc] == m:
                        sink.append(loc)
                        i = loc[0]
                        j = loc[1]
                        found_new_loc = True
                    else:
                        catchment.append(loc)
                    visited.append(loc)
                idx_locs += 1
            else:
                found_new_loc = True
                trapped = True
    return sink, catchment, visited

def sink_catchment(np.ndarray[FLOAT_t, ndim=2] arr not None, float m,
                   list catchment, list visited, float fill_val,
                   ):
    """
    Returns the catchment area for the given sink in the given array.
    """
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef list i_range, j_range
    cdef int ii, jj
    cdef Py_ssize_t idx_locs, k, u, v, i, j
    trapped = False
    cdef int nc = 0
    idx_locs = 0
    cdef tuple loc
    cdef list nbrs
    cdef tuple nbr

#    # first identify those points which are greater than m
#    cdef Py_ssize_t i, j
#    cdef list checklist = []
#    for i in range(1, nx - 1):
#        for j in range(1, ny - 1):
#            if arr[i, j] is not fill_val:
#                x = min(# minimum of neighbours
#                        arr[i-1, j], arr[i+1, j],       # north & south
#                        arr[i, j-1], arr[i, j+1],       # east & west
#                        arr[i-1, j-1], arr[i-1, j+1],   # nw & ne
#                        arr[i+1, j-1], arr[i+1, j+1]    # sw & se
#                        )
#                if x is not fill_val:
#                    if (i, j) not in visited:
#                        if arr[i, j] > m:
#                            checklist.append((i, j))
#    import sys
#    print len(checklist)
#    print nx * ny
#    sys.exit()
    cdef float x
    cdef int stop = 0
    while not trapped:
        loc = catchment[idx_locs]
        x = arr[loc]
        # get neighbour list
        # writing the entire function here to save some overhead
        i = loc[0]
        j = loc[1]
        if i == 0:
            i_range = range(i, i + 2)
        elif i == nx - 1:
            i_range = range(nx - 2, nx)
        else:
            i_range = range(i - 1, i + 2)
        if j == 0:
            j_range = range(j, j + 2)
        elif j == ny - 1:
            j_range = range(ny - 2, ny)
        else:
            j_range = range(j - 1, j + 2)
        k = 0
        nbrs = [0] * 9
        for u in i_range:
            for v in j_range:
                nbrs[k] = (u, v)
                k += 1
        nbrs = nbrs[:k]
        # go through the nbr list and add to catchment if necessary
        for nbr in nbrs:
            if nbr not in catchment:
                if arr[nbr] >= x:
                    catchment.append(nbr)
                    stop += 1
        idx_locs += 1
        nc = len(catchment)
        if idx_locs == nc:
            trapped = True
        if stop % 1000 == 0:
            print stop
        if stop in [5000, 50000, 51000, 52000, 53000, 54000, 55000]:
            print "plotting..."
            import matplotlib.pyplot as pl
            arr[arr==fill_val] = np.nan
            pl.imshow(arr)
            for i in range(nc):
                pl.plot(catchment[i][1], catchment[i][0],
                        "ks", mec="none", ms=0.75, alpha=0.5)
            pl.show()
            import sys
            sys.exit()
    visited.extend(catchment)
    return catchment, visited

def fill_value_multi_pixel(np.ndarray[FLOAT_t, ndim=2] arr not None,
                           list catchment, list sink):
    """
    Returns fill value for pre-identified multi-pixel depressions.
    """
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef float m = 999999.
    cdef float fval = 999999.
    cdef list minlocs = []
    cdef tuple outlet
    cdef list border_out = []
    cdef list border_in = []
    for loc in catchment:
        i = loc[0]
        j = loc[1]
        if i == 0:
            i_range = range(i, i + 2)
        elif i == nx - 1:
            i_range = range(nx - 2, nx)
        else:
            i_range = range(i - 1, i + 2)
        if j == 0:
            j_range = range(j, j + 2)
        elif j == ny - 1:
            j_range = range(ny - 2, ny)
        else:
            j_range = range(j - 1, j + 2)
        for ii in i_range:
            for jj in j_range:
                if (ii, jj) not in sink:
                    if (ii, jj) not in catchment:
                        if (ii, jj) not in border_out:
                            border_out.append((ii, jj))
                        if (i, j) not in border_in:
                            border_in.append((i, j))
    for loc in border_out:
        if arr[loc] <= m:
            m = arr[loc]
            minlocs.append(loc)
    for loc in minlocs:
        i = loc[0]
        j = loc[1]
        if i == 0:
            i_range = range(i, i + 2)
        elif i == nx - 1:
            i_range = range(nx - 2, nx)
        else:
            i_range = range(i - 1, i + 2)
        if j == 0:
            j_range = range(j, j + 2)
        elif j == ny - 1:
            j_range = range(ny - 2, ny)
        else:
            j_range = range(j - 1, j + 2)
        for ii in i_range:
            for jj in j_range:
                if (ii, jj) in border_in:
                    if arr[ii, jj] <= fval:
                        fval = arr[ii, jj]
                        outlet = (ii, jj)
    return fval, outlet

def fill_multi_pixel(np.ndarray[FLOAT_t, ndim=2] arr not None,
                     list sinks,
                     np.ndarray[FLOAT_t, ndim=1] fval
                     ):
    """
    Fills multi-pixel depressions in given DEM array.
    """
    cdef list sink
    cdef tuple loc
    cdef int k
    cdef float f
    k = 0
    for sink in sinks:
        f = fval[k]
        for loc in sink:
            arr[loc] = f
        k += 1
    return arr

def neighbour_indices(int i, int j, int nx, int ny):
    """
    Returns indices of neighbours depending on the location of (i, j).
    """
    cdef list i_range
    cdef list j_range
    if i == 0:
        i_range = range(i, i + 2)
    elif i == nx - 1:
        i_range = range(nx - 2, nx)
    else:
        i_range = range(i - 1, i + 2)
    if j == 0:
        j_range = range(j, j + 2)
    elif j == ny - 1:
        j_range = range(ny - 2, ny)
    else:
        j_range = range(j - 1, j + 2)
    return i_range, j_range



