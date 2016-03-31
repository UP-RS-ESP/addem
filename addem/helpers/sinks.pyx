"""
SINKS.PYX

Created: Wed Mar 16, 2016  02:41PM
Last modified: Thu Mar 31, 2016  12:08PM

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
    cdef Py_ssize_t i, j, k
    cdef float m, f
    cdef np.ndarray[FLOAT_t, ndim=1] fval = np.zeros([max_ns], dtype=FLOAT)
    cdef np.ndarray[INT_t, ndim=2] ca = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] va = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ca_idx
    cdef list sinks = [0] * max_ns
    cdef list catchments = [0] * max_ns
    cdef list visited = []
    cdef list si = []
    k = 0
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            if arr[i, j] is not fill_val:
                if va[i, j] == 0:
                    m = min(# minimum of neighbours
                            arr[i-1, j], arr[i+1, j],       # north & south
                            arr[i, j-1], arr[i, j+1],       # east & west
                            arr[i-1, j-1], arr[i-1, j+1],   # nw & ne
                            arr[i+1, j-1], arr[i+1, j+1]    # sw & se
                            )
                    if m is not fill_val:
                        if arr[i, j] == m:
                            si = []
                            ca_idx = np.zeros([2, max_ca], dtype=INT)
                            print "plateau at (%d, %d) ..."%(i, j)
#                             si, ca, va, ca_idx, nc = plateau(arr, ca, va,
#                                                              ca_idx,
#                                                              i, j, nx, ny, m,
#                                                              si)
                            import sys
                            print "catchment at (%d, %d)..."%(i, j)
                            ca, ca_idx = catchment(arr, 
                                                           #ca, va, 
                                                           #ca_idx,
                                                           i, j,
                                                           nx, ny,# nc, 
                                                           m, fill_val)
                            #print "estimated catchment with %d pixels"%nc
                            print "estimating sinks within catchment..."
                            si, va = sinks_in_catchment(ca, ca_idx)
                            f, outlet = fill_value_multi_pixel(arr,
                                                               catchment,
                                                               si)
                            sinks[k] = si
                            catchments[k] = catchment
                            fval[k] = f
                            k += 1
    sinks = sinks[:k]
    catchments = catchments[:k]
    fval = fval[:k]
    return sinks, catchments, fval, outlet

def sinks_in_catchment(np.ndarray[INT_t, ndim=2] ca not None,
                       np.ndarray[INT_t, ndim=2] ca_idx not None,
                       ):
    """
    Returns multi-pixel sinks within a given catchment.
    """
    cdef np.ndarray[INT_t, ndim=1] row = ca_idx[0, :]
    cdef np.ndarray[INT_t, ndim=1] col = ca_idx[1, :]
    cdef int min_row = min(row)
    cdef int max_row = max(row)
    cdef int min_col = min(col)
    cdef int max_col = max(col)
    cdef int nr = max_row - min_row
    cdef int nc = max_col - min_col
    cdef int nx = ca.shape[0]
    cdef int ny = ca.shape[1]
    cdef np.ndarray[INT_t, ndim=2] s = np.zeros([nr, nc], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] v = np.zeros([nr, nc], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] c = np.zeros([nr, nc], dtype=INT)
    c = ca[min_row:max_row, min_col:max_col]
    cdef Py_ssize_t i, j, k
    cdef list sinks = []
    cdef ns = 0
    cdef nv = 0
    for i in range(nr):
        for j in range(nc):
            if c[i, j] == 0 and v[i, j] == 0:
                s, v, s_id, nv = plateau_int(c, s, v, i, j, 
                                             nr, nc, nv, ca[i, j])
                sinks.append(s_id)
                ns += 1
    print "filtering sinks which extend to borders..."
    cdef list si = [0] * ns
    k = 0
    for i in range(ns):
        s_id = sinks[i]
        rlo = min(s_id[0])
        rhi = max(s_id[0])
        clo = min(s_id[0])
        chi = max(s_id[0])
        if not (rlo == 0 or clo == 0 or rhi == nr - 1 or chi == nc - 1):
            si[k] = s_id
            k += 1
            c[s_id[0], s_id[1]] = 1
    print "done. storing results in original ca array and plotting..."
    ca[min_row:max_row, min_col:max_col] = c
    import sys
    import matplotlib.pyplot as pl
    pl.imshow(c, cmap="gray_r", interpolation="none")
    pl.show()
    sys.exit()
    si = si[:k]
    return si, ca

def plateau_int(np.ndarray[INT_t, ndim=2] arr not None, 
                np.ndarray[INT_t, ndim=2] sa not None,
                np.ndarray[INT_t, ndim=2] va not None,
                int i, int j, int nr, int nc, int nv, int m,
                ):
    """
    Returns the plateau region in given array starting at position i, j.
    """
    cdef list i_range, j_range
    cdef list nbrs
    cdef tuple pos
    cdef Py_ssize_t u, v, k
    cdef max_sa = (nr * nc) - nv
    cdef np.ndarray[INT_t, ndim=2] s_id = np.zeros([2, max_sa], dtype=INT)
    cdef num = 1
    trapped = False
    cdef Py_ssize_t idx = 0
    s_id[0, idx] = i
    s_id[1, idx] = j
    while not trapped:
        i = s_id[0, idx]
        j = s_id[1, idx]
        x = arr[i, j]
        # get neighbour list
        i_range, j_range = neighbour_indices(i, j, nr, nc)
        k = 0
        nbrs = [0] * 9
        for u in i_range:
            for v in j_range:
                nbrs[k] = (u, v)
                k += 1
        nbrs = nbrs[:k]
        # go through the nbr list and add to sink if necessary
        for pos in nbrs:
            if va[pos] == 0:   # i.e. pixel is not yet in sink
                if arr[pos] == m:
                    sa[pos] = 2
                    s_id[0, num] = pos[0]
                    s_id[1, num] = pos[1]
                    num += 1
                    va[pos] = 1
                    nv += 1
        idx += 1
        if idx == num:
            trapped = True
    s_id = s_id[:, :num]
    return sa, va, s_id, nv

def plateau(np.ndarray[FLOAT_t, ndim=2] arr not None, 
            np.ndarray[INT_t, ndim=2] ca not None,
            np.ndarray[INT_t, ndim=2] va not None,
            np.ndarray[INT_t, ndim=2] ca_idx not None,
            int i, int j, int nx, int ny, float m, 
            list sink,
            ):
    """
    Returns the plateau region in given array starting at position i, j.
    """
    cdef list i_range, j_range
    cdef list nbrs
    cdef tuple pos
    cdef Py_ssize_t idx, u, v, k
    cdef int nc = 0
    trapped = False
    while not trapped:
        k = 0
        nbrs = [0] * 9
        i_range, j_range = neighbour_indices(i, j, nx, ny)
        for u in i_range:
            for v in j_range:
                nbrs[k] = (u, v)
                k += 1
        nbrs = nbrs[:k]
        idx = 0
        shift_pos = False
        while not shift_pos:
            if idx < k:
                pos = nbrs[idx]
                if va[pos] == 0:    # i.e., this pixel is not visited yet
                    if arr[pos] == m:
                        sink.append(pos)
                        i = pos[0]
                        j = pos[1]
                        shift_pos = True
                    else:
                        ca[pos] = 1
                        ca_idx[0, nc] = pos[0]
                        ca_idx[1, nc] = pos[1]
                        nc += 1
                    va[pos] = 1
                idx += 1
            else:
                shift_pos = True
                trapped = True
    return sink, ca, va, ca_idx, nc

def catchment(np.ndarray[FLOAT_t, ndim=2] arr not None,
              #np.ndarray[INT_t, ndim=2] ca not None,
              #np.ndarray[INT_t, ndim=2] va not None,
              #np.ndarray[INT_t, ndim=2] ca_idx not None,
              int i, int j,
              int nx, int ny, 
              #int nc, 
              float m, float fill_val
              ):
    """
    Returns the catchment area for the given sink in the given array.
    """
    ########################################################################
    ########################################################################
    ############
    ## Find initial pixel nbr to (i, j) which is part of catchment
    ############
    cdef Py_ssize_t k, u, v
    cdef list nbrs, i_range, j_range
    cdef tuple pos
    cdef int num_nbrs
    ############
    #### 1. Find nbrs
    ############
    i_range, j_range = neighbour_indices(i, j, nx, ny)
    k = 0
    nbrs = [0] * 9
    for u in i_range:
        for v in j_range:
            nbrs[k] = (u, v)
            k += 1
    nbrs = nbrs[:k]
    ############
    #### 2. Go through nbrs and note down 1st nbr > minimum
    ############
    num_nbrs = len(nbrs)
    notFound = True
    k = 0
    while notFound and k < num_nbrs:
        pos = nbrs[k]
        if arr[pos] > m:
            i = pos[0]
            j = pos[1]
            notFound = False
        else:
            k += 1
    ############
    ########################################################################
    ########################################################################
    ############
    ## Using above (i, j) as IC, walk through DEM and identify catchment
    ############
    cdef float x
    cdef int max_ca = nx * ny
    cdef np.ndarray[INT_t, ndim=2] ca_idx = np.zeros([2, max_ca], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ca = np.zeros([nx, ny], dtype=INT)
    cdef nc = 1
    cdef Py_ssize_t idx = 0
    ca_idx[0, idx] = i
    ca_idx[1, idx] = j
    trapped = False
    while not trapped:
        i = ca_idx[0, idx]
        j = ca_idx[1, idx]
        x = arr[i, j]
        # get neighbour list
        i_range, j_range = neighbour_indices(i, j, nx, ny)
        k = 0
        nbrs = [0] * 9
        for u in i_range:
            for v in j_range:
                nbrs[k] = (u, v)
                k += 1
        nbrs = nbrs[:k]
        # go through the nbr list and add to catchment if necessary
        for pos in nbrs:
            if ca[pos] == 0:   # i.e. pixel is not yet in catchment
                if arr[pos] >= x:
                    ca[pos] = 1
                    ca_idx[0, nc] = pos[0]
                    ca_idx[1, nc] = pos[1]
                    nc += 1
#                     va[pos] = 1
        idx += 1
        if idx == nc:
            trapped = True
    ca_idx = ca_idx[:, :nc]
    return ca, ca_idx#, nc

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



