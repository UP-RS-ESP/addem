"""
SINKS.PYX

Created: Wed Mar 16, 2016  02:41PM
Last modified: Fri Apr 15, 2016  06:43PM

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
    print "Hello! I am the 'addem.helpers.sinks' module!"
    return None

@cython.boundscheck(False) # turn off bounds-checking for this function
def fill(np.ndarray[FLOAT_t, ndim=2] arr not None,
         float missing_value, int verbose=1):
    """
    Finds and fills depressions in given DEM array.
    """
    print "running sinks.fill..."
    cdef np.ndarray[FLOAT_t, ndim=2] arr_out = arr
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef Py_ssize_t i, j, k
    cdef float x, m
    ####
    ####
    # get data pixels
    cdef int max_ndp = nx * ny
    cdef np.ndarray[INT_t, ndim=2] dp_ind = np.ndarray([2, max_ndp], dtype=INT)
    cdef int dp_num
    dp_ind, dp_num = get_data_pixels(arr_out[1:nx-1, 1:ny-1],
                                     nx-2, ny-2, missing_value)
    dp_ind = dp_ind + 1
    ####
    ####
    # fill single pixel sinks...
    print "fill single-pixel sinks.."
    cdef np.ndarray[INT_t, ndim=1] row = dp_ind[0]
    cdef np.ndarray[INT_t, ndim=1] col = dp_ind[1]
    cdef float nn, ss, ee, ww, ne, nw, se, sw
    for k in range(dp_num):
        i = <Py_ssize_t>row[k]
        j = <Py_ssize_t>col[k]
        x = arr_out[i, j]
        m = min(
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
    print "done."
    ####
    ####
    # fill multi-pixel sinks iteratively
    cdef np.ndarray[INT_t, ndim=2] visited = np.ndarray([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ind
    cdef bint is_sink
    print "fill multi-pixel sinks..."
    for k in range(dp_num):
        i = <Py_ssize_t>row[k]
        j = <Py_ssize_t>col[k]
        x = arr_out[i, j]
        m = min(
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
                    is_sink = floodplain(arr_out,
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
    print "done."
    return arr_out

def get_data_pixels(np.ndarray[FLOAT_t, ndim=2] arr not None, 
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
cdef bint floodplain(np.ndarray[FLOAT_t, ndim=2] arr, 
                   Py_ssize_t i, Py_ssize_t j, 
                   int nx, int ny, int dp_num, 
                   float missing_value):
    """
    Returns flatland at (i, j); min. of bordering pixels; and if way down.
    """
    cdef float y
    cdef float x
    cdef int num = 1
    cdef int nbr_i_min = 0
    cdef int nbr_i_max = 0
    cdef int nbr_j_min = 0
    cdef int nbr_j_max = 0
    cdef int max_nbrs = 8
    cdef bint way_down = 0
    cdef bint at_border = 0
    cdef bint trapped = 0
    cdef Py_ssize_t ii, jj, kk
    cdef Py_ssize_t k
    cdef Py_ssize_t idx = 0
    cdef np.ndarray[INT_t, ndim=2] va = np.ndarray([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ind = np.ndarray([2, dp_num], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] nbr = np.ndarray([2, max_nbrs], dtype=INT)
    ####
    ind[0, idx] = i
    ind[1, idx] = j
    x = arr[i, j]
    va[i, j] = 1
    while not trapped:
        i = <Py_ssize_t>ind[0, idx]
        j = <Py_ssize_t>ind[1, idx]
        ###
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
        ###
        ### Now, go through the range of i, j values and create a nbr list.
        k = 0
        for ii in range(nbr_i_min, nbr_i_max):
            for jj in range(nbr_j_min, nbr_j_max):
                if not (ii == i and jj == j):
                    nbr[0, k] = ii
                    nbr[1, k] = jj
                    k += 1
        # go through the nbr list and add to flatland if necessary
        for kk in range(k):
            ii = <Py_ssize_t>nbr[0, kk]
            jj = <Py_ssize_t>nbr[1, kk]
            y = arr[ii, jj]
            if va[ii, jj] == 0:   # i.e. pixel is not yet visited
                if y == x:
                    ind[0, num] = ii
                    ind[1, num] = jj
                    num += 1
                    va[ii, jj] = 1
                elif y < x:
                    way_down = 1
                if y == missing_value:
                    at_border = 1
        idx += 1
        if idx == num:
            trapped = 1
    # decide if floodplain is sink or not
    cdef bint is_sink = 0
    if not at_border:
        if not way_down:
            is_sink = 1
    return is_sink

@cython.boundscheck(False) # turn off bounds-checking for this function
def floodfill(np.ndarray[FLOAT_t, ndim=2] arr not None, 
              np.ndarray[INT_t, ndim=2] visited not None,
              Py_ssize_t i, Py_ssize_t j, int nx, int ny, int dp_num):
    """
    blabla
    """
    cdef float y
    cdef float x
    cdef float m
    cdef int num = 1
    cdef int max_nbrs = 8
    cdef int nbr_i_min = 0
    cdef int nbr_i_max = 0
    cdef int nbr_j_min = 0
    cdef int nbr_j_max = 0
    cdef Py_ssize_t ii, jj, kk
    cdef Py_ssize_t k
    cdef Py_ssize_t idx = 0
    cdef np.ndarray[INT_t, ndim=2] va = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] ind = np.ndarray([2, dp_num], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] nbrs = np.ndarray([2, max_nbrs], dtype=INT)
    cdef bint way_down = 0
    cdef bint trapped = 0
    ####
    ####
    ind[0, idx] = i
    ind[1, idx] = j
    va[i, j] = 1
    while not way_down:
        trapped = 0
        m = 100000
        x = arr[i, j]
        idx = <Py_ssize_t>0
        while not trapped:
            i = <Py_ssize_t>ind[0, idx]
            j = <Py_ssize_t>ind[1, idx]
            ###
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
            ###
            ### Now, go through the range of i, j values and create a nbr list.
            k = 0
            for ii in range(nbr_i_min, nbr_i_max):
                for jj in range(nbr_j_min, nbr_j_max):
                    if not (ii == i and jj == j):
                        nbrs[0, k] = ii
                        nbrs[1, k] = jj
                        k += 1
            # go through the nbr list and add to flatland if necessary
            for kk in range(k):
                ii = <Py_ssize_t>nbrs[0, kk]
                jj = <Py_ssize_t>nbrs[1, kk]
                y = arr[ii, jj]
                if va[ii, jj] == 0:   # i.e. pixel is not yet visited
                    if y == x:
                        ind[0, num] = ii
                        ind[1, num] = jj
                        num += 1
                        va[ii, jj] = 1
                    elif y < x:
                        way_down = 1
                    elif y > x:
                        if y < m:
                            m = y
            idx += 1
            if idx == num:
                trapped = 1
        if not way_down:
            for k in range(num):
                ii = <Py_ssize_t>ind[0, k]
                jj = <Py_ssize_t>ind[1, k]
                arr[ii, jj] = <float>m
    for k in range(num):
        ii = <Py_ssize_t>ind[0, k]
        jj = <Py_ssize_t>ind[1, k]
        visited[ii, jj] = 1
    return arr, visited

@cython.boundscheck(False) # turn off bounds-checking for this function
def find(np.ndarray[FLOAT_t, ndim=2] arr not None, float missing_value):
    """
    Finds multi-pixel depressions in given DEM array.
    """
    ####################################################
    ####################################################
    #### A. FIRST SWEEP:
    ####    1. Check if pixel has missing data
    ####    2. Check if pixel is at border of landscape
    ####    3. Check if there's a way down from pixel
    ####################################################
    ####################################################
    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]
    cdef int max_sa = (nx * ny)
    cdef int ndp = 0                        # no. of data pixels
    cdef int ns = 0                         # no. of sink pixels
    cdef Py_ssize_t i, j
    cdef Py_ssize_t k = 0
    cdef float x, m
    cdef float nn, ss, ee, ww, ne, nw, se, sw
    cdef int way_down = 0
    cdef np.ndarray[FLOAT_t, ndim=1] nbrs = np.zeros([8], dtype=FLOAT)
    cdef np.ndarray[FLOAT_t, ndim=2] mi = np.zeros([nx, ny], dtype=FLOAT)
    cdef np.ndarray[INT_t, ndim=2] si = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] si_id = np.zeros([2, max_sa], dtype=INT)
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            ####
            # 1. Check if pixel is missing data
            ####
            x = arr[i, j]
            mi[i, j] = missing_value
            if x is not missing_value:
                ndp += 1
                ####
                # 2. Check if pixel is at border of landscape
                ####
                # 2.a. Define neighbours of (i, j)
                nn = arr[i - 1, j]          # north-north
                ss = arr[i + 1, j]          # south-south
                ee = arr[i, j - 1]          # east-east
                ww = arr[i, j + 1]          # west-west
                ne = arr[i - 1, j - 1]      # north-east
                nw = arr[i - 1, j + 1]      # north-west
                se = arr[i + 1, j - 1]      # south-east
                sw = arr[i + 1, j + 1]      # south-west
                # 2.b. Find minimum of neighbours
                m = min(
                        nn, ss,
                        ee, ww,
                        ne, nw,
                        se, sw
                        )
                mi[i, j] = m
                # 2.c. m equals missing value => pixel is at border
                if m is not missing_value:
                    ####
                    # 3. Check if there's a way down from the pixel
                    ####
                    # 3.a. Assign elements to nbrs array
                    nbrs[0] = nn
                    nbrs[1] = ss
                    nbrs[2] = ee
                    nbrs[3] = ww
                    nbrs[4] = ne
                    nbrs[5] = nw
                    nbrs[6] = se
                    nbrs[7] = sw
                    # 3.b. Check among neighbours for a way down
                    way_down = 0
                    k = 0
                    while way_down == 0 and k < 8:
                        if nbrs[k] < x:
                            way_down = 1
                        else:
                            k += 1
                    # 3.c. No way down => pixel could be part of sink
                    if way_down == 0:
                        si[i, j] = 1
                        si_id[0, ns] = i
                        si_id[1, ns] = j
                        ns += 1
    print "%d possible sink pixels found..."%ns
    print "\t...out of a total of %d pixels in data array"%ndp
    ####################################################
    ####################################################
    #### B. SECOND SWEEP
    ####    1. Identify flat "trapping regions" in sinks
    ####    2. Check if:
    ####        2.a. "edge sinks" which exit at borders
    ####        2.b. "hill tops" which have a way down
    ####    3. Store the trapping region only if:
    ####        3.a. it's not at the border
    ####        3.b. there's no way down
    ####################################################
    ####################################################
    print "getting trapping regions..."
    cdef int nv = 0
    cdef int ntr = 0
    cdef int np_tr_old = 0
    cdef int np_tr_new = 0
    cdef Py_ssize_t ii, jj, kk
    cdef np.ndarray[INT_t, ndim=2] tr = np.zeros([2, ndp], dtype=INT)
    cdef np.ndarray[INT_t, ndim=1] np_tr = np.zeros([ns], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] sa = np.zeros([2, ns], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] va = np.zeros([nx, ny], dtype=INT)
    si = np.zeros([nx, ny], dtype=INT)
    for k in range(ns/100):
        i = si_id[0, k]
        j = si_id[1, k]
        ####
        # 1. Identify "trapping regions" in sink pixels from above
        ####
        if va[i, j] == 0:
            sa, np_tr_new, va, nv = trapping_region(arr, va, 
                                                    i, j, nx, ny, nv, ndp)
            si[sa[0], sa[1]] = 1
            ####
            # 2.a. Check if trapping region is at border
            ####
            is_at_border = 0
            kk = 0
            while kk < np_tr_new and is_at_border == 0:
                ii = sa[0, kk]
                jj = sa[1, kk]
                m = mi[ii, jj]
                if m is missing_value:
                    is_at_border = 1
                else:
                    kk += 1
            if is_at_border == 0:
                ####
                # 2.b. Check if trapping region has a way down
                ####
                way_down = 0
                kk = 0
                while kk < np_tr_new and way_down == 0:
                    ii = sa[0, kk]
                    jj = sa[1, kk]
                    m = mi[ii, jj]
                    if m < arr[ii, jj]:
                        way_down = 1
                    else:
                        kk += 1
                if way_down == 0:
                    ####
                    # 3. Accept only if not at border and no way down
                    ####
                    si[sa[0], sa[1]] = 1
                    tr[0, np_tr_old:np_tr_old + np_tr_new] = sa[0]
                    tr[1, np_tr_old:np_tr_old + np_tr_new] = sa[1]
                    np_tr[ntr] = np_tr_new
                    np_tr_old += np_tr_new
                    ntr += 1
    tr = tr[:, :ntr]
    np_tr = np_tr[:ntr]
    return tr, np_tr, ntr

@cython.boundscheck(False) # turn off bounds-checking for this function
def trapping_region(np.ndarray[FLOAT_t, ndim=2] arr not None, 
                    np.ndarray[INT_t, ndim=2] va not None,
                    int i, int j, int nx, int ny, int nv, int ndp,
                    ):
    """
    Returns the plateau region in given array starting at position i, j.
    """
    cdef float x
    x = arr[i, j]
    cdef Py_ssize_t ii, jj, k, m
    cdef max_sa = ndp - nv
    cdef np.ndarray[INT_t, ndim=2] s_id = np.ndarray([2, max_sa], dtype=INT)
    cdef num = 1
    cdef Py_ssize_t idx = 0
    #### for neighbour indices
    cdef int nbr_i_min = 0
    cdef int nbr_i_max = 0
    cdef int nbr_j_min = 0
    cdef int nbr_j_max = 0
    cdef int max_nbrs = 8
    cdef np.ndarray[INT_t, ndim=2] nbrs = np.ndarray([2, max_nbrs], dtype=INT)
    ####
    s_id[0, idx] = i
    s_id[1, idx] = j
    va[i, j] = 1
    trapped = False
    while not trapped:
        i = s_id[0, idx]
        j = s_id[1, idx]
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
        ###
        ### Now, go through the range of i, j values and create a nbr list.
        k = 0
        for ii in range(nbr_i_min, nbr_i_max):
            for jj in range(nbr_j_min, nbr_j_max):
                if not (ii == i and jj == j):
                    nbrs[0, k] = ii
                    nbrs[1, k] = jj
                    k += 1
        # go through the nbr list and add to sink if necessary
        for m in range(k):
            ii = nbrs[0, m]
            jj = nbrs[1, m]
            if va[ii, jj] == 0:   # i.e. pixel is not yet in sink
                if arr[ii, jj] == x:
                    s_id[0, num] = ii
                    s_id[1, num] = jj
                    num += 1
                    va[ii, jj] = 1
                    nv += 1
        idx += 1
        if idx == num:
            trapped = True
    s_id = s_id[:, :num]
    return s_id, num, va, nv

