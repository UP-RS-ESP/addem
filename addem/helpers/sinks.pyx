"""
SINKS.PYX

Created: Wed Mar 16, 2016  02:41PM
Last modified: Fri Apr 01, 2016  06:38PM

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
                      float missing_value):
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
    k = 0
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            if arr[i, j] is not missing_value:
                m = min(# minimum of neighbours
                        arr[i-1, j], arr[i+1, j],       # north & south
                        arr[i, j-1], arr[i, j+1],       # east & west
                        arr[i-1, j-1], arr[i-1, j+1],   # nw & ne
                        arr[i+1, j-1], arr[i+1, j+1]    # sw & se
                        )
                if m is not missing_value: # i.e. (i, j) not at border
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
                     float missing_value):
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
    cdef int n_data_pixels = 0
    cdef int ns = 0
    cdef Py_ssize_t i, j
    cdef Py_ssize_t k = 0
    cdef float x, m
    cdef float nn, ss, ee, ww, ne, nw, se, sw
    cdef int way_down = 0
    cdef np.ndarray[FLOAT_t, ndim=1] nbrs = np.zeros([8], dtype=FLOAT)
    cdef np.ndarray[INT_t, ndim=2] si = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] si_id = np.zeros([2, max_sa], dtype=INT)
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            ####
            # 1. Check if pixel is missing data
            ####
            x = arr[i, j]
            if x is not missing_value:
                n_data_pixels += 1
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
    print "\t...out of a total of %d pixels in data array"%n_data_pixels
    ####################################################
    ####################################################
    #### B. SECOND SWEEP
    ####    1. Walk through potential sink pixels from A
    ####    2. Identify flat "trapping" regions 
    ####    3. Discard:
    ####        3.a. "hill tops" which have a way down
    ####        3.b. "edge sinks" which exit at borders
    ####################################################
    ####################################################
    cdef nv = 0
    cdef np.ndarray[INT_t, ndim=2] sink = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] va = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] sa = np.zeros([nx, ny], dtype=INT)
    cdef np.ndarray[INT_t, ndim=2] s = np.zeros([2, ns], dtype=INT)
    for k in range(ns):
        i = si_id[0, k]
        j = si_id[1, k]
        if va[i, j] == 0:
            s, va, nv = plateau_float(arr, va, 
                                      i, j, nx, ny, nv, n_data_pixels)
            sink[s[0], s[1]] = 1
    # TODO: 
    #   1. Store the above identified sinks efficiently in a ndarray
    #       a. This array can have a max row length
    #       b. Another array can store the no. of elements in each row,
    #          beyond which all remaining tail-end elements are zero.
    #   2. Go through the above identified 'sinks' and discard the 
    #       spurious sinks which are either hill tops, slopes, or edges
    return None

def plateau_float(np.ndarray[FLOAT_t, ndim=2] arr not None, 
                  np.ndarray[INT_t, ndim=2] va not None,
                  int i, int j, int nx, int ny, int nv, int n_data_pixels,
                  ):
    """
    Returns the plateau region in given array starting at position i, j.
    """
    cdef np.ndarray[INT_t, ndim=2] nbrs
    cdef float x
    x = arr[i, j]
    cdef Py_ssize_t ii, jj, k, m
    cdef max_sa = n_data_pixels - nv
    cdef np.ndarray[INT_t, ndim=2] s_id = np.zeros([2, max_sa], dtype=INT)
    cdef num = 1
    cdef Py_ssize_t idx = 0
    s_id[0, idx] = i
    s_id[1, idx] = j
    trapped = False
    while not trapped:
        i = s_id[0, idx]
        j = s_id[1, idx]
        # get neighbour list
        nbrs, k = neighbour_indices(i, j, nx, ny)
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
    return s_id, va, nv

def sinks_in_catchment(np.ndarray[INT_t, ndim=2] ca not None,
                       np.ndarray[INT_t, ndim=2] ca_idx not None,
                       np.ndarray[INT_t, ndim=2] va not None
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
#     import matplotlib.pyplot as pl
#     pl.imshow(c, cmap="gray_r", interpolation="none")
#     pl.show()
#     import sys
#     sys.exit()
    cdef int ii, jj
    ii = 4 - min_row
    jj = 4133 - min_col
    cdef Py_ssize_t i, j, k
    cdef list sinks = []
    cdef ns = 0
    cdef nv = 0
    for i in range(nr):
        for j in range(nc):
            if c[i, j] == 0 and s[i, j] == 0:
                s, s_id, nv = plateau_int(c, s, i, j, 
                                          nr, nc, nv, c[i, j])
                sinks.append(s_id)
                ns += 1
#     cdef np.ndarray[INT_t, ndim=2] S = np.zeros([nr, nc], dtype=INT)
#     for i in range(ns):
#         s_id = sinks[i]
#         S[s_id[0], s_id[1]] = 2
#     import matplotlib.pyplot as pl
#     pl.imshow(c+S, cmap="jet", interpolation="none")
#     pl.colorbar()
#     pl.show()
#     import sys
#     sys.exit()
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
            ca[s_id[0] + min_row, s_id[1] + min_col] = 1
            va[s_id[0] + min_row, s_id[1] + min_col] = 1
    va[ca_idx[0], ca_idx[1]] = 1
    si = si[:k]
#     import matplotlib.pyplot as pl
#     pl.imshow(c, cmap="gray_r", interpolation="none")
#     pl.show()
#     import sys
#     sys.exit()
    return si, ca, va, nv

def plateau_int(np.ndarray[INT_t, ndim=2] arr not None, 
                np.ndarray[INT_t, ndim=2] sa not None,
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
    if nr == 2 and nc == 1:
        print "(i, j) ", i, j
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
            if nr == 2 and nc == 1:
                print "i, j", i, j
                print "pos", pos
                print "num, nv ", num, nv
#                 import sys
#                 sys.exit()
            if sa[pos] == 0:   # i.e. pixel is not yet in sink
                if arr[pos] == m:
                    sa[pos] = 1
                    try:
                        s_id[0, num] = pos[0]
                    except IndexError:
                        print "(nr, nc)", nr, nc
                        print "nbrs", nbrs
                        print "pos", pos
                        print "sa", sa
                        raise IndexError, "something's wrong!"
                    s_id[0, num] = pos[0]
                    s_id[1, num] = pos[1]
                    num += 1
                    nv += 1
        idx += 1
        if idx == num:
            trapped = True
    s_id = s_id[:, :num]
    return sa, s_id, nv

def catchment(np.ndarray[FLOAT_t, ndim=2] arr not None,
              np.ndarray[INT_t, ndim=2] va not None,
              int i, int j, int nx, int ny, int nv,
              float m, float fill_val
              ):
    """
    Returns the catchment area for the given sink in the given array.
    """
    ########################################################################
    ########################################################################
    ############
    ## Set nbrs of (i, j) as initial condition (IC) for catchment
    ############
    cdef Py_ssize_t k, u, v
    cdef list nbrs, i_range, j_range
    cdef tuple pos
    cdef int max_ca = (nx * ny) - nv
    cdef np.ndarray[INT_t, ndim=2] ca_idx = np.zeros([2, max_ca], dtype=INT)
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
    #### 2. Go through nbrs and note down those nbrs > minimum
    ############
    k = 0
    for pos in nbrs:
        if arr[pos] > m:
            ca_idx[0, k] = pos[0]
            ca_idx[1, k] = pos[1]
            k += 1
    ############
    ########################################################################
    ########################################################################
    ############
    ## Using above (i, j) as IC, walk through DEM and identify catchment
    ############
    cdef float x
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
            if ca[pos] == 0 and va[pos] == 0:   # pixel not yet in catchment
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
    Returns indices & no. of neighbours depending on location of (i, j).
    """
    ### From the point of view of neighbourhood indices, 
    ### for each of nx or ny, we find 3 possibilities:
    ###     A. n >= 3
    ###     B. n = 2
    ###     C. n = 1
    ###
    ### For nx pr ny falling in each of the 3 cases above, we find
    ### 9 different possibilities:
    ###      Case       nx      ny
    ###         1       A       A
    ###         2       A       B
    ###         3       B       A
    ###         4       B       B
    ###         5       B       C
    ###         6       C       B
    ###         7       C       C
    ###         8       A       C
    ###         9       C       A
    ###
    ### Below, we outline the different possible neighbour indices
    ### for each of the above 9 cases.
    ###
    ###
    cdef list i_range
    cdef list j_range
    ###
    ###
    if nx >= 3 and ny >= 3:   # a min. of 3x3 array
        ###
        ### Case 1. nx >= 3 and ny >=3
        ###
        ###### set i range
        if i == 0:
            i_range = [0, 1]
        elif i == nx - 1:
            i_range = [nx - 2, nx - 1]
        else:
            i_range = range(i - 1, i + 2)
        ###### set j range
        if j == 0:
            j_range = [0, 1]
        elif j == ny - 1:
            j_range = [ny - 2, ny - 1]
        else:
            j_range = range(j - 1, j + 2)
        ###
        ###
    elif nx >= 3 and ny == 2:
        ###
        ### Case 2. nx >= 3 and ny = 2
        ###
        ###### set i range
        if i == 0:
            i_range = [0, 1]
        elif i == nx - 1:
            i_range = [nx - 2, nx - 1]
        else:
            i_range = range(i - 1, i + 2)
        ###### set j range
        j_range = [0, 1]
        ###
        ###
    elif nx == 2 and ny >= 3:
        ###
        ### Case 3. nx = 2 and ny >= 3
        ###
        ###### set i range
        i_range = [0, 1]
        ###### set j range
        if j == 0:
            j_range = [0, 1]
        elif j == ny - 1:
            j_range = [ny - 2, ny - 1]
        else:
            j_range = range(j - 1, j + 2)
        ###
        ###
    elif nx == 2 and ny == 2:
        ###
        ### Case 4. nx = 2 and ny = 2
        ###
        ###### set i range
        i_range = [0, 1]
        ###### set j range
        j_range = [0, 1]
        ###
        ###
    elif nx == 2 and ny == 1:
        ###
        ### Case 5. nx = 2 and ny = 1
        ###
        ###### set i range
        i_range = [0, 1]
        ###### set j range
        j_range = [0]
        ###
        ###
    elif nx == 1 and ny == 2:
        ###
        ### Case 6. nx = 1 and ny = 2
        ###
        ###### set i range
        i_range = [0]
        ###### set j range
        j_range = [0, 1]
        pass
        ###
        ###
    elif nx == 1 and ny == 1:
        ###
        ### Case 7. nx = 1 and ny = 1
        ###
        ###### set i range
        i_range = [0]
        ###### set j range
        j_range = [0]
        ###
        ###
    elif nx >= 3 and ny == 1:
        ###
        ### Case 8. nx >= 3 and ny = 1
        ###
        ###### set i range
        if i == 0:
            i_range = [0, 1]
        elif i == nx - 1:
            i_range = [nx - 2, nx - 1]
        else:
            i_range = range(i - 1, i + 2)
        ###### set j range
        j_range = [0]
        ###
        ###
    elif nx == 1 and ny >= 3:
        ###
        ### Case 9. nx = 2 and ny >= 3
        ###
        ###### set i range
        i_range = [0]
        ###### set j range
        if j == 0:
            j_range = [0, 1]
        elif j == ny - 1:
            j_range = [ny - 2, ny - 1]
        else:
            j_range = range(j - 1, j + 2)
        ###
        ###
    ###
    ###
    ### Now, go through the range of i, j values and create a nbr list.
    cdef int max_nbrs = 8
    cdef np.ndarray[INT_t, ndim=2] nbrs = np.zeros([2, max_nbrs], dtype=INT)
    cdef Py_ssize_t ii, jj
    cdef Py_ssize_t k = 0
    for ii in i_range:
        for jj in j_range:
            if not (ii == i and jj == j):
                nbrs[0, k] = ii
                nbrs[1, k] = jj
                k += 1
    nbrs = nbrs[:k]
    return nbrs, k



#     import matplotlib.pyplot as pl
#     arr[arr==missing_value] = np.nan
#     pl.imshow(arr, cmap="jet")
#     pl.colorbar()
#     pl.imshow(si, cmap="gray_r", alpha=0.4)
#     pl.show()
#     import sys
#     print "starting 'trapping' region estimation..."
#     sys.stdout.flush()
#     cdef float pc = 0
#     for k in range(ns):
#         pc = float(k * 100) / float(ns)
#         print "\b\b\b\b\b\b\b\b%.1f %%"%pc,
#         sys.stdout.flush()
#     print "\ndone!"
