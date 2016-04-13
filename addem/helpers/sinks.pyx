"""
SINKS.PYX

Created: Wed Mar 16, 2016  02:41PM
Last modified: Wed Apr 13, 2016  04:07PM

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
            sa, np_tr_new, va, nv = plateau_float(arr, va, 
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
    print "done."
    print "estimate fill value routine..."
    fv = fill_value_multi_pixel(arr, si, tr, np_tr, ntr, nx, ny,
                                missing_value)
    ################
    ################ Temporary plotting for visualizing the sinks
#     import matplotlib.pyplot as pl
#     arr[arr==missing_value] = np.nan
#     pl.imshow(arr, cmap="jet")
#     pl.colorbar()
#     pl.imshow(si, cmap="gray_r", alpha=0.4)
#     pl.show()
#     import sys
#     sys.exit()
    ################
    ################
    # TODO: 
    #   1. Use fill_multi_pixel to estimate fill values and fill sink
    #       1.a. Go through sink indices
    #       1.b. Estimate catchment for sink
    #       1.c. Get the border of catchment and find the outlet
    #       1.d. estimate the fill value for all sinks in catchment
    #       1.c. Get all the sinks that fall in that catchment from the
    #       results of find_multi_pixel which are now passed to
    #       fill_multi_pixel
    #       1.d. Fill all sinks in the catchment with estimated fill value
    #       1.e. Mark all sinks in that catchment as filled
    ################
    ################
    return None

def plateau_float(np.ndarray[FLOAT_t, ndim=2] arr not None, 
                  np.ndarray[INT_t, ndim=2] va not None,
                  int i, int j, int nx, int ny, int nv, int ndp,
                  ):
    """
    Returns the plateau region in given array starting at position i, j.
    """
    cdef np.ndarray[INT_t, ndim=2] nbrs
    cdef float x
    x = arr[i, j]
    cdef Py_ssize_t ii, jj, k, m
    cdef max_sa = ndp - nv
    cdef np.ndarray[INT_t, ndim=2] s_id = np.zeros([2, max_sa], dtype=INT)
    cdef num = 1
    cdef Py_ssize_t idx = 0
    s_id[0, idx] = i
    s_id[1, idx] = j
    va[i, j] = 1
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
    return s_id, num, va, nv

def fill_value_multi_pixel(np.ndarray[FLOAT_t, ndim=2] arr not None,
                           np.ndarray[INT_t, ndim=2] si not None,
                           np.ndarray[INT_t, ndim=2] si_id not None,
                           np.ndarray[INT_t, ndim=1] si_np not None,
                           int ns, int nx, int ny,
                           float missing_value,
                           ):
    """
    Returns fill value for pre-identified multi-pixel depressions.
    """
    ################
    ################
    #   1. Go through sink indices
    #   2. Estimate catchment for sink
    #   3. Get the border of catchment and find the outlet
    #   4. Estimate the fill value for all sinks in catchment
    #   5. Get all the sinks that fall in that catchment from the
    #      results of find_multi_pixel which are now passed to
    #      fill_multi_pixel
    #   6. Fill all sinks in the catchment with estimated fill value
    #   7. Mark all sinks in that catchment as filled
    ################
    ################
    print "this"
    print "%d sinks to fill!"%ns
    cdef Py_ssize_t k
    cdef np.ndarray[INT_t, ndim=2] si_curr_id
    cdef Py_ssize_t start, stop
    cdef Py_ssize_t i, j
    cdef int nv = ns
    cdef np.ndarray[INT_t, ndim=2] ca
    cdef np.ndarray[INT_t, ndim=2] ca_idx
    cdef int nc
    start = 0
    ####
    #### Temp. plotting commands
    cdef int mr, Mr, mc, Mc
    cdef np.ndarray[INT_t, ndim=2] tmp_ca
    cdef np.ndarray[INT_t, ndim=2] tmp_si
    cdef np.ndarray[FLOAT_t, ndim=2] tmp_ar
    import matplotlib.pyplot as pl
    ####
    ####
    for k in range(ns):
        stop = start + si_np[k]
        si_curr_id = si_id[:, start:stop]
        i = si_curr_id[0, 0]
        j = si_curr_id[1, 0]
        print "catchment initial condition = ", i, j
        ca, ca_idx, nc = catchment(arr, i, j, nx, ny, nv, missing_value)
        print nc
        nv =+ nc
        print "(i, j) = (%d, %d)"%(i, j)
        print "(start, stop, si_np) = (%d, %d, %d)"%(start, stop, si_np[k])
        mr = min(ca_idx[0])
        Mr = max(ca_idx[0])
        mc = min(ca_idx[1])
        Mc = max(ca_idx[1])
        tmp_ca = ca[:10, 4130:4150]
        tmp_si = si[:10, 4130:4150]
        tmp_ar = arr[:10, 4130:4150]
#         tmp_ca = ca[mr:Mr, mc:Mc]
#         tmp_ar = arr[mr:Mr, mc:Mc]
        tmp_ar[tmp_ar == -9999] = np.nan
        pl.imshow(tmp_ar,
                  cmap="jet",
                  interpolation="none",
                  #extent=(mc, Mc, mr, Mr)
                  extent=(4130, 4150, 10, 0)
                  )
        pl.colorbar()
        pl.imshow(tmp_ca,
                  cmap="gray_r",
                  interpolation="none",
                  alpha=0.25,
#                   extent=(mc, Mc, mr, Mr)
                  extent=(4130, 4150, 10, 0)
                  )
        pl.plot(si_curr_id[1] + 0.5, si_curr_id[0] + 0.5, "ks")
#         pl.imshow(tmp_si,
#                   cmap="gray_r",
#                   interpolation="none",
#                   alpha=0.25,
# #                   extent=(mc, Mc, mr, Mr)
#                   extent=(4130, 4150, 10, 0)
#                   )
        pl.title("Catchment for sink at (%d, %d)"%(i, j))
        pl.show()
        if k == 1:
            import sys
            sys.exit()
        start = stop
    import sys
    sys.exit()
    return None

def catchment(np.ndarray[FLOAT_t, ndim=2] arr not None,
              int i, int j, int nx, int ny, int nv,
              float missing_value,
              ):
    """
    Returns the catchment area for the given sink in the given array.
    """
    cdef Py_ssize_t k
    cdef np.ndarray[INT_t, ndim=2] nbrs = np.zeros([2, 8], dtype=INT)
    cdef Py_ssize_t ii, jj, kk
    cdef int max_ca = (nx * ny) - nv
    cdef np.ndarray[INT_t, ndim=2] ca_idx = np.zeros([2, max_ca], dtype=INT)
    cdef float x
    cdef np.ndarray[INT_t, ndim=2] ca = np.zeros([nx, ny], dtype=INT)
    cdef ca_np = 1
    cdef Py_ssize_t idx = 0
    ca_idx[0, idx] = i
    ca_idx[1, idx] = j
    trapped = False
    while not trapped:
        i = ca_idx[0, idx]
        j = ca_idx[1, idx]
        x = arr[i, j]
        # get neighbour list
        nbrs, k = neighbour_indices(i, j, nx, ny)
        # go through the nbr list and add to catchment if necessary
        for kk in range(k):
            ii = nbrs[0, kk]
            jj = nbrs[1, kk]
            if arr[ii, jj] is not missing_value:
                if ca[ii, jj] == 0:# and va[pos] == 0:   # pixel not yet in catchment
                    if arr[ii, jj] >= x:
                        ca[ii, jj] = 1
                        ca_idx[0, ca_np] = ii
                        ca_idx[1, ca_np] = jj
                        ca_np += 1
        idx += 1
        if idx == ca_np:
            trapped = True
    ca_idx = ca_idx[:, :ca_np]
    return ca, ca_idx, ca_np

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

