"""
TEST_HELPERS.PY

Created: Fri Mar 11, 2016  02:10PM
Last modified: Tue Mar 22, 2016  11:05AM

"""

import numpy as np

def landscape_with_sinks():
    """
    Creates a square landscape with sinks and depressions.

    This example is taken from the Fill sinks algorithm given at:
    http://spatial-analyst.net/ILWIS/htm/ilwisapp/fill_sinks_algorithm.htm
    """
    height = np.array([
                       [2, 2, 4, 5, 6, 6, 8, 6, 5],
                       [3, 3, 5, 5, 5, 5, 7, 6, 4],
                       [5, 5, 4, 5, 4, 4, 7, 8, 6],
                       [4, 4, 4, 3, 3, 3, 5, 7, 8],
                       [2, 3, 3, 3, 1, 3, 4, 5, 6],
                       [4, 4, 3, 2, 2, 3, 4, 7, 7],
                       [3, 4, 4, 2, 2, 3, 6, 8, 9],
                       [2, 3, 4, 4, 4, 4, 5, 9, 7],
                       [2, 2, 3, 4, 4, 5, 6, 8, 7],
                       ],
                       dtype="float64"
            )
    return height

def __landscape_with_sinks(size=1000, num_sinks=500):
    """
    Creates a flat square landscape with sinks.
    """
    block = np.array([
                     [10, 12, 13, 14, 15, 17],
                     [23, 25, 26, 27, 28, 29],
                     [31, 34, 36, 37, 38, 39],
                     [41, 43, 45, 46, 47, 49],
                     [51, 52, 53, 55, 56, 58],
                     [60, 61, 63, 64, 66, 67],
                     ]
                    )
    rows, cols = [1, 2, 4], [3, 1, 4]
    block[rows, cols] = -5
    block = np.c_[block, block[:, ::-1]]
    landscape = np.r_[block, block[::-1, :]]
    rows, cols = np.where(landscape == -5.)
    return landscape, (rows, cols)

def _randint_unique(lo=0, hi=10, size=5):
    """
    Array of random integers without replacement between given bounds.
    """
    arr = range(lo, hi)
    res = np.random.permutation(arr)[:size]
    return res

def indices_1d_to_2d(indices1d, shape):
    """
    Change array indices from 1D to 2D based on given shape of array.
    """
    r, c = shape[0], shape[1]
    assert max(indices1d) < (r * c), "Index is greater than array size!"
    indices1d = np.array(indices1d).astype("int")
    i = indices1d / int(c)
    j = indices1d % int(c)
    return i, j

def minimum_of_nbrs(arr, r, c):
    """
    Minimum of all the neighbours of the entry located at [r,c] in arr.
    """
    return arr[neighbours(r, c)].min()

def neighbours(r, c):
    """
    Returns the indices of eight neighbours for array entry at [r, c].
    """
    rows = [r - 1, 
            r + 1, 
            r, 
            r, 
            r - 1, 
            r - 1,
            r + 1, 
            r + 1]
    cols = [c, 
            c, 
            c + 1, 
            c - 1, 
            c + 1, 
            c - 1, 
            c + 1, 
            c - 1]
    return rows, cols
