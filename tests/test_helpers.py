"""
TEST_HELPERS.PY

Created: Fri Mar 11, 2016  02:10PM
Last modified: Mon Mar 14, 2016  04:57PM

"""

import numpy as np

def flatland_with_sinks(size=1000, num_sinks=5000):
    """
    Creates a flat square landscape with sinks.
    """
    sink_offset = 0.5
    flatland = np.ones((size, size))
    fl_size = flatland.size
    fl_inner_size = fl_size - 2 * (size + (size - 2))
    sink_indices = _randint_unique(lo=0, hi=fl_inner_size, size=num_sinks)
    r, c = indices_1d_to_2d(sink_indices, shape=(size - 2, size - 2))
    for k in range(num_sinks):
        i, j = r[k], c[k]
        m = minimum_of_nbrs(flatland, i, j)
        flatland[i, j] = sink_offset * m
    return flatland, (r, c)

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
    m = min(
                arr[ r + 1, c ],
                arr[ r - 1, c ],
                arr[ r, c + 1 ],
                arr[ r, c - 1 ],
                arr[ r + 1 , c + 1 ],
                arr[ r + 1 , c - 1 ],
                arr[ r - 1 , c - 1 ],
                arr[ r - 1 , c + 1 ],
            )
    return m
