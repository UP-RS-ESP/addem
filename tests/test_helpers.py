"""
TEST_HELPERS.PY

Created: Fri Mar 11, 2016  02:10PM
Last modified: Fri Mar 11, 2016  03:53PM

"""

import numpy as np
#import matplotlib.pyplot as pl

def flatland_with_sinks(size=1000, num_sinks=5000):
    """
    Creates a flat landscape with sinks.
    """
    flatland = np.ones( ( size, size ) )
    sink_indices = np.random.randint( low=0, high=size, size=( num_sinks, 2 ))
    row, col = [], []
    sink_offset = - 0.5
    for i in range(num_sinks):
        r, c = sink_indices[ i, 0 ], sink_indices[ i, 1]
        if ( r not in row ) and ( c not in col ):
            try:
                m = minimum_of_nbrs( flatland, r, c )
                row.append(r)
                col.append(c)
                flatland[ r, c ] = m - sink_offset
            except IndexError:
                pass
    sink_indices = (row, col)
    return flatland, sink_indices

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
