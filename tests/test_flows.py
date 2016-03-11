"""
TEST_FLOWS.PY

Created: Thu Mar 10, 2016  01:57PM
Last modified: Fri Mar 11, 2016  03:56PM

"""
    
import numpy as np
from addem import flows
import test_helpers as thp

def sinks():
    """Test sink identification in the addem.flows module"""
    arr, sink_idx = thp.flatland_with_sinks( size=1000, num_sinks=5000 )
    num_sinks = flows.sinks(arr, pbar=True)
    print num_sinks
    print len(sink_idx[0])
    return None

