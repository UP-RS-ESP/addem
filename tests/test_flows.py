"""
TEST_FLOWS.PY

Created: Thu Mar 10, 2016  01:57PM
Last modified: Tue Mar 15, 2016  02:06PM

"""
    
import numpy as np
from addem import flows
import test_helpers as thp

def sinks():
    """Test sink identification in the addem.flows module"""
    arr, sink_idx = thp.landscape_with_sinks(size=1000, num_sinks=500)
    sinks_detected = flows.sinks(arr, pbar=True)
    sinks_given = len(sink_idx[0])
    assert sinks_detected == sinks_given, "Sink detection test failed!"
    return None

