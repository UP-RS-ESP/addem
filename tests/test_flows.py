"""
TEST_FLOWS.PY

Created: Thu Mar 10, 2016  01:57PM
Last modified: Tue Mar 29, 2016  01:54PM

"""
    
from addem import flows
import test_helpers as thp

def sink_filling():
    """Test sink identification in the addem.flows module"""
    arr = thp.landscape_with_sinks()
    print "original array"
    print arr
    filled_arr, ns, row, col, fval = flows.sink_filling(arr, True)
    print "%d sinks detected at:"%ns
    print zip(row, col)
    print "which were filled with values:"
    print fval
    print "filled array"
    print filled_arr
    return None

