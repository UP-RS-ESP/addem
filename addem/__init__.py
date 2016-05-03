#! /usr/bin/env python
# __INIT__.PY

try:
    from addem import flows
    from addem import distributions
    from addem import sinks
except ImportError:
    pass

__all__ = ["flows", "distributions", "sinks"]
__version__ = '0.0.1.dev1'
