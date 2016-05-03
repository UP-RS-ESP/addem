Analysis of Distributions from Digital Elevation Models (ADDEM)
---------------------------------------------------------------

This repository contains a Python implementation of ongoing efforts to
analyse the distributions of slopes and curvatures of Digital Elevation
Models (DEMs) in the hope that these would be able to provide a more
wholesome perspective that can supplement the traditional approaches which
involve the analysis of the slope and curvatures of the only the river flows
obtained from the DEMs. 

**NOTE: This is a work in progress! There might be serious bugs in the code!**

Objectives
----------

The suite of Python codes in this repository has two main objectives: (i) a
high-speed, parallelised Python implementation of traditional DEM analyses
such as sink filling, flow direction estimation, flow accumulation
estimation, etc., and (ii) to create a suite of new functions that enable
the end-user to analyse the distributions of various geophysically relevant
quantities that are estimable from DEMs, such as the slop, curvature and
drainage areas. The idea is that the end-user should be able to work all of
the routines available in this repository with the full utilisation of
typically available multi-core, multi-processor architectures of present day
computers.

In this regard, we have currently set the following two
categories of *goals* for the project:

### Short-term goals
* Sink filling
* Flow direction
* Flow accumulation
* River flow networks
* Drainage areas
* Catchment areas

### Long-term goals
* Representation, abstraction and analysis of probability distributions
* Release of the code as a cross-platform software package

Installation
------------

1. Clone the GitHub repository:
```bash
git clone https://github.com/bedartha/addem.git
cd addem/
```

2. Use the Python setup file to install. You may also consider the use of
   [Virtual Envinonments](http://docs.python-guide.org/en/latest/dev/virtualenvs/):
```bash
python setup.py build
python setup.py install
```

Usage
-----

After including the proper path in the Python paths list, you can simply
``import addem`` as a module and test out the dummy ``hello``
functions given in the submodules ``addem.distributions`` and
``addem.flows``.

```python
import addem
addem.distributions.hello()
# Hello! I am the 'addem.distributions' module!
addem.flows.hello()
#Hello! I am the 'addem.flows' module!
```

You can query the help documentation about the entire module or its
submodules as well.

```python
help(addem) # shows help documentation on the entire module
help(addem.flows) # shows help documentation on the flows submodule
```

Code Example
------------

*coming soon*


License
-------

Copyright (c) 2016 Bedartha Goswami 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
MA  02110-1301,USA.

Authors
-------

* Bedartha Goswami <goswami@uni-potsdam.de>
* Aljoscha Rheinwalt <aljoscha.rheinwalt@uni-potsdam.de>
* Bodo Bookhagen <bodo.bookhagen@uni-potsdam.de>

About this file
---------------

Created: Wed Mar 09, 2016  03:23PM
Last modified: Tue May 03, 2016  02:19PM

