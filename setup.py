#! /usr/bin/env python
"""
Setup file for ADDEM

Created: Wed Mar 16, 2016  02:41PM
Last modified: Tue May 03, 2016  03:48PM

"""
import os
from distutils.core import setup, Extension
import numpy as np
import addem


# Utility function to read the README file.
# Source: http://pythonhosted.org/an_example_pypi_project/setuptools.html
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


# define extension module objects to pass as args to setup()
ext_mod_sinks = Extension('sinks',
                          sources=['src/sinks.c'],
                          include_dirs=[np.get_include()],
                          )
ext_mod_flows = Extension('flows',
                          sources=['src/flows.c'],
                          include_dirs=[np.get_include()],
                          )


setup(name='addem',
      version=addem.__version__,
      description='Analysis of Distributions from Digital Elevation Models',
      long_description=read('README.md'),
      author='Bedartha Goswami',
      author_email='goswami@uni-potsdam.de',
      url='https://github.com/bedartha/addem',
      license='GNU',
      packages=['addem'],
      py_modules=['addem/distributions'],
      ext_package='addem',
      ext_modules=[
                    ext_mod_sinks,
                    ext_mod_flows,
                    ],
      install_requires=['numpy'],
      classifiers=[
                   'Development Status :: 2 - Pre-Alpha',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Physics',
                   ],
      keywords='research climate geoscience'
      )
