#!/usr/bin/env python
"""
setup.py - script for building MyApplication

Usage:
      % python setup.py py2app
"""
#import sys
#sys.path.append("/Users/marzio/Trilinos/G4_SERIAL/lib/python2.3/site-packages/")
from distutils.core import setup
import py2app

setup(
    app=['Trilinos.py'],
)
