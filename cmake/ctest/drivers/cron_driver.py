#!/usr/bin/env python

# Locate the TriBITS system level cron_driver.
import os
import sys

def join_paths(*args):
    """ Put together a list of paths. """
    return reduce(os.path.join, args)
    
this_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

# Get the path to the TriBITS build system.
tribits_path = os.path.abspath(join_paths(this_path, '..', '..', 'tribits'))
tdd_path = join_paths(tribits_path, 'dashboard_driver')
sys.path.append(tdd_path)
print "sys.path =", sys.path

# Get the path to the repository root.
repo_path = os.path.abspath(join_paths(this_path, '..', '..', '..'))

import tdd_driver
tdd_driver.run_driver(this_path, repo_path)
