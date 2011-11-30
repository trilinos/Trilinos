#!/usr/bin/env python

# This file allows project-level configuration of the checkin-test
# system. If a project does not require any special options, then this
# file can be replaced by a simple symlink.

# This is a dictionary with key-value pairs that map to the script's
# command line arguments. These are used to add project-specific
# defaults for the arguments that all developers should use.
configuration = {
    'extra-cmake-options': '-DTPL_ENABLE_Pthread:BOOL=OFF -DTPL_ENABLE_BinUtils:BOOL=OFF'
    'send-email-to-on-push': 'trilinos-checkin-tests@software.sandia.gov',
}

# Load the main implementation and execute.
import sys
from os.path import dirname, abspath, realpath, join
sys.path.append(join(dirname(abspath(realpath(__file__))), 'cmake', 'tribits', 'python'))
import CheckinTestImpl
CheckinTestImpl.main(configuration)
