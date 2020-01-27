#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
# Change above to '/usr/bin/python -3' for python 3.x porting warnings
"""
This script will clean the current workspace IFF
   It is called post-build
  or
   The workspace has not been cleaned since the date in the sentinel file
"""
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import pickle

from clean_sentinel import set_reference_date


class CleanReference(object):
    """Sets the current time as the reference clean date"""

    def run(self):
        """Run the actual routine in clean_sentinel"""
        try:
            set_reference_date()
        except BaseException as e:
            print(e.args, file=sys.stderr)
            sys.exit()
        except pickle.PicklingError:
            print('Pickle Failure', file=sys.stderr)
            sys.exit()

if __name__ == '__main__':
    cleanRef = CleanReference()
    cleanRef.run()
