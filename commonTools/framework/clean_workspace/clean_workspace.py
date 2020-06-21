#!/usr/bin/env python3
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

import os

import argparse
import shutil

from clean_sentinel import clean_reference_date
from clean_sentinel import last_clean_date
from clean_sentinel import update_last_clean_date


class Cleaner(object):
    """Class to hold the capability
         Parse if we are at start or end
         Call appropriate function
    """

    def __init__(self):
        self.args = None

    def run(self):
        self.args = self.parse_args()
        if self.args.dir is None:
            raise SystemExit("No directory passed - exiting!")
        if self.args.force_clean:
            print("Cleaning directory {clean_dir} due to command line option".format(clean_dir=self.args.dir))
            self.force_clean_space()
        else:
            self.clean_space_by_date()

    def parse_args(self):
        """
        Get the options that were passed in and return them in a namespace
        """
        parser = argparse.ArgumentParser(description='Parser for the clean_workspace class')
        parser.add_argument('dir',
                            help="directory to clean",
                            action="store")
        parser.add_argument('--force-clean',
                            help="force a local cleanup",
                            action="store_true",
                            default=False)
        options = parser.parse_args()
        return options

    def force_clean_space(self):
        """Do the actual cleanup
             basically just os.unlink()
        """
        if os.path.isdir(self.args.dir):
            shutil.rmtree(self.args.dir)

    def clean_space_by_date(self):
        if last_clean_date() < clean_reference_date():
            print("Cleaning directory {clean_dir} due to newer reference date".format(clean_dir=self.args.dir))
            self.force_clean_space()
            update_last_clean_date()
        
if __name__ == '__main__':
    clean = Cleaner()
    clean.run()
