#!/bin/env python3
"""This is superclass for all Log packages"""

class Log(object):
    def states(self):
        raise

    def transitions(self):
        raise

    def run(self, filename):
        raise
