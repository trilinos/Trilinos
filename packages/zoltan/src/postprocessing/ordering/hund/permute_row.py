#!/usr/bin/env python

## Script which compares orderings provided by Scotch and ParMetis by using
## zdrive on MatrixMarket files.

import re, sys, getopt, commands, time, os
import numpy as np

# Paths to executables and names
# User can modify this

# Displays help

def usage():
    """Print the usage information of the command"""
    print "Usage: ./ordering [options] inputfile"
    print "Options:"
    print "\t--help (-h):    This page"


#    displayresults(result, ('fflop', 'ftime', 'nzL', 'nzU', 'memory'))
#    print result


# The main routine.

def main(argv=None):
    if argv is None:
        argv = sys.argv

    filename=argv[1]
    permname=argv[2]

    infile = open(filename)
    inperm = np.loadtxt(permname)
    working = False
    for line in infile:
        if line[0] == '%':
            print line,
            continue
        if not working:
            working = True
            print line,
            continue
        values=line.split()
        print "%d %s %s" % (inperm[int(values[0])-1]+1, values[1], values[2])


if __name__ == "__main__":
    main()

