#! /usr/bin/env python

# tree.py
#
# Written by Doug Dahms
#
# Modified by Ross Bartlett
#
# Prints the tree structure for the path specified on the command line

from os import listdir, sep
from os.path import abspath, basename, isdir
from sys import argv, exit
from optparse import OptionParser

def tree(dir, padding, print_files=False, print_compact=False, top_level=False):
  init_prefix = padding[:-1]
  if top_level: init_prefix += '  '
  else: init_prefix += '+-'
  print init_prefix + basename(abspath(dir)) + '/'
  padding = padding + ' '
  files = []
  if print_files:
    files = listdir(dir)
  else:
    files = [x for x in listdir(dir) if isdir(dir + sep + x)]
  count = 0
  files.sort()
  for file in files:
    count += 1
    if not print_compact:
      print padding + '|'
    path = dir + sep + file
    if isdir(path):
      if count == len(files):
        tree(path, padding + ' ', print_files, print_compact)
      else:
        tree(path, padding + '|', print_files, print_compact)
    else:
      print padding + '+-' + file

usageHelp = r""" tree.py [-f] [-c] <PATH>
Print tree structure of specified <PATH>.
"""

def main():

  # Get command-line arguments

  clp = OptionParser(usage=usageHelp)

  clp.add_option(
    "-f", dest="printFiles", action="store_true",
    help="Show files in addition to just directoires.",
    default=False )

  clp.add_option(
    "-c", dest="printCompact", action="store_true",
    help="Make output more compact.",
    default=False )
  
  (options, args) = clp.parse_args()

  path = args[0]

  if not isdir(path):
    print "ERROR: \'"+path+"\' is not a directory!"
    print "See --help!"
    exit(1)

  tree(path, ' ', options.printFiles, options.printCompact, True)


if __name__ == '__main__':
  main()
