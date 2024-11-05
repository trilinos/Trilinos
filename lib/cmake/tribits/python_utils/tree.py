#! /usr/bin/env python3

# tree.py
#
# Written by Doug Dahms
#
# Modified by Ross Bartlett
#
# Prints the tree structure for the path specified on the command line

import sys
from os import listdir, sep
from os.path import abspath, basename, isdir
from sys import argv, exit
from optparse import OptionParser

def tree(dir, padding, options, depth, top_level=False):

  if depth != None:
    if depth == 0:
      return
    else:
      depth = depth - 1

  print_files = options.printFiles
  print_compact = options.printCompact

  if options.noDirectorySep:
    verticalLineChar = ' '
    fileDirPrefix = '  '
  else:
    verticalLineChar = '|'
    fileDirPrefix = '+-'
    
  init_prefix = padding[:-1]
  if top_level: init_prefix += '  '
  else: init_prefix += fileDirPrefix

  print(init_prefix + basename(abspath(dir)) + '/')

  padding = padding + ' '

  raw_listdir = listdir(dir)
  files = []
  if print_files:
    files = raw_listdir
  else:
    files = [x for x in raw_listdir if isdir(dir + sep + x)]
  files.sort()
  for file in files:
    if (file[0] == '.') and (options.hideHiddenFiles):
      continue
    if not print_compact:
      print(padding + verticalLineChar)
    path = dir + sep + file
    if isdir(path):
      if (depth == None) or (depth > 0):
        tree(path, padding + verticalLineChar, options, depth)
      else:
        print(padding + fileDirPrefix + file + "/")
    else:
      print(padding + fileDirPrefix + file)


usageHelp = r""" tree.py [-f] [-c] [-x] [-n] <PATH>
Print tree structure of specified <PATH> to given depth.
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

  clp.add_option(
    "-x", dest="noDirectorySep", action="store_true",
    help="Remove the directory separators and continuation lines.",
    default=False )

  clp.add_option(
    "-n", dest="hideHiddenFiles", action="store_true",
    help="Hide hidden files.",
    default=False )

  clp.add_option(
    "--depth", dest="depth", type="int", default=None,
    help="Depth (integer) to recurse into.  Default = '' or unbounded.")
  
  (options, args) = clp.parse_args()

  if len(args) != 1:
    print("Error: Need to specify path a single path argument.  See, see --help")
    sys.exit(1)
  path = args[0]

  if not isdir(path):
    print("ERROR: \'" + path + "\' is not a directory!")
    print("See --help!")
    exit(1)

  depth = options.depth
  tree(path, ' ', options, depth, True)


if __name__ == '__main__':
  main()
