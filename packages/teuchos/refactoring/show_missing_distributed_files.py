#!/bin/env python
#

# Used to find out files that have not been moved into a subpackage yet.  Run
# under teuchos soruce directory.
#
#

import sys
import os
import commands


from optparse import OptionParser

#
# Functions
#


def convertToSortedLines(aSet):
  arrayOfLines = ("\n".join(aSet)).split("\n")
  arrayOfLines.sort()
  return "\n".join(arrayOfLines)


def getFilesListSet(opts):
  rawFilesList = commands.getoutput("ls").split("\n")
  if opts.excludeFileNames:
    excludeFilesList = opts.excludeFileNames.split(":")
    for excludeFile in excludeFilesList:
      try:
        # Ignore if you can't remove it
        rawFilesList.remove(excludeFile)
      except Exception, e:
        None
  return set(rawFilesList)


#
# Executable statments
#

# A) Get commandline options

clp = OptionParser()
clp.add_option("--ref-dir", dest="refDir", type="string", default="")
clp.add_option("--other-dirs", dest="otherDirs", type="string", default="", \
  help="path1:path2:...")
clp.add_option("--exclude-file-names", dest="excludeFileNames", type="string", default="", \
  help="file1:file2:...")
(options, args) = clp.parse_args()

baseDir = os.getcwd()

print "Differencing files in reference directory'to other directories ...\n"

os.chdir(options.refDir)
refFilesListSet = getFilesListSet(options)
os.chdir(baseDir)
print "Reference directory = '"+options.refDir+"', num files = "+str(len(refFilesListSet))
#print "\nrefDir="+options.refDir+", refFilesListSet:\n", convertToSortedLines(refFilesListSet)

otherDirsList = options.otherDirs.split(":")
otherDirsFilesListSet = set()

for otherDir in otherDirsList:
  os.chdir(otherDir)
  otherDirFilesSet = getFilesListSet(options)
  os.chdir(baseDir)
  print "Other directory = '"+otherDir+"', num files = "+str(len(otherDirFilesSet))
  #print "\notherDir="+otherDir+", otherDirFilesSet:\n", convertToSortedLines(otherDirFilesSet)
  duplicatedOtherFilesSet = otherDirsFilesListSet.intersection(otherDirFilesSet)
  if (len(duplicatedOtherFilesSet) > 0):
    print "Error!, the following files in this directory are duplicateded in" + \
      " other dirs not including --ref-dir:\n", convertToSortedLines(duplicatedOtherFilesSet)
    sys.exit(1)
  otherDirsFilesListSet = otherDirsFilesListSet.union(otherDirFilesSet)

nonConvertedSourceFiles = refFilesListSet.difference(otherDirsFilesListSet)

print "\nFiles in reference directory not in other directories"+ \
  " (num = "+str(len(nonConvertedSourceFiles))+"):\n\n", \
  convertToSortedLines(nonConvertedSourceFiles)
