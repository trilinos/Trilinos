#!/usr/bin/env python

"""
Python module/program that compares to Thyra SPMD multi-vector files
and checks that the relative error is between a given set of bounds.
"""


#
# 2007/08/06: rabartl: ToDo:
#
# * Add the ability to read multiple files, one for each process,
#   when computing the diffs and norms.  This would take some work
#   but would not be too hard to do.


#
# Import commands
#

import sys
import os
import re
import math


#
# Functions
#


# Read in a multi-vector from file
def readInMultiVector( multiVectorFile ):

  # Read the header
  line = multiVectorFile.readline();

  (numRows,numCols) = line.split();
  numRows = int(numRows)
  numCols = int(numCols)
  #print "\nnumRows = ", numRows
  #print "\nnumCols = ", numCols

  # Read the rows of the multi-vector!

  multiVec = []
  for i in range(numRows):
    line = multiVectorFile.readline();
    indexAndVals = line.split()
    assert len(indexAndVals) == (numCols+1)
    vals = indexAndVals[1:]
    assert len(vals) == numCols
    for val in vals:
      multiVec.append(float(val));

  return multiVec


# Compute the difference
def diffMultiVectors(
  multiVector1, multiVector2
  ):
  assert len(multiVector1) == len(multiVector2)
  multiVectorDiff = []
  for i in range(len(multiVector1)):
    multiVectorDiff.append( multiVector1[i] - multiVector2[i] )
  return multiVectorDiff


# Compute the norm of the difference
def computeNormDiff( multiVectorDiff, normType ):

  normVal = 0.0

  if normType == "inf":

    for x in multiVectorDiff:
      normVal = max(normVal,math.fabs(x))

  elif normType == "one":

    for x in multiVectorDiff:
      normVal += math.fabs(x)

  else:
    
    raise "\nError, invalid norm type \"" + normType + "\"!"
  
  return normVal


# Compare to multi-vector files
def compareThyraSpmdMultiVectorFiles(
  multiVectorFileName1,
  multiVectorFileName2,
  normType,
  minNormError,
  maxNormError,
  verbose,
  dumpAll
  ):

  # Assert the input arguments
  assert multiVectorFileName1 != ""
  assert multiVectorFileName2 != ""
  assert normType != ""
  assert minNormError >= 0.0
  assert maxNormError >= minNormError

  if verbose:
    print "Comparing multi-vectors in files ", \
          multiVectorFileName1, " and ", \
          multiVectorFileName2, " ..."
  
  # Open the files
  multiVectorFile1 = open(multiVectorFileName1,'r')
  multiVectorFile2 = open(multiVectorFileName2,'r')

  # Read the multi-vectors into memory
  multiVector1 = readInMultiVector(multiVectorFile1)
  multiVector2 = readInMultiVector(multiVectorFile2)
  if dumpAll:
    print "multiVector1 = ", multiVector1
    print "multiVector2 = ", multiVector2
  assert len(multiVector1) == len(multiVector2)
  
  # Compute the difference between the multi-vectors and norm of the
  # difference.
  multiVectorDiff = diffMultiVectors(multiVector1, multiVector2)
  normDiff = computeNormDiff(multiVectorDiff,normType)
  if dumpAll:
    print "multiVectorDiff = ", multiVectorDiff

  if verbose:
    print "norm(diff) = ", normDiff
  
  # Check the error

  testPassed = True

  if normDiff < minNormError:
    if verbose:
      print "Error, normDiff = ", normDiff, \
            " < minNormDiff = ", minNormError, "!"
      testPassed = False

  if normDiff > maxNormError:
    if verbose:
      print "Error, normDiff = ", normDiff, \
            " > maxNormDiff = ", maxNormError, "!"
      testPassed = False

  return testPassed


#
# Main program
#

from optparse import OptionParser

if __name__ == '__main__':

  #
  # Read in the commandline arguments
  #

  clp = OptionParser()

  clp.add_option(
    "--file1", dest="multiVectorFileName1", type="string",
    help="File name for first mulit-vector to compare" );

  clp.add_option(
    "--file2", dest="multiVectorFileName2", type="string",
    help="File name for second mulit-vector to compare" );

  clp.add_option(
    "--norm-type", dest="normType", type="string", default="inf",
    help="The type of norm to compute: one, two, or inf" );

  clp.add_option(
    "--min-norm-err", dest="minNormError", type="float", default=0.0,
    help="The minimum the norm of the diff can be" );

  clp.add_option(
    "--max-norm-err", dest="maxNormError", type="float", default=0.0,
    help="The minimum the norm of the diff can be" );

  clp.add_option(
    "-v", "--verbose", action="store_true", dest="verbose",default=True);
  clp.add_option("-q", "--quiet", action="store_false", dest="verbose")

  clp.add_option(
    "--dump", action="store_true", dest="dumpAll" ,default=False );
    
  (options,args) = clp.parse_args()

  #print "\nmultiVectorFileName1 = ", options.multiVectorFileName1


  #
  # Compare the multi-vector files
  #

  success = compareThyraSpmdMultiVectorFiles(
    options.multiVectorFileName1,
    options.multiVectorFileName2,
    options.normType,
    options.minNormError,
    options.maxNormError,
    options.verbose,
    options.dumpAll
    )

  if success:
    sys.exit(0)
    
  sys.exit(1) # Failure!
