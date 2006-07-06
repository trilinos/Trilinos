#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#              PyTrilinos.Epetra: Python Interface to Epetra
#                   Copyright (2005) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

# ----------------------------------------------------------------------------
# This example shows how to work with FECrsMatrix to set non-local elements.
#
# This example should be run with more than one processor. A diagonal matrix
# is created, and all elements are set on processor zero only. Then the
# matrix is distributed, and the local entries are modified.
#
# \author Marzio Sala, SNL 9214
#
# \date Last updated on 02-Aug-05
# ----------------------------------------------------------------------------

# Standard procedure to load PyTrilinos modules.
from optparse import *

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
    import setpath
    import Epetra
else:
  try:
    import setpath
    import Epetra
  except ImportError:
    from PyTrilinos import Epetra
    print "Using system-installed Epetra"

def main():

  # Defines a communicator.
  Comm = Epetra.PyComm()
  NumGlobalRows = 10
  Map = Epetra.Map(NumGlobalRows, 0, Comm)
  Matrix = Epetra.FECrsMatrix(Epetra.Copy, Map, 0)

  # Fills the entire matrix on processor 0, setting both local
  # and non-local entries. Communication among processors will occur when
  # GlobalAssemble() is called. After GlobalAssemble(), elements cannot be
  # added to matrix; although the value local elements can be modified.
  #
  # Method Matrix[GlobalRow, GlobalCol] requires global IDs for rows
  # and columns.
  if Comm.MyPID() == 0:
    for i in xrange(NumGlobalRows):
      # Here i is the global ID or a given row
      Matrix[i, i] = 1.0
  
  Matrix.GlobalAssemble()
  print Matrix
  
  Comm.Barrier()
  
  # Gets a list containing the global ID or each local row 
  MyGlobalElements = Map.MyGlobalElements()
  
  # We can now extract local rows of the matrix. Method Matrix[i] returns
  # the nonzeros indices and values for the locally hosted global row `i'.
  for i in MyGlobalElements:
    print Matrix[i]
  
  # new reset the local values.
  for i in MyGlobalElements:
    Indices, Values = Matrix[i]
    for j in xrange(len(Indices)):
      Matrix[i, Indices[j]] = 10 * Values[j]
  
  print Matrix

  # synchronize processors
  Comm.Barrier()

  if Comm.MyPID() == 0: print "End Result: TEST PASSED"

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
