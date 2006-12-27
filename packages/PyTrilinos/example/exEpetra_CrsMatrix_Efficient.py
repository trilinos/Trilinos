#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

# Shows how to build an Epetra_CrsMatrix in python
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
  except:
    from PyTrilinos import Epetra
    print "Using system-installed Epetra"

def main():
    # Creates an Epetra.SerialComm in serial mode, or an Epetra.MpiComm
    # if configured with MPI support
    Comm  = Epetra.PyComm()
    NumGlobalRows = 5
    Map   = Epetra.Map(NumGlobalRows, 0, Comm)
    A     = Epetra.CrsMatrix(Epetra.Copy, Map, 0)
    NumMyRows = Map.NumMyElements()
    # Loop over all local rows to create a tridiagonal matrix
    for ii in xrange(NumMyRows):
      # `i' is the global ID of local ID `ii'
      i = Map.GID(ii)
      if i != NumGlobalRows - 1:
        Indices = [i, i + 1]
        Values  = [2.0, -1.0]
      else:
        Indices = [i]
        Values  = [2.0]
      A.InsertGlobalValues(i, Values, Indices)
    # transform the matrix into local representation -- no entries
    # can be added after this point. However, existing entries can be
    # modified.
    ierr = A.FillComplete()

    print A
    NormInf = A.NormInf()
    if Comm.MyPID() == 0:
      print "inf norm of A =", NormInf
 
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
