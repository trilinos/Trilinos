#! /usr/bin/env python

# @header
# ************************************************************************
#
#           PyTrilinos.EpetraExt: Python interface to EpetraExt
#                   Copyright (2005) Sandia Corporation
#
# Under terms of contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# license, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @header

from   optparse import *
import sys

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
    import EpetraExt
else:
    try:
        import setpath
        import Epetra, EpetraExt
    except ImportError:
        from PyTrilinos import Epetra, EpetraExt
        print >>sys.stderr, "Using system-installed Epetra, EpetraExt"

# Create a global communicator
comm    = Epetra.PyComm()
numProc = comm.NumProc()
iAmRoot = comm.MyPID() == 0

def main():
    n   = 10 * numProc
    map = Epetra.Map(n, 0, comm)

    # ================================================================= #
    # Creates two matrices, one is diagonal (A), the other contains the #
    # first sub- and super-diagonal (B), then sum them (B = A + B).     #
    # Note that B cannot be FillComplete()'d before calling Add()       #
    # unless it already contains the structure of A + B.                #
    # ================================================================= #

    A    = Epetra.CrsMatrix(Epetra.Copy, map, 0)
    rows = map.MyGlobalElements()
    for i in rows:
        err = A.InsertGlobalValues(i, 2.0, i)
        if err < 0:
            raise RunTimeError, "Processor %d, global row %d of A, error code %d" \
                  % (comm.MyPID(), i, err)
    A.FillComplete()

    B = Epetra.CrsMatrix(Epetra.Copy, map, 0)
    for i in rows:
        indices = [ ]
        values  = [ ]
        if i > 0:
            indices.append(i-1)
            values.append(-1)
        elif i < n - 1:
            indices.append(i+1)
            values.append(-1)
        err = B.InsertGlobalValues(i, values, indices)
        if err < 0:
            raise RunTimeError, "Processor %d, global row %d of B, error code %d" \
                  % (comm.MyPID(), i, err)

    EpetraExt.Add(A, False, 1.0, B, 1.0)
    print B

    return 0

if __name__ == "__main__":
    failures = main()
    failures = comm.SumAll(failures)
    if failures == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(failures)
