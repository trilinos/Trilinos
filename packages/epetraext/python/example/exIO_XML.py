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
        print >>sys.stderr, "Using system-installed Epetra, EptraExt"

# Build a global communicator
comm    = Epetra.PyComm()
numProc = comm.NumProc()
myPID   = comm.MyPID()
iAmRoot = myPID == 0

def main():

    tolerance = 1e-5

    # Construct a vector x and populate with random values.
    # Then construct a diagonal Epetra.CrsMatrix
    n   = 10 * numProc
    map = Epetra.Map(n, 0, comm)
    x   = Epetra.Vector(map)
    x.Random()

    matrix = Epetra.CrsMatrix(Epetra.Copy, map, 0)

    for lrid in range(matrix.NumMyRows()):
        grid = matrix.GRID(lrid)
        matrix[grid,grid] = grid

    matrix.FillComplete()

    # -------------------------------- #
    # Part I: Write elements to a file #
    # -------------------------------- #

    XMLWriter = EpetraExt.XMLWriter(comm, "data.xml")

    XMLWriter.Create("test xml file")
    if iAmRoot: print "Writing map ...",
    XMLWriter.Write("map"   , map   )
    if iAmRoot: print "ok"
    if iAmRoot: print "Writing vector ...",
    XMLWriter.Write("x"     , x     )
    if iAmRoot: print "ok"
    if iAmRoot: print "Writing matrix ...",
    XMLWriter.Write("matrix", matrix)
    if iAmRoot: print "ok"
    XMLWriter.Close()

    # --------------------------------- #
    # Part II:Read elements from a file #
    # --------------------------------- #

    XMLReader = EpetraExt.XMLReader(comm, "data.xml")

    if iAmRoot: print "Reading map ...",
    map2 = XMLReader.ReadMap("map")
    if map2 is None:
        print "FAIL:  map2 is None on processor", myPID
        return 1
    if iAmRoot: print "ok"

    if iAmRoot: print "Reading vector ...",
    x2 = XMLReader.ReadMultiVector("x")
    if x2 is None:
        print "FAIL:  x2 is None on processor", myPID
        return 1
    if iAmRoot: print "ok"

    if iAmRoot: print "Reading matrix ...",
    matrix2 = XMLReader.ReadCrsMatrix("matrix")
    if matrix2 is None:
        print "FAIL: matrix2 is None on processor", myPID
        return 1
    if iAmRoot: print "ok"

    x2.Update(1.0, x, -1.0)
    norm = x2.Norm2()

    if iAmRoot: print "Checking tolerance ...",
    if abs(norm) < tolerance:
        if iAmRoot: print "ok"
    else:
        if iAmRoot:
            print "FAILED:"
            print "  abs(norm) =", abs(norm)
            print "  tolerance =", tolerance
        return 1

    return 0

################################################################

if __name__ == "__main__":
    failures = main()
    failures = comm.SumAll(failures)
    if failures == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(failures)
