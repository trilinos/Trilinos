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

import sys

try:
  import setpath
  import Epetra, EpetraExt
except ImportError:
  from PyTrilinos import Epetra, EpetraExt
  print >>sys.stderr, "Using system-installed Epetra, EptraExt"

# Build a global communicator
comm    = Epetra.PyComm()
numProc = comm.NumProc()
iAmRoot = comm.MyPID() == 0

def main():

  failures = 0
  tolerance = 1e-5

  # Construct a vector x and populate with random values.
  # Then constructs a diagonal Epetra_CrsMatrix.
  n       = 10 * numProc
  map     = Epetra.Map(n, 0, comm)
  x       = Epetra.Vector(map)
  x.Random()

  matrix = Epetra.CrsMatrix(Epetra.Copy, map, 0);

  for lrid in range(matrix.NumMyRows()):
      grid = matrix.GRID(lrid)
      matrix[grid,grid] = grid

  matrix.FillComplete()

  # -------------------------------- #
  # Part I: Write elements to a file #
  # -------------------------------- #

  XMLWriter = EpetraExt.XMLWriter(comm, "data.xml");

  XMLWriter.Create("test xml file")
  XMLWriter.Write("map",    map)
  XMLWriter.Write("x",      x)
  XMLWriter.Write("matrix", matrix)
  XMLWriter.Close()

  # --------------------------------- #
  # Part II:Read elements from a file #
  # --------------------------------- #

  XMLReader = EpetraExt.XMLReader(comm, "data.xml");

  map2    = XMLReader.ReadMap("map")
  if map2 == None:
    failures += 1

  x2      = XMLReader.ReadMultiVector("x")
  if x2 == None:
    failures += 1

  matrix2 = XMLReader.ReadCrsMatrix("matrix")
  if matrix2 == None:
    failures += 1

  x2.Update(1.0, x, -1.0)
  norm = x2.Norm2()

  if abs(norm) < tolerance:
      if iAmRoot: print "ok"
  else:
      if iAmRoot: print "FAILED"
      failures += 1

  return failures

################################################################

if __name__ == "__main__":
  failures = main()
  failures = comm.SumAll(failures)
  if failures == 0 and iAmRoot: print "End Result: TEST PASSED"
  sys.exit(failures)
