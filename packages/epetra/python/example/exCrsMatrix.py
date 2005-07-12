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

# Shows how to build an Epetra_CrsMatrix in python

try:
  import setpath
  import Epetra
except:
  try:
    from PyTrilinos import Epetra
  except ImportError:
    raise ImportError, "error w/ Epetra"

def main():
    Comm  = Epetra.SerialComm()
    NumGlobalRows = 5
    Map   = Epetra.Map(NumGlobalRows, 0, Comm)
    A     = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
    for i in range(0, NumGlobalRows):
      if i != NumGlobalRows - 1:
        Indices = [i, i + 1]
        Values = [2.0, -1.0]
      else:
        Indices = [i]
        Values = [2.0];
    A.InsertGlobalValues(i, Values, Indices);
    ierr = A.FillComplete();

    print A
    print "inf norm of A =", A.NormInf()

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main()
