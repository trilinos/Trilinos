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

# Imports
import setpath
import unittest
from   Numeric    import *
from   PyTrilinos import Epetra

##########################################################################

class EpetraCrsGraphTestCase(unittest.TestCase):
    "TestCase class for Epetra CrsGraphs"

    def setUp(self):
        self.size = 11
        self.comm = Epetra.SerialComm()
        self.map  = Epetra.Map(self.size, 0, self.comm)

    def testConstructor1(self):
        "Test Epetra.CrsGraph constructor with fixed number of indices per row"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.map, 3)

    def testInsertGlobalIndices(self):
        "Test Epetra.CrsGraph InsertGlobalIndices method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.map, 3)
        crsg.InsertGlobalIndices(0,2,array([0,1]))
        for i in range(1,self.size-1):
            crsg.InsertGlobalIndices(i,3,array([i-1,i,i+1]))
        crsg.InsertGlobalIndices(self.size-1,2,array([self.size-2,self.size-1]))

    def testInsertMyIndices(self):
        "Test Epetra.CrsGraph InsertMyIndices method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.map, 3)
        crsg.InsertMyIndices(0,2,array([0,1]))
        for i in range(1,self.size-1):
            crsg.InsertMyIndices(i,3,array([i-1,i,i+1]))
        crsg.InsertMyIndices(self.size-1,2,array([self.size-2,self.size-1]))

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraCrsGraphTestCase))

    # Run the test suite
    print "\n***********************\nTesting Epetra.CrsGraph\n***********************\n"
    unittest.TextTestRunner(verbosity=2).run(suite)
