#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#         PyTrilinos.New_Package: Python Interface to New_Package
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

# Imports.  Users importing an installed version of PyTrilinos should use the
# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
import setpath
import Epetra
import New_Package
import unittest

####################################################################

# Swahili support
haveSwahili = "Newp_Jambo" in dir(New_Package)

####################################################################

class NewPackageTestCase(unittest.TestCase):
    "TestCase class for New_Package module"

    def testVersion(self):
        "Test New_Package New_Package_Version function"
        front   = "New_Package Version "
        version = New_Package.New_Package_Version()
        self.assertEquals(version[:len(front)], front)

    def setUp(self):
        self.comm = Epetra.SerialComm()

    def testHelloConstructor(self):
        "Test New_Package Hello constructor"
        hello = New_Package.Newp_Hello(self.comm)
        # No assert here; just executing without an exception is enough

    def testHelloPrint(self):
        "Test New_Package Hello print operator"
        hello = New_Package.Newp_Hello(self.comm)
        string = hello.__str__()
        self.assertEquals(string, "This will print out one line for each" +
                          " of the 1 processes \n\nHello.  I am process 0")

    if haveSwahili:
        def testJamboConstructor(self):
            "Test New_Package Jambo constructor"
            jambo = New_Package.Newp_Jambo(self.comm)
            # No assert here; just executing without an exception is enough

        def testJamboPrint(self):
            "Test New_Package Jambo print operator"
            jambo = New_Package.Newp_Jambo(self.comm)
            string = jambo.__str__()
            self.assertEquals(string, "This will print out one line for each" +
                              " of the 1 processes \n\nJambo.  I am process 0")

####################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(NewPackageTestCase))

    # Run the test suite
    unittest.TextTestRunner(verbosity=2).run(suite)
