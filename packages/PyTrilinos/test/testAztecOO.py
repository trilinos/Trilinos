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
from   optparse import *
import sys
import unittest

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
    import AztecOO
else:
    try:
        import setpath
        import Epetra
        import AztecOO
    except ImportError:
        from PyTrilinos import Epetra
        from PyTrilinos import AztecOO
        print >>sys.stderr, "Using system-installed Epetra, AztecOO"

####################################################################

class TestOperator(Epetra.Operator):
    """
    This is a dummy Epetra.Operator subclass.  Only a bare minimum of member
    functions will be implemented.
    """
    def __init__(self, map):
        Epetra.Operator.__init__(self)
        self.__comm  = map.Comm()
        self.__map   = map
        self.__label = "Test Operator"
    def Label(self):             return self.__label
    def Comm(self):              return self.__comm
    def OperatorDomainMap(self): return self.__map
    def OperatorRangeMap(self):  return self.__map
    def HasNormInf(self):        return False

##########################################################################

class TestRowMatrix(Epetra.RowMatrix):
    """
    This is a dummy Epetra.RowMatrix subclass.  Only a bare minimum of member
    functions will be implemented.
    """
    def __init__(self, map):
        Epetra.RowMatrix.__init__(self)
        self.__comm = map.Comm()
        self.__map  = map
        self.__label = "Test Operator"
    def Label(self):             return self.__label
    def Comm(self):              return self.__comm
    def OperatorDomainMap(self): return self.__map
    def OperatorRangeMap(self):  return self.__map
    def HasNormInf(self):        return False
    def RowMatrixImporter(self): return Epetra.Import(self.__map,self.__map)

####################################################################

class AztecOOTestCase(unittest.TestCase):
    "TestCase class for AztecOO module"

    def setUp(self):
        self.comm = Epetra.PyComm()
        self.comm.Barrier()
        self.size = self.comm.NumProc() * 8
        self.map  = Epetra.Map(self.size, 0, self.comm)

    def tearDown(self):
        # This will help tame the printing
        self.comm.Barrier()

    def testVersion(self):
        "Test AztecOO AztecOO_Version function"
        front   = "AztecOO Version "
        version = AztecOO.AztecOO_Version()
        self.assertEquals(version[:len(front)], front)

    def testConstructor0(self):
        "Test AztecOO default constructor"
        solver = AztecOO.AztecOO()
        self.failUnless(isinstance(solver, AztecOO.AztecOO))
        #print solver.CheckInput()

    def testConstructor1(self):
        "Test AztecOO (Epetra.Operator, Epetra.MultiVector, Epetra.MultiVector) constructor"
        op = TestOperator(self.map)
        x  = Epetra.Vector(self.map)
        b  = Epetra.Vector(self.map)
        x.Random()
        b.Random()
        solver = AztecOO.AztecOO(op, x, b)
        self.failUnless((solver.GetLHS() == x).all())
        self.failUnless((solver.GetRHS() == b).all())

    def testConstructor2(self):
        "Test AztecOO (Epetra.RowMatrix, Epetra.MultiVector, Epetra.MultiVector) constructor"
        rm = TestRowMatrix(self.map)
        x  = Epetra.Vector(self.map)
        b  = Epetra.Vector(self.map)
        x.Random()
        b.Random()
        solver = AztecOO.AztecOO(rm, x, b)
        self.failUnless((solver.GetLHS() == x).all())
        self.failUnless((solver.GetRHS() == b).all())

    def testConstructor3(self):
        "Test AztecOO (Epetra.LinearProblem) constructor"
        rm = TestRowMatrix(self.map)
        x  = Epetra.Vector(self.map)
        b  = Epetra.Vector(self.map)
        x.Random()
        b.Random()
        lp = Epetra.LinearProblem(rm, x, b)
        solver = AztecOO.AztecOO(lp)
        self.failUnless((solver.GetLHS() == x).all())
        self.failUnless((solver.GetRHS() == b).all())

    # This currently triggers a bus error
    #def testConstructor4(self):
    #    "Test AztecOO copy constructor, with non-Epetra.LinearOperator source"
    #    op = TestOperator(self.map)
    #    x  = Epetra.Vector(self.map)
    #    b  = Epetra.Vector(self.map)
    #    x.Random()
    #    b.Random()
    #    solver1 = AztecOO.AztecOO(op, x, b)
    #    solver2 = AztecOO.AztecOO(solver1)
    #    self.failUnless((solver2.GetLHS() == x).all())
    #    self.failUnless((solver2.GetRHS() == b).all())

    def testConstructor5(self):
        "Test AztecOO copy constructor with Epetra.LinearOperator source"
        rm = TestRowMatrix(self.map)
        x  = Epetra.Vector(self.map)
        b  = Epetra.Vector(self.map)
        x.Random()
        b.Random()
        lp = Epetra.LinearProblem(rm,x,b)
        solver1 = AztecOO.AztecOO(lp)
        solver2 = AztecOO.AztecOO(solver1)
        self.failUnless((solver2.GetLHS() == x).all())
        self.failUnless((solver2.GetRHS() == b).all())

    def testSetProblem(self):
        "Test AztecOO SetProblem method"
        rm = TestRowMatrix(self.map)
        x  = Epetra.Vector(self.map)
        b  = Epetra.Vector(self.map)
        x.Random()
        b.Random()
        lp = Epetra.LinearProblem(rm, x, b)
        solver = AztecOO.AztecOO()
        self.failUnless(solver.GetProblem() is None)
        solver.SetProblem(lp)
        self.failUnless(isinstance(solver.GetProblem()     , Epetra.LinearProblem))
        self.failUnless(isinstance(solver.GetUserOperator(), Epetra.Operator     ))
        self.failUnless(isinstance(solver.GetUserMatrix()  , Epetra.RowMatrix    ))
        self.failUnless((solver.GetLHS() == x).all())
        self.failUnless((solver.GetRHS() == b).all())

    def testSetUserOperator(self):
        "Test AztecOO SetUserOperator method"
        op = TestRowMatrix(self.map)
        solver = AztecOO.AztecOO()
        self.failUnless(solver.GetUserOperator() is None)
        solver.SetUserOperator(op)
        self.failUnless(isinstance(solver.GetUserOperator(), Epetra.Operator))

####################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(AztecOOTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n*******************\nTesting AztecOO\n*******************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
