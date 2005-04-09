#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                  PyTrilinos: Rapid Prototyping Package
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
import NOX
import unittest

class ParameterTestCase(unittest.TestCase):
    "TestCase for NOX.Parameter.Lists"

    def setUp(self):
        self.pl = NOX.Parameter.List()
        self.pl.setParameter("Integer Parameter", 2003    )
        self.pl.setParameter("Float Parameter",   3.14    )
        self.pl.setParameter("String Parameter",  "Sandia")

    def testIntExistence(self):
        "Test NOX.Parameter.List integer existence"
        self.assertEqual(self.pl.isParameter("Integer Parameter" ), True )

    def testFloatExistence(self):
        "Test NOX.Parameter.List float existence"
        self.assertEqual(self.pl.isParameter("Float Parameter"   ), True )

    def testStringExistence(self):
        "Test NOX.Parameter.List string existence"
        self.assertEqual(self.pl.isParameter("String Parameter"  ), True )

    def testFakeExistence(self):
        "Test NOX.Parameter.List fake existence"
        self.assertEqual(self.pl.isParameter("Fake Parameter"    ), False)

    def testIntValues(self):
        "Test NOX.Parameter.List integer value"
        self.assertEqual(self.pl.getParameter("Integer Parameter"), 2003)

    def testFloatValues(self):
        "Test NOX.Parameter.List float value"
        self.assertEqual(self.pl.getParameter("Float Parameter"  ), 3.14)

    def testStringValues(self):
        "Test NOX.Parameter.List string value"
        self.assertEqual(self.pl.getParameter("String Parameter" ), "Sandia")

    def testSubList(self):
        "Test NOX.Parameter.List sublists"
        sl1 = self.pl.sublist("Sublist 1")
        self.assertEquals(self.pl.isParameter("Sublist 1"), True)
        sl2 = sl1.sublist("Sublist 2")
        self.assertEquals(sl1.isParameter("Sublist 2"), True)
        outputInfo = NOX.Parameter.Utils.Warning        + \
                     NOX.Parameter.Utils.OuterIteration + \
                     NOX.Parameter.Utils.InnerIteration + \
                     NOX.Parameter.Utils.Parameters     + \
                     NOX.Parameter.Utils.Details        + \
                     NOX.Parameter.Utils.OuterIterationStatusTest
        sl2.setParameter("New Parameter", outputInfo)
        self.assertEquals(sl2.getParameter("New Parameter"), outputInfo)

class StatusTestTestCase(unittest.TestCase):
    "TestCase class for NOX.StatusTest objects"

    def setUp(self):
        self.absResid  = NOX.StatusTest.NormF(1.0e-6)
        self.relResid  = NOX.StatusTest.NormF(1.0e-6,NOX.StatusTest.NormF.Unscaled)
        self.update    = NOX.StatusTest.NormUpdate(1.0e-5)
        self.wrms      = NOX.StatusTest.NormWRMS(1.0e-2, 1.0e-8)
        self.maxIters  = NOX.StatusTest.MaxIters(20)
        self.stagnate  = NOX.StatusTest.Stagnation()
        self.finiteVal = NOX.StatusTest.FiniteValue()

    def testAbsNormFGetStatus(self):
        "Test NOX.StatusTest.NormF getStatus method"
        self.assertEqual(self.absResid.getStatus(),  NOX.StatusTest.Unevaluated)

    def testRelNormFGetStatus(self):
        "Test NOX.StatusTest.NormF getStatus method"
        self.assertEqual(self.relResid.getStatus(),  NOX.StatusTest.Unevaluated)

    def testNormUpdateGetStatus(self):
        "Test NOX.StatusTest.NormUpdate getStatus method"
        self.assertEqual(self.update.getStatus(),    NOX.StatusTest.Unevaluated)

    def testNormWRMSGetStatus(self):
        "Test NOX.StatusTest.NormWRMS getStatus method"
        self.assertEqual(self.wrms.getStatus(),      NOX.StatusTest.Unconverged)

    def testMaxItersGetStatus(self):
        "Test NOX.StatusTest.MaxIters getStatus method"
        self.assertEqual(self.maxIters.getStatus(),  NOX.StatusTest.Unevaluated)

    def testStagnateGetStatus(self):
        "Test NOX.StatusTest.Stagnate getStatus method"
        self.assertEqual(self.stagnate.getStatus(),  NOX.StatusTest.Unevaluated)

    def testFiniteValueGetStatus(self):
        "Test NOX.StatusTest.FiniteValue getStatus method"
        self.assertEqual(self.finiteVal.getStatus(), NOX.StatusTest.Unevaluated)

    def testComboAnd(self):
        "Test NOX.StatusTest.Combo constructor with AND option"
        converged = NOX.StatusTest.Combo(NOX.StatusTest.Combo.AND)
        converged.addStatusTest(self.absResid)
        converged.addStatusTest(self.relResid)
        converged.addStatusTest(self.wrms    )
        converged.addStatusTest(self.update  )

    def testComboOr(self):
        "Test NOX.StatusTest.Combo constructor with OR option"
        combo = NOX.StatusTest.Combo(NOX.StatusTest.Combo.OR)
        combo.addStatusTest(self.maxIters)
        combo.addStatusTest(self.stagnate)
        combo.addStatusTest(self.finiteVal)


if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(ParameterTestCase ))
    suite.addTest(unittest.makeSuite(StatusTestTestCase))

    # Run the test suite
    unittest.TextTestRunner(verbosity=2).run(suite)
