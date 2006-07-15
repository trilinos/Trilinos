#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#             PyTrilinos.Teuchos: Python Interface to Teuchos
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
    import Teuchos
else:
    try:
        import setpath
        import Teuchos
    except ImportError:
        from PyTrilinos import Teuchos
        print >>sys.stderr, "Using system-installed Teuchos"

####################################################################

class ParameterListTestCase(unittest.TestCase):
    "TestCase class for Teuchos.ParameterList"

    def setUp(self):
        self.plist = Teuchos.ParameterList()
        self.name  = "Solver Params"

    def testConstructor0(self):
        "Test Teuchos.ParameterList default constructor"
        self.assertEqual(isinstance(self.plist,Teuchos.ParameterList), True)

    def testConstructor1(self):
        "Test Teuchos.ParameterList string constructor"
        plist = Teuchos.ParameterList(self.name)
        self.assertEqual(isinstance(plist,Teuchos.ParameterList), True)
        self.assertEqual(plist.name(), self.name)

    def testConstructor2(self):
        "Test Teuchos.ParameterList copy constructor"
        plist_copy = Teuchos.ParameterList(self.plist)
        self.assertEqual(isinstance(plist_copy,Teuchos.ParameterList), True)

    def testSetName(self):
        "Test Teuchos.ParameterList name and setName methods"
        self.assertEqual(self.plist.name(), "ANONYMOUS")
        self.plist.setName(self.name)
        self.assertEqual(self.plist.name(), self.name)


    def testSetParameters(self):
        "Test Teuchos.ParameterList setParameters method"
        intName    = "int parameter"
        intValue   = 8
        floatName  = "float parameter"
        floatValue = 3.14
        self.plist.set(intName,  intValue  )
        self.plist.set(floatName,floatValue)
        newList = Teuchos.ParameterList()
        newList.setParameters(self.plist)
        self.assertEqual(newList.get(intName  ), intValue  )
        self.assertEqual(newList.get(floatName), floatValue)

    def testSetInt(self):
        "Test Teuchos.ParameterList set and get methods for an integer"
        name  = "int parameter"
        value = 12
        self.plist.set(name, value)
        self.assertEqual(self.plist.get(name), value)

    def testSetFloat(self):
        "Test Teuchos.ParameterList set and get methods for a float"
        name  = "float parameter"
        value = 12.0
        self.plist.set(name, value)
        self.assertEqual(self.plist.get(name), value)

    def testSetString(self):
        "Test Teuchos.ParameterList set and get methods for a string"
        name  = "string parameter"
        value = "12"
        self.plist.set(name, value)
        self.assertEqual(self.plist.get(name), value)

#     def testSetDictionary(self):
#         "Test Teuchos.ParameterList set method for a dictionary"
#         name = "dict parameter"
#         d    = {"a" : 1,
#                 "b" : 2,
#                 "c" : 3 }
#         print "\nrefcount(d) =", sys.getrefcount(d)
#         self.plist.set(name, d)
#         print "refcount(d) =", sys.getrefcount(d)
#         print "list ="
#         self.plist._print()
#         sublist = self.plist.get(name)
#         print "sublist ="
#         sublist._print()
#         self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
#         for key in d.keys():
#             self.assertEqual(sublist.get(key), d[key])

    def testSetParameterList(self):
        "Test Teuchos.ParameterList set and get methods for a ParameterList"
        name    = "sublist"
        sublist = Teuchos.ParameterList()
        self.plist.set(name, sublist)
        self.assertEqual(isinstance(self.plist.get(name),
                                    Teuchos.ParameterList), True)
        self.assertEqual(self.plist.isSublist(name), True)

    def testSetNone(self):
        "Test Teuchos.ParameterList set method for None"
        self.assertRaises(TypeError, self.plist.set, "bad parameter", None)

    def testSetBadType(self):
        "Test Teuchos.ParameterList set method for an unsupported type"
        self.assertRaises(TypeError, self.plist.set, "bad parameter", (1,2,3))

    def testGetDefault(self):
        "Test Teuchos.ParameterList get method using default value"
        default = None
        self.assertEqual(self.plist.get("junk",None), None)

    def testSublistNew(self):
        "Test Teuchos.ParameterList sublist method for new sublist"
        sublist = self.plist.sublist("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
        sublist = self.plist.get("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)

    def testSublistOld(self):
        "Test Teuchos.ParameterList sublist method for existing sublist"
        sublist = self.plist.sublist("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
        sublist = self.plist.sublist("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)

    def testSublistBad(self):
        "Test Teuchos.ParameterList sublist method for non-sublist"
        self.plist.set("new", 1)
        self.assertRaises(RuntimeError, self.plist.sublist, "new")

    def testIsParameterTrue(self):
        "Test Teuchos.ParameterList isParameter method existing parameter"
        name = "string parameter"
        self.plist.set(name,"Hello")
        self.assertEqual(self.plist.isParameter(name), True)

    def testIsParameterFalse(self):
        "Test Teuchos.ParameterList isParameter method nonexisting parameter"
        name = "parameter"
        self.assertEqual(self.plist.isParameter(name), False)

    def testIsSublistTrue(self):
        "Test Teuchos.ParameterList isSublist method for existing sublist"
        name = "sublist"
        self.plist.sublist(name)
        self.assertEqual(self.plist.isSublist(name), True)

    def testIsSublistFalse1(self):
        "Test Teuchos.ParameterList isSublist method for existing non-sublist parameter"
        name = "string parameter"
        self.plist.set(name,"Hello")
        self.assertEqual(self.plist.isSublist(name), False)

    def testIsSublistFalse2(self):
        "Test Teuchos.ParameterList isSublist method for nonexisting parameter"
        name = "parameter"
        self.assertEqual(self.plist.isSublist(name), False)

    def testPrint0(self):
        "Test Teuchos.ParameterList _print method for empty list"
        fName = "testParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f)
        f.close()
        self.assertEqual(open(fName,"r").read(), "[empty list]\n")

    def testPrint1(self):
        "Test Teuchos.ParameterList _print method for non-empty list"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        fName = "testParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i], "%s = %g   [unused]\n" % (names[i], values[i]))

    def testPrint2(self):
        "Test Teuchos.ParameterList _print method for non-empty list and indentation"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        fName = "testParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f,2)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i], "  %s = %g   [unused]\n" % (names[i], values[i]))

    def testPrint3(self):
        "Test Teuchos.ParameterList _print method for non-empty list, indentation and types"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        types  = ["int"    ,"double"   ]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        fName = "testParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f,2,True)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i], "  %s : %s = %g   [unused]\n" %
                             (names[i], types[i], values[i]))

    def testUnused(self):
        "Test Teuchos.ParameterList unused method"
        names  = ["s1"   , "s2"   , "s3"         ]
        values = ["Hello", "World", "Albuquerque"]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        fName = "testParameterList.dat"
        f = open(fName,"w")
        self.plist.unused(f)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i],
                             'WARNING: Parameter "%s" %s   [unused] is unused\n' %
                             (names[i], values[i]))

    def testCurrentParametersString(self):
        "Test Teuchos.ParameterList currentParametersString method"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        types  = ["int"    ,"double"   ]
        result = "{"
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
            result += '"%s":%s=%g,' % (names[i], types[i], values[i])
        result = result[:-1] + "}"
        self.assertEqual(self.plist.currentParametersString(), result)

    def testType(self):
        "Test Teuchos.ParameterList type method"
        sublist = Teuchos.ParameterList()
        names  = ["iParm", "fParm", "sParm", "lParm"              ]
        values = [2006   , 2.71828, "Hello", sublist              ]
        types  = [int    , float  , str    , Teuchos.ParameterList]
        for i in range(len(names)):
            self.plist.set(names[i],values[i])
        for i in range(len(names)):
            self.assertEqual(self.plist.type(names[i]), types[i])

####################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(ParameterListTestCase))

    # Create a communicator
    #comm    = Epetra.PyComm()
    #iAmRoot = comm.MyPID() == 0
    iAmRoot = True

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n*****************************\n" + \
       "Testing Teuchos.ParameterList\n" + \
       "*****************************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    #errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    errsPlusFails = len(result.errors) + len(result.failures)
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
