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
    import Teuchos
else:
    try:
        import setpath
        import Teuchos
    except ImportError:
        from PyTrilinos import Teuchos
        print >>sys.stderr, "Using system-installed Teuchos"

####################################################################

class TeuchosTestCase(unittest.TestCase):
    "TestCase class for Teuchos module"

    def testVersion(self):
        "Test Teuchos Teuchos_Version function"
        front   = "Teuchos Version "
        version = Teuchos.Teuchos_Version()
        self.assertEquals(version[:len(front)], front)

####################################################################

class TeuchosXMLTestCase(unittest.TestCase):
    "TestCase class for Teuchos module XML classes"

    def testXMLObjectConstructor0(self):
        "Test Teuchos XMLObject empty constructor"
        xmlObj = Teuchos.XMLObject()
        self.failUnless(xmlObj.isEmpty())

    def testXMLObjectConstructor1(self):
        "Test Teuchos XMLObject string constructor"
        tag    = "Parameters"
        xmlObj = Teuchos.XMLObject(tag)
        self.failIf(xmlObj.isEmpty())

    def testXMLObjectDeepCopy(self):
        "Test Teuchos XMLObject deepCopy method"
        tag    = "Parameters"
        xmlObj = Teuchos.XMLObject(tag)
        xmlCpy = xmlObj.deepCopy()
        self.assertEquals(xmlCpy.getTag(), xmlObj.getTag())

    def testXMLObjectDeepCopyEmpty(self):
        "Test Teuchos XMLObject deepCopy method, empty"
        xmlObj = Teuchos.XMLObject()
        self.assertRaises(RuntimeError, xmlObj.deepCopy)

    def testXMLObjectGetTag(self):
        "Test Teuchos XMLObject getTag method"
        tag    = "Parameters"
        xmlObj = Teuchos.XMLObject(tag)
        self.assertEquals(xmlObj.getTag(), tag)

    def testXMLObjectGetTagEmpty(self):
        "Test Teuchos XMLObject getTag method, empty"
        xmlObj = Teuchos.XMLObject()
        self.assertRaises(RuntimeError, xmlObj.getTag)

    def testXMLObjectAddAttributeBad(self):
        "Test Teuchos XMLObject addAttribute method, no tag"
        xmlObj = Teuchos.XMLObject()
        self.assertRaises(RuntimeError, xmlObj.addAttribute, "date", "today")

    def testXMLObjectAttribute(self):
        "Test Teuchos XMLObject add/get attribute methods"
        xmlObj = Teuchos.XMLObject("tag")
        xmlObj.addAttribute("date","today")
        self.assertEquals(xmlObj.getAttribute("date"), "today")

    def testXMLObjectBadAttribute(self):
        "Test Teuchos XMLObject getAttribute method, nonexistent"
        xmlObj = Teuchos.XMLObject()
        self.assertRaises(RuntimeError, xmlObj.getAttribute, "date")

    def testXMLObjectHasAttribute(self):
        "Test Teuchos XMLObject hasAttribute methods"
        xmlObj = Teuchos.XMLObject("tag")
        self.failIf(xmlObj.hasAttribute("date"))
        xmlObj.addAttribute("date","today")
        self.failUnless(xmlObj.hasAttribute("date"))

    def testXMLObjectGetRequired(self):
        "Test Teuchos XMLObject getRequired method"
        xmlObj = Teuchos.XMLObject("tag")
        self.assertRaises(RuntimeError, xmlObj.getRequired, "date")
        xmlObj.addAttribute("date","today")
        self.assertEquals(xmlObj.getRequired("date"), "today")

    def testXMLObjectGetRequiredDouble(self):
        "Test Teuchos XMLObject getRequiredDouble method"
        xmlObj = Teuchos.XMLObject("tag")
        self.assertRaises(RuntimeError, xmlObj.getRequiredDouble, "pi")
        xmlObj.addDouble("pi",3.14)
        self.assertEquals(xmlObj.getRequiredDouble("pi"), 3.14)

    def testXMLObjectGetRequiredInt(self):
        "Test Teuchos XMLObject getRequiredInt method"
        xmlObj = Teuchos.XMLObject("tag")
        self.assertRaises(RuntimeError, xmlObj.getRequiredInt, "year")
        xmlObj.addInt("year",2007)
        self.assertEquals(xmlObj.getRequiredInt("year"), 2007)

    def testXMLObjectGetRequiredBool(self):
        "Test Teuchos XMLObject getRequiredBool method"
        xmlObj = Teuchos.XMLObject("tag")
        self.assertEquals(xmlObj.getRequiredBool("flag"), False)
        xmlObj.addBool("flag",True)
        self.assertEquals(xmlObj.getRequiredBool("flag"), True)

    def testXMLObjectGetWithDefault(self):
        "Test Teuchos XMLObject getWithDefault method"
        xmlObj = Teuchos.XMLObject("tag")
        self.assertEquals(xmlObj.getWithDefault("country","USA"), "USA")
        xmlObj.addAttribute("country", "Canada")
        self.assertEquals(xmlObj.getWithDefault("country","USA"), "Canada")

    def testXMLObjectChild(self):
        "Test Teuchos XMLObject add/get child methods"
        p = Teuchos.XMLObject("parent")
        self.assertEquals(p.numChildren(), 0)
        #self.assertRaises(RuntimeError, p.getChild, 0)
        c = Teuchos.XMLObject("child")
        p.addChild(c)
        self.assertEquals(p.numChildren(), 1)
        r = p.getChild(0)
        self.assertEquals(r.getTag(), c.getTag())

    def testXMLObjectContent(self):
        "Test Teuchos XMLObject add/get content methods"
        poem = ["Twas brillig, and the slithy toves",
                "Did gyre and gimble in the wabe:"  ,
                "All mimsy were the borogoves,"     ,
                "And the mome raths outgrabe."       ]
        p = Teuchos.XMLObject("Jabberwocky")
        self.assertEquals(p.numContentLines(), 0)
        #self.assertRaises(RuntimeError, p.getContentLine, 0)
        for i in range(len(poem)):
            p.addContent(poem[i])
            self.assertEquals(p.numContentLines(), i+1)
        for i in range(len(poem)):
            self.assertEquals(p.getContentLine(i), poem[i])

    def testXMLObjectToStringEmpty(self):
        "Test Teuchos XMLObject toString method, empty"
        xmlObj = Teuchos.XMLObject()
        self.assertRaises(RuntimeError, xmlObj.toString)

    def testXMLObjectHeader(self):
        "Test Teuchos XMLObject header method"
        xmlObj = Teuchos.XMLObject("header")
        self.assertEquals(xmlObj.header(), "<header>")

    def testXMLObjectTerminatedHeader(self):
        "Test Teuchos XMLObject terminatedHeader method"
        xmlObj = Teuchos.XMLObject("header")
        self.assertEquals(xmlObj.terminatedHeader(), "<header/>")

    def testXMLObjectFooter(self):
        "Test Teuchos XMLObject footer method"
        xmlObj = Teuchos.XMLObject("footer")
        self.assertEquals(xmlObj.footer(), "</footer>")

    def testXMLObjectCheckTag(self):
        "Test Teuchos XMLObject checkTag method"
        xmlObj = Teuchos.XMLObject("tag")
        xmlObj.checkTag("tag")
        self.assertRaises(RuntimeError, xmlObj.checkTag, "junk")

    def testXMLParameterListReaderWriter(self):
        "Test Teuchos XMLParameterList Reader/Writer classes"
        dict = { "maxiters"       : 100,
                 "tol"            : 1.0e-6,
                 "preconditioner" : "IC" }
        writer = Teuchos.XMLParameterListWriter()
        reader = Teuchos.XMLParameterListReader()
        # Transform: dict -> XMLObject -> ParameterList
        xmlObj = writer.toXML(dict)
        pList  = reader.toParameterList(xmlObj)
        for parameter in dict:
            value = dict[parameter]
            pType = type(value)
            self.assertEquals(pType(pList[parameter]), value)

####################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(TeuchosTestCase))
    suite.addTest(unittest.makeSuite(TeuchosXMLTestCase))

    iAmRoot = True

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n***************\nTesting Teuchos\n***************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = len(result.errors) + len(result.failures)
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
