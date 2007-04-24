#! /usr/bin/env python

# @header
# ************************************************************************
#
#                PyTrilinos: Python interface to Trilinos
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
# Questions? contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @header

# This example script builds Teuchos.ParameterList objects and functionally
# equivalent python dictionaries, demonstrating various capabilities of the
# Teuchos.ParameterList class.  These are then stored and read as XML files.

# Author: Bill Spotz, Sandia National Laboratories
# Date: 24 Apr 2007

from   optparse import *
import sys
from   types    import *

################################################################################

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

################################################################################

def main():

    # Create a ParameterList in the C++ way:
    pList = Teuchos.ParameterList()
    pList.set("maxiters", 100)
    pList.set("tol", 1.0e-6)
    pList.set("precond", "ILUT")
    print
    print pList

    # The isType() method is disabled.  In its place is a type() method that
    # returns a python type object:
    assert(pList.type("maxiters") is IntType   )
    assert(pList.type("tol"     ) is FloatType )
    assert(pList.type("precond" ) is StringType)

    # The python ParameterList has been augmented to support python dictionary
    # methods:
    keys = pList.keys()  # Returns a list of parameter names
    assert("maxiters" in keys)
    assert("tol"      in keys)
    assert("precond"  in keys)

    # Actually, the "in" operator works on ParameterLists, too:
    assert("maxiters" in pList)
    assert("tol"      in pList)
    assert("precond"  in pList)

    # These additional methods include special methods that enable iteration:
    print
    for name in pList:
        print name, ":", pList[name]

    # PyTrilinos as a whole has been designed to accept python dictionaries
    # where ParameterLists are expected.  For example, there is an overloaded
    # version of "set" that accepts a ParameterList.  Here we give it a
    # dictionary:
    pList.set("coefs", {"a":1.0, "b":-0.5, "c":0.0})
    assert(pList.type("coefs") is Teuchos.ParameterList)

    # To retrieve data from a ParameterList, you can use the get() method, or
    # the more dictionary-like square-brackets:
    m1 = pList.get("maxiters")
    m2 = pList["maxiters"]
    assert(m1 == m2)

    # An asDict() method has been added, to convert a ParameterList to a
    # dictionary:
    d = pList.asDict()
    assert(isinstance(d, dict))
    del d["maxiters"]    # Deletion is not supported by ParameterList
    pList2 = Teuchos.ParameterList(d)  # New constructor that accepts a dictionary
    print
    print pList2

    # The len() operator returns the number of top-level parameters:
    assert(len(pList) == len(pList2)+1)

    # ParameterList also supports comparison, and the square brackets can be
    # used for assignment, too:
    assert(pList != pList2)
    pList2["maxiters"] = 100
    assert(pList == pList2)

    # You can convert ParameterList objects and dictionaries to XMLObjects:
    writer = Teuchos.XMLParameterListWriter()
    xmlObj1 = writer.toXML(d)       # d is a dictionary
    xmlObj2 = writer.toXML(pList2)  # ParameterList argument

    # XMLObjects support the __str__() method
    assert(isinstance(xmlObj1, Teuchos.XMLObject))
    assert(isinstance(xmlObj2, Teuchos.XMLObject))
    print
    print xmlObj2

    # Write the XMLObject to disk, read it back and check it:
    open("params.xml","w").write(xmlObj2.toString())
    source  = Teuchos.FileInputSource("params.xml")
    xmlObj3 = source.getObject()
    reader  = Teuchos.XMLParameterListReader()
    pList3  = reader.toParameterList(xmlObj3)
    assert(pList3 == pList2)

################################################################################

if __name__ == "__main__":

    print "******************"
    print "** exTeuchos.py **"
    print "******************"

    try:
        main()
        print "End Result: TEST PASSED"

    except Exception, e:
        print e
