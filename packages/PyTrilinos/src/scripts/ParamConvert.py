#! /usr/bin/env python
# -*- python -*-
"""
ParamConvert.py - A script for converting between XML and python representations
                  of Teuchos ParameterList objects.  The "from" and "to" formats
                  are determined by the extensions of the input and output
                  files, respectively.  Valid extensions are ".py" for python
                  and ".xml" for XML.  If the output file name is omitted, it is
                  built from the base of the input file name and the opposite
                  extension.
"""

__version__ = "1.0"
__author__  = "Bill Spotz"
__date__    = "Apr 23 2007"

# System imports
from   optparse import *
import os.path
from   pprint   import pprint
import sys

# Trilinos import
from PyTrilinos import Teuchos

# Valid file extensions
PYTHON_EXT = ".py"
XML_EXT    = ".xml"

################################################################################

def validExtension(ext):
    return (ext in (PYTHON_EXT, XML_EXT))

################################################################################

def outputExtension(ext):
    if ext == PYTHON_EXT: return XML_EXT
    if ext == XML_EXT   : return PYTHON_EXT
    raise ValueError, "Unsupported file extension"

################################################################################

def main():

    # Initialization
    prog = os.path.split(sys.argv[0])[1]

    # Set up the command-line parser object
    usage   = prog + " [options] input_file [output_file]\n" + __doc__
    version = prog + " " + __version__ + " " + __date__
    parser  = OptionParser(usage=usage,version=version)
    parser.set_defaults(verbose=False)
    parser.add_option("-v", "--verbose", action="store_true" , dest="verbose",
                      help="run in verbose mode")
    parser.add_option("-q", "--quiet"  , action="store_false", dest="verbose",
                      help="run in quiet mode [default]")

    # Get the options and arguments from the command line and check for valid
    # number of arguments
    (options,args) = parser.parse_args()
    if len(args) == 0 or len(args) > 2:
        print "usage:", usage
        sys.exit(1)

    # Validate the input file name
    inputFile = args[0]
    (inputRoot,inputExt) = os.path.splitext(inputFile)
    if not validExtension(inputExt):
        print "Unrecognized input file extension:", inputExt
        sys.exit(2)

    # Build and/or validate the output file name
    if len(args) == 1:
        outputRoot = inputRoot
        outputExt  = outputExtension(inputExt)
        outputFile = outputRoot + outputExt
    else:
        outputFile = args[1]
        (outputRoot,outputExt) = os.path.splitext(outputFile)
        if not validateExtension(outputExt):
            print "Unrecognized output file extension:", outputExt
            sys.exit(3)

    # Check that a conversion is requested
    if inputExt == outputExt:
        print "Input and output files have same extension:", inputExt
        sys.exit(4)

    # Convert python input to XML
    if inputExt == PYTHON_EXT:
        try:
            pyObj  = eval(open(inputFile).read().strip())
            xmlObj = Teuchos.XMLParameterListWriter().toXML(pyObj)
            open(outputFile,"w").write(xmlObj.toString())
            if options.verbose: print "Wrote '%s'" % outputFile
        except Exception, e:
            print "Could not convert python file to XML file:"
            print e
            sys.exit(5)

    # Convert XML input to python
    if inputExt == XML_EXT:
        try:
            xmlObj = Teuchos.FileInputSource(inputFile).getObject()
            pList  = Teuchos.XMLParameterListReader().toParameterList(xmlObj)
            stream = open(outputFile,"w")
            pprint(pList.asDict(), stream)   # Pretty-print the dictionary
            stream.close()
            if options.verbose: print "Wrote '%s'" % outputFile
        except Exception, e:
            print "Could not convert XML file to python file:"
            print e
            sys.exit(6)

################################################################################

if __name__ == "__main__":
    main()
