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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

"""
pkg_info [options] filename

This pkg_info script takes a SWIG interface filename as its command-line
argument, and outputs requested data about the package name.  It does this by
searching the interface file for a %module directive.  If that directive is
found then the name of the package is extracted, including anything specified by
a "package=..." option.  Output depends on the option given.

This script is useful for the SWIG "-outdir" option.  For example, if swig is
invoked on "file.i" via

    swig ... -outdir `pkg_info -o file.i` ... file.i

then the output will be placed in the directory specified by its package,
assuming that directory exists.
"""

__version__ = "1.1"
__author__  = "Bill Spotz"
__date__    = "Dec 9 2006"

# System imports
from   optparse import *
import os.path
import re
import sys

# Local imports
from MakefileVariables import *

################################################################################

def extractPackageNames(filename):

    """
    extractPackageNames(str) -> list of strings
    
    Read the contents of the given filename for a SWIG interface file.  If it
    has a %module directive, parse it to obtain the module name, including
    anything specified by the "package=..." option.  Return the results as a
    list of strings.
    """

    # Regular expression object
    moduleRE = re.compile(r"%module(\([^\)]*\))\s*(\S+)",re.DOTALL)

    # Search for first occurence of the regular expression within the given
    # filename given on the command line
    match = moduleRE.search(open(filename,"r").read())

    # Default package is empty list
    package = [ ]
    # Check for match
    if match:

        # Check for module options
        if match.group(1):
            modOptions = match.group(1)[1:-1].split(",")
            for opt in modOptions:
                name,val = opt.split("=")
                name = name.strip()
                val  = val.strip()[1:-1]
                if name == "package": package = val.split(".")

        # Append module name
        package.append(match.group(2))

    return package

################################################################################

def main():
    usage   = __doc__
    version = "%prog " + __version__ + " " + __date__
    parser  = OptionParser(usage=usage,version=version)
    parser.set_defaults(verbose = False,
                        outdir  = False,
                        name    = False,
                        upper   = False,
                        lower   = False,
                        include = False)
    #parser.add_option("-v", "--verbose", action="store_true" , dest="verbose",
    #                  help="run in verbose mode")
    #parser.add_option("-q", "--quiet"  , action="store_false", dest="verbose",
    #                  help="run in quiet mode [default]")
    parser.add_option("-o", "--outdir" , action="store_true" , dest="outdir" ,
                      help="Output the package in the form of a directory")
    parser.add_option("-n", "--name"   , action="store_true" , dest="name"   ,
                      help="Output the package capitalized name")
    parser.add_option("-u", "--upper"  , action="store_true" , dest="upper"  ,
                      help="Output the package uppercase name")
    parser.add_option("-l", "--lower"  , action="store_true" , dest="lower"  ,
                      help="Output the package lowercase name")
    parser.add_option("-i", "--include", action="store_true" , dest="include",
                      help="Output the package include files list")

    # Get the options and arguments from the command line
    (options,args) = parser.parse_args()

    # Initialization
    pyt = "PyTrilinos"

    # Search for first occurence of the regular expression within the given
    # filename given on the command line
    filename = args[0]
    package = extractPackageNames(filename)

    # Output for -o, --outdir option
    if options.outdir:
        if package:
            print os.path.join(*package[:-1])
        else:
            print "."

    # Output for -n, --name option
    elif options.name:
        if package:
            if pyt in package: package.remove(pyt)
            print package[0]

    # Output for -u, --upper option
    elif options.name:
        if package:
            if pyt in package: package.remove(pyt)
            print package[0].upper()

    # Output for -l, --lower option
    elif options.name:
        if package:
            if pyt in package: package.remove(pyt)
            print package[0].lower()

    # Output for -i, --include
    elif options.include:
        if package:
            if pyt in package: package.remove(pyt)
            pkgName = package[0].lower()
            exportMakefile = os.path.join("..","..",pkgName,
                                          "Makefile.export."+pkgName)
            if os.path.isfile(exportMakefile):
                vars = processMakefile(exportMakefile)
                name = pkgName.upper()+"_INCLUDES"
                if vars.has_key(name):
                    print vars[name]

    # Default output
    else:
        print filename,
        if package:
            if package[-1] == "__init__":
                print "defines module", ".".join(package[:-1])
            else:
                print "defines module", ".".join(package)
        else:
            print "does not define a python module"

################################################################################

if __name__ == "__main__":
    main()
