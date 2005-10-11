#! /usr/bin/env python

"""
CreatePythonDir.py is used to create a directory suitable for maintaining the
wrappers necessary to generate a python interface to your Trilinos package.
Correct usage is

    CreatePythonDir.py [options] [your_package]

where the options are

    -c Capitalized_Name
    --capital=Capitalized_Name
            Specify the capitalized version of your package name

    -l lower_name
    --lower=lower_name
            Specify the lower case version of your package name

    -u UPPER_NAME
    --upper=UPPER_NAME
            Specify the upper case version of your package name

    -h
    --help
            Print this message and exit

    -v
    --verbose
            Print verbose messages while executing (default False)

    -V
    --version
            Print the version number and exit

The argument "your_package" can also be used to specify your package name.  The
script will determine if this name is best described as lower case, upper case,
or capitalized, and determine the other name versions appropriately, unless
these other name versions are specified via the command line options.

CreatePythonDir.py will create a directory named <lower_name>_python that mimics
the new_package/python directory with the exception that all of the references
to "New_Package" in its various forms will have been translated to
"Your_Package".  You can then move <lower_name>_python into your source
directory, renaming it "python" as the first step in adding a python interface
to your package.
"""

__version__ = "1.0"
__author__  = "Bill Spotz"
__date__    = "Oct 10 2005"

# System imports
from   getopt  import *
import os
import os.path
import sys

# Data
filesToProcess = ["Makefile.am",
                  os.path.join("example", "Makefile.am"         ),
                  os.path.join("example", "setpath.py"          ),
                  os.path.join("src",     "Makefile.am"         ),
                  os.path.join("src",     "New_Package.i"       ),
                  os.path.join("src",     "TRILINOS_HOME_DIR.in"),
                  os.path.join("src",     "__init__.py.in"      ),
                  os.path.join("src",     "setup.py"            ),
                  os.path.join("test",    "Makefile.am"         ),
                  os.path.join("test",    "setpath.py"          ),
                  os.path.join("test",    "testNew_Package.py"  ) ]
oldCapital     = "New_Package"
oldLower       = "new_package"
oldUpper       = "NEW_PACKAGE"

def substitute(lines, old, new):
    for i in range(len(lines)):
        lines[i] = lines[i].replace(old,new)

def processFile(fileName,oldDir,newDir,newCapital,newLower,newUpper):
    oldPath = os.path.join(oldDir,fileName)
    newPath = os.path.join(newDir,fileName.replace(oldCapital,newCapital))
    lines   = file(oldPath,"r").readlines()
    substitute(lines, oldCapital, newCapital)
    substitute(lines, oldLower,   newLower  )
    substitute(lines, oldUpper,   newUpper  )
    file(newPath,"w").writelines(lines)

def main():

    # Initialization
    (progDir,progName) = os.path.split(sys.argv[0])
    options      = "c:hl:u:vV"
    long_options = ["capital=", "help", "lower=", "upper=", "verbose", "version"]
    newCapital   = None
    newLower     = None
    newUpper     = None
    verbose      = False

    # Get the options and arguemnts from the command line
    (opts,args) = getopt(sys.argv[1:], options, long_options)

    # Loop over options and implement
    for flag in opts:
        if flag[0] in ("-h","--help"):
            print __doc__
            sys.exit()
        elif flag[0] in ("-V", "--version"):
            print progName, __version__, __date__
            sys.exit()
        elif flag[0] in ("-v", "--verbose"):
            verbose = True
        elif flag[0] in ("-c", "--capital"):
            newCapital = flag[1]
        elif flag[0] in ("-l", "--lower"):
            newLower = flag[1]
        elif flag[0] in ("-u", "--upper"):
            newUpper = flag[1]
        else:
            print "Unrecognized flag:", flag[0]
            print __doc__
            sys.exit()

    # Check the arguments
    if len(args) > 1:
        print "Only one package name, please."
        sys.exit()
    if len(args) == 1:
        if   args[0].islower(): newLower   = args[0]
        elif args[0].isupper(): newUpper   = args[0]
        else:                   newCapital = args[0]

    # Fill in the blanks
    if newCapital:
        if not newLower:   newLower   = newCapital.lower()
        if not newUpper:   newUpper   = newCapital.upper()
    elif newLower:
        if not newCapital: newCapital = newLower.capitalize()
        if not newUpper:   newUpper   = newLower.upper()
    elif newUpper:
        if not newLower:   newLower   = newUpper.lower()
        if not newCapital: newCapital = newUpper.capitalize()
    else:
        print __doc__
        sys.exit()

    # Print package names
    if verbose:
        print "Capitalized package name =", newCapital
        print "Lower case  package name =", newLower
        print "Upper case  package name =", newUpper

    # Directories
    oldDir  = "python"
    newDir  = newLower + "_" + oldDir
    newDirs = [ ]
    for file in filesToProcess:
        dirname = os.path.join(newDir, os.path.dirname(file))
        if not dirname in newDirs: newDirs.append(dirname)
    for dir in newDirs:
        if not os.path.isdir(dir):
            if verbose: print "Making", dir
            os.makedirs(dir)

    # Process the files
    for file in filesToProcess:
        if verbose: print "Processing", file
        processFile(file, oldDir, newDir, newCapital, newLower, newUpper)

if __name__ == "__main__":
    main()
