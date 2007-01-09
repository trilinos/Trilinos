# @HEADER
# ************************************************************************
# 
#              PyTrilinos: Python Interface to Trilinos
#                 Copyright (2004) Sandia Corporation
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

# System imports
import commands
import glob
import os
import sys

# Local import
from MakefileVariables import *

######################################################################

def buildSharedLibraries():
    """
    Return boolean indicating whether the shared libraries need to be built.
    Returns True if the operating system is either Darwin or Linux and this is
    an MPI build.
    """
    #if os.uname()[0] not in ("Darwin", "Linux"): return False
    #makefileVars = processMakefile("Makefile")
    #return bool(makefileVars["HAVE_MPI"])
    return os.uname()[0] in ("Darwin","Linux")

######################################################################

def needsToBeBuilt(target,dependencies):
    """
    Return True if a target file does not exist or is older than any of the
    files in the list of dependencies.
    """
    if not os.path.isfile(target): return True
    targetTime = os.path.getmtime(target)
    for dep in dependencies:
        if os.path.getmtime(dep) > targetTime: return True
    return False

######################################################################

def changeDirectory(path):
    """
    Change directory verbosely
    """
    print "changing directory to", path
    os.chdir(path)

######################################################################

def runCommand(cmd):
    """
    Print a Unix command, run it, print its output and raise an exception if its
    return status is non-zero.
    """
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if output: print output
    if status: raise RuntimeError, "Exit status = %d" % status

######################################################################

def deleteFile(path):
    """
    If a file exists, delete it verbosely; if not, do nothing verbosely.
    """
    if os.path.isfile(path):
        print "deleting", path
        os.remove(path)
    else:
        print "nothing needs to be done for", path

######################################################################

def removeExtensionModules(path):
    """
    Remove every file whose name matches the pattern _*.so under the given
    path.
    """
    if True: return
    print "removing extension modules"
    runCommand("find %s -name _*.so -print -delete" % path)

######################################################################

class SharedTrilinosBuilder:
    """
    Class for performing various build processes (build, install, clean,
    uninstall) related to converting a Trilinos library from static to dynamic.
    """

    def __init__(self,package):
        """
        Initialize a SharedTrilinosBuilder object for the specified Trilinos
        package.  The package name should be in its capitalized form.
        """
        self.__package      = package
        self.__packageLower = package.lower()
        self.__packageUpper = package.upper()
        self.__sysName      = os.uname()[0]
        self.__libPathVar   = self.getLibPathVarName()
        self.__libOption    = "-l"  + self.__packageLower
        self.__dylibName    = self.getDylibName()
        if self.__packageLower == "loca":
            self.__topBuildDir  = os.path.join("..", "..", "nox")
            self.__buildDir     = os.path.join(self.__topBuildDir, "src-loca",
                                               "src")
        elif self.__packageLower == "noxepetra":
            self.__topBuildDir  = os.path.join("..", "..", "nox")
            self.__buildDir     = os.path.join(self.__topBuildDir, "src-epetra")
        else:
            self.__topBuildDir  = os.path.join("..", "..", self.__packageLower)
            self.__buildDir     = os.path.join(self.__topBuildDir, "src")
        self.__pythonDir    = os.path.join(self.__topBuildDir, "python", "src")
        self.__thisDir      = os.getcwd()
        self.__absDylibName = os.path.join(self.__thisDir, self.__dylibName)
        self.__makeMacros   = self.getMakeMacros()
        self.__linkCmd      = self.getLinkCmd()
        self.__supported    = buildSharedLibraries()

    ######################################################################

    def getLibPathVarName(self):
        """
        Determine the name of environment variable used to specify the (dynamic)
        library load path.
        """
        if self.__sysName in ("Darwin",):
            return "DYLD_LIBRARY_PATH"
        else:
            return "LD_LIBRARY_PATH"

    ######################################################################

    def getDylibName(self):
        """
        Determine what the dynamic library name should be, based on the type of
        system we are on.
        """
        if self.__sysName in ("Darwin",):
            return "lib" + self.__packageLower + ".dylib"
        elif self.__sysName in ("Linux",):
            return "lib" + self.__packageLower + ".so"
        else:
            return "lib" + self.__packageLower + ".?"

    ######################################################################

    def getMakeMacros(self):
        """
        Obtain the dictionary of Makefile macros appropriate for PyTrilinos.
        """
        return processMakefile("Makefile")

    ######################################################################

    def getLinkCmd(self):
        """
        Determine the appropriate dynamic link command for this Trilinos
        package.
        """
        cxx = self.__makeMacros["CXX"]
        if self.__packageLower == "noxepetra":
            libVar = "NOX_LIBS"
        else:
            libVar = self.__packageUpper + "_LIBS"
        default = self.__makeMacros[libVar]
        ldFlags = self.__makeMacros.get(self.__packageUpper + "_PYTHON_LIBS",
                                        default)
        ldFlags = ldFlags.replace(self.__libOption+" ", "")  # Remove -lpackage
        ldFlags = "-L" + self.__thisDir + " " + ldFlags
        if self.__sysName == "Darwin":
            linkCmd = "%s -dynamiclib -o %s *.o -single_module %s" % (cxx,
                                                                      self.__dylibName,
                                                                      ldFlags)
        elif self.__sysName == "Linux":
            linkCmd = "%s -shared -Wl,-soname,%s -o %s *.o %s" % (cxx,
                                                                  self.__dylibName,
                                                                  self.__dylibName,
                                                                  ldFlags)
        else:
            linkCmd = "echo %s not supported" % self.__sysName
        return linkCmd

    ######################################################################

    def buildShared(self):
        """
        Change directory to the Trilinos package and, if necessary, re-link the
        package library as a dynamic/shared library.
        """
        if not self.__supported: return

        # Build a new shared library
        if needsToBeBuilt(self.__absDylibName, glob.glob(os.path.join(self.__buildDir,"*.o"))):
            print "\nConverting", self.__package, "to a shared library"
            changeDirectory(self.__buildDir)
            runCommand(self.__linkCmd)
            print "Moving", self.__dylibName, "to", self.__absDylibName
            os.rename(self.__dylibName, self.__absDylibName)
            # Restore the current working directory
            changeDirectory(self.__thisDir)

    ######################################################################

    def clean(self):
        """
        Remove all files created by the buildShared() method.
        """
        if not self.__supported: return

        # Remove the shared library
        deleteFile(self.__dylibName)

    ######################################################################

    def install(self):
        """
        Install the new shared Trilinos library and the newly linked python
        extension module in the appropriate directories.
        """
        if not self.__supported: return

        # Initialize
        install    = self.__makeMacros["INSTALL"]
        prefix     = self.__makeMacros["prefix" ]
        installDir = os.path.join(prefix, "lib")

        # Install the shared library and python module
        runCommand(" ".join([install, self.__dylibName, installDir]))

    ######################################################################

    def uninstall(self):
        """
        Uninstall the shared Trilinos library.
        """
        if not self.__supported: return

        # Initialize
        prefix     = self.__makeMacros["prefix"]
        installDir = os.path.join(prefix, "lib")

        # Uninstall the shared library and python module
        deleteFile(os.path.join(installDir, self.__dylibName))
