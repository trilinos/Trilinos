#! /usr/bin/env python
# @HEADER
# *****************************************************************************
#                 Belos: Block Linear Solvers Package
#
# Copyright 2004-2016 NTESS and the Belos contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import os
import shlex
import subprocess
import re
import string
import shutil
import sys
import glob

################################################################################
# Run a command, given as a single string argument.  If the verbose flag is set,
# print the command.  Run the command and wait for the process to end.  The
# **kwargs argument can be used to pass any keyword arguments that are
# appropriate for the subprocess.Popen constructor.  If the command returns with
# a non-zero code, print a warning message.  Return the subprocess return code.
################################################################################
def runCommand(cmd, verbose=False, **kwargs):
    if verbose: print cmd
    args = shlex.split(cmd)
    p = subprocess.Popen(args, **kwargs)
    returncode = p.wait()
    if returncode:
        print "Command '%s' exited with code %s" % (args[0], str(returncode))
    return returncode

##############################

def main(buildDir):

    myFile = os.path.abspath(__file__)
    myPath = os.path.split(myFile)[0]

    # Get the Trilinos base directory that this script is in
    trilinosBasePath = os.path.dirname(
        os.path.dirname(os.path.dirname(myPath)))

    # Get all necessary paths
    package = os.path.basename(os.path.dirname(myPath))
    packageSrcDir = os.path.join(trilinosBasePath, 'packages', package)
    packageBuildDir = os.path.join(buildDir, 'packages', package)

    # Assumptions:
    # - source/packages/<pkg>/doc/html created by build_docs script
    # - XML files created in build/packages/<pkg>/doc/parameterList
    if os.path.isdir(os.path.join(packageSrcDir, 'doc', 'html')):
        xmlFilePath = os.path.join(packageBuildDir, 'doc', 'parameterList')
        if os.path.isdir(xmlFilePath):
            xmlFiles = glob.glob(os.path.join(xmlFilePath, '*.xml'))
            if not xmlFiles:
                print "ERROR: XML files not found in" + xmlFilePath
            else:
                for file in xmlFiles:
                    print "  Copying file " + file
                    shutil.copy(file, os.path.join(packageSrcDir, 'doc', 'html'))
        else:
            print "ERROR: " + xmlFilePath + " does not exist"
    else:
        print "ERROR: html directory not found in " + packageSrcDir + "/doc"
        print '  run the build_doc script first'

    # For debugging
    print "**********************************************************"
    print
    print "trilinosBasePath: " + trilinosBasePath
    print 'package source directory: ' + packageSrcDir + '/' + package
    print 'package build directory: ' + packageBuildDir + '/' + package
    print
    print "*********************************************************"


##############################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage:", sys.argv[0], "build_dir"
        sys.exit(1)
    build_dir = sys.argv[1]
    main(build_dir)
