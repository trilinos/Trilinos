#! /usr/bin/env python

import os
import shlex
import subprocess
import re
import string
import shutil

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

def main():

    myFile = os.path.abspath(__file__)
    myPath = os.path.split(myFile)[0]

    # Get the Trilinos base directory that this script is in
    trilinosBasePath = os.path.dirname(os.path.dirname(os.path.dirname(myPath)))

    # Get the package name
    package = os.path.basename(os.path.dirname(myPath))
    print "package: " + package

    # Detemine if in source or build tree
    inBuild = 'CMakeCache.txt' in os.listdir(trilinosBasePath)

    # Set TRILINOS_HOME
    if os.environ.has_key('TRILINOS_HOME'):
        print "TRILINOS_HOME has already been set!"
    else:
        if inBuild:
            cacheFile = open(os.path.join(trilinosBasePath, 'CMakeCache.txt'), 'r')
            for line in cacheFile:
                if re.match('Trilinos_SOURCE_DIR', line):
                    trilinosSrcDir = line.split("=")[1]
                    trilinosSrcDir = trilinosSrcDir.replace('\n', '')
                    os.environ['TRILINOS_HOME'] = trilinosSrcDir
                    break
        else:
            print "TRILINOS_HOME has not been set.  Setting it!"
            os.environ['TRILINOS_HOME'] = trilinosBasePath
    trilinosHome = os.environ['TRILINOS_HOME']

    # Prepare to run Doxygen in Trilinos source directory
    #TODO: Devise long term strategy for building Doxygen files in the
    # build tree. For now we will build in the source directory and
    # copy doxygen html subtree to the build directory if building 
    # docs in the build directory.
    newDir = os.path.join(trilinosHome, 'packages', package, 'doc')
    os.chdir(newDir)

    # Remove html directory from source/<package>/doc directory, if exists
    htmlPackage = os.path.join(newDir, package)
    if os.path.isdir(htmlPackage):
        shutil.rmtree(htmlPackage)

    # Build doxygen browser documentation
    print
    print "Building doxygen browser documentation for all of stratimikos " + \
          "as a single doxygen collection ..."
    print
    doxyfile = os.path.join('..','browser','doc','Doxyfile')
    runCommand('doxygen ' + doxyfile, True)

    # Build main doxygen documentation
    print
    print "Building main doxygen documentation page for stratimikos in " + \
        newDir + "..."
    print
    runCommand('doxygen Doxyfile', True)

    # Copy generated html directory to Trlinos_BINARY_DIR/html/<pkg>
    if inBuild:
        #os.rename('html', package)
        htmlBuildDir = os.path.join(trilinosBasePath, 'html') 
        htmlPackageDir = os.path.join(htmlBuildDir, package) 
        if os.path.isdir(htmlPackageDir):
            shutil.rmtree(htmlPackageDir)
        shutil.move(package, htmlBuildDir)
        packageXml = os.path.join(myPath, "parameterList", package + ".xml")
        shutil.copy(packageXml, htmlPackageDir)

    # For debugging
    print "**********************************************************"
    print
    print "trilinosBasePath: " + trilinosBasePath
    if inBuild:
        print "inBuild: True"
    else:
        print "inBuild: False"

    print "TRILINOS_HOME: " + trilinosHome
    print "newDir: " + newDir
    print
    print "*********************************************************"


##############################

if __name__ == "__main__":
    main()
