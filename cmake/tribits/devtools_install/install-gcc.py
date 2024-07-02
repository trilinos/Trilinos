#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


#
# Defaults
#

gccBaseName = "gcc"
gccDefaultVersion = "4.8.3"
gccSupportedVersions = ["4.8.3"]

#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *


class GccInstall:

  def __init__(self):
    self.dummy = None

  #
  # Called before even knowing the product version
  #

  def getScriptName(self):
    return "install-gcc.py"

  def getProductBaseName(self):
    return gccBaseName

  def getProductDefaultVersion(self):
    return gccDefaultVersion

  def getProductSupportedVersions(self):
    return gccSupportedVersions

  #
  # Called after knowing the product version but before parsing the
  # command-line.
  #

  def getProductName(self, version):
    return gccBaseName+"-"+version

  def getBaseDirName(self, version):
    return gccBaseName+"-"+version+"-base"

  def getExtraHelpStr(self, version):
    return """
This script builds """+self.getProductName(version)+""" from source compiled with the
configured C compilers in your path.

NOTE: The assumed directory structure of the download source provided by the
command --download-cmnd=<download-cmnd> is:

   gcc-<version>-base/
     gcc-<full-version>.tar.gz
"""

  def injectExtraCmndLineOptions(self, clp, version):
    setStdDownloadCmndOption(self, clp, version)
    clp.add_option(
      "--extra-configure-options", dest="extraConfigureOptions", type="string", \
      default="", \
      help="Extra options to add to the 'configure' command for " \
        + self.getProductName(version)+"." \
        +"  Note: This does not override the hard-coded configure options." )

  def echoExtraCmndLineOptions(self, inOptions):
    cmndLine = ""
    cmndLine += "  --download-cmnd='"+inOptions.downloadCmnd+"' \\\n"
    cmndLine += "  --extra-configure-options='"+inOptions.extraConfigureOptions+"' \\\n"
    return cmndLine

  #
  # Called after parsing the command-line
  #
    
  def setup(self, inOptions):
    self.inOptions = inOptions
    self.baseDir = os.getcwd()
    self.gccBaseDir = self.baseDir+"/"+self.getBaseDirName(self.inOptions.version)
    self.gccSrcDir = "gcc-"+self.inOptions.version
    self.gccBuildBaseDir = self.gccBaseDir+"/gcc-build"
    self.scriptBaseDir = getScriptBaseDir()

  #
  # Called after setup()
  #

  def doDownload(self):
    removeDirIfExists(self.gccBaseDir, True)
    echoRunSysCmnd(self.inOptions.downloadCmnd)

  def doUntar(self):
    print "Nothing to untar!"

  def doConfigure(self):
    createDir(self.gccBuildBaseDir)
    echoRunSysCmnd(
      "../"+self.gccSrcDir+"/configure --disable-multilib --enable-languages='c,c++,fortran'"+\
      " "+self.inOptions.extraConfigureOptions+\
      " --prefix="+self.inOptions.installDir,
      workingDir=self.gccBuildBaseDir,
      extraEnv={"CFLAGS":"-O3"},
      )

  def doBuild(self):
    echoChDir(self.gccBuildBaseDir)
    echoRunSysCmnd("make " + getParallelOpt(self.inOptions, "-j") \
      + self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.gccBuildBaseDir)
    echoRunSysCmnd("make " + getParallelOpt(self.inOptions, "-j") \
      + self.inOptions.makeOptions + " install")

  def getFinalInstructions(self):
    return """
To use the installed version of gcc-"""+self.inOptions.version+""" add the path:

  """+self.inOptions.installDir+"""/bin

to your path and that should be it!

Also, when you link shared libs or executables, pass in:

   -Wl,-rpath,"""+self.inOptions.installDir+"""/lib[64]

That will make it so that you don't need to add this GCC libs to your
LD_LIBRARY_PATH.
"""


#
# Executable statements
#

gccInstaller = InstallProgramDriver(GccInstall())
gccInstaller.runDriver()
