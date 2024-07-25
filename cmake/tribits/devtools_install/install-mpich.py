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

mpichBaseName = "mpich"
mpichDefaultVersion = "3.1.3"
mpichSupportedVersions = ["3.1.3"]
mpichTarballVersions = {
  "3.1.3" : "3.1.3"
  }

#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *


class MpichInstall:

  def __init__(self):
    self.dummy = None

  #
  # Called before even knowing the product version
  #

  def getScriptName(self):
    return "install-mpich.py"

  def getProductBaseName(self):
    return mpichBaseName

  def getProductDefaultVersion(self):
    return mpichDefaultVersion

  def getProductSupportedVersions(self):
    return mpichSupportedVersions

  #
  # Called after knowing the product version but before parsing the
  # command-line.
  #

  def getProductName(self, version):
    return mpichBaseName+"-"+version

  def getBaseDirName(self, version):
    return mpichBaseName+"-"+version+"-base"

  def getExtraHelpStr(self, version):
    return """
This script builds """+self.getProductName(version)+""" from source compiled
with the configured C/C++/Fortran compilers in your path or set by the env
vars CC, CXX, and FC.

NOTE: The assumed directory structure of the download source provided by the
command --download-cmnd=<download-cmnd> is:

   mpich-<version>-base/
     mpich-<version>.tar.gz
"""

  def injectExtraCmndLineOptions(self, clp, version):
    setStdDownloadCmndOption(self, clp, version)
    clp.add_option(
      "--extra-configure-options", dest="extraConfigureOptions", type="string", \
      default="", \
      help="Extra options to add to the 'configure' command for "+self.getProductName(version)+"." \
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
    self.mpichBaseDir = self.baseDir+"/"+self.getBaseDirName(self.inOptions.version)
    mpichVersionFull = mpichTarballVersions[self.inOptions.version]
    self.mpichTarball = "mpich-"+mpichVersionFull+".tar.gz"
    self.mpichSrcDir = "mpich-"+mpichVersionFull
    self.mpichBuildBaseDir = self.mpichBaseDir+"/mpich-build"
    self.scriptBaseDir = getScriptBaseDir()

  #
  # Called after setup()
  #

  def doDownload(self):
    removeDirIfExists(self.mpichBaseDir, True)
    echoRunSysCmnd(self.inOptions.downloadCmnd)

  def doUntar(self):
    # Find the full name of the source tarball
    echoChDir(self.mpichBaseDir)
    echoRunSysCmnd("tar -xzf "+self.mpichTarball)
    # NOTE: I found that you have to untar the tarball and can't store the
    # open source in the git repo.  Otherwise the timestamps are messed up and
    # it 'make' tries to recreate some generated files.

  def doConfigure(self):
    createDir(self.mpichBuildBaseDir, True, True)
    echoRunSysCmnd(
      "../"+self.mpichSrcDir+"/configure "+\
      " "+self.inOptions.extraConfigureOptions+\
      " --prefix="+self.inOptions.installDir,
      extraEnv={"CFLAGS":"-O3", "CXXFLAGS":"-O3", "FFLAGS":"-O3"},
      )

  def doBuild(self):
    echoChDir(self.mpichBuildBaseDir)
    echoRunSysCmnd("make " + getParallelOpt(self.inOptions, "-j") \
      + self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.mpichBuildBaseDir)
    echoRunSysCmnd("make " + getParallelOpt(self.inOptions, "-j") \
      + self.inOptions.makeOptions + " install")

  def getFinalInstructions(self):
    return """
To use the installed version of mpich-"""+self.inOptions.version+""" add the path:

  """+self.inOptions.installDir+"""/bin

to your path and that should be it!

Also, when you link shared libs or executables, pass in:

   -Wl,-rpath,"""+self.inOptions.installDir+"""/lib[64]

That will make it so that you don't need to add this MPICH libs to your
LD_LIBRARY_PATH.
"""


#
# Executable statements
#

mpichInstaller = InstallProgramDriver(MpichInstall())
mpichInstaller.runDriver()
