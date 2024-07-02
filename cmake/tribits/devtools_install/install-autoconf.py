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

autoconfBaseName = "autoconf"
autoconfDefaultVersion = "2.69"
autoconfSupportedVersions = ["2.69"]
autoconfTarballVersions = {
  "2.69" : "2.69"
  }

#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *


class AutoconfInstall:

  def __init__(self):
    self.dummy = None

  #
  # Called before even knowing the product version
  #

  def getScriptName(self):
    return "install-autoconf.py"

  def getProductBaseName(self):
    return autoconfBaseName

  def getProductDefaultVersion(self):
    return autoconfDefaultVersion

  def getProductSupportedVersions(self):
    return autoconfSupportedVersions

  #
  # Called after knowing the product version but before parsing the
  # command-line.
  #

  def getProductName(self, version):
    return autoconfBaseName+"-"+version

  def getBaseDirName(self, version):
    return autoconfBaseName+"-"+version+"-base"

  def getExtraHelpStr(self, version):
    return """
This script builds Autoconf"""+self.getProductName(version)+""" from source compiled with the
configured C compilers in your path.

NOTE: The assumed directory structure of the download source provided by the
command --download-cmnd=<download-cmnd> is:

   autoconf-<version>-base/
     autoconf-<full-version>.tar.gz
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
    self.autoconfBaseDir = self.baseDir+"/"+self.getBaseDirName(self.inOptions.version)
    autoconfVersionFull = autoconfTarballVersions[self.inOptions.version]
    self.autoconfTarball = "autoconf-"+autoconfVersionFull+".tar.gz"
    self.autoconfSrcDir = "autoconf-"+autoconfVersionFull
    self.autoconfBuildBaseDir = self.autoconfBaseDir+"/autoconf-build"
    self.scriptBaseDir = getScriptBaseDir()

  #
  # Called after setup()
  #

  def doDownload(self):
    removeDirIfExists(self.autoconfBaseDir, True)
    echoRunSysCmnd(self.inOptions.downloadCmnd)

  def doUntar(self):
    # Find the full name of the source tarball
    echoChDir(self.autoconfBaseDir)
    echoRunSysCmnd("tar -xzf "+self.autoconfTarball)

  def doConfigure(self):
    createDir(self.autoconfBuildBaseDir, True, True)
    echoRunSysCmnd(
      "../"+self.autoconfSrcDir+"/configure "+\
      " "+self.inOptions.extraConfigureOptions+\
      " --prefix="+self.inOptions.installDir,
      extraEnv={"CFLAGS":"-O3"},
      )

  def doBuild(self):
    echoChDir(self.autoconfBuildBaseDir)
    echoRunSysCmnd("make " + getParallelOpt(self.inOptions, "-j") \
      + self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.autoconfBuildBaseDir)
    echoRunSysCmnd("make " + getParallelOpt(self.inOptions, "-j") \
      + self.inOptions.makeOptions + " install")

  def getFinalInstructions(self):
    return """
To use the installed version of autoconf-"""+self.inOptions.version+""" add the path:

  """+self.inOptions.installDir+"""/bin

to your path and that should be it!
"""


#
# Executable statements
#

autoconfInstaller = InstallProgramDriver(AutoconfInstall())
autoconfInstaller.runDriver()
