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

gitBaseName = "git"
gitDefaultVersion = "2.6.4"
gitSupportedVersions = ["2.6.4"]
gitTarballVersions = {
  "2.6.4" : "2.6.4",
  }


#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *


class GitInstall:

  def __init__(self):
    self.dummy = None

  #
  # Called before even knowing the product version
  #

  def getScriptName(self):
    return "install-git.py"

  def getProductBaseName(self):
    return gitBaseName

  def getProductDefaultVersion(self):
    return gitDefaultVersion

  def getProductSupportedVersions(self):
    return gitSupportedVersions

  #
  # Called after knowing the product version but before parsing the
  # command-line.
  #

  def getProductName(self, version):
    return gitBaseName+"-"+version

  def getBaseDirName(self, version):
    return gitBaseName+"-"+version+"-base"

  def getExtraHelpStr(self, version):
    return """
This script builds """+self.getProductName(version)+""" from source compiled with the
configured C compiler in your path.   This also installs the git-subtree provided in the
contributed folder on install.

To also build and install the documentation, additionally, pass in:

  --with-doc --with-info

(but note the extra packages that must be installed on the system).

For more details on installing from source and required system dependencies
before attempting this build, see:

  https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

In particular, to build everything you will need to install the documentation
(--with-doc and --with-info), you will need to install the packages
'asciidoc', 'xmlto', and 'docbook2X'.  Again, see the above webpage for
details (it is not 100% pretty).

NOTE: The assumed directory structure of the download source provided by the
command --download-cmnd=<download-cmnd> is:

   git-<version>-base/
     git-<full-version>.tar.gz
"""

  def injectExtraCmndLineOptions(self, clp, version):
    setStdDownloadCmndOption(self, clp, version)
    clp.add_option(
      "--extra-configure-options", dest="extraConfigureOptions", type="string", \
      default="", \
      help="Extra options to add to the 'configure' command for "+self.getProductName(version)+"." \
        +"  Note: This does not override the hard-coded configure options." )
    clp.add_option(
      "--with-doc", dest="withDoc", action="store_true", default=False,
      help="Build and install manpage documentation (requires asciidoc and xmlto)." )
    clp.add_option(
      "--with-info", dest="withInfo", action="store_true", default=False,
      help="Build and install info documentation (requires docbook2x)." )

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
    self.gitBaseDir = self.baseDir+"/"+self.getBaseDirName(self.inOptions.version)
    gitVersionFull = gitTarballVersions[self.inOptions.version]
    self.gitTarball = "git-"+gitVersionFull+".tar.gz"
    self.gitSrcDir = "git-"+gitVersionFull
    self.gitSrcBuildDir = self.gitBaseDir+"/"+self.gitSrcDir
    self.gitSubtreeSrcBuildDir = self.gitSrcBuildDir+"/contrib/subtree" 
    self.scriptBaseDir = getScriptBaseDir()

  #
  # Called after setup()
  #

  def doDownload(self):
    removeDirIfExists(self.gitBaseDir, True)
    echoRunSysCmnd(self.inOptions.downloadCmnd)

  def doUntar(self):
    # Find the full name of the source tarball
    echoChDir(self.gitBaseDir)
    echoRunSysCmnd("tar -xzf "+self.gitTarball)

  def doConfigure(self):
    echoChDir(self.gitSrcBuildDir)
    echoRunSysCmnd("make configure")
    echoRunSysCmnd(
      "./configure "+\
      " "+self.inOptions.extraConfigureOptions+\
      " --prefix="+self.inOptions.installDir
      )
    # NOTE: Git appears to only allow an in-source build :-(

  def doBuild(self):
    echoChDir(self.gitSrcBuildDir)
    echoRunSysCmnd("make "+getParallelOpt(self.inOptions, "-j")+\
      self.inOptions.makeOptions+" all")
    if self.inOptions.withDoc:
      echoRunSysCmnd("make "+getParallelOpt(self.inOptions, "-j")+\
        self.inOptions.makeOptions+" doc")
    if self.inOptions.withInfo:
      echoRunSysCmnd("make "+getParallelOpt(self.inOptions, "-j")+\
        self.inOptions.makeOptions+" info")
    # Build git-subtree to get ready to install
    echoChDir(self.gitSubtreeSrcBuildDir)


  def doInstall(self):
    echoChDir(self.gitSrcBuildDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions+" install")
    if self.inOptions.withDoc:
      echoRunSysCmnd("make "+self.inOptions.makeOptions+" install-doc")
      echoRunSysCmnd("make "+self.inOptions.makeOptions+" install-html")
    if self.inOptions.withInfo:
      echoRunSysCmnd("make "+self.inOptions.makeOptions+" install-info")
    # Install git-subtree and documentation
    echoChDir(self.gitSubtreeSrcBuildDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions+" install")
    if self.inOptions.withDoc:
      echoRunSysCmnd("make "+self.inOptions.makeOptions+" install-doc")

  def getFinalInstructions(self):
    return """
To use the installed version of git-"""+self.inOptions.version+""" add the path:

  """+self.inOptions.installDir+"""/bin

to your path and that should be it!
"""

#
# Executable statements
#

gitInstaller = InstallProgramDriver(GitInstall())
gitInstaller.runDriver()
