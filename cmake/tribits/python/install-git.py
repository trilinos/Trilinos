#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

gitVersion = "1.7.0.4"
egVersion = "1.7.0.4"


from InstallProgramDriver import *
from GeneralScriptSupport import *


class GitInstall:

  def __init__(self):
    self.dummy = None

  def getProductName(self):
    return "git-"+gitVersion

  def getScriptName(self):
    return "install-git.py"

  def getExtraHelpStr(self):
    return """
This script installs both git-"""+gitVersion+""" and eg-"""+egVersion+""".
"""

  def getBaseDirName(self):
    return "Git.base"

  def injectExtraCmndLineOptions(self, clp):
    clp.add_option(
      "--checkout-cmnd", dest="checkoutCmnd", type="string",
      default="cvs -d :ext:software.sandia.gov:/space/CVS co -d "+self.getBaseDirName()+\
      " Trilinos3PL/Git.base",
      help="Command used to check out "+self.getProductName()+" tarball(s)." )

  def echoExtraCmndLineOptions(self, inOptions):
    cmndLine = ""
    cmndLine += "  --checkout-cmnd='" + inOptions.checkoutCmnd + "' \\\n"
    return cmndLine
    
  def setup(self, inOptions):
    self.inOptions = inOptions
    self.baseDir = os.getcwd()
    self.gitBaseDir = self.baseDir+"/"+self.getBaseDirName()
    self.gitSrcBuildDir = self.gitBaseDir+"/"+self.getProductName()

  def doCheckout(self):
    echoRunSysCmnd(self.inOptions.checkoutCmnd)

  def doUntar(self):
    echoChDir(self.gitBaseDir)
    echoRunSysCmnd("tar -xzvf git-"+gitVersion+".tar.gz")
    echoRunSysCmnd("tar -xzvf git-subtree.tar.gz")

  def doConfigure(self):
    echoChDir(self.gitSrcBuildDir)
    echoRunSysCmnd("./configure --prefix="+self.inOptions.installDir)

  def doBuild(self):
    echoChDir(self.gitSrcBuildDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.gitSrcBuildDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions+" install")
    echoRunSysCmnd("cd "+self.inOptions.installDir+"/share/man" \
      +" && tar -xzvf "+self.gitBaseDir+"/git-manpages-"+gitVersion+".tar.gz")
    echoChDir(self.gitBaseDir)
    gitSubtreeInstallFile = self.inOptions.installDir+"/bin/git-subtree"
    echoRunSysCmnd("cp git-subtree/git-subtree.sh "+gitSubtreeInstallFile)
    echoRunSysCmnd("chmod a+rx "+gitSubtreeInstallFile)
    egInstallFile = self.inOptions.installDir+"/bin/eg"
    echoRunSysCmnd("cp eg."+egVersion+".pl "+egInstallFile)
    echoRunSysCmnd("chmod a+rx "+egInstallFile)

  def getFinalInstructions(self):
    return """
In order to use """+self.getProductName()+""" append

   """+self.inOptions.installDir+"""/bin

to your PATH!
"""


#
# Executable statements
#

gitInstaller = InstallProgramDriver(GitInstall())

gitInstaller.runDriver()
