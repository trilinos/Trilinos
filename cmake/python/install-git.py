#!/usr/bin/env python

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
    return "install-git"

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
