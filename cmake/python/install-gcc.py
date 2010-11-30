#!/usr/bin/env python

#
# Version info that will change with new versions
#

gccVersion = "4.5.1"
gccTarball = "gcc-"+gccVersion+".tar.gz"
gccSrcDir = "gcc-"+gccVersion

gmpVersion = "4.3.2"
gmpTarball = "gmp-"+gmpVersion+".tar.gz"
gmpSrcDir = "gmp-"+gmpVersion

mpfrVersion = "2.4.2"
mpfrTarball = "mpfr-"+mpfrVersion+".tar.gz"
mpfrSrcDir = "mpfr-"+mpfrVersion

mpcVersion = "0.8.1"
mpcTarball = "mpc-"+mpcVersion+".tar.gz"
mpcSrcDir = "mpc-"+mpcVersion

# ToDo: Move js and dehydra to their own install-js.py and install-dehydra.py
# scripts!

jsVersion = "20101126"
jsTarball = "js-20101126-autoconfed.tar.gz"
jsSrcDir = "js"

dehydraVersion = "20101125"
dehydraTarball = "dehydra-20101125.tar.gz"
dehydraSrcDir = "dehydra"


#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *


class GccInstall:

  def __init__(self):
    self.dummy = None

  def getProductName(self):
    return "gcc-"+gccVersion

  def getScriptName(self):
    return "install-gcc.py"

  def getExtraHelpStr(self):
    return """
This script builds gcc-"""+gccVersion+""" from source which includes the
sources for gmp-"""+gmpVersion+""", mpfr-"""+mpfrVersion+""", and mpc-"""+mpcVersion+"""
which are built and installed along with GCC.
"""

  def getBaseDirName(self):
    return "gcc.BASE"

  def injectExtraCmndLineOptions(self, clp):
    clp.add_option(
      "--checkout-cmnd", dest="checkoutCmnd", type="string",
      default="git clone software.sandia.gov:/space/git/TrilinosToolset/gcc.BASE "+\
      self.getBaseDirName(),
      help="Command used to check out "+self.getProductName()+" and dependent source tarball(s)." )

  def echoExtraCmndLineOptions(self, inOptions):
    cmndLine = ""
    cmndLine += "  --checkout-cmnd='" + inOptions.checkoutCmnd + "' \\\n"
    return cmndLine
    
  def setup(self, inOptions):
    self.inOptions = inOptions
    self.baseDir = os.getcwd()
    self.gccBaseDir = self.baseDir+"/"+self.getBaseDirName()
    self.gccSrcBaseDir = self.gccBaseDir+"/"+gccSrcDir
    self.gccBuildBaseDir = self.gccBaseDir+"/gcc-build"

  def doCheckout(self):
    echoRunSysCmnd(self.inOptions.checkoutCmnd)

  def doUntar(self):
    echoChDir(self.gccBaseDir)
    echoRunSysCmnd("tar -xzvf "+gccTarball)
    echoRunSysCmnd("tar -xzvf "+gmpTarball)
    echoRunSysCmnd("tar -xzvf "+mpfrTarball)
    echoRunSysCmnd("tar -xzvf "+mpcTarball)

  def doConfigure(self):
    echoChDir(self.gccSrcBaseDir)
    echoRunSysCmnd("ln -sf ../"+gmpSrcDir+" gmp")
    echoRunSysCmnd("ln -sf ../"+mpfrSrcDir+" mpfr")
    echoRunSysCmnd("ln -sf ../"+mpcSrcDir+" mpc")
    createDir(self.gccBuildBaseDir, True, True)
    echoRunSysCmnd("../"+gccSrcDir+"/configure --enable-languages='c,c++,fortran'"+\
      " --disable-gnu-unique-object --prefix="+self.inOptions.installDir)

  def doBuild(self):
    echoChDir(self.gccBuildBaseDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.gccBuildBaseDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions + " install")

  def getFinalInstructions(self):
    return """
In order to use """+self.getProductName()+""" prepend

   """+self.inOptions.installDir+"""/bin

to your PATH env variable.

Also, you must prepend

   """+self.inOptions.installDir+"""/lib[64]

to your LD_LIBRARY_PATH env variable.
"""


#
# Executable statements
#

gitInstaller = InstallProgramDriver(GccInstall())

gitInstaller.runDriver()
