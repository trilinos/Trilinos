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

#
# Version info that will change with new versions
#

openmpiVersion = "1.4.3"
openmpiTarball = "openmpi-"+openmpiVersion+".tar.gz"
openmpiSrcDir = "openmpi-"+openmpiVersion


#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *


class OpenMpiInstall:

  def __init__(self):
    self.dummy = None

  def getProductName(self):
    return "openmpi-"+openmpiVersion

  def getScriptName(self):
    return "install-openmpi.py"

  def getExtraHelpStr(self):
    return """
This script builds openmpi-"""+openmpiVersion+""" from source compiled with
the configured GCC compiler (see options --with-gcc --with-lib-dir).
"""

  def getBaseDirName(self):
    return "openmpi.BASE"

  def injectExtraCmndLineOptions(self, clp):
    clp.add_option(
      "--checkout-cmnd", dest="checkoutCmnd", type="string",
      default="git clone software.sandia.gov:/space/git/TrilinosToolset/openmpi.BASE "+\
      self.getBaseDirName(),
      help="Command used to check out "+self.getProductName()+" and dependent source tarball(s)." )
    clp.add_option(
      "--extra-configure-options", dest="extraConfigureOptions", type="string", default="",
      help="Extra options to add to the 'configure' cmmand for "+self.getProductName()+"." \
      +"  Note: This does not override the hard-coded configure options." )
    clp.add_option(
      "--with-path", dest="withPath", type="string", default="",
      help="Prepended path to find executables." )
    clp.add_option(
      "--with-lib-dir", dest="withLibDir", type="string", default="",
      help="Libraries to add to LD_LIBRARY_PATH when configuring, building, and installing." )
      
  def echoExtraCmndLineOptions(self, inOptions):
    cmndLine = ""
    cmndLine += "  --checkout-cmnd='"+inOptions.checkoutCmnd+"' \\\n"
    if inOptions.extraConfigureOptions:
      cmndLine += "  --extra-configure-options='"+inOptions.extraConfigureOptions+"' \\\n"
    cmndLine += "  --with-path='"+inOptions.withPath+"' \\\n"
    cmndLine += "  --with-lib-dir='"+inOptions.withLibDir+"' \\\n"
    return cmndLine
    
  def setup(self, inOptions):
    self.inOptions = inOptions
    self.baseDir = os.getcwd()
    self.openmpiBaseDir = self.baseDir+"/"+self.getBaseDirName()
    self.openmpiSrcBaseDir = self.openmpiBaseDir+"/"+openmpiSrcDir
    self.openmpiBuildBaseDir = self.openmpiBaseDir+"/openmpi-build"
    self.scriptBaseDir = getScriptBaseDir()

  def doCheckout(self):
    echoRunSysCmnd(self.inOptions.checkoutCmnd)

  def doUntar(self):
    echoChDir(self.openmpiBaseDir)
    echoRunSysCmnd("tar -xzvf "+openmpiTarball)

  def doConfigure(self):
    createDir(self.openmpiBuildBaseDir, True, True)
    echoRunSysCmnd(
      self.getEnvCmnd()+"../"+openmpiSrcDir+"/configure "+\
      " "+self.inOptions.extraConfigureOptions+\
      " --prefix="+self.inOptions.installDir)

  def doBuild(self):
    echoChDir(self.openmpiBuildBaseDir)
    echoRunSysCmnd(self.getEnvCmnd()+"make "+self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.openmpiBuildBaseDir)
    echoRunSysCmnd(self.getEnvCmnd()+"make "+self.inOptions.makeOptions+" install")

  def getFinalInstructions(self):
    return """
To use the installed version of openmpi-"""+openmpiVersion+""" access the include
headers in the directory:

  """+self.inOptions.installDir+"""/include

the libraries in the directory:

  """+self.inOptions.installDir+"""/lib

and the executables and compiler wrappers (e.g. mpicc, mpic++, mpirun, etc.) in:

  """+self.inOptions.installDir+"""/bin

Actually, using the MPI compiler wrappers mpicc, mpic++, etc. is prefered.
"""

  def getEnvCmnd(self):
    envCmnd = "env "
    if self.inOptions.withPath:
      envCmnd = "PATH="+self.inOptions.withPath+":$PATH "
    if self.inOptions.withLibDir:
      envCmnd += "LD_LIBRARY_PATH="+self.inOptions.withLibDir+":$LD_LIBRARY_PATH "
    envCmnd += "echo 'PATH='$PATH && echo 'LD_LIBRARY_PATH='$LD_LIBRARY_PATH && which gcc && gcc --version && "
    return envCmnd


#
# Executable statements
#

openmpiInstaller = InstallProgramDriver(OpenMpiInstall())

openmpiInstaller.runDriver()
