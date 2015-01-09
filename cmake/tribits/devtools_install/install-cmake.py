#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
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

cmakeVersion = "2.8.11"
cmakeVersionFull = "2.8.11.2"
# ToDo: Allow user to select a new version from the env
# (e.g. TRIBITS_INSTALL_CMAKE_VERSION=2.8.12,
# TRIBITS_INSTALL_CMAKE_VERSION_FULL=2.8.12.1).

cmakeTarball = "cmake-"+cmakeVersionFull+".tar.gz"
cmakeSrcDir = "cmake-"+cmakeVersionFull
cmakeDefaultCheckoutCmnd = \
  "git clone https://github.com/TriBITSPub/devtools-cmake-"+cmakeVersion+"-base" \
  + " cmake-"+cmakeVersion+"-base" 

#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *


class CMakeInstall:

  def __init__(self):
    self.dummy = None

  def getProductName(self):
    return "cmake-"+cmakeVersion

  def getBaseDirName(self):
    return "cmake-"+cmakeVersion+"-base"

  def getScriptName(self):
    return "install-cmake.py"

  def getExtraHelpStr(self):
    return """
This script builds CMake"""+self.getProductName()+""" from source compiled with the
configured C/C++ compilers in your path.  Note that the provided CMake
configure script actually builds a local bootstrap copy of CMake first, before
building the final version of CMake and the rest of the tools that gets
installed.

WARNING: By default CMake builds an unoptimized version!  To build an optimized
version you must set CXXFLAGS and CFLAGS with:

  env CXXFLAGS=-O3 CFLAGS=-O3 install-cmake.py [other options]

There does not appear to be any other way to set compiler flags than
on the env.

NOTE: The assumed directory structure of the download source provided by the
command --checkout-cmnd=<checkout-cmnd> is:

   cmake-<version>-base/
     cmake-<version>.tar.gz

ToDo: Allow user to select different cmake versions.
"""

  def injectExtraCmndLineOptions(self, clp):
    clp.add_option(
      "--checkout-cmnd", dest="checkoutCmnd", type="string",
      default=cmakeDefaultCheckoutCmnd,
      help="Command used to check out CMake and dependent source tarball(s)." \
        +"  (Default ='"+cmakeDefaultCheckoutCmnd+"')" )
    clp.add_option(
      "--extra-configure-options", dest="extraConfigureOptions", type="string", \
      default="", \
      help="Extra options to add to the 'configure' command for "+self.getProductName()+"." \
      +"  Note: This does not override the hard-coded configure options." )
    clp.add_option(
      "--parallel", dest="parallel", type="int", \
      default=0, \
      help="Uses parallelism in build if set to > 0." )

  def echoExtraCmndLineOptions(self, inOptions):
    cmndLine = ""
    cmndLine += "  --checkout-cmnd='"+inOptions.checkoutCmnd+"' \\\n"
    return cmndLine
    
  def setup(self, inOptions):
    self.inOptions = inOptions
    self.baseDir = os.getcwd()
    self.cmakeBaseDir = self.baseDir+"/"+self.getBaseDirName()
    self.cmakeSrcBaseDir = self.cmakeBaseDir+"/"+cmakeSrcDir
    self.cmakeBuildBaseDir = self.cmakeBaseDir+"/cmake-build"
    self.scriptBaseDir = getScriptBaseDir()

  def getParallelOpt(self, optName):
    if self.inOptions.parallel > 0:
      return " "+optName+str(self.inOptions.parallel)
    return " "

  def doCheckout(self):
    echoRunSysCmnd(self.inOptions.checkoutCmnd)

  def doUntar(self):
    echoChDir(self.cmakeBaseDir)
    echoRunSysCmnd("tar -xzvf "+cmakeTarball)

  def doConfigure(self):
    createDir(self.cmakeBuildBaseDir, True, True)
    echoRunSysCmnd(
      "../"+cmakeSrcDir+"/configure "+\
      " "+self.inOptions.extraConfigureOptions+\
      self.getParallelOpt("--parallel=")+\
      " --prefix="+self.inOptions.installDir)

  def doBuild(self):
    echoChDir(self.cmakeBuildBaseDir)
    echoRunSysCmnd("make "+self.getParallelOpt("-j")+self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.cmakeBuildBaseDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions+" install")

  def getFinalInstructions(self):
    return """
To use the installed version of cmake-"""+cmakeVersion+""" add the path:

  """+self.inOptions.installDir+"""/bin

to your path and that should be it!
"""


#
# Executable statements
#

cmakeInstaller = InstallProgramDriver(CMakeInstall())
cmakeInstaller.runDriver()
