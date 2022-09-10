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
# Defaults
#

cmakeBaseName = "cmake"
cmakeDefaultVersion = "3.17.4"


#
# Script code
#


from InstallProgramDriver import *
from GeneralScriptSupport import *

import os.path

devtools_install_dir = os.path.dirname(os.path.abspath(__file__))

class CMakeInstall:

  def __init__(self):
    self.dummy = None

  #
  # Called before even knowing the product version
  #

  def getScriptName(self):
    return "install-cmake.py"

  def getProductBaseName(self):
    return cmakeBaseName

  def getProductDefaultVersion(self):
    return cmakeDefaultVersion

  def getProductSupportedVersions(self):
    return []  # Support all versions!

  #
  # Called after knowing the product version but before parsing the
  # command-line.
  #

  def getProductName(self, version):
    return cmakeBaseName+"-"+version

  def getBaseDirName(self, version):
    return cmakeBaseName+"-"+version

  def getExtraHelpStr(self, version):
    return """
This script builds """+self.getProductName(version)+""" from source compiled with the
configured C/C++ compilers in your path.  To select a C and C++ compiler which is not
the default in your path, run with:

  env CC=<path-to-c-compiler> CXX=<path-to-cxx-compiler> install-cmake.py [other options]

This downloads tarballs from GitHub by default for the given cmake version.

This build script sets the environment vars CXXFLAGS=-O3 AND CFLAGS=-O3
when doing the configure.  Therefore, this builds and installs an optimized
version of CMake by default.

If CMake 3.17 is selected, a patch is applied which adds the CTEST_RESOURCE_SPEC_FILE
variable.  (For versions 3.18+ this is not needed.)

NOTE: To install CMake from the tip of a branch such as 'master', one must
override the 'download' command and eliminate the 'untar' command.  One must
stick with the directory structure that is assumed by the underlying code
which creates temp directories based in the CMake version.  For example, to
cone and install the tip of the CMake 'master' branch one can run:

  install-cmake.py \\
    --install-dir=<install-prefix./cmake-master-YYYYMMDD \\
    --cmake-version=master \\
    --download-cmnd="git clone git@github.com:kitware/cmake.git cmake-master" \\
    --parallel=15 --download --configure --build --install

What this does is to combine the 'download' and 'untar' commands together to
produce the source dir 'cmake-master' in one shot.  Note that the 'master' in
'cmake-master' in the git clone command must match the 'master' passed in the
argument '--cmake-version=master'.
 """

  def injectExtraCmndLineOptions(self, clp, version):
    setStdGithubDownloadCmndOption(self, "kitware", "cmake", clp, version)
    clp.add_option(
      "--extra-configure-options", dest="extraConfigureOptions", type="string", \
      default="", \
      help="Extra options to add to the 'configure' command for "\
        +self.getProductName(version)+"." \
        +"  Note: This does not override the hard-coded configure options." )
    addOptionParserChoiceOption(
      "--use-native-cmake-config", "cmakeNativeConfig",
      ["on", "off", "auto"], 1,
      "Use an already installed version of CMake for configuring " \
      +self.getProductName(version)+".", clp)
    addOptionParserChoiceOption(
      "--use-openssl", "useOpenSSL",
      ["on", "off"], 1,
      "Build in support for OpenSSL to submit to CDash via HTTPS for " \
      +self.getProductName(version)+".", clp)

  def echoExtraCmndLineOptions(self, inOptions):
    cmndLine = ""
    cmndLine += "  --download-cmnd='"+inOptions.downloadCmnd+"' \\\n"
    cmndLine += "  --extra-configure-options='"+inOptions.extraConfigureOptions+"' \\\n"
    cmndLine += "  --use-native-cmake-config='"+inOptions.cmakeNativeConfig+"' \\\n"
    cmndLine += "  --use-openssl='"+inOptions.useOpenSSL+"' \\\n"
    return cmndLine

  #
  # Called after parsing the command-line
  #
    
  def setup(self, inOptions):
    self.inOptions = inOptions
    self.baseDir = os.getcwd()
    self.cmakeBaseDir = self.baseDir+"/"+self.getBaseDirName(self.inOptions.version)
    cmakeVersionFull = self.inOptions.version
    self.cmakeTarball = "cmake-"+cmakeVersionFull+".tar.gz"
    self.cmakeSrcDir = "cmake-"+cmakeVersionFull
    self.cmakeBuildBaseDir = self.cmakeBaseDir+"/cmake-build"
    self.scriptBaseDir = getScriptBaseDir()

  #
  # Called after setup()
  #

  def doDownload(self):
    removeDirIfExists(self.cmakeBaseDir, True)
    createDir(self.cmakeBaseDir, cdIntoDir=True, verbose=True)
    echoRunSysCmnd(self.inOptions.downloadCmnd)

  def doUntar(self):
    # Find the full name of the source tarball
    echoChDir(self.cmakeBaseDir)
    createDir(self.cmakeSrcDir, verbose=True)
    echoRunSysCmnd("tar -xzf "+self.cmakeTarball \
     +" -C "+self.cmakeSrcDir+" --strip-components 1")
    if self.inOptions.version.startswith("3.17"):
      echoRunSysCmnd("patch -d "+self.cmakeSrcDir+" -p1 -i " \
       +os.path.join(devtools_install_dir, "0001-CTest-Add-CTEST_RESOURCE_SPEC_FILE-variable.patch"))

  def doConfigure(self):
    createDir(self.cmakeBuildBaseDir, True, True)

    # Look for an already installed CMake
    if self.inOptions.cmakeNativeConfig is not "off":
      cmakeCmd = self.which("cmake")
    else:
      # Configuring NOT with CMake was explicitly requested
      cmakeCmd = None

    if self.inOptions.useOpenSSL == "on":
      openSSLCachVarStr=" -DCMAKE_USE_OPENSSL=ON"
    else:
      openSSLCachVarStr=" -DCMAKE_USE_OPENSSL=OFF"

    if cmakeCmd is not None:
      # CMake is installed, configure CMake with CMake
      echoRunSysCmnd(
        cmakeCmd+" ../"+self.cmakeSrcDir+\
        " "+self.inOptions.extraConfigureOptions+\
        " -DCMAKE_BUILD_TYPE=Release"+\
        openSSLCachVarStr+\
        " -DCMAKE_INSTALL_PREFIX="+self.inOptions.installDir
        )

    elif self.inOptions.cmakeNativeConfig is "on":
      # Configuring CMake with CMake was explicitly requested but
      # the CMake executable was not found in PATH
      raise Exception("Could not find 'cmake' in PATH and --use-native-cmake-config=on!")
    else:
      # CMake is not already installed on system, so use configure
      echoRunSysCmnd(
        "../"+self.cmakeSrcDir+"/bootstrap "+\
        " "+self.inOptions.extraConfigureOptions+\
        getParallelOpt(self.inOptions, "--parallel=")+\
        " --prefix="+self.inOptions.installDir+\
        " -- "+openSSLCachVarStr,
        extraEnv={"CXXFLAGS":"-O3", "CFLAGS":"-O3"},
        )

  def doBuild(self):
    echoChDir(self.cmakeBuildBaseDir)
    echoRunSysCmnd("make "+getParallelOpt(self.inOptions, "-j")+self.inOptions.makeOptions)

  def doInstall(self):
    echoChDir(self.cmakeBuildBaseDir)
    echoRunSysCmnd("make "+self.inOptions.makeOptions+" install")

  def getFinalInstructions(self):
    return """
To use the installed version of cmake-"""+self.inOptions.version+""" add the path:

  """+self.inOptions.installDir+"""/bin

to your path and that should be it!
"""
  @staticmethod
  def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

#
# Executable statements
#

cmakeInstaller = InstallProgramDriver(CMakeInstall())
cmakeInstaller.runDriver()
