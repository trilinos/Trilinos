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

cmakeBaseName = "cmake"
cmakeDefaultVersion = "3.23.4"


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
