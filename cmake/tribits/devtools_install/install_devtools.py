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
# Imports
#

from FindGeneralScriptSupport import *
import InstallProgramDriver


#
# Defaults and constants
#

devtools_install_dir = os.path.dirname(os.path.abspath(__file__))

scratch_dir = os.getcwd()

sourceGitUrlBase_default = "https://github.com/tribitsdevtools/"

# tool default versions
autoconf_version_default = "2.69"
cmake_version_default = "3.3.2"
gcc_version_default = "4.8.3"
mpich_version_default = "3.1.3"

# Common (compile independent) tools
commonToolsArray = [ "gitdist", "autoconf", "cmake" ]
commonToolsChoices = (["all"] + commonToolsArray + [""])

# Compiler toolset
compilerToolsetArray = [ "gcc", "mpich" ];
compilerToolsetChoices = (["all"] + compilerToolsetArray + [""])


#
# Utility functions
#

usageHelp = r"""install-devtools.py [OPTIONS]

This script drives the installation of a number of tools needed by many
TriBITS-based projects.  The most typical usage is to first create a scratch
directory with::

  mkdir scratch
  cd scratch

and then run:

  install-devtools.py --install-dir=<dev_env_base> \
   --parallel=<num-procs> --do-all

By default, this installs the following tools in the dev env install
directory:

  <dev_env_base>/
    common_tools/
      autoconf-<autoconf-version>/
      cmake-<cmake-version>/
      gitdist
    gcc-<gcc-version>/
      load_dev_env.[sh,csh]
      toolset/
        gcc-<gcc-version>/
        mpich-<mpich-version>/

The default versions of the tools installed are:

* autoconf-"""+autoconf_version_default+"""
* cmake-"""+cmake_version_default+"""
* gcc-"""+gcc_version_default+"""
* mpich-"""+mpich_version_default+"""

The tools installed under common_tools/ only need to be installed once
independent of any compilers that may be used to build TriBITS-based projects.

The tools installed under gcc-<gcc-version>/ are specific to a GCC compiler
and MPICH configuration and build.

The download and install of each of these tools is drive by its own
install-<toolname>.py script in the same directory as install-devtools.py.

Before running this script, some version of a C and C++ compiler must already
be installed on the system.  

At a high level, this script performs the following actions.

1) Create the base directories (if they don't already exist) and install
   load_dev_env.sh (csh).  (if --initial-setup is passed in.)

2) Download the sources for all of the requested common tools and compiler
   toolset.  (if --download is passed in.)

3) Configure, build, and install the requested common tools under
   common_tools/. (if --install is passed in.)

4) Configure, build, and install the downloaded GCC and MPICH tools.  First
   install GCC then MPICH using the installed GCC and install under
   gcc-<gcc-version>/.  (if --install is passed in.)

The informational arguments to this function are:

  --install-dir=<dev_env_base>

    The base directory that will be used for the install.  There is not
    default.  If this is not specified then it will abort.

  --source-git-url-base=<url_base>
  
    Gives the base URL for to get the tool sources from.  The default is:
  
      """+sourceGitUrlBase_default+"""
  
    This is used to build the full git URL as:
  
      <url_base><tool_name>-<tool_version>-base
  
    This can also accommodate gitolite repos and other directory structures,
    for example, with:
  
      git@<host-name>:prerequisites/
  
  --common-tools=all
  
    Specifies the tools to download and install under common_tools/.  One can
    pick specific tools with:
  
      --common-tools=autoconf,cmake,...
  
    This will download and install the default versions of these tools.  To
    select specific versions, use:
  
      --common-tools=autoconf:"""+autoconf_version_default+""",cmake:"""+cmake_version_default+""",...

    The default is 'all'.  To install none of these, pass in empty:

      --common-tools=''

    (NOTE: A version of 'git' is *not* installed using this script but can be
    installed using the script install-git.py.  But note the extra packages
    that must be installed on a system in order to fully install git and its
    documentation.  All of the git-related TriBITS tools can use any recent
    version of git and most systems will already have a current-enough version
    of git so there is no need to install one to be effective doing
    development.)
  
  --compiler-toolset=all
  
    Specifies GCC and MPICH (and other compiler-specific tools) to download
    and install under gcc-<gcc-version>/toolset/.  One can pick specific
    components with:
  
      --compiler-toolset=gcc,mpich
  
    or specific versions with:
  
      --compiler-toolset=gcc:"""+gcc_version_default+""",mpich:"""+mpich_version_default+"""
  
    Of course if one is only installing GCC with an existing installed MPICH,
    one will need to also reinstall MPICH as well.

    The default is 'all'.  To install none of these, pass in empty:

      --compiler-toolset=''

The action arguments are:

  --initial-setup: Create <dev_env_base>/ directories and install
    load_dev_env.sh
  
  --download: Download all of the requested tools
  
  --install: Configure, build, and install all of the requested tools
  
  --do-all: Do everything.  Implies --initial-setup --downlaod --install

To change modify the permissions of the installed files, see the options
--install-owner, --install-group, and --install-for-all.

Note that the user can see what operations and command would be run without
actually running them by passing in --no-op.  This can be used to show how to
run each of the individual install command so that the user can run it for
him/her-self and customize it as needed.

If the user needs more customization, then they can just run with --do-all
--no-op and see what commands are run to install things and then they can run
the commands themselves manually and make whatever modifications they need.

NOTE: The actual tool installs are performed using the scripts:

* install-autoconf.py
* install-cmake.py
* install-gcc.py
* install-git.py
* install-mpich.py
* install-openmpi.py

More information about what versions are installed, how they are installed,
etc. is found in these scripts.  Note that some of these scripts apply patches
for certain versions.  For details, look at the --help output from these
scripts and look at the implementation of these scripts.
"""        


# Get and process command-line arguments
def getCmndLineOptions(cmndLineArgs, skipEchoCmndLine=False):

  from optparse import OptionParser
  
  clp = OptionParser(usage=usageHelp)

  clp.add_option(
    "--install-dir", dest="installDir", type="string", default="",
    help="The base directory <dev_env_base> that will be used for the install." \
      +"  There is not default.  If this is not specified then will abort.")

  InstallProgramDriver.insertInstallPermissionsOptions(clp)

  clp.add_option(
    "--source-git-url-base", dest="sourceGitUrlBase", type="string",
    default=sourceGitUrlBase_default,
    help="Gives the base URL <url_base> for the git repos to object the source from.")

  clp.add_option(
    "--load-dev-env-file-base-name", dest="loadDevEnvFileBaseName",
    type="string", default="load_dev_env",
    help="Base name of the load dev env script that will be installed." \
      "  (Default = 'load_dev_env')" )

  clp.add_option(
    "--common-tools", dest="commonTools", type="string", default="all",
    help="Specifies the common tools to download and install under common_tools/." \
      "  Can be 'all', or empty '', or any combination of" \
      " '"+(",".join(commonToolsArray))+"' (separated by commas, no spaces).")

  clp.add_option(
    "--compiler-toolset", dest="compilerToolset", type="string", default="all",
    help="Specifies GCC and MPICH and other compiler-specific tools to" \
      " download and install under gcc-<gcc-version>/toolset/." \
      "  Can be 'all', or empty '', or any combination of" \
      " '"+(",".join(compilerToolsetArray))+"' (separated by commas, no spaces).")

  clp.add_option(
    "--parallel", dest="parallelLevel", type="string", default="1",
    help="Number of parallel processes to use in the build.  The default is" \
      " just '1'.  Use something like '8' to get faster parallel builds." )

  clp.add_option(
    "--do-op", dest="skipOp", action="store_false",
    help="Do all of the requested actions [default].")
  clp.add_option(
    "--no-op", dest="skipOp", action="store_true", default=False,
    help="Skip all of the requested actions and just print what would be done.")
    
  clp.add_option(
    "--show-defaults", dest="showDefaults", action="store_true", default=False,
    help="[ACTION] Show the defaults and exit." )

  clp.add_option(
    "--initial-setup", dest="doInitialSetup", action="store_true", default=False,
    help="[ACTION] Create base directories under <dev_env_base>/ and install" \
      " load_dev_env.[sh,csh].")

  clp.add_option(
    "--download", dest="doDownload", action="store_true", default=False,
    help="[ACTION] Download all of the tools specified by --common-tools" \
      " and --compiler-toolset.  WARNING:  If the source for a tool has" \
      " already been downloaded, it will be deleted (along with the build" \
      " directory) and downloaded from scratch!")

  clp.add_option(
    "--install", dest="doInstall", action="store_true", default=False,
    help="[ACTION] Configure, build, and install all of the tools specified by" \
      " --common-tools and --compiler-toolset.")
    
  clp.add_option(
    "--show-final-instructions", dest="showFinalInstructions", action="store_true",
    default=False,
    help="[ACTION] Show final instructions for using the installed dev env." )

  clp.add_option(
    "--do-all", dest="doAll", action="store_true", default=False,
    help="[AGGR ACTION] Do everything.  Implies --initial-setup --downlaod" \
      +" --install --show-final-instructions")

  (options, args) = clp.parse_args(args=cmndLineArgs)

  # NOTE: Above, in the pairs of boolean options, the *last* add_option(...) 
  # takes effect!  That is why the commands are ordered the way they are!

  #
  # Echo the command-line
  #

  if not skipEchoCmndLine:

    cmndLine = "**************************************************************************\n"
    cmndLine +=  "Script: install-devtools.py \\\n"
    cmndLine +=  "  --install-dir='"+options.installDir+"' \\\n"
    cmndLine += InstallProgramDriver.echoInsertPermissionsOptions(options)
    cmndLine +=  "  --source-git-url-base='"+options.sourceGitUrlBase+"' \\\n"
    cmndLine +=  "  --load-dev-env-file-base-name='"+options.loadDevEnvFileBaseName+"' \\\n"
    cmndLine +=  "  --common-tools='"+options.commonTools+"' \\\n"
    cmndLine +=  "  --compiler-toolset='"+options.compilerToolset+"' \\\n"
    cmndLine +=  "  --parallel='"+options.parallelLevel+"' \\\n"
    if not options.skipOp:
      cmndLine +=  "  --do-op \\\n"
    else:
      cmndLine +=  "  --no-op \\\n"
    if options.doInitialSetup:
      cmndLine +=  "  --initial-setup \\\n"
    if options.doDownload:
      cmndLine +=  "  --download \\\n"
    if options.doInstall:
      cmndLine +=  "  --install \\\n"
    if options.showFinalInstructions:
      cmndLine +=  "  --show-final-instructions \\\n"
    if options.doAll:
      cmndLine +=  "  --do-all \\\n"

    print(cmndLine)

    if options.showDefaults:
      sys.exit(0);

  #
  # Check the input arguments
  #

  if options.installDir == "":
    print("\nError, you must set --install-dir=<dev_env_base>!")
    raise Exception("Bad input option --install-dir")
  options.installDir = os.path.abspath(os.path.expanduser(options.installDir))

  if options.commonTools == "all":
    options.commonTools = ",".join(commonToolsArray)
  #print("options.commonTools = '"+options.commonTools+"'")

  if options.compilerToolset == "all":
    options.compilerToolset = ",".join(compilerToolsetArray)
  #print("options.compilerToolset = '"+options.compilerToolset+"'")

  if options.doAll:
    options.doInitialSetup = True
    options.doDownload = True
    options.doInstall = True
    options.showFinalInstructions = True

  #
  # Return the options
  #

  return options


#
# Get array of selected tools (can be empty)
#
def getToolsSelectedArray(toolsSelectedStr, validToolsArray):
  validToolsArraySet = set(validToolsArray)
  if toolsSelectedStr == "":
    return []
  toolsArray = []
  for toolName in toolsSelectedStr.split(","):
    if not toolName in validToolsArraySet:
      raise Exception("Error, '"+toolName+"' is not one of" \
        " '"+(",".join(validToolsArray))+"'")
    toolsArray.append(toolName)
  return toolsArray


#
# Do substututions in a string given replacements
#
def substituteStrings(inputStr, subPairArray):
  outputStr = ""
  inputStrArray = inputStr.splitlines()
  if inputStrArray[-1] == "": inputStrArray = inputStrArray[0:-1]
  for line in inputStrArray:
    #print("line = '"+line+"'")
    for (str1, str2) in subPairArray:
      #print("(str1, str2) =", (str1, str2))
      line = line.replace(str1, str2)
    outputStr += (line + "\n")
  return outputStr


#
# Configure a file by substituting strings
#
def configureFile(fileInPath, subPairArray, fileOutPath):
  fileInStr = open(fileInPath, 'r').read()
  fileOutStr = substituteStrings(fileInStr, subPairArray)
  open(fileOutPath, 'w').write(fileOutStr)

#
# Assert an install directory exists
#
def assertInstallDirExists(dirPath, inOptions):
  if not os.path.exists(dirPath) and not inOptions.skipOp:
    raise Exception(
      "Error, the install directory '"+dirPath+"'" \
       " does not exist!")


#
# Write the files load_dev_env.[sh,csh]
#
def writeLoadDevEnvFiles(devEnvBaseDir, compilersToolsetBaseDir, inOptions):

  subPairArray = [
    ("@DEV_ENV_BASE@", devEnvBaseDir),
    ("@CMAKE_VERSION@", cmake_version_default),
    ("@AUTOCONF_VERSION@", autoconf_version_default),
    ("@GCC_VERSION@", gcc_version_default),
    ("@MPICH_VERSION@", mpich_version_default)
    ]

  load_dev_env_base = inOptions.loadDevEnvFileBaseName

  configureFile(
    os.path.join(devtools_install_dir, "load_dev_env.sh.in"),
    subPairArray,
    os.path.join(compilersToolsetBaseDir, load_dev_env_base+".sh")
    )

  configureFile(
    os.path.join(devtools_install_dir, "load_dev_env.csh.in"),
    subPairArray,
    os.path.join(compilersToolsetBaseDir, load_dev_env_base+".csh")
    )


#
# Download the source for tool
#
def downloadToolSource(toolName, toolVer, gitUrlBase, inOptions):

  toolDir = toolName+"-"+toolVer

  print("\nDownloading the source for " + toolDir + " ...")

  outFile = toolDir+"-download.log"
  workingDir=scratch_dir
  toolSrcBaseDir = toolDir+"-base"
  targetToolSrcDir = workingDir+"/"+toolSrcBaseDir

  if os.path.exists(targetToolSrcDir):
    print("\nRemoving existing directory '" + targetToolSrcDir + "' ...")
    cmnd = "rm -rf "+targetToolSrcDir
    if not inOptions.skipOp:
      echoRunSysCmnd(cmnd)
    else:
      print("\nRunning: " + cmnd)

  cmnd = "git clone "+gitUrlBase+toolSrcBaseDir+" "+targetToolSrcDir
  if not inOptions.skipOp:
    echoRunSysCmnd(cmnd, workingDir=workingDir, outFile=outFile, timeCmnd=True)
  else:
    print("\nRunning: " + cmnd)
    print("\n  Running in working directory: " + workingDir)
    print("\n   Writing console output to file " + outFile)


#
# Install downloaded tool from source
#
def installToolFromSource(toolName, toolVer, installBaseDir,
  extraEnv, inOptions \
  ):

  toolDir = toolName+"-"+toolVer

  print("\nInstalling " + toolDir + " ...")

  outFile = toolDir+"-install.log"
  workingDir=scratch_dir
  toolInstallDir = installBaseDir+"/"+toolDir

  cmnd = devtools_install_dir+"/install-"+toolName+".py" \
    +" --"+toolName+"-version="+toolVer \
    +" --untar --configure --build --install --show-final-instructions" \
    +" --install-dir="+toolInstallDir \
    +" --parallel="+inOptions.parallelLevel \
    +" --install-owner="+inOptions.installOwner \
    +" --install-group="+inOptions.installGroup
  if inOptions.installForAll:
    cmnd += "  --install-for-all"
  if not inOptions.skipOp:
    echoRunSysCmnd(cmnd, workingDir=workingDir, outFile=outFile, timeCmnd=True,
      extraEnv=extraEnv)
  else:
    print("\nRunning: " + cmnd)
    print("\n  Running in working directory: " + workingDir)
    print("\n  Appending environment: " + str(extraEnv))
    print("\n  Writing console output to file " + outFile)


#
# Main
#

def main(cmndLineArgs):

  #
  # Get the command-line options
  #

  inOptions = getCmndLineOptions(cmndLineArgs)

  if inOptions.skipOp:
    print("\n***")
    print("*** NOTE: --no-op provided, will only trace actions and not touch the filesystem!")
    print("***\n")

  commonToolsSelected = \
    getToolsSelectedArray(inOptions.commonTools, commonToolsArray)
  print("\nSelected common tools = " + str(commonToolsSelected))
  commonToolsSelectedSet = set(commonToolsSelected)

  compilerToolsetSelected = \
    getToolsSelectedArray(inOptions.compilerToolset, compilerToolsetArray)
  print("\nSelected compiler toolset = " + str(compilerToolsetSelected))
  compilerToolsetSelectedSet = set(compilerToolsetSelected)

  gccVersion=gcc_version_default #ToDo: Make variable!

  dev_env_base_dir = inOptions.installDir

  ###
  print("\n\nA) Setup the install directory <dev_env_base> ='" +
        dev_env_base_dir + "':\n")
  ###

  dev_env_base_exists = os.path.exists(dev_env_base_dir)

  common_tools_dir = os.path.join(dev_env_base_dir, "common_tools")
  common_tools_exists = os.path.exists(common_tools_dir)

  compiler_toolset_base_dir = os.path.join(dev_env_base_dir, "gcc-"+gccVersion)
  compiler_toolset_base_exists = os.path.exists(compiler_toolset_base_dir)

  compiler_toolset_dir = os.path.join(compiler_toolset_base_dir, "toolset")
  compiler_toolset_exists = os.path.exists(compiler_toolset_dir)

  if inOptions.doInitialSetup:

    if not dev_env_base_exists:
      print("Creating directory '" + dev_env_base_dir + "' ...")
      if not inOptions.skipOp:
        os.makedirs(dev_env_base_dir)

    if not common_tools_exists:
      print("Creating directory '" + common_tools_dir + "' ...")
      if not inOptions.skipOp:
        os.makedirs(common_tools_dir)

    # Always create this directory so we can write the load_dev_env.sh script!
    if not compiler_toolset_base_exists:
      print("Creating directory '" + compiler_toolset_base_dir + "' ...")
      if not inOptions.skipOp:
        os.makedirs(compiler_toolset_base_dir)

    if not compiler_toolset_exists:
      print("Creating directory '" + compiler_toolset_dir + "' ...")
      if not inOptions.skipOp:
        os.makedirs(compiler_toolset_dir)

    print("Writing new files " + inOptions.loadDevEnvFileBaseName +
          ".[sh,csh] ...")
    if not inOptions.skipOp:
      writeLoadDevEnvFiles(dev_env_base_dir, compiler_toolset_base_dir, inOptions)

  else:

    print("Skipping setup of the install directory by request!")

    assertInstallDirExists(dev_env_base_dir, inOptions)
    assertInstallDirExists(common_tools_dir, inOptions)
    assertInstallDirExists(compiler_toolset_base_dir, inOptions)
    assertInstallDirExists(compiler_toolset_dir, inOptions)

  ###
  print("\n\nB) Download all sources for each selected tool:\n")
  ###

  if inOptions.doDownload:

    if "cmake" in commonToolsSelectedSet:
      downloadToolSource("cmake", cmake_version_default,
        inOptions.sourceGitUrlBase, inOptions)

    if "autoconf" in commonToolsSelectedSet:
      downloadToolSource("autoconf", autoconf_version_default,
        inOptions.sourceGitUrlBase, inOptions)

    if "gcc" in compilerToolsetSelectedSet:
      downloadToolSource("gcc", gcc_version_default,
        inOptions.sourceGitUrlBase, inOptions)

    if "mpich" in compilerToolsetSelectedSet:
      downloadToolSource("mpich", mpich_version_default,
        inOptions.sourceGitUrlBase, inOptions)

  else:

    print("Skipping download of the source for the tools on request!")
    if inOptions.doInstall:
      print("NOTE: The downloads had better be there for the install!")

  ###
  print("\n\nC) Untar, configure, build and install each selected tool:\n")
  ###

  if inOptions.doInstall:

    if "gitdist" in commonToolsSelectedSet:
      print("\nInstalling gitdist ...")
      echoRunSysCmnd("cp "+pythonUtilsDir+"/gitdist "+common_tools_dir+"/")
      InstallProgramDriver.fixupInstallPermissions(inOptions, common_tools_dir)

    if "cmake" in commonToolsSelectedSet:
      installToolFromSource("cmake", cmake_version_default,
        common_tools_dir, None, inOptions )

    if "autoconf" in commonToolsSelectedSet:
      installToolFromSource("autoconf", autoconf_version_default,
        common_tools_dir, None, inOptions )

    if "gcc" in compilerToolsetSelectedSet:
      installToolFromSource("gcc", gcc_version_default,
        compiler_toolset_dir, None, inOptions )

    if "mpich" in compilerToolsetSelectedSet:
      gccInstallDir = compiler_toolset_dir+"/gcc-"+gcc_version_default
      if not os.path.exists(gccInstallDir) and not inOptions.skipOp:
        raise Exception("Error, gcc has not been installed yet." \
          "  Missing directory '"+gccInstallDir+"'") 
      LD_LIBRARY_PATH = os.environ.get("LD_LIBRARY_PATH", "")
      installToolFromSource(
        "mpich",
        mpich_version_default,
        compiler_toolset_dir,
        {
          "CC" : gccInstallDir+"/bin/gcc",
          "CXX" : gccInstallDir+"/bin/g++",
          "FC" : gccInstallDir+"/bin/gfortran",
          "LD_LIBRARY_PATH" : gccInstallDir+"/lib64:"+LD_LIBRARY_PATH
         },
        inOptions
        )

  else:

    print("Skipping install of the tools on request!")

  ###
  print("\n\nD) Final instructions for using installed dev env:")
  ###
    
  if inOptions.showFinalInstructions:
    print("\nTo use the new dev env, just source the file:\n")
    print("  source " + compiler_toolset_base_dir + "/load_dev_env.sh\n")
    print("for sh or bash shells (or load_dev_env.csh for csh shell).\n")
    print("TIP: Add this source to your ~/.bash_profile!")
  else:
    print("Skipping on request ...")

  if inOptions.skipOp:
    print("\n***")
    print("*** NOTE: --no-op provided, only traced actions that would have been taken!")
    print("***")
  
  print("\n[End]")


#
# Script driver
#

if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    print(e)
    print()
    printStackTrace()
    sys.exit(1)
