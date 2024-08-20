# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

from FindGeneralScriptSupport import *
from gitdist import addOptionParserChoiceOption

from optparse import OptionParser


#
# Implement behavior of an tool install
#
class InstallProgramDriver:


  def __init__(self, installObj):
    self.installObj = installObj
    self.productVersion = self.installObj.getProductDefaultVersion()


  def getProductVersion(self):
    return self.productVersion


  def runDriver(self):

    # Get the production version out of command-line

    productBaseName = self.installObj.getProductBaseName()

    versionCmndArgName = "--"+productBaseName+"-version"    
    for arg in sys.argv[1:]:
      #print("arg = '"+arg+"'")
      arg_and_value = arg.split("=")
      if len(arg_and_value) and arg_and_value[0] == versionCmndArgName:
        self.productVersion = arg_and_value[1].strip()

    # Get basic info after knowing the version

    productName = self.installObj.getProductName(self.productVersion)
    baseDirName = self.installObj.getBaseDirName(self.productVersion)

    scriptName = self.installObj.getScriptName();


    #
    # 1) Set up the help text
    #

    productSupportedVersions = self.installObj.getProductSupportedVersions()
    if productSupportedVersions:
      supportedVersionsTxt = "Versions supported include:\n   "+\
        str(self.installObj.getProductSupportedVersions())
    else:
      supportedVersionsTxt = "Arbitrary versions are supported."
      
    usageHelp = scriptName+\
r""" [OPTIONS] [--install-dir=<install-dir> ...]

This script checks out source, untars, configures, builds, and installs
"""+productName+r""" in one shot.  """+supportedVersionsTxt+r"""

The version to install is set with --"""+productBaseName+r"""-version=<version>.

By default, if you just type:

   $ """+scriptName+r""" --install-dir=<install-dir> --parallel=<num-procs> --do-all

then the directory """+baseDirName+r""" will be created in the local working directory
and it will contain the source for """+productName+r""" and the build files. NOTE: This
requires that you not run as root or your userid on the download computer
will not be correct.  If you want to install as root, see below.

You can control various parts of the process with various options (see below).

If you do not install as root then you must override the option --install-dir
which is set to /usr/local/bin by default.  For example, you might just type:

  $ """+scriptName+r""" --install-dir=$HOME/install --parallel=8 --do-all

and then it would install """+productName+r""" and the other executables in $HOME/install/bin.
NOTE: You will have to update your PATH variable to include whatever directory
you choose to install """+productName+r""" in.

NOTE: If you need to use sudo to install in /usr/local/bin or some other place
that needs root privileges, do:

  $ """+scriptName+r""" --install-dir=$HOME/install --parallel=8 \
     --download --untar --configure --build
  $ sudo """+scriptName+r""" --install-dir=$HOME/install --parallel=8 \
     --install-owner=<owner> --install-group=<group> [--install-for-all] \
     --install

This appears to work on most systems.

After you have done a successful install, you might want to do:

  $ rm -r """+baseDirName+r"""

in order to remove the intermediate source and build files.
""" + self.installObj.getExtraHelpStr(self.productVersion)

    #
    # 2) Parse the command-line
    #

    clp = OptionParser(usage=usageHelp)

    supportedVersions = self.installObj.getProductSupportedVersions()
    defaultVersion = self.installObj.getProductDefaultVersion()
    defaultVersionIdx = findInSequence(supportedVersions, defaultVersion)

    if supportedVersions:
      addOptionParserChoiceOption(
        versionCmndArgName, "version", supportedVersions, defaultVersionIdx,
        "Version to install for "+productName+".", clp)
    else:
      clp.add_option(
        versionCmndArgName, dest="version", type="string",
        default=defaultVersion,
        help="Version to install for "+productName+" (list of versions open-ended).")

    clp.add_option(
      "--install-dir", dest="installDir", type="string",
      default="",
      help="The install directory <install-dir> for "+productName+ \
        " (default = '').  This can be a relative or absolute path, it can" \
        " start with ~/, etc." )

    clp.add_option(
      "--install-dir-base", dest="installDirBase", type="string",
      default="",
      help="The base install directory <install-dir> for "+productName+ \
        " (default = '').  In this case the subdir "+productName+\
        " will be created under this directory for the install prefix." )

    insertInstallPermissionsOptions(clp)

    clp.add_option(
      "--parallel", dest="parallel", type="int", \
      default=0, \
      help="Uses parallelism in build if set to > 0." )
    
    clp.add_option(
      "--make-options", dest="makeOptions", type="string",
      default="",
      help="The options to pass to make for "+productName+"." )

    self.installObj.injectExtraCmndLineOptions(clp, self.productVersion)
    
    clp.add_option(
      "--show-defaults", dest="showDefaults", action="store_true", default=False,
      help="[ACTION] Show the defaults and exit." )
    
    clp.add_option(
      "--download", dest="download", action="store_true", default=False,
      help="[ACTION] Do the download of the tarball" )
    
    clp.add_option(
      "--untar", dest="untar", action="store_true", default=False,
      help="[ACTION] Do the untar of the "+productName+" sources" )
    
    clp.add_option(
      "--configure", dest="configure", action="store_true", default=False,
      help="[ACTION] Configure "+productName+" to build" )
    
    clp.add_option(
      "--build", dest="build", action="store_true", default=False,
      help="[Action] Build "+productName+" and related executables" )
    
    clp.add_option(
      "--install", dest="install", action="store_true", default=False,
      help="[ACTION] Install "+productName )
    
    clp.add_option(
      "--show-final-instructions", dest="showFinalInstructions", action="store_true",
      default=False,
      help="[ACTION] Show final instructions for using "+productName )
    
    clp.add_option(
      "--do-all", dest="doAll", action="store_true", default=False,
      help="[AGGR ACTION] Same as --download --untar --configure --build --install" \
      +" --show-final-instructions")
    
    (options, args) = clp.parse_args()
     

    #
    # 3) Echo the command-line options
    #

    cmndLine = \
      "******************************************************************************\n"
    cmndLine += scriptName + " \\\n"
    cmndLine += "  "+versionCmndArgName + "='"+options.version+"' \\\n"
    cmndLine += "  --install-dir='" + options.installDir + "' \\\n"
    cmndLine += "  --install-dir-base='" + options.installDirBase + "' \\\n"
    cmndLine += echoInsertPermissionsOptions(options)
    cmndLine += "  --parallel='" + str(options.parallel) + "' \\\n"
    cmndLine += "  --make-options='" + options.makeOptions + "'\\\n"
    cmndLine += self.installObj.echoExtraCmndLineOptions(options)
    if options.download:
      cmndLine += "  --download \\\n"
    if options.untar:
      cmndLine += "  --untar \\\n"
    if options.configure:
      cmndLine += "  --configure \\\n"
    if options.build:
      cmndLine += "  --build \\\n"
    if options.install:
      cmndLine += "  --install \\\n"
    if options.showFinalInstructions:
      cmndLine += "  --show-final-instructions \\\n"
    if options.doAll:
      cmndLine += "  --do-all \\\n"

    print(cmndLine)

    if options.showDefaults:
      return 0;

    # Check the options

    if options.installDirBase != "" and options.installDir == "":
      options.installDir=options.installDirBase+"/"+productName

    elif options.installDir == "":
      raise Exception("Error, --install-dir=<install-dir> can't be empty!")
    options.installDir = os.path.abspath(os.path.expanduser(options.installDir))

    #
    # 4) Execute the commands
    #

    if options.doAll:
      options.download = True
      options.untar = True
      options.configure = True
      options.build = True
      options.install = True
      options.showFinalInstructions = True
    
    baseDir = os.getcwd()
    
    productBaseDir = baseDir+"/"+baseDirName

    self.installObj.setup(options)

    print("")
    print("A) Download the source for "+productName+" ...")
    print("")
    
    if options.download:
      self.installObj.doDownload()
    else:
      print("Skipping on request ...")
    
    print("")
    print("B) Untar the tarball(s) and set up ready to configure ...")
    print("")
    
    if options.untar:
      self.installObj.doUntar()
    else:
      print("Skipping on request ...")
    
    print("")
    print("C) Configure "+productName+" ...")
    print("")
    
    if options.configure:
      self.installObj.doConfigure()
    else:
      print("Skipping on request ...")
    
    print("")
    print("D) Build "+productName+" ...")
    print("")
    
    if options.build:
      self.installObj.doBuild()
    else:
      print("Skipping on request ...")
    
    print("")
    print("E) Install "+productName+" ...")
    print("")
    
    if options.install:
      self.installObj.doInstall()
      fixupInstallPermissions(options, options.installDir)
    else:
      print("Skipping on request ...")
    
    
    print("")
    print("D) Final instructions for using "+productName+" ...")
    print("")
    
    if options.showFinalInstructions:
      print(self.installObj.getFinalInstructions())
    else:
      print("Skipping on request ...")
    
    print("\n[End]")


#
# Insert the standard --download-cmnd option
#

def setStdDownloadCmndOption(installObj, clp, version):
  productName = installObj.getProductBaseName()+"-"+version
  productBaseDirName = productName+"-base"
  defaultDownloadCmnd = \
    "git clone https://github.com/tribitsdevtools/"+productBaseDirName
  clp.add_option(
    "--download-cmnd", dest="downloadCmnd", type="string",
    default=defaultDownloadCmnd,
    help="Command used to download source for "+productName+"." \
      +"  (Default ='"+defaultDownloadCmnd+"')  WARNING: This will delete" \
      +" an existing directory '"+productBaseDirName+"' if it already exists!")

def setStdGithubDownloadCmndOption(installObj, githubOrg, githubRepo, clp, version):
  productName = installObj.getProductBaseName()+"-"+version
  productBaseDirName = productName+"-base"
  defaultDownloadCmnd = \
    "wget -O "+productName+".tar.gz https://github.com/"+githubOrg+"/"+githubRepo+"/tarball/v"+version
  clp.add_option(
    "--download-cmnd", dest="downloadCmnd", type="string",
    default=defaultDownloadCmnd,
    help="Command used to download source tarball for "+productName+"." \
      +"  (Default ='"+defaultDownloadCmnd+"')  WARNING: This will delete" \
      +" an existing directory '"+productBaseDirName+"' if it already exists!")

#
# Get the parallel option
#

def getParallelOpt(inOptions, optName):
  if inOptions.parallel > 0:
    return " "+optName+str(inOptions.parallel)
  return " "


#
# Insert options for fixing owner, group, permissions
#

def insertInstallPermissionsOptions(clp):
  
  clp.add_option(
    "--install-owner", dest="installOwner", type="string", default="",
    help="If set, then 'chown -R <install-owner> <install-dir>' will be run after install." \
      "  Note that you can only change the owner when running this script as sudo." )
  
  clp.add_option(
    "--install-group", dest="installGroup", type="string", default="",
    help="If set, then 'chgrp -R <install-group> <install-dir>' and " \
      "'chmod -R g+rX <install-dir> will be run after install." \
      "  Note that you can only change a to a group that the owner is" \
      " a member of." )

  clp.add_option(
    "--install-for-all", dest="installForAll", action="store_true",
    help="If set, then 'chmod -R a+rX <install-dir>' will be run after install.")
  clp.add_option(
    "--no-install-for-all", dest="installForAll", action="store_false", default=False,
    help="If set, then <install-dir> is not opened up to everyone." )


def echoInsertPermissionsOptions(inOptions):
  cmndLine = ""
  cmndLine += "  --install-owner='" + inOptions.installOwner + "' \\\n"
  cmndLine += "  --install-group='" + inOptions.installGroup + "' \\\n"
  if inOptions.installForAll:
    cmndLine += "  --install-for-all \\\n"
  else:
    cmndLine += "  --no-install-for-all \\\n"
  return cmndLine


def fixupInstallPermissions(inOptions, installDir):
  if inOptions.installOwner:
    print("\n*** Changing owner to '" + inOptions.installOwner + "':")
    echoRunSysCmnd("chown -R " + inOptions.installOwner + " " + installDir)
  if inOptions.installGroup:
    print("\n*** Changing group to '" + inOptions.installGroup + "' and giving "
          "group read/execute:")
    echoRunSysCmnd("chgrp -R " + inOptions.installGroup + " " + installDir)
    echoRunSysCmnd("chmod -R g+rX " + installDir)
  if inOptions.installForAll:
    print("\n*** Allowing everyone to read/execute:")
    echoRunSysCmnd("chmod -R a+rX " + installDir)
