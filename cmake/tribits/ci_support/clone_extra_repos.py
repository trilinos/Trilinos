#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

from CheckinTestConstants import *
from FindGeneralScriptSupport import *
from GeneralScriptSupport import *
from gitdist import addOptionParserChoiceOption
import gitdist


#
# Constants
#

g_extraReposTypes = g_knownTribitsTestRepoTypes
g_extraReposTypeDefault = "Nightly"
g_extraReposTypesDefaulIdx = findInSequence(g_extraReposTypes, g_extraReposTypeDefault)

g_extraRerposFileDefault = "cmake/ExtraRepositoriesList.cmake"

g_verbosityLevels = [ "none", "minimal",  "more", "most" ]
g_verbosityLevelDefault = "more"
g_verbosityLevelDefaultIdx = findInSequence(g_verbosityLevels, g_verbosityLevelDefault)


#
# Help message
#

# This part can be reused in other scripts that are project-specific
genericUsageHelp = \
r"""
By default, this will clone all the 'Nightly' extra repos that are listed in
the file:

  <projectDir>/"""+g_extraRerposFileDefault+r"""

(other repo types can be selected using --extra-repos-type).

The list of which repos to clone can be "white-list" selected with the option
--extra-repos (see options below for details).  Extra repos can in addition be
"back-listed" using the option --not-extra-repos.

To see the full list of repos that can be cloned, pass in just:

  --skip-clone --verbosity=more

That will print out a table like:

  ------------------------------------------------------------------------------
  | ID | Repo Name  | Repo Dir   | VC  | Repo URL                 | Category   |
  |----|------------|------------|-----|--------------------------|------------|
  |  1 | ExtraRepo1 | ExtraRepo1 | GIT | someurl.com.ExtraRepo1   | Continuous |
  |  2 | ExtraRepo3 | ExtraRepo3 | GIT | someurl3.com:/ExtraRepo3 | Continuous |
  ------------------------------------------------------------------------------

If the git repo server is using gitolite, one can set
--gitolite-root=<gitolite-root> and that will result in git repos being
selected only if the selected repos are listed in 'ssh <gitolite-root> info'.
This allows one to automatically exclude repos from being cloned that the user
has no permissions to clone.  NOTE: See warning about the --gitolite-root option below!

TIP: After cloning the set of repos, a nice way to interact with the repos is
to use the tool 'gitdist'.  If your project does not have a version controlled
.gitdist.default file, you can generate one using the
--create-gitdist-file=<gitdist-file> argument, for example with:

  --create-gitdist-file=.gitdist

This will restrict the list of repos processed by gitdist to just the repos cloned.
"""


usageHelp = r"""clone_extra_repos.py [options]

This script clones one more extra repos listed in a TriBITS
ExtraRepositoriesList.cmake file.  The standard usage is:

  $ cd base <projectDir>
  $ ./cmake/tribits/ci_support/clone-extra-repos.py

where <projectDir> is the base TriBITS project dir and base git repo.
""" + \
genericUsageHelp


#
# Helper functions
#


def injectCmndLineOptionsInParser(clp, gitoliteRootDefault=""):
  
  clp.add_option(
    "--extra-repos", dest="extraRepos", type="string", default="",
    help="List of names of extra repos to be cloned <extra-repos>" \
      " (i.e. \"repo0,repo1,,...\").  When set to empty '' (the default value)" \
      " then all repos that match <extra-repos-type> listed in <extra-repos-file>" \
      " will be selected.  But the repos listed in <extra-repos> must always" \
      " be a subset of the repos of type <extra-repos-type> selected from" \
      " <extra-repos-file>.   (Default '')" )
  
  clp.add_option(
    "--not-extra-repos", dest="notExtraRepos", type="string", default="",
    help="List of names of extra repos *NOT* to clone (i.e. \"repo0,repo1,...\")." \
      "  (Default '')" )
  
  clp.add_option(
    "--extra-repos-file", dest="extraReposFile", type="string",
     default=g_extraRerposFileDefault,
    help="The file path <extra-repos-file> for the ExtraRepositoriesList.cmake file." \
      "  This can be an absolute or relative path." \
      "  (Default = '"+g_extraRerposFileDefault+"')")

  addOptionParserChoiceOption(
    "--extra-repos-type", "extraReposType",
    g_extraReposTypes , g_extraReposTypesDefaulIdx,
    "Type of extra repositories <extra-repos-type> to select from " \
      "<extra-repos-file>.  When --extra-repos is set, then this argument" \
      " is ignored.",
    clp )
  
  clp.add_option(
    "--gitolite-root", dest="gitoliteRoot", type="string", default=gitoliteRootDefault,
    help="Gives the root for a gitolite repos <gitolite-root> (e.g. git@<some-url>)." \
      "  If specified, then any git repos with the <gitolite-root> listed as their" \
      " root will only be selected if they are listed with 'R' permissions returned" \
      " from 'ssh <gitolite-root> info'.  WARNING: Make sure that you have your" \
      " gitoliote SSH registered correctly before using this option by typing" \
      " the command 'ssh <gitlite-root> info' and make sure that it does *not*"
      " ask for a password! (Default = '"+gitoliteRootDefault+"')" )

  clp.add_option(
    "--with-cmake", dest="withCmake", type="string", default="cmake",
    help="CMake executable to use with cmake -P scripts internally (only set" \
    +" by unit testing code).  (Default = 'cmake')")

  addOptionParserChoiceOption(
    "--verbosity", "verbLevel", g_verbosityLevels, g_verbosityLevelDefaultIdx,
    "Verbosity of the script (levels are cumulative):" \
    "  none = no output at all (except for commands with --no-op). " \
    "  minimal = print script args echo and clone commands." \
    "  more = print basic repo include/exclude logic and print repo table." \
    "  most = print output from cmake script called, the output from gitolite," \
    " and other detailed info." \
    ,
    clp )

  clp.add_option(
    "--do-clone", dest="doClone", action="store_true",
    help="Do the clone of the selected repos. [default]")
  clp.add_option(
    "--skip-clone", dest="doClone", action="store_false",
    help="Skip the clone of the repos and just show what would be done.",
    default=True )

  clp.add_option(
    "--do-op", dest="doOp", action="store_true",
    help="Do the clone of the selected repos. [default]" )
  clp.add_option(
    "--no-op", dest="doOp", action="store_false",
    help="Skip cloning the repos and just show the clone commands.",
    default=True )
  
  clp.add_option(
    "--create-gitdist-file", dest="createGitdistFile", type="string", default="",
    help="If specified, the file <gitdist-file> will get generated with the list" \
      " of git repos (the same list that is cloned with --do-clone)." \
      "  (Default = '')")
  
  clp.add_option(
    "--show-defaults", dest="showDefaults", action="store_true",
    help="Show the default option values and do nothing at all.",
    default=False )


def getCmndLineOptions():
  from optparse import OptionParser
  clp = OptionParser(usage=usageHelp)
  injectCmndLineOptionsInParser(clp)
  (options, args) = clp.parse_args()
  return options


def isVerbosityLevel(inOptions, testVerbLevel):
  requestedVerbLevelInt = findInSequence(g_verbosityLevels, inOptions.verbLevel)
  testVerbLevelInt = findInSequence(g_verbosityLevels, testVerbLevel)
  if testVerbLevelInt <= requestedVerbLevelInt:
    return True
  return False


def fwdCmndLineOptions(inOptions, terminator=""):
  cmndLineOpts = \
    "  --extra-repos='"+inOptions.extraRepos+"'"+terminator +  \
    "  --not-extra-repos='"+inOptions.notExtraRepos+"'"+terminator +  \
    "  --extra-repos-file='"+inOptions.extraReposFile+"'"+terminator +  \
    "  --extra-repos-type='"+inOptions.extraReposType+"'"+terminator +  \
    "  --gitolite-root='"+inOptions.gitoliteRoot+"'"+terminator +  \
    "  --with-cmake='"+inOptions.withCmake+"'"+terminator +  \
    "  --verbosity='"+inOptions.verbLevel+"'"+terminator
  if inOptions.doClone:
    cmndLineOpts += "  --do-clone" + terminator
  else:
    cmndLineOpts +="  --skip-clone" + terminator
  if inOptions.doOp:
    cmndLineOpts += "  --do-op" + terminator
  else:
    cmndLineOpts += "  --no-op" + terminator
  cmndLineOpts += \
    "  --create-gitdist-file='"+inOptions.createGitdistFile+"'"+terminator
  return cmndLineOpts 


def echoCmndLineOptions(inOptions):
  print(fwdCmndLineOptions(inOptions, " \\\n"))


def echoCmndLine(inOptions):

  print("")
  print("**************************************************************************")
  print("Script: clone_extra_repos.py \\")

  echoCmndLineOptions(inOptions)


def getHeaderOutputAndExtraReposDictList(rawOutputFromCmakefile):
  rawOutputFromCmakefile = s(rawOutputFromCmakefile)
  headerOuput = ""
  pythonDictListStr = ""
  processingPythonDict = False
  for line in rawOutputFromCmakefile.splitlines():
    if line == "*** Extra Repositories Python Dictionary":
      processingPythonDict=True
      continue
    if processingPythonDict:
      pythonDictListStr += (line + "\n")
    else:
      headerOuput += (line + "\n")
  #print("\nheaderOuput:\n\n", headerOuput)
  #print("\npythonDictListStr = '"+pythonDictListStr+"'")
  pythonDictList = eval(pythonDictListStr)
  return (headerOuput, pythonDictList)


def getExtraReposDictListFromCmakefile(projectDir, extraReposFile, withCmake,
  extraReposType=g_extraReposTypeDefault, extraRepos="",
  tribitsDir=tribitsDir, verbose=True,
  ):
  cmnd = "\""+withCmake+"\""+ \
    " -DPROJECT_SOURCE_DIR="+projectDir+ \
    " -DTRIBITS_BASE_DIR="+tribitsDir+ \
    " -DEXTRA_REPOS_FILE="+extraReposFile+ \
    " -DENABLE_KNOWN_EXTERNAL_REPOS_TYPE="+extraReposType+ \
    " -DEXTRA_REPOS="+extraRepos+\
    " -DCHECK_EXTRAREPOS_EXIST=FALSE"
#  " -DTRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG=TRUE"
  cmnd += \
    " -P "+ciSupportDir+"/TribitsGetExtraReposForCheckinTest.cmake"
  rawOutput = getCmndOutput(cmnd, throwOnError=True, getStdErr=True)
  (headerOutput, extraReposPytonDictList) = getHeaderOutputAndExtraReposDictList(rawOutput)
  if verbose:
    print("\n" + headerOutput)
  return extraReposPytonDictList


# Parse raw 'ssh <gitolite-root> info' output and return list of gitolite
# repos.
def parseRawSshGitoliteRootInfoOutput(rawSshGitoliteRootInfoOutput, verbose=False):

  rawSshGitoliteRootInfoOutputList = rawSshGitoliteRootInfoOutput.splitlines()

  gitoliteSshHeader = rawSshGitoliteRootInfoOutputList[0]
  if verbose:
    print(gitoliteSshHeader)

  gitolioteReposList = []
  parsingRepos = False
  for line in rawSshGitoliteRootInfoOutputList:
    #print("line = '"+line+"'")
    # Look for the first blank line and then the next line will be the first
    # listed repo.
    if line == "":
      parsingRepos = True
      continue
    # Parse the repos
    if parsingRepos:
      # The permissions and the repo name are split by a tab such as:
      #  R W	ExtraRepo1
      repoSplit = line.split("\t")
      #print("repoSplit =", repoSplit)
      if len(repoSplit) == 2:
        gitolioteReposList.append(repoSplit[1])
  return gitolioteReposList


# Get the list of gitolite repos from 'ssh <gitolite-root> info'.
def getGitoliteReposList(inOptions, trace=False, dumpGitoliteOutput=False):
  cmnd = "ssh "+inOptions.gitoliteRoot+" info"
  if trace:
    "Running: "+cmnd
  rawSshGitoliteRootInfoOutput = getCmndOutput(cmnd)
  if dumpGitoliteOutput:
    print(rawSshGitoliteRootInfoOutput)
  if (trace and not dumpGitoliteOutput): verbose=True
  else: verbose=False
  return parseRawSshGitoliteRootInfoOutput(rawSshGitoliteRootInfoOutput, verbose=verbose)


def filterOutNotExtraRepos(extraRepoDictList_in, notExtraRepos, verbose=False):
  notExtraReposSet = set(notExtraRepos)
  extraRepoDictList = []
  for extraReposDict in extraRepoDictList_in:
    extraRepoName = extraReposDict["NAME"]
    if extraRepoName in notExtraReposSet:
      if verbose:
        print("Excluding extra repo '"+extraRepoName+"'!")
    else:
      extraRepoDictList.append(extraReposDict)
  return extraRepoDictList


def filterOutMissingGitoliteRepos(extraRepoDictList_in,
  gitoliteRepos, verbose=False):
  gitoliteReposSet = set(gitoliteRepos)
  extraRepoDictList = []
  for extraReposDict in extraRepoDictList_in:
    extraRepoName = extraReposDict["NAME"]
    if not extraRepoName in gitoliteReposSet:
      if verbose:
        print("Excluding extra repo '" + extraRepoName +
              "' not listed by gitolite!")
    else:
      extraRepoDictList.append(extraReposDict)
  return extraRepoDictList


def getExtraReposTable(extraRepoDictList):

  # Get the lists for each column in the table
  repoIdList = []
  repoNameList = []
  repoDirList = []
  repoVcTypeList = []
  repoUrlList = []
  repoCategoryList = []
  repoId = 1 # Use one-base indexing to match gitdist IDs!
  for extraRepoDict in extraRepoDictList:
    repoIdList.append(str(repoId))
    repoNameList.append(extraRepoDict["NAME"])
    repoDirList.append(extraRepoDict["DIR"])
    repoVcTypeList.append(extraRepoDict["REPOTYPE"])
    repoUrlList.append(extraRepoDict["REPOURL"])
    repoCategoryList.append(extraRepoDict["CATEGORY"])
    repoId += 1

  # Create the table
  extraReposTableDictList = [
    { "label":"ID", "align":"R", "fields":repoIdList },
    { "label":"Repo Name", "align":"L", "fields":repoNameList },
    { "label":"Repo Dir", "align":"L", "fields":repoDirList },
    { "label":"VC", "align":"L", "fields":repoVcTypeList },
    { "label":"Repo URL", "align":"L", "fields":repoUrlList },
    { "label":"Category", "align":"L", "fields":repoCategoryList },
    ]
  #print("extraReposTableDictList =", extraReposTableDictList)
  
  extraReposTable = gitdist.createTable(extraReposTableDictList)

  # Return the table
  return extraReposTable


def createGitDistFile(extraRepoDictList, gitdistFile, verbose=False):
  gitdistFileStr = ""
  for extraRepoDict in extraRepoDictList:
    gitdistFileStr += extraRepoDict["DIR"]+"\n"
  if verbose:
    print("Writing the file '" + gitdistFile + "' ...")
  open(gitdistFile, 'w').write(gitdistFileStr)


def cloneExtraRepo(inOptions, extraRepoDict):
  repoName = extraRepoDict["NAME"]
  repoDir = extraRepoDict["DIR"]
  repoUrl = extraRepoDict["REPOURL"]
  repoVcType = extraRepoDict["REPOTYPE"]
  verbLevelIsMinimum = isVerbosityLevel(inOptions, "minimal")
  if verbLevelIsMinimum:
    print("\nCloning repo " + repoName + " ...")
  if os.path.exists(repoDir):
    if verbLevelIsMinimum:
      print("\n  ==> Repo dir = '" + repoDir +
            "' already exists.  Skipping clone!")
    return
  if repoVcType != "GIT":
    print("\n  ==> ERROR: Repo type '"+repoVcType+"' not supported!")
    sys.exit(1)
  cmnd = "git clone "+repoUrl+" "+repoDir
  if inOptions.doOp:
    echoRunSysCmnd(cmnd, timeCmnd=True, verbose=verbLevelIsMinimum)
  elif verbLevelIsMinimum:
    print("\nRunning: " + cmnd)


def cloneExtraRepos(inOptions):

  verbLevelIsMinimal = isVerbosityLevel(inOptions, "minimal")

  #
  # A) Get the list of extra repos
  #

  extraRepoDictList = getExtraReposDictListFromCmakefile(
    projectDir=os.getcwd(),
    extraReposFile=inOptions.extraReposFile,
    extraReposType=inOptions.extraReposType,
    extraRepos=inOptions.extraRepos,
    withCmake=inOptions.withCmake,
    verbose=isVerbosityLevel(inOptions, "most")
    )
  #print("extraRepoDictList =" + str(extraRepoDictList))

  #
  # B) Get the list of gitolite repos
  #
  
  if inOptions.gitoliteRoot:
    if verbLevelIsMinimal:
      print("\n***")
      print("*** Get the list of repos that can be cloned from gitolite server:")
      print("***\n")
    gitoliteRepos = getGitoliteReposList(inOptions,
      trace=verbLevelIsMinimal,
      dumpGitoliteOutput=isVerbosityLevel(inOptions, "most"))
  else:
    gitoliteRepos = []
  #print("gitoliteRepos =", gitoliteRepos)

  #
  # C) Filter the list of extra repos
  #

  if gitoliteRepos:
    if verbLevelIsMinimal:
      print("\n***")
      print("*** Filtering the set of extra repos based on gitolite repos:")
      print("***\n")
    extraRepoDictList = filterOutMissingGitoliteRepos(
      extraRepoDictList, gitoliteRepos, verbose=verbLevelIsMinimal)

  if inOptions.notExtraRepos:
    if verbLevelIsMinimal:
      print("\n***")
      print("*** Filtering the set of extra repos based on --not-extra-repos:")
      print("***\n")
    extraRepoDictList = filterOutNotExtraRepos(
      extraRepoDictList,
      inOptions.notExtraRepos.split(","),
      verbose=verbLevelIsMinimal)

  #
  # D) print out table of repos
  #

  if isVerbosityLevel(inOptions, "more"):

    print("\n***")
    print("*** List of selected extra repos to clone:")
    print("***\n")
  
    extraReposTable = getExtraReposTable(extraRepoDictList)
    print(extraReposTable)

  #
  # E) Clone the repos
  #

  if inOptions.doClone:

    if verbLevelIsMinimal:
      print("\n***")
      print("*** Clone the selected extra repos:")
      print("***\n")

    for extraRepoDict in extraRepoDictList:
      cloneExtraRepo(inOptions, extraRepoDict)

  #
  # F) Create the gitdist file
  #

  if inOptions.createGitdistFile:

    if verbLevelIsMinimal:
      print("\n***")
      print("*** Create the gitdist file:")
      print("***\n")

    createGitDistFile(extraRepoDictList, inOptions.createGitdistFile,
      verbose=verbLevelIsMinimal)


#
# Run the script
#

if __name__ == '__main__':

  inOptions = getCmndLineOptions()
  if isVerbosityLevel(inOptions, "minimal"):
    echoCmndLine(inOptions)
  if inOptions.showDefaults:
    sys.exit(0)

  cloneExtraRepos(inOptions)
