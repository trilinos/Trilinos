#!/usr/bin/env python


distRepoStatusLegend = r"""Legend:
* ID: Repository ID, zero based (order git commands are run)
* Repo Dir: Relative to base repo (base repo shown first with '(Base)')
* Branch: Current branch (or detached HEAD)
* Tracking Branch: Tracking branch (or empty if no tracking branch exists)
* C: Number local commits w.r.t. tracking branch (empty if zero or no TB)
* M: Number of tracked modified (uncommitted) files (empty if zero)
* ?: Number of untracked, non-ignored files (empty if zero)
"""

usageHelp = r"""gitdist [gitdist arguments] [git arguments]
       gitdist [gitdist arguments] dist-repo-status

Run git recursively over extra repos multi-repository git projects.  Also
includes other tools like printing a repo status table and tracking versions
through RepoVersion.txt files.

Instead of typing:

  $ git [git arguments]

type:

  $ gitdist [gitdist options] [git arguments]

This will distribute git commands across all the listed git repos, including
the base git repo.  The options in [gitdist options] are prefixed with
'--dist-' and are pulled out before passing the remaining arguments in [git
arguments] to all the git repos.  See --help to see the [gitdist options].

For example, consider the following git repos with other git repos cloned
under it:

  BaseRepo/
    .git/
    ExtraRepo1/
      .git/
    ExtraRep2/
      .git/

The gitdist command is run in the base git repo 'BaseRepo' and the git
commands are run on it as well as the git repos 'ExtraRepo1' and 'ExtraRepo2'.

The set of repos processed is determined by the argument --dist-extra-repos or
the files .gitdist or .gitdist.default.  If --dist-extra-repos="", then the
list of extra repos will be read from the file '.gitdist' in the current
working directory.  If the file '.gitdist' does not exist, then the list of
extra repos will be read from the file '.gitdist.default' in the current
working directory.  The format of this files '.gitdist' and '.gitdist.default'
is to have one repo directory/name per line as in:

  ExtraRepo1
  ExtraRepo2
  ...

where each repo is the relative path to the repo under the base git repo
(e.g. under 'BaseRepo/').  The file .gitdist.default is meant to be committed
to the base git repo so that gitdist is ready to use right away after the
extra repos are cloned (see the tool clone_extra_repos.py).

If an extra repository directory listed in the .gitdist file (or specified
using --dist-extra-repos) does not exist, then it will be ignored by the
script.  Therefore, be careful to manually verify that the script recognizes
the repositories that you list.  The best way to do that is to run 'gitdist
dist-repo-status' and see which repos are listed.

This script is self-contained and has no dependencies other than standard
python 2.4 packages so it can be copied to anywhere and used.

TIPS:

 - To see the status of all repos in a compact way, use the special
   'dist-repo-status' command (see below).

 - To process only repos that are changed w.r.t. their tracking branch, run
   'gitdist --dist-mod-only [git arguments]'.  For example, to see the status
   of only changed repos use 'gitdist --dist-mod-only status'.

 - 'gitdist --help' will run gitdist help, not git help.  If you want raw git
   help, run 'git --help'.

 - By default, gitdist will use 'git' in the environment.  If it can't find
   'git' in the environment, it will require that the user specify the git
   command to run with --dist-use-git=<the command>.

 - Pass in --dist-not-base-repo and/or
   --dist-not-extra-repos=RepoX,RepoZ,... to exclude the processing of either
   the base git repo and/or other git repos listed in .gitdist, respectively.

SUMMARY OF REPO STATUS:

This script supports the special command 'dist-repo-status' which prints a
table showing the current status of all the repos.  For the repos shown above:

  $ gitdist dist-repo-status

prints a table like:

  ----------------------------------------------------------------
  | ID | Repo Dir        | Branch | Tracking Branch | C | M  | ? |
  |----|-----------------|--------|-----------------|---|----|---|
  |  0 | BaseRepo (Base) | dummy  |                 |   |    |   |
  |  1 | ExtraRepo1      | master | origin/master   | 1 |  2 |   |
  |  2 | ExtraRepo2      | HEAD   |                 |   | 25 | 4 |
  ----------------------------------------------------------------

If the option --dist-legend is passed in, it will print the legend:

"""+distRepoStatusLegend+\
r"""
REPO VERSION FILES:

This script supports the options --dist-version-file=<versionFile> and
--dist-version-file2=<versionFile2> which are used to provide different SHA1
versions for each repo.  Each of these version files is expected to represent
a compatible set of versions of the repos.

The format of these repo version files is shown in the following example:

-----------------------------------------------------
*** Base Git Repo: BaseRepo
e102e27 [Mon Sep 23 11:34:59 2013 -0400] <author1@someurl.com>
First summary message
*** Git Repo: ExtraRepo1
b894b9c [Fri Aug 30 09:55:07 2013 -0400] <author2@someurl.com>
Second summary message
*** Git Repo: ExtraRepo2
97cf1ac [Thu Dec 1 23:34:06 2011 -0500] <author3@someurl.com>
Third summary message
-----------------------------------------------------

(the lines '---------' are *not* included in the file.)

Each repository entry can have a summary message or not (i.e. use two or three
lines per repo in the file).  A compatible repo version file can be generated
with this script using, for example, using:

  $ gitdist --dist-no-color log -1 --pretty=format:"%h [%ad] <%ae>%n%s" \
    | grep -v "^$" &> RepoVersion.txt

using three lines per repo, or just:

  $ gitdist --dist-no-color log -1 --pretty=format:"%h [%ad] <%ae>" \
    | grep -v "^$" &> RepoVersion.txt

using two lines per repo in the output file.

This allows checking out consistent versions of the repos, diffing two
consistent versions of the repos, etc.

To checkout an older set of consistent versions of the set of repos
represented by the set of versions given in the file RepoVersion.<date>.txt,
use:

  $ gitdist fetch origin
  $ gitdist --dist-version-file=RepoVersion.<date>.txt checkout _VERSION_

The '_VERSION_' string will be replaced with the SHA1 for each of the repos.

To tag and branch the set of repos using a consistent set of versions, use:

  $ gitdist --dist-version-file=RepoVersion.<date>.txt \
      tag -a <some_tag> _VERSION_

To diff two sets of versions of the repos, for example, use:

  $ gitdist --dist-version-file=RepoVersion.<new-date>.txt \
      --dist-version-file2=RepoVersion.<old-date>.txt \
      diff _VERSION_ ^_VERSION2_

Here, _VERSION_ is replaced by the SHA1s listed in RepoVersion.<new-date>.txt
and _VERSION2_ is replaced by the SHA1s listed in RepoVersion.<old-date>.txt.

One can construct any git commit taking one or two different repo version
arguments (SHA1s) using this approach.

Note that the set of repos listed in the RepoVersion.txt file must be a
super-set of those processed by this script or an error will occur and the
script will stop.

USEFUL ALIASES:

Two very useful Linux/Unix aliases are:

  $ alias gitdist-status="gitdist dist-repo-status"
  $ alias gitdist-mod="gitdist --dist-mod-only"

This avoids lot of typing as these arguments are used a lot when working with
many git repos.  For example, to see the status of all your repos, do:

  $ gitdist-status | less

To see the changes in the repos use:

  $ gitdist-mod log --name-status HEAD ^origin/master

or

  $ gitdist-mod local-stat

(where 'local-stat' is a useful git alias defined in the script
'git-config-alias.sh').

SCRIPT DEPENDENCIES:

This Python script only depends on the Python 2.4+ standard modules 'sys',
'os', 'subprocess', and 're'. Also, of course, it requires some compatible
version of 'git' in your path.

"""

import sys
import os
import subprocess
import re

from optparse import OptionParser


#
# Format an ASCII table
#


# Fill in a field
def getTableField(field, width, just):
  if just == "R":
    return " "+field.rjust(width)+" |"
  return " "+field.ljust(width)+" |"


# Format an ASCII table from a set of fields
#
# The format is of tableData input is:
#
#   [ { "label":"<label0>:, "align":"<align0>, "fields":[<fld00>, ... ]}, 
#     { "label":"<label1>:, "align":"<align1>, "fields":[<fld10>, ... ]},
#     ...
#     ]
#
# The "algin" field is either "R" for right, or "L" for left.
#
def createAsciiTable(tableData):

  asciiTable = ""

  # Table size
  numFields = len(tableData)
  numRows = len(tableData[0]["fields"])

  # a) Get the max field width for each column.
  fullTableWidth = 1  # The left '|'
  tableFieldWidth = []
  for fieldDict in tableData:
    label = fieldDict["label"]
    maxFieldWidth = len(label)
    if len(fieldDict["fields"]) != numRows:
      raise Exception("Error: column '"+label+"' numfields = " + \
        str(len(fieldDict["fields"])) + " != numRows = "+str(numRows)+"\n" )
    for field in fieldDict["fields"]:
      fieldWidth = len(field)
      if fieldWidth > maxFieldWidth: maxFieldWidth = fieldWidth 
    fullTableWidth += (maxFieldWidth + 3) # begin " ", end " ", '|'
    tableFieldWidth.append(maxFieldWidth)

  # b) Write the header of the table (always left-align the colume labels)
  asciiTable += ('-'*fullTableWidth)+"\n"
  asciiTable += "|"
  fieldIdx = 0
  for fieldDict in tableData:
    asciiTable += getTableField(fieldDict["label"], tableFieldWidth[fieldIdx], "L")
    fieldIdx += 1
  asciiTable += "\n"
  asciiTable += "|"
  for field_i in range(numFields):
    asciiTable += ('-'*(tableFieldWidth[field_i]+2))+"|"
    fieldIdx += 1
  asciiTable += "\n"

  # c) Write each row of the table
  for row_i in range(numRows):
    asciiTable += "|"
    field_i = 0
    for fieldDict in tableData:
      asciiTable += getTableField(fieldDict["fields"][row_i],
        tableFieldWidth[field_i], fieldDict["align"] )
      field_i += 1
    asciiTable += "\n"
  asciiTable += ('-'*fullTableWidth)+"\n"
  
  return asciiTable


#
# Helper functions for gitdist
#


# Get output from command
def getCmndOutput(cmnd, rtnCode=False):
  child = subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE,
    stderr = subprocess.STDOUT)
  output = child.stdout.read()
  child.wait()
  if rtnCode:
    return (output, child.returncode)
  return output


# Run a command and syncronize the output
def runCmnd(options, cmnd):
  if options.debug:
    print "*** Running command:", cmnd
  if options.noOpt:
    print cmnd
  else:
    child = subprocess.Popen(cmnd, stdout=subprocess.PIPE).stdout
    output = child.read()
    sys.stdout.flush()
    print output
    sys.stdout.flush()


# Determine if a command exists:
def commandExists(cmnd):
  whichCmnd = getCmndOutput("which "+cmnd).strip()
  #print "whichCmnd =", whichCmnd
  if os.path.exists(whichCmnd):
    return True
  return False


# Get the terminal colors
txtbld=getCmndOutput(r"tput bold")       # Bold
txtblu=getCmndOutput(r"tput setaf 4")    # Blue
txtred=getCmndOutput(r"tput setaf 1")    # Red
txtrst=getCmndOutput(r"tput sgr0")       # Text reset


# Add color to the repo dirs printed out
def addColorToRepoDir(useColor, strIn):
  if useColor:
    return txtbld+txtblu+strIn+txtrst
  return strIn


# Add color to the error messages printed out
def addColorToErrorMsg(useColor, strIn):
  if useColor:
    return txtred+strIn+txtrst
  return strIn


# Get the commandline options
def getCommandlineOps():

  #
  # A) Define the native gitdist command-line arguments
  #

  clp = OptionParser(usage=usageHelp)

  helpArgName = "--help"
  withGitArgName = "--dist-use-git"
  extraRepoArgName = "--dist-extra-repos"
  notExtraRepoArgName = "--dist-not-extra-repos"
  notBaseRepoArgName = "--dist-not-base-repo"
  versionFileName = "--dist-version-file"
  versionFile2Name = "--dist-version-file2"
  noColorArgName = "--dist-no-color"
  debugArgName = "--dist-debug"
  noOptName = "--dist-no-opt"
  modifiedOnlyName = "--dist-mod-only"
  legendName = "--dist-legend"

  nativeArgNames = [ helpArgName, withGitArgName, \
    extraRepoArgName, notExtraRepoArgName, notBaseRepoArgName, \
    versionFileName, versionFile2Name, noColorArgName, debugArgName, noOptName, \
    modifiedOnlyName, legendName ]

  distRepoStatus = "dist-repo-status"
  nativeCmndNames = [ distRepoStatus ]

  # Find a default git to use

  # Select a version of git (see above help documentation)

  defaultGit = "git" # Try system git
  if not commandExists(defaultGit):
    defaultGit = "" # Give up and make the user specify

  clp.add_option(
    withGitArgName, dest="useGit", type="string",
    default=defaultGit,
    help="The (path) to the git executable to use for each git repo command (default='"+defaultGit+"')"
    )

  clp.add_option(
    extraRepoArgName, dest="extraRepos", type="string",
    default="",
    help="Comma-separated list of extra repos to forward git commands to."
    +"  If the list is empty, it will look for a file called .gitdist to"
    +" get the list of extra repos separated by newlines."
    )

  clp.add_option(
    notExtraRepoArgName, dest="notExtraRepos", type="string",
    default="",
    help="Comma separated list of extra repos to *not* forward git commands to."
    +"  This removes any repos from being processed that would otherwise be."
    )

  clp.add_option(
    notBaseRepoArgName, dest="processBaseRepo", action="store_false",
    help="If set, don't pass the git command on to the base git repo.",
    default=True )

  clp.add_option(
    versionFileName, dest="versionFile", type="string",
    default="",
    help="Path to a file contains a list of extra repo directories and git versions (replaces _VERSION_)."
    )

  clp.add_option(
    versionFile2Name, dest="versionFile2", type="string",
    default="",
    help="Path to a second file contains a list of extra repo directories and git versions (replaces _VERSION2_)."
    )

  clp.add_option(
    noColorArgName, dest="useColor", action="store_false",
    help="If set, don't use color in the output for gitdist (better for output to a file).",
    default=True )

  clp.add_option(
    debugArgName, dest="debug", action="store_true",
    help="If set, then debugging info is printed.",
    default=False )

  clp.add_option(
    noOptName, dest="noOpt", action="store_true",
    help="If set, then no git commands will be run but instead will just be printed.",
    default=False )

  clp.add_option(
    modifiedOnlyName, dest="modifiedOnly", action="store_true",
    help="If set, then the listed git command will be only be run and output for the" \
      " repo will only be produced if the command 'git diff --name-only ^<tracking-branch>'"
      " returns non-empty output where <tracking-branch> is returned" \
      " from 'rev-parse --abbrev-ref --symbolic-full-name @{u}'.  In order words," \
      " if a git repo is unchanged w.r.t. its tracking branch, then the git command is" \
      " skipped for that repo.  If a repo does not have a tracking branch, then the repo will" \
      " be skipped as well.  Therefore, be careful to first run with dist-local-stat to see the" \
      " status of each local repo to know which repos don't have tracking branches.",
    default=False )

  clp.add_option(
    legendName, dest="printLegend", action="store_true",
    help="If set, then a legend will be printed below the repo summary table"\
      " for the special dist-repo-status command.  Only applicable with dist-repo-status.",
    default=False )

  #
  # B) Pull the native commandline arguments out of the commandline
  #

  argv = sys.argv[1:]
  nativeArgs = []
  nativeCmnds = []
  otherArgs = []

  for arg in argv:
    #print "\narg = '"+arg+"'"
    matchedNativeArg = False
    for nativeArgName in nativeArgNames:
      #print "\nnativeArgName ='"+nativeArgName+"'"
      currentArgName = arg[0:len(nativeArgName)]
      #print "currentArgName = '"+currentArgName+"'"
      if currentArgName == nativeArgName:
        #print "\nMatches native arg!"
        nativeArgs.append(arg)
        matchedNativeArg = True
        break
    matchedNativeCmnd = False
    for nativeCmndName in nativeCmndNames:
      if arg == nativeCmndName:
        #print "\nMatches native cmnd!"
        nativeCmnds.append(nativeCmndName)
        matchedNativeCmnd = True
        break
    if not (matchedNativeArg or matchedNativeCmnd):
      #print "\nDoes *not* match native arg!"
      otherArgs.append(arg)
    #print "\nnativeArgs =", nativeArgs
    #print "otherArgs =", otherArgs

  #print "\nnativeArgs =", nativeArgs
  #print "nativeCmnds =", nativeCmnds
  #print "otherArgs =", otherArgs

  if len(nativeCmnds) == 0:
    nativeCmnd = None
  elif len(nativeCmnds) == 1:
    nativeCmnd = nativeCmnds[0]
  elif len(nativeCmnds) > 1:
    raise Exception("Error: Can't have more than one dist-xxx command "+\
      " but was passed in "+str(nativeCmnds))

  (options, args) = clp.parse_args(nativeArgs)

  debugFromEnv = os.environ.get("GITDIST_DEBUG_OVERRIDE")
  if debugFromEnv:
    options.debug = True

  # Check for valid usage

  if not nativeCmnd and len(otherArgs) == 0:
    print addColorToErrorMsg(options.useColor,
      "Must specify git command. See 'git --help' for options.")
    sys.exit(1)

  if not options.useGit:
    print addColorToErrorMsg(options.useColor,
      "Can't find git, please set --dist-use-git")
    sys.exit(1)

  # Get the list of extra repos

  if options.extraRepos:
    extraReposFullList = options.extraRepos.split(",")
  else:
    if os.path.exists(".gitdist"):
      gitdistfile = ".gitdist"
    elif os.path.exists(".gitdist.default"):
      gitdistfile = ".gitdist.default"
    else:
      gitdistfile = None
    if gitdistfile:
      extraReposFullList = open(gitdistfile, 'r').read().split()
    else:
      extraReposFullList = []

  # Get list of not extra repos

  if options.notExtraRepos:
    notExtraReposFullList = options.notExtraRepos.split(",")
  else:
    notExtraReposFullList = []

  return (options, nativeCmnd, otherArgs, extraReposFullList,
    notExtraReposFullList)


# Requote commandline arguments into an array
def requoteCmndLineArgsIntoArray(inArgs):
  argsArray = []
  for arg in inArgs:
    splitArg = arg.split("=")
    newArg = None
    if len(splitArg) == 1:
      newArg = arg
    else:
      newArg = splitArg[0]+"="+'='.join(splitArg[1:])
    #print "\nnewArg =", newArg
    argsArray.append(newArg)
  return argsArray


# Get a data-structure for a set of repos from a string
def getRepoVersionDictFromRepoVersionFileString(repoVersionFileStr):
  repoVersionFileStrList = repoVersionFileStr.split("\n")
  repoVersionDict = {}
  i = 0
  while i < len(repoVersionFileStrList):
    #print "i = ", i
    repoDirLine = repoVersionFileStrList[i]
    #print "repoDirLine = '"+repoDirLine+"'"
    if repoDirLine[0:3] == "***":
      repoDir = repoDirLine.split(":")[1].strip()
      #print "repoDir = '"+repoDir+"'"
      repoVersionLine = repoVersionFileStrList[i+1]
      #print "repoVersionLine = '"+repoVersionLine+"'"
      repoSha1 = repoVersionLine.split(" ")[0].strip()
      #print "repoSha1 = '"+repoSha1+"'"
      repoVersionDict.update({repoDir : repoSha1})
    else:
      break
    if repoVersionFileStrList[i+2][0:3] == "***":
      # Has no summary line
      i = i + 2
    else:
      # Has a summary line
      i = i + 3
  return repoVersionDict


# Get a data-structure for a set of repos from a file
def getRepoVersionDictFromRepoVersionFile(repoVersionFileName):
  if repoVersionFileName:
    repoVersionFileStr = open(repoVersionFileName, 'r').read()
    return getRepoVersionDictFromRepoVersionFileString(repoVersionFileStr)
  else:
    None


def assertAndGetRepoVersionFromDict(repoDirName, repoVersionDict):
  if repoVersionDict:
    repoSha1 = repoVersionDict.get(repoDirName, "")
    if not repoSha1:
      print addColorToErrorMsg(options.useColor,
        "Extra repo '"+repoDirName+"' is not in the list of extra repos "+\
        str(repoVersionDict.keys()[1:])+" read in from version file.")
      sys.exit(3)
    return repoSha1
  else:
    return ""


def replaceRepoVersionInCmndLineArg(cmndLineArg, verToken, repoDirName, repoSha1):
  if repoSha1:
    newCmndLineArg = re.sub(verToken, repoSha1, cmndLineArg)
    return newCmndLineArg
  return cmndLineArg


def replaceRepoVersionInCmndLineArgs(cmndLineArgsArray, repoDirName, \
  repoVersionDict, repoVersionDict2 \
  ):
  #print "repoDirName =", repoDirName
  repoSha1 = assertAndGetRepoVersionFromDict(repoDirName, repoVersionDict)
  repoSha1_2 = assertAndGetRepoVersionFromDict(repoDirName, repoVersionDict2)
  #print "repoSha1 =", repoSha1
  #print "repoSha1_2 =", repoSha1_2
  cmndLineArgsArrayRepo = []
  for cmndLineArg in cmndLineArgsArray:
    #print "cmndLineArg =", cmndLineArg
    newCmndLineArg = replaceRepoVersionInCmndLineArg(cmndLineArg, \
      "_VERSION_", repoDirName, repoSha1)
    #print "newCmndLineArg =", newCmndLineArg
    newCmndLineArg = replaceRepoVersionInCmndLineArg(newCmndLineArg, \
      "_VERSION2_", repoDirName, repoSha1_2)
    #print "newCmndLineArg =", newCmndLineArg
    cmndLineArgsArrayRepo.append(newCmndLineArg)
  return cmndLineArgsArrayRepo


# Generate the command line arguments
def runRepoCmnd(options, cmndLineArgsArray, repoDirName, baseDir, \
  repoVersionDict, repoVersionDict2 \
  ):
  cmndLineArgsArryRepo = replaceRepoVersionInCmndLineArgs(cmndLineArgsArray, \
    repoDirName, repoVersionDict, repoVersionDict2)
  egCmndArray = [ options.useGit ] + cmndLineArgsArryRepo
  runCmnd(options, egCmndArray)


# Determine if the extra repo should be processed or not
def repoExistsAndNotExcluded(options, extraRepo, notExtraReposList):
  if not os.path.isdir(extraRepo): return False
  if extraRepo in notExtraReposList: return False
  return True


# Get the tracking branch for a repo
def getLocalBranch(options):
  (branch, rtnCode) = getCmndOutput(
    options.useGit + " rev-parse --abbrev-ref HEAD",
    rtnCode=True )
  if rtnCode == 0:
    return branch.strip()
  return ""


# Get the tracking branch for a repo
def getTrackingBranch(options):
  (trackingBranch, rtnCode) = getCmndOutput(
    options.useGit + " rev-parse --abbrev-ref --symbolic-full-name @{u}",
    rtnCode=True )
  if rtnCode == 0:
    return trackingBranch.strip()
  return ""
  # Above, if the command failed, there is likely no tracking branch.
  # However, this could fail for other reasons so it is a little dangerous to
  # just fail and return "" but I don't know of another way to do this.


# Get number of commits as a str wr.t.t tracking branch
def getNumCommitsWrtTrackingBranch(options, trackingBranch):
  if trackingBranch == "":
    return ""
  (summaryLines, rtnCode) = getCmndOutput(
    options.useGit + " shortlog -s HEAD ^"+trackingBranch, rtnCode=True )
  if rtnCode != 0:
    raise Exception(summaryLines)
  numCommits = 0
  summaryLines = summaryLines.strip()
  if summaryLines:
    for summaryLine in summaryLines.split("\n"):
      #print "summaryLine = '"+summaryLine+"'"
      numAuthorCommits = int(summaryLine.strip().split()[0].strip())
      #print "numAuthorCommits =", numAuthorCommits
      numCommits += numAuthorCommits
  return str(numCommits)
  # NOTE: Above, we would like to use 'git ref-list --count' but that is not
  # supported in older versions of git (at least not in 1.7.0.4).  Using 'git
  # shortlog -s' will return just one line per author so this is not likley to
  # return a lot of data and the cost of the python code to process this
  # should be insignificant compared to the process execution command.


# Get the number of modified
def getNumModifiedAndUntracked(options):
  (rawStatusOutput, rtnCode) = getCmndOutput(
    options.useGit + " status --porcelain", rtnCode=True )
  if rtnCode == 0:
    numModified = 0
    numUntracked = 0
    for line in rawStatusOutput.split("\n"):
      if line.find(" M ") == 0 or line.find("M  ") == 0:
        numModified += 1
      elif line.find(" D ") == 0 or line.find("D  ") == 0:
        numModified += 1
      elif line.find(" T ") == 0 or line.find("T  ") == 0:
        numModified += 1
      elif line.find("??") == 0:
        numUntracked += 1
    return (str(numModified), str(numUntracked))
  return ("", "")


#
# Get the repo statistics
#

class RepoStatsStruct:

  def __init__(self, branch, trackingBranch, numCommits, numModified, numUntracked):
    self.branch = branch
    self.trackingBranch = trackingBranch
    self.numCommits = numCommits
    self.numModified = numModified
    self.numUntracked = numUntracked

  def __str__(self):
    return "{" \
     "branch='"+self.branch+"'," \
     " trackingBranch='"+self.trackingBranch+"'," \
     " numCommits='"+self.numCommits+"'," \
     " numModified='"+self.numModified+"'," \
     " numUntracked='"+self.numUntracked+"'" \
     "}"

  def numCommitsInt(self):
    if self.numCommits == '': return 0
    return int(self.numCommits)

  def numModifiedInt(self):
    if self.numModified == '': return 0
    return int(self.numModified)

  def numUntrackedInt(self):
    if self.numUntracked == '': return 0
    return int(self.numUntracked)

  def hasLocalChanges(self):
    if self.numCommitsInt()+self.numModifiedInt()+self.numUntrackedInt() > 0:
      return True
    return False


def getRepoStats(options):
  branch = getLocalBranch(options)
  trackingBranch = getTrackingBranch(options)
  numCommits = getNumCommitsWrtTrackingBranch(options, trackingBranch)
  (numModified, numUntracked) = getNumModifiedAndUntracked(options)
  return RepoStatsStruct(branch, trackingBranch,
   numCommits, numModified, numUntracked)


def convertZeroStrToEmpty(strIn):
  if strIn == "0":
    return ""
  return strIn


class RepoStatTable:
  
  def __init__(self):
    self.tableData = [
      { "label" : "ID", "align" : "R", "fields" : [] },
      { "label" : "Repo Dir", "align" : "L", "fields" : [] },
      { "label" : "Branch", "align":"L", "fields" : [] },
      { "label" : "Tracking Branch", "align":"L", "fields" : [] },
      { "label" : "C", "align":"R", "fields" : [] },
      { "label" : "M", "align":"R", "fields" : [] },
      { "label" : "?", "align":"R", "fields" : [] },
      ]
    self.ID = 0

  def insertRepoStat(self, repoDir, repoStat):
    self.tableData[0]["fields"].append(str(self.ID))
    self.tableData[1]["fields"].append(repoDir)
    self.tableData[2]["fields"].append(repoStat.branch)
    self.tableData[3]["fields"].append(repoStat.trackingBranch)
    self.tableData[4]["fields"].append(convertZeroStrToEmpty(repoStat.numCommits))
    self.tableData[5]["fields"].append(convertZeroStrToEmpty(repoStat.numModified))
    self.tableData[6]["fields"].append(convertZeroStrToEmpty(repoStat.numUntracked))
    self.ID += 1

  def getTableData(self):
    return self.tableData


#
# Run the script
#

if __name__ == '__main__':

  (options, nativeCmnd, otherArgs, extraReposFullList, notExtraReposList) = \
    getCommandlineOps()

  if nativeCmnd == "dist-repo-status":
    distRepoStatus = True
    if len(otherArgs) > 0:
      print "Error, passing in extra git commands/args =" \
        "'"+" ".join(otherArgs)+"' with special comamnd 'dist-repo-status" \
        " is not allowed!"
      sys.exit(1)
  else:
    distRepoStatus = False

  # Get the repo version files
  repoVersionDict = getRepoVersionDictFromRepoVersionFile(options.versionFile)
  repoVersionDict2 = getRepoVersionDictFromRepoVersionFile(options.versionFile2)

  # Reform the commandline arguments correctly
  #print "otherArgs =", otherArgs
  cmndLineArgsArray = requoteCmndLineArgsIntoArray(otherArgs)

  # Get the reference base directory
  baseDir = os.getcwd()

  if options.debug:
    print "*** Using git:", options.useGit

  # Get the name of the base repo
  baseDirArray = baseDir.split("/")
  baseRepoName = baseDirArray[-1]

  repoStatTable = RepoStatTable()

  # Compute base repo stats
  if options.modifiedOnly or distRepoStatus:
    baseRepoStats = getRepoStats(options)
  else:
    baseRepoStats = None

  # See if we should process the base repo or not
  processBaseRepo = True
  if not options.processBaseRepo:
    processBaseRepo = False
  elif options.modifiedOnly and not baseRepoStats.hasLocalChanges():
    processBaseRepo = False

  # Process the base git repo
  if processBaseRepo:
    if distRepoStatus:
      repoStatTable.insertRepoStat(baseRepoName+" (Base)", baseRepoStats)
    else:
      print ""
      print "*** Base Git Repo: "+addColorToRepoDir(options.useColor, baseRepoName)
      if options.debug:
        print "*** Tracking branch for git repo" \
          " '"+baseRepoName+"' = '"+baseRepoStats.trackingBranch+"'"
      sys.stdout.flush()
      runRepoCmnd(options, cmndLineArgsArray, baseRepoName, baseDir,
        repoVersionDict, repoVersionDict2)

  for extraRepo in extraReposFullList:

    # Determine if we should process this extra repo
    processThisExtraRepo = True
    if not repoExistsAndNotExcluded(options, extraRepo, notExtraReposList):
      processThisExtraRepo = False
    if processThisExtraRepo:
      # cd into extrarepo dir
      if options.debug:
        print "\n*** Changing to directory "+extraRepo,
      os.chdir(extraRepo)
      # Get repo stats
      if options.modifiedOnly or distRepoStatus:
        extraRepoStats = getRepoStats(options)
      else:
        extraRepoStats = None
      # See if we should process based on --dist-mod-only
      if options.modifiedOnly and not extraRepoStats.hasLocalChanges():
         processThisExtraRepo = False

    # Process the extra repo
    if processThisExtraRepo:
      if distRepoStatus:
        repoStatTable.insertRepoStat(extraRepo, extraRepoStats)
        processThisExtraRepo = False
      else:
        print ""
        print "*** Git Repo: "+addColorToRepoDir(options.useColor, extraRepo)
        sys.stdout.flush()
        if options.debug:
          print "*** Tracking branch for git repo" \
           " '"+extraRepo+"' = '"+extraRepoStats.trackingBranch+"'"
        runRepoCmnd(options, cmndLineArgsArray, extraRepo, baseDir, \
          repoVersionDict, repoVersionDict2)
        if options.debug:
          print "*** Changing to directory "+baseDir

    os.chdir(baseDir)

  if distRepoStatus:
    print createAsciiTable(repoStatTable.getTableData())
    if options.printLegend:
      print distRepoStatusLegend
    else:
      print "(tip: to see a legend, pass in --dist-legend.)" 
  else:
    print ""

  sys.stdout.flush()

