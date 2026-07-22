# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER
#


#
# Imports
# 

import sys
import time
import pprint

pp = pprint.PrettyPrinter(indent=4)

from GeneralScriptSupport import *
from Python2and3 import b, s


#
# Class for defining default options
#

class DefaultOptions:

  def __init__(self):
    self.origDir = ""
    self.destDir = ""

  def setDefaultDefaults(self):
    self.origDir = addTrailingSlashToPath(os.path.abspath(os.path.dirname(sys.argv[0])))
    self.destDir = addTrailingSlashToPath(os.getcwd())

  def setDefaultOrigDir(self, origDirIn):
    self.origDir = origDirIn

  def getDefaultOrigDir(self):
    return self.origDir

  def setDefaultDestDir(self, destDirIn):
    self.destDir = destDirIn

  def getDefaultDestDir(self):
    return self.destDir


#
# Main help string
#

usageHelp = r"""

This tool uses rsync and some git commands to snapshot the contents of an
origin directory ('orig-dir') in one Git repo to destination directory
('dest-dir') in another Git repo and logs version information in a new Git
commit message.

WARNING: Because this tool uses 'rsync --delete', it can be quite destructive
to the contents of dest-dir if one is not careful and selects the wrong
orig-dir or wrong dest-dir!  Therefore, please carefully read this entire help
message and especially the warnings given below and always start by running
with --no-op!

To sync between any two arbitrary directories between two different git repos
invoking this script from any directory location, one can do:

  $ <some-base-dir>/snapshot-dir.py \
    --orig-dir=<dest-dir> \
    --dest-dir=<dest-dir> \
    [--no-op]

(where --no-op should be used on the initial trial run to check that the
correct rsync and other commands will be run, see below).

To describe how this script is used, consider the desire to snapshot the
directory tree from one git repo:

  <some-orig-base-dir>/orig-dir/

and exactly duplicate it in another git repo under:

  <some-dest-base-dir>/dest-dir/

Here, the directories can be any two directories from local git repos with any
names. (Note if the paths don't end with '/', then '/' will be added.
Otherwise, rsync will copy the contents from 'orig-dir' into a subdir of
'dest-dir' which is usually not what you want.)

A typical case is to have snapshot-dir.py soft linked into orig-dir/ to allow
a simple sync process.  The linked-in location of snapshot-dir.py gives the
default 'orig-dir' directory automatically (but can be overridden with
--orig-dir=<orig-dir> option).

When snapshot-dir.py is soft-linked into the 'orig-dir' directory base, the
way to run this script would be:

  $ cd <some-dest-base-dir>/dest-dir/
  $ <some-orig-base-dir>/orig-dir/snapshot-dir.py

By default, this assumes that git repos are used for both the 'orig-dir' and
'dest-dir' locations.  The info about the origin of the snapshot from
'orig-dir' is recorded in the commit message of the 'dest-dir' git repo to
provide tractability for the versions (see below).

By default, this script does the following:

1) Assert that 'orig-dir' in its git repo is clean (i.e. no uncommitted
   files).  (Can be disabled by passing in --allow-dirty-orig-dir.)

2) Assert that <dest-dir>/ in its git repo is clean (same checks as for
   'orig-dir' above).  (Can be disabled by passing in --allow-dirty-dest-dir.
   Also, this must be skipped on the initial snapshot where <some-dist-dir>/
   does not exist yet.)

3) Clean out the ignored files from <some-source-dir>/orig-dir using 'git
   clean -xdf' run in that directory.  (Only if --clean-ignored-files-orig-dir
   is passed.)

4) Run 'rsync -cav --delete [other options] <orig-dir>/ <dest-dir>/' to copy
   the contents from 'orig-dir' to 'dest-dir', excluding the '.git/' directory
   if it exists in either git repo dir.  After this runs, <dest-dir>/ should
   be an exact duplicate of <orig-dir>/ (except for otherwise noted excluded
   files).  This rsync will delete any files in 'dest-dir' that are not in
   'orig-dir'.  Note that if there are any ignored untracked files in
   'orig-dir' that get copied over, then the copied .gitignore files should
   avoid treating them as tracked files in the 'dest-dir' git repo.  (The
   rsync command is skipped if the argument --no-op is passed.)

5) Run 'git add .' in <dest-dir>/ to stage any files copied over.  (Note that
   git will automatically stage deletes for any files removed by the 'rsync
   -cav --delete' command.  Also note that any untracked, unknown or ignored
   files in 'orig-dir' that get copied over and are not ignored in 'dist-dir'
   in the copied-over '.gitignore' files, or other git ignore files in
   'dist-dir', then they will be added to the new commit in 'dist-dir'.)

6) Get the git remote URL from the orig-dir git repo, and the git log for the
   last commit for the directory orig-dir from its git repo.  (If the orig-dir
   repo is on a tracking branch, then this remote will be guaranteed to be
   correct.  However, if the orig-dir repo not on a tracking branch, then the
   first remote returned from 'git remote -v' will be used. This information
   is used to provide tractability of the version info back to the originating
   git repo in the commit created in the dest-dir repo.)

7) Commit the updated dest-dir directory using a commit message with the
   orig-dir snapshot version info.  (This will only commit files in 'dest-dir'
   and not in other directories in the destination git repo!)  (The 'git
   commit' will be skipped if the options --skip-commit or --no-op are
   passed.)

NOTES:

* When first running this tool, use the --no-op option to see what commands
  would be run without actually performing any mutable operations.  Especially
  pay attention to the rsync command that would be run and make sure it is
  operating on the desired directories.

* On the first creation of <dest-dir>/, one must pass in
  --allow-dirty-dest-dir to avoid checks of <dest-dir>/ (because it does not yet
  exist).

* This script allows the syncing between base git repos or subdirs within git
  repos.  This is allowed because the rsync command is told to ignore the
  .git/ directory when syncing.

* The cleaning of orig-dir/ using 'git clean -xdf' may be somewhat dangerous
  but it is recommended by passing in --clean-ignored-files-orig-dir to avoid
  copying locally-ignored files in orig-dir/ (e.g. ignored in
  .git/info/excludes and not ignored in a committed .gitignore file in
  orig-dir/) that get copied to and then committed in the dest-dir/ repo.
  Therefore, be sure you don't have any of these type of ignored files in
  orig-dir/ that you want to keep before you run this tool with the option
  --clean-ignored-files-orig-dir!

* Snapshotting with this script will create an exact duplicate of 'orig-dir'
  in 'dest-dir' and therefore if there are any local changes to the files or
  changes after the last snapshot, they will get wiped out.  To avoid this,
  one can the snapshot on a branch in the 'dest-dir' git repo, then merge that
  branch into the main branch (e.g. 'master') in 'dest-dir' repo.  As long as
  there are no merge conflicts, this will preserve local changes for the
  mirrored directories and files.  This strategy can work well as a way to
  allow for local modifications but still do snapshotting.

WARNINGS:

* Make sure that orig-dir is a non-ignored subdir of the origin git repo and
  that it does not contain any ignored subdirs and files that are ignored by
  listing them in the .git/info/exclude file instead of .gitignore files that
  would get copied over to dist-dir.  (If ignored files are ignored using
  .gitignore files within orig-dir, then those files will be copied by the
  rsync command along with the rest of the files from orig-dir and those files
  will also be ignored after being copied to dest-dir.  However, if these
  files and dirs are ignored because of being listed in the .git/info/exclude
  file in orig-dir, then those files will not be ignored when copied over to
  dest-dir and will therefore be added to the git commit in the dist-dir git
  repo.  That is why it is recommended to run with
  --clean-ignored-files-orig-dir to ensure that all ignore files are removed
  from orig-dir before doing the rsync.)

* Make sure that the top commit in orig-dir has been pushed (or will be
  pushed) to the listed remote git repo showed in the output line 'origin
  remote name' or 'origin remote URL'.  Otherwise, others will not be able to
  trace that exact version by cloning that repo and traceability is lost.

* Make sure that dest-dir is a non-ignored subdir of the destination git repo
  and does not contain any ignored subdirs or files that you don't care if
  they are deleted.  (Otherwise, the 'rsync --delete' command will delete any
  files in dest-dir that are not also in orig-dir and since these non-ignored
  files in dest-dir are not under version control, they will not be
  recoverable.)

* As a corollary to the above warnings about ignored files and directories, do
  not snapshot from an orig-dir or to a dest-dir that contain build
  directories or other generated files that you want to keep!  Even if those
  files and directories are ignored in copied-over .gitignore files, the copy
  of those ignored files will still occur. (Just another reason to keep your
  build directories completely outside of the source tree!)

"""


# Direct script driver (taking in command-line arguments)
#
# This runs given a set of command-line arguments and a stream to do
# outputting to.  This allows running an a unit testing environment very
# efficiently.
#
def snapshotDirMainDriver(cmndLineArgs, defaultOptionsIn = None, stdout = None):

  oldstdout = sys.stdout

  try:
  
    if stdout:
      sys.stdout = stdout

    if defaultOptionsIn:
      defaultOptions = defaultOptionsIn
    else:
      defaultOptions = DefaultOptions()
      defaultOptions.setDefaultDefaults()

    #print("cmndLineArgs = " + str(cmndLineArgs))
  
    #
    # A) Get the command-line options
    #
  
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
  
    clp = ArgumentParser(
      description=usageHelp,
      formatter_class=RawDescriptionHelpFormatter)
  
    clp.add_argument(
      "--show-defaults", dest="showDefaults", action="store_true",
      help="Show the default option values and do nothing at all.",
      default=False )

    clp.add_argument(
      "--orig-dir", dest="origDir",
      default=defaultOptions.getDefaultOrigDir(),
      help="Original directory that is the source for the snapshotted directory." \
      +"  If a trailing '/' is missing then it will be added." \
      +"  The default is the directory where this script lives (or is soft-linked)." \
      +"  [default: '"+defaultOptions.getDefaultOrigDir()+"']")

    clp.add_argument(
      "--dest-dir", dest="destDir",
      default=defaultOptions.getDefaultDestDir(),
      help="Destination directory that is the target for the snapshoted directory." \
      +"  If a trailing '/' is missing then it will be added." \
      +"  The default dest-dir is current working directory." \
      +"  [default: '"+defaultOptions.getDefaultDestDir()+"']" \
      )

    clp.add_argument(
      "--exclude", dest="exclude", default=None, nargs="*",
      help="List of files/directories/globs to exclude from orig-dir when " \
      +"snapshotting.")

    clp.add_argument(
      "--assert-clean-orig-dir", dest="assertCleanOrigDir", action="store_true",
      default=True,
      help="Check that orig-dir is committed and clean. [default]" )
    clp.add_argument(
      "--allow-dirty-orig-dir", dest="assertCleanOrigDir", action="store_false",
      help="Skip clean check of orig-dir." )

    clp.add_argument(
      "--assert-clean-dest-dir", dest="assertCleanDestDir", action="store_true",
      default=True,
      help="Check that dest-dir is committed and clean. [default]" )
    clp.add_argument(
      "--allow-dirty-dest-dir", dest="assertCleanDestDir", action="store_false",
      help="Skip clean check of dest-dir." )

    clp.add_argument(
      "--clean-ignored-files-orig-dir", dest="cleanIgnoredFilesOrigDir", action="store_true",
      default=False,
      help="Clean out the ignored files from orig-dir/ before snapshotting." )
    clp.add_argument(
      "--no-clean-ignored-files-orig-dir", dest="cleanIgnoredFilesOrigDir", action="store_false",
      help="Do not clean out orig-dir/ ignored files before snapshotting. [default]" )

    clp.add_argument(
      "--do-commit", dest="doCommit", action="store_true",
      default=True,
      help="Actually do the commit. [default]" )
    clp.add_argument(
      "--skip-commit", dest="doCommit", action="store_false",
      help="Skip the commit." )

    clp.add_argument(
      "--verify-commit", dest="noVerifyCommit", action="store_false",
      default=False,
      help="Do not pass --no-verify to git commit.  [default]" )
    clp.add_argument(
      "--no-verify-commit", dest="noVerifyCommit", action="store_true",
      help="Pass --no-verify to git commit." )

    clp.add_argument(
      "--no-op", dest="noOp", action="store_true",
      default=False,
      help="Don't actually run any commands that would change the state other"
      +" (other than the natural side-effects of running git query commands)" )

    options = clp.parse_args(cmndLineArgs)
  
    #
    # B) Echo the command-line
    #
  
    print("")
    print("**************************************************************************")
    print("Script: snapshot-dir.py \\")

    print("  --orig-dir='" + options.origDir + "' \\")
    print("  --dest-dir='" + options.destDir + "' \\")
    if options.exclude:
      print("  --exclude " + " ".join(options.exclude) + " \\")
    if options.assertCleanOrigDir:
      print("  --assert-clean-orig-dir \\")
    else:
      print("  --allow-dirty-orig-dir \\")
    if options.assertCleanDestDir:
      print("  --assert-clean-dest-dir \\")
    else:
      print("  --allow-dirty-dest-dir \\")
    if options.cleanIgnoredFilesOrigDir:
      print("  --clean-ignored-files-orig-dir \\")
    else:
      print("  --no-clean-ignored-files-orig-dir \\")
    if options.doCommit:
      print("  --do-commit \\")
    else:
      print("  --skip-commit \\")
    if options.noVerifyCommit:
      print("  --no-verify-commit \\")
    else:
      print("  --verify-commit \\")
    if options.noOp:
      print("  --no-op \\")
  
    if options.showDefaults:
      return  # All done!
  
    #
    # C) Execute the 
    #
  
    snapshotDir(options)
  
  finally:
    sys.stdout = oldstdout



#
# Implement the guts of snapshoting after reading in options
#

def snapshotDir(inOptions):

  addTrailingSlashToPaths(inOptions)

  #
  print("\nA) Assert that orig-dir is 100% clean with all changes committed\n")
  #

  if inOptions.assertCleanOrigDir:
   assertCleanGitDir(inOptions.origDir, "origin",
      "The created snapshot commit would not have the correct origin commit info!" )
  else:
    print("Skipping on request!")

  #
  print("\nB) Assert that dest-dir is 100% clean with all changes committed\n")
  #

  if inOptions.assertCleanDestDir:
    assertCleanGitDir(inOptions.destDir, "destination",
      "Location changes in the destination directory would be overritten and lost!")
  else:
    print("Skipping on request!")

  #
  print("\nC) Cleaning out ignored files in orig-dir\n")
  #

  if inOptions.cleanIgnoredFilesOrigDir:
    cleanIgnoredFilesFromGitDir(inOptions.origDir, inOptions.noOp, "origin")
  else:
    print("Skipping on request!")

  #
  print("\nD) Get info for git commit from orig-dir [optional]\n")
  #

  # Get the repo for origin
  (remoteRepoName, remoteBranch, remoteRepoUrl) = \
     getGitRepoRemoteNameBranchAndUrl(inOptions.origDir)
  print("origin remote name = '" + remoteRepoName + "'")
  print("origin remote branch = '" + remoteBranch + "'")
  print("origin remote URL = '" + remoteRepoUrl + "'")
  gitDescribe = getGitDescribe(inOptions.origDir)
  print("Git describe = '" + gitDescribe + "'")

  # Get the last commit message
  originLastCommitMsg = getLastCommitMsg(inOptions.origDir)
  print("\norigin commit message:")
  print("---------------------------------------")
  print(originLastCommitMsg)
  print("---------------------------------------")

  #
  print("\nE) Run rsync to add and remove files and dirs between two directories\n")
  #

  excludes = r"""--exclude=\.git"""
  if inOptions.exclude:
      excludes += " " + " ".join(map(lambda ex: "--exclude="+ex,
                                     inOptions.exclude))
      print("Excluding files/directories/globs: " +
            " ".join(inOptions.exclude))
  # Note that when syncing one git repo to another, we want to sync the
  # .gitingore and other hidden files as well.

  # When we support syncing from hg repos, add these excludes as well:
  #    --exclude=\.hg --exclude=.hgignore --exclude=.hgtags

  rsyncCmnd = \
    r"rsync -cav --delete "+excludes+" "+inOptions.origDir+" "+inOptions.destDir

  if not inOptions.noOp:
    rtn = echoRunSysCmnd(rsyncCmnd,  throwExcept=False,  timeCmnd=True)
    if rtn != 0:
      print("Rsync failed, aborting!")
      return False
  else:
    print("Would be running: "+rsyncCmnd)

  #
  print("\nE) Create a new commit in dest-dir [optional]")
  #

  origDirLast = inOptions.origDir.split("/")[-2]
  origSha1 = getCommitSha1(inOptions.origDir)

  commitMessage = \
    "Automatic snapshot commit from "+origDirLast+" at "+origSha1+"\n"+\
    "\n"

  if remoteBranch:
    commitMessage += \
      "Origin repo remote tracking branch: '"+remoteRepoName+"/"+remoteBranch+"'\n"+\
      "Origin repo remote repo URL: '"+remoteRepoName+" = "+remoteRepoUrl+"'\n"
  else:
    commitMessage += \
      "Origin repo remote repo URL: '"+remoteRepoName+" = "+remoteRepoUrl+"'\n"

  commitMessage += \
    "Git describe: "+gitDescribe+"\n" +\
    "\n"+\
    "At commit:\n"+\
    "\n"+\
    originLastCommitMsg

  print("\nGenerating commit in dest-dir with commit message:\n")
  print("---------------------------------------")
  print(commitMessage)
  print("---------------------------------------")

  if inOptions.doCommit:

    gitAddCmnd = "git add ."
    if not inOptions.noOp:
      echoRunSysCmnd(gitAddCmnd, workingDir=inOptions.destDir)
    else:
      print("\nWould be running: "+gitAddCmnd+"\n" \
        +"\n    in directory '"+inOptions.destDir+"'" )

    if inOptions.noVerifyCommit:
      noVerifyCommitArgStr = " --no-verify"
    else:
      noVerifyCommitArgStr = ""

    gitCommitCmndBegin = "git commit"+noVerifyCommitArgStr+" -m "
    if not inOptions.noOp:
      echoRunSysCmnd(gitCommitCmndBegin+"\""+commitMessage+"\"",
        workingDir=inOptions.destDir)
    else:
      print("\nWould be running: "+gitCommitCmndBegin+"\"<commit-msg>\"\n"
        +"\n    in directory '"+inOptions.destDir+"'" )

  else:

    print("\nSkipping commit on request!\n")

  if inOptions.noOp:

    print(
      "\n***\n"
      "*** NOTE: No modifying operations were performed!\n"
      "***\n"
      "*** Run again removing the option --no-op to make modifying\n"
      "*** changes.\n"
      "***\n"
      "*** But first, carefully look at the orig-dir, dest-dir and the\n"
      "*** various operations performed above to make sure that\n"
      "*** everything is as it should be before removing the option --no-op.\n"
      "***\n"
      "*** In particular,  look carefully at the 'git clean' and  'rsync' commands\n"
      "*** on the lines that begin with 'Would be running:'\n"
      "***\n"
      "***\n"
      "***\n"
      "***\n"
      )

  #
  # F) Success! (if you get this far)
  #

  return True


#
# Helper functions
#


def addTrailingSlashToPath(path):
  if path[-1] != "/":
    return path + "/"
  return path


def addTrailingSlashToPaths(options):
  options.origDir = addTrailingSlashToPath(options.origDir)
  options.destDir = addTrailingSlashToPath(options.destDir)


def assertCleanGitDir(dirPath, dirName, explanation):

  changedFiles = getCmndOutput(
    "git diff --name-status HEAD -- .",
    stripTrailingSpaces = True,
    workingDir = dirPath )

  if changedFiles:
    raise Exception(
      "Error, aborting snapshot!\n" \
      "The "+dirName+" git directory '"+dirPath+"' is not clean and" \
      +" has the changed files:\n"+changedFiles+"\n" \
      +explanation
      )
  else:
    print("The " + dirName + " git directory '" + dirPath + "' is clean!")

  # NOTE: The above git diff command will not catch unknown files but that is
  # not a huge risk for the use cases that I am concerned with.


def cleanIgnoredFilesFromGitDir(dirPath, noOp, dirName):
  gitCleanCmnd = r"git clean -xdf"
  if not noOp:
    rtn = echoRunSysCmnd(gitCleanCmnd,  workingDir=dirPath,
      throwExcept=False, timeCmnd=True)
    if rtn != 0:
      raise Exception(
        "Error, cleaning of origin `"+dirPath+"` failed!")
  else:
    print("Would be running: "+gitCleanCmnd+"\n" \
      +"\n    in directory '"+dirPath+"'" )


def getCommitSha1(gitDir):
  return getCmndOutput("git log -1 --pretty=format:'%h' -- .", workingDir=gitDir).strip()


def getGitRepoRemoteNameBranchAndUrl(gitDir):

  remoteRepoName = ""
  remoteBranch = ""
  remoteRepoUrl = ""

  # Get the remote tracking branch
  (trackingBranchStr, trackingBranchErrCode) = getCmndOutput(
     "git rev-parse --abbrev-ref --symbolic-full-name @{u}", workingDir=gitDir,
     throwOnError=False, rtnCode=True)

  if trackingBranchErrCode == 0:
    (remoteRepoName, remoteBranch) = trackingBranchStr.strip().split("/")
  else:
    remoteRepoName = ""
    remoteBranch = ""

  # Get the list of remote repos
  remoteReposListStr = getCmndOutput("git remote -v", workingDir=gitDir)
  #print("remoteReposListStr = " + remoteReposListStr)

  # Loop through looking for remoteRepoName
  for remoteRepo in remoteReposListStr.splitlines():

    #print("remoteRepo = '" + remoteRepo + "'")
    if remoteRepo == "":
      continue
    
    remoteRepoList = remoteRepo.split(" ")
    #print("remoteRepoList = " + str(remoteRepoList))

    # Remove empty items
    k = 0
    while k < len(remoteRepoList):
      if remoteRepoList[k] == "":
        del remoteRepoList[k]
      k += 1
    #print("remoteRepoList = " + str(remoteRepoList))

    # Get the remote name and URL
    (repoName, repoUrl) = remoteRepoList[0].split("\t")
    #print("repoName = '" + repoName + "'")
    #print("repoUrl  = '" + repoUrl  + "'")

    if remoteRepoName:
      # Grab the URL if the remote name matches
      if repoName == remoteRepoName:
        remoteRepoUrl = repoUrl
        break
    else:
      # Just grab the first remote name you find if there is no tracking branch
      remoteRepoName = repoName
      remoteRepoUrl = repoUrl
      break

  # end for

  return (remoteRepoName, remoteBranch, remoteRepoUrl)


def getGitDescribe(gitDir):

  gitDescribe = getCmndOutput( "git describe", workingDir=gitDir, stripTrailingSpaces=True)

  return gitDescribe


def getLastCommitMsg(gitDir):
  return getCmndOutput(
    "git log " \
    +" --pretty=format:'commit %H%nAuthor:  %an <%ae>%nDate:    %ad%nSummary: %s%n'" \
    +" -1 -- .",
    workingDir=gitDir
    )

#  LocalWords:  traceability TriBITS Snapshotting snapshotting
