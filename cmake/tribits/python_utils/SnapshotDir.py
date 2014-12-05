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


#
# Imports
# 

import sys
import time
import pprint

pp = pprint.PrettyPrinter(indent=4)

from GeneralScriptSupport import *


#
# Class for defining default options
#

class DefaultOptions:

  def __init__(self):
    self.origDir = ""
    self.destDir = ""

  def setDefaultDefaults(self):
    self.origDir = os.path.abspath(os.path.dirname(sys.argv[0])) + "/"
    self.destDir = os.getcwd() + "/"

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

usageHelp = r"""snapshot-dir.py [OPTIONS]

This tool snapshots the contents of an origin directory ('orig-dir') to
destination directory ('dest-dir') and creates linkages between the two git
repos in the commit message in the 'dest-dir' git branch.  The command 'git'
must be in the path for this script to be used.

To demonstrate how this script is used, consider the desire to snapshot the
directory tree:

  <some-orig-base-dir>/orig-dir/

and duplicate it in the directory tree

  <some-dest-base-dir>/dest-dir/

Here, the directories can be any two directories from local git repos with any
names as long as they are given a final '/' at the end.  Otherwise, if you are
missing the final '/', then rsync will copy the contents from 'orig-dir' into
a subdir of 'dest-dir' which is usually not what you want.

A typical case is to have snapshot-dir.py soft linked into orig-dir/ to allow
a simple sync process.  This is the case, for example, with the 'tribits'
source tree.  The linked-in location of snapshot-dir.py gives the default
'orig-dir' directory automatically (but can be overridden with --orig-dir
option).

When snapshot-dir.py is soft-linked into the 'orig-dir' directory base, the
way to run this script would be:

  $ cd <some-dest-base-dir>/dest-dir/
  $ <some-orig-base-dir>/orig-dir/snapshot-dir.py

By default, this assumes that git repos are used for both the 'orig-dir' and
'dest-dir' locations.  The info about the origin of the snapshot from
'orig-dir' is recorded in the commit message of the 'dest-dir' git repo to
provide tractability for the versions (see below).

To sync between any two arbitrary directories invoking this script from any
directory location, one can do:

  $ <some-base-dir>/snapshot-dir.py \
    --orig-dir=<some-orig-dir>/ \
    --dest-dir=<some-dest-dir>/

Note the trailing '/' is critical for the correct functioning of rsync.

By default, this script does the following:

1) Asserts that the git repo for 'orig-dir' is clean (i.e. no uncommitted
   files, no unknown files, etc.).  This can be disabled by passing in
   --allow-dirty-orig-dir.

2) Asserts that the git repo for <some-dest-dir>/ is clean (see above).  This
   can be disabled by passing in --allow-dirty-dest-dir.

3) Run 'rsync -av --delete' to copy the contents from 'orig-dir' to
   'dest-dir', excluding the '.git/' directory if it exists in either git repo
   dir.  After this runs, <some-dest-dir>/ should be an exact duplicate of
   <some-orig-dir>/ (except for otherwise noted excluded files).  This rsync
   will delete any files in 'dest-dir' that are not in 'orig-dir'.  Note that
   if there are ignored untracked files, then the copied .gitignore files
   should avoid showing them as tracked or unknown files in the 'dest-dir' git
   repo as well.

4) Run 'git add .' in <some-dest-dir>/ to stage any new files.  Note that git
   will automatically stage deletes for any files removed by the 'rsync -av
   --delete' command!

5) Get the git remote URL from the orig-dir git repo, and the git log for the
   last commit for the directory from orig-dir.  This information is used to
   define perfect tracing of the version info when doing the snapshot.

6) Commit the updated dest-dir directory using a commit message with the
   orig-dir repo URL and log info.  This will only commit files in 'dest-dir'
   and not in other directories in the destination git repo!

NOTES:

* This script allows the syncing between base git repos or subdirs within git
  repos.  This is allowed because the rsync command is told to ignore the
  .git/ directory when syncing.

* Snapshotting with this script will create an exact duplicate of 'orig-dir'
  in 'dest-dir' and therefore if there are any local changes to the files or
  chagnes after the last snapshot, they will get wiped out.  To avoid this, do
  the snapshot on a branch in the 'dest-dir' git repo, then merge that branch
  into the master branch in 'dest-dir' repo that has the local changes.  As
  long as there are no merge conflicts, this will preserve local changes for
  the mirrored directories and files.  This works amazingly well in most
  cases.
"""


#
# Direct script driver (taking in command-line arguments)
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

    #print "cmndLineArgs =", cmndLineArgs
  
    #
    # A) Get the command-line options
    #
  
    from optparse import OptionParser
  
    clp = OptionParser(usage=usageHelp)
  
    clp.add_option(
      "--show-defaults", dest="showDefaults", action="store_true",
      help="Show the default option values and do nothing at all.",
      default=False )

    clp.add_option(
      "--orig-dir", dest="origDir", type="string",
      default=defaultOptions.getDefaultOrigDir(),
      help="Original directory that is the source for the snapshotted directory." \
      +"  Note that it is important to add a final /' to the directory name." \
      +"  The default is the directory where this script lives (or is soft-linked)." \
      +"  [default: '"+defaultOptions.getDefaultOrigDir()+"']")

    clp.add_option(
      "--dest-dir", dest="destDir", type="string",
      default=defaultOptions.getDefaultDestDir(),
      help="Destination directory that is the target for the snapshoted directory." \
      +"  Note that a final '/' just be added or the origin will be added as subdir." \
      +"  The default dest-dir is current working directory." \
      +"  [default: '"+defaultOptions.getDefaultDestDir()+"']" \
      )

    clp.add_option(
      "--assert-clean-orig-dir", dest="assertCleanOrigDir", action="store_true",
      help="Check that orig-dir is committed and clean. [default]" )
    clp.add_option(
      "--allow-dirty-orig-dir", dest="assertCleanOrigDir", action="store_false",
      default=True,
      help="Skip clean check of orig-dir." )

    clp.add_option(
      "--assert-clean-dest-dir", dest="assertCleanDestDir", action="store_true",
      help="Check that dest-dir is committed and clean. [default]" )
    clp.add_option(
      "--allow-dirty-dest-dir", dest="assertCleanDestDir", action="store_false",
      default=True,
      help="Skip clean check of dest-dir." )

    clp.add_option(
      "--do-rsync", dest="doRsync", action="store_true",
      help="Actually do the rsync. [default]" )
    clp.add_option(
      "--skip-rsync", dest="doRsync", action="store_false",
      default=True,
      help="Skip the rsync (testing only?)." )

    clp.add_option(
      "--do-commit", dest="doCommit", action="store_true",
      help="Actually do the commit. [default]" )
    clp.add_option(
      "--skip-commit", dest="doCommit", action="store_false",
      default=True,
      help="Skip the commit." )
    
    (options, args) = clp.parse_args(cmndLineArgs)
  
    #
    # B) Echo the command-line
    #
  
    print ""
    print "**************************************************************************"
    print "Script: snapshot-dir.py \\"

    print "  --orig-dir='"+options.origDir+"' \\"
    print "  --dest-dir='"+options.destDir+"' \\"
    if options.assertCleanOrigDir:
      print "  --assert-clean-orig-dir \\"
    else:
      print "  --allow-dirty-orig-dir \\"
    if options.assertCleanDestDir:
      print "  --assert-clean-dest-dir \\"
    else:
      print "  --allow-dirty-dest-dir \\"
    if options.doCommit:
      print "  --do-commit \\"
    else:
      print "  --skip-commit \\"
    if options.doRsync:
      print "  --do-rsync \\"
    else:
      print "  --skip-rsync \\"
  
    if options.showDefaults:
      return  # All done!
  
    #
    # C) Exectute the 
    #
  
    snapshotDir(options)
  
  finally:
    sys.stdout = oldstdout



#
# Implement the guts of snapshoting after reading in options
#

def snapshotDir(inOptions):

  #
  print "\nA) Assert that orig-dir is 100% clean with all changes committed\n"
  #

  if inOptions.assertCleanOrigDir:
   assertCleanGitDir(inOptions.origDir, "origin",
      "The created snapshot commit would not have the correct origin commit info!" )
  else:
    print "Skipping on request!"

  #
  print "\nB) Assert that dest-dir is 100% clean with all changes committed\n"
  #

  if inOptions.assertCleanDestDir:
    assertCleanGitDir(inOptions.destDir, "destination",
      "Location changes in the destination directory would be overritten and lost!")
  else:
    print "Skipping on request!"

  #
  print "\nC) Get info from git commit from origDir [optional]\n"
  #

  # Get the repo for origin
  (remoteRepoName, remoteBranch, remoteRepoUrl) = \
     getGitRepoUrl(inOptions.origDir)
  print "origin remote name = '"+remoteRepoName+"'"
  print "origin remote branch = '"+remoteBranch+"'"
  print "origin remote URL = '"+remoteRepoUrl+"'"

  # Get the last commit message
  originLastCommitMsg = getLastCommitMsg(inOptions.origDir)
  print "\norigin commit message:\n\n" + originLastCommitMsg + "\n"

  #
  print "\nD) Run rsync to add and remove files and dirs between two directories\n"
  #

  if inOptions.doRsync:

    excludes = r"""--exclude=\.git"""
    # Note that when syncing one git repo to another, we want to sync the
    # .gitingore and other hidden files as well.
  
    # When we support syncing from hg repos, add these excludes as well:
    #    --exclude=\.hg --exclude=.hgignore --exclude=.hgtags
  
    rtn = echoRunSysCmnd(
      r"rsync -av --delete "+excludes+" "+inOptions.origDir+" "+inOptions.destDir,
      throwExcept=False,
      timeCmnd=True
      )
  
    if rtn != 0:
      "Rsync failed, aborting!"
      return False

  else:

    print "\nSkipping rsync on request!"

  #
  print "\nE) Create a new commit in destination directory [optional]"
  #

  origDirLast = inOptions.origDir.split("/")[-2]
  origSha1 = getCommitSha1(inOptions.origDir)

  commitMessage = \
    "Automatic snapshot commit from "+origDirLast+" at "+origSha1+"\n"+\
    "\n"+\
    "Origin repo remote tracking branch: '"+remoteRepoName+"/"+remoteBranch+"'\n"+\
    "Origin repo remote repo URL: '"+remoteRepoName+" = "+remoteRepoUrl+"'\n"+\
    "\n"+\
    "At commit:\n"+\
    "\n"+\
    originLastCommitMsg

  print "\nGeneratting commit with commit message:\n"
  print "---------------------------------------"
  print commitMessage
  print "---------------------------------------"

  if inOptions.doCommit:

    echoRunSysCmnd(
      "git add .",
      workingDir=inOptions.destDir
      )
  
    echoRunSysCmnd(
      "git commit -m \""+commitMessage+"\" -- .",
      workingDir=inOptions.destDir
      )

  else:

    print "\nSkipping commit on request!\n"

  #
  # F) Success! (if you get this far)
  #

  return True


#
# Helper functions
#


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
    print "The "+dirName+" git directory '"+dirPath+"' is clean!"

  # NOTE: The above git diff command will not catch unknown files but that is
  # not a huge risk for the use cases that I am concerned with.


def getCommitSha1(gitDir):
  return getCmndOutput("git log -1 --pretty=format:'%h'", workingDir=gitDir).strip()


def getGitRepoUrl(gitDir):

  remoteRepoName = ""
  remoteBranch = ""
  remoteRepoUrl = ""

  # Get the remote tracking branch
  trackingBranchStr = getCmndOutput(
     "git rev-parse --abbrev-ref --symbolic-full-name @{u}", workingDir=gitDir)

  (remoteRepoName, remoteBranch) = trackingBranchStr.strip().split("/")

  # Get the list of remote repos
  remoteReposListStr = getCmndOutput("git remote -v", workingDir=gitDir)
  #print "remoteReposListStr =", remoteReposListStr

  # Loop through looking for remoteRepoName
  for remoteRepo in remoteReposListStr.split("\n"):

    #print "remoteRepo = '"+remoteRepo+"'"
    if remoteRepo == "":
      continue
    
    remoteRepoList = remoteRepo.split(" ")
    #print "remoteRepoList =", remoteRepoList

    # Remove empty items
    k = 0
    while k < len(remoteRepoList):
      if remoteRepoList[k] == "":
        del remoteRepoList[k]
      k += 1
    #print "remoteRepoList =", remoteRepoList

    # Get the remote name and URL
    (repoName, repoUrl) = remoteRepoList[0].split("\t")
    #print "repoName = '"+repoName+"'"
    #print "repoUrl = '"+repoUrl+"'"

    # Grab the URL if the remote name matches
    if repoName == remoteRepoName:
      remoteRepoUrl = repoUrl
      break

  # end for

  return (remoteRepoName, remoteBranch, remoteRepoUrl)


def getLastCommitMsg(gitDir):
  return getCmndOutput(
    "git log " \
    +" --pretty=format:'%h \'%s\'%nAuthor: %an <%ae>%nDate: %ad%n'" \
    +" -1 -- .",
    workingDir=gitDir
    )
