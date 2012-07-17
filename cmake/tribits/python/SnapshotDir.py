# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
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

This tool snapshots the contents of an origin directory (origDir) to
destination directory (destDir) and creates linkages between the two git repos
in the commit messages.

To demonstrate how this script is used, consider the desire to snapshot from:

  orig-dir/

to    

  dest-dir/

Here, the directories can be any two directories with any names as long as
they are given a final '/' at the end.  Otherwise, if you are missing the
final '/', then rsync will copy the contents from orig-dir into a subdir of
dest-dir which is usually not what you want.

The typical case is to have snapshot-dir.py soft linked into orig-dir/ to
allow a simple sync prodess.  This will be the assumption in the short text
below.  The linked-in location of snapshot-dir.py gives the default orig-dir
(but can be overridden with --orig-dir option).

The way to run this script would be:

  $ cd dest-dir
  $ orig-dir/snapshot-dir.py

By default, this assumes that git repos are used on for both the orig-dir and
dest-dir locations.  The info about the origin of the snapshot is recorded to
provide tracability for the versions.

NOTES:

* This script allows the syncing between base git repos or subdirs within git
  repos.  This is allows as the rsync command is told to ignore the .git/
  directory.
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
      help="Original directory that is the source for the snapshoted directory." \
      +"  Note that it is important to add a final /' to the directory name." \
      +"  The default is the directory where this script lives (or is softlinked)." \
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
  originRepoUrl = getGitRepoUrl(inOptions.origDir, "origin")
  print "origin URL = '"+originRepoUrl+"'"

  # Get the last commit message
  originLastCommitMsg = getLastCommitMsg(inOptions.origDir)
  print "\norigin commit message:\n\n" + originLastCommitMsg + "\n"

  #
  print "\nD) Run rsync to add and remove files and dirs between two directories\n"
  #

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

  #
  print "\nE) Create a new commit in destination directory [optional]"
  #

  echoRunSysCmnd(
    "git add .",
    workingDir=inOptions.destDir
    )

  commitMessage = \
    "Automatic snapshot commit\n"+\
    "\n"+\
    "origin: '"+originRepoUrl+"'\n"+\
    "\n"+\
    "At commit:\n"+\
    "\n"+\
    originLastCommitMsg

  echoRunSysCmnd(
    "git commit -m \""+commitMessage+"\" -- .",
    workingDir=inOptions.destDir
    )

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



def getGitRepoUrl(gitDir, repoRemoteName):

  remoteRepoUrl = ""

  # Get the list of remote repos
  remoteReposListStr = getCmndOutput("git remote -v", workingDir=gitDir)
  #print "remoteReposListStr =", remoteReposListStr

  # Loop through looking for repoRemoteName
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
    if repoName == 'origin':
      remoteRepoUrl = repoUrl
      break

  # end for

  return remoteRepoUrl


def getLastCommitMsg(gitDir):
  return getCmndOutput(
    "git log " \
    +" --pretty=format:'%h \'%s\'%nAuthor: %an <%ae>%nDate: %ad%n'" \
    +" -1 -- .",
    workingDir=gitDir
    )
