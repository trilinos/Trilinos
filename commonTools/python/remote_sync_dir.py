#!/usr/bin/env python
#
# Usage: remote_sync_repo.py [options] to|from [subdir]

import sys, os
from optparse import OptionParser

#
# Parse the commandline
#

description = r"""
This program is used to sync back and forth between two directory trees
using rsync on Unix/Linux systems.

The options --dir-name=DIRNAME and --remote-base-dir=REMOTEBASEDIR
specify the base directories on the local and remote machine to sync
between.  This uses rsync (use --no-opt to see the options) to sync
the directory trees preserving important properites between each set
of source trees.

The way this script is be used is to run it from the local machine in
the base directory of DIRNAME (i.e. where './DIRNAME' exists) and then
use this to move changes back and forth.

To sync the entire local directory tree to the remote machine do:

  $ cd LOCALBASEDIR
  $ remote_sync_dir --dir-name=DIRNAME --remote-base-dir=REMOTEBASEDIR to

You can then modify the files on the remote location.  Just be careful
not to touch the local files while you are making modifications on the
remote copy.

Then to sync back to the local directory tree do:

  $ cd LOCALBASEDIR
  $ remote_sync_dir --dir-name=DIRNAME --remote-base-dir=REMOTEBASEDIR from

After you have synced changed back to the local directory tree, you
have to be careful not to change things on the remote machine.

Note that you can sync back changes to just a single file or
subdirectory of DIRNAME by adding a 'subdir' option such as with:

  $ cd LOCALBASEDIR
  $ remote_sync_dir ... to someSubDir

and 

  $ cd LOCALBASEDIR
  $ remote_sync_dir ... from someSubDir

NOTE: This script was written to avoid problems with running git on
slow machines (e.g. cygwin laptops).  This is not ideal but it works
if you are very careful.
"""

parser = OptionParser(usage="usage: %prog [options] to|from [subdir]\n"+description,
 version="1.0" )
parser.add_option("--dir-name", dest="dirName", type="string",
  help="Name of the directory (or file) to sync back and forth (e.g. Trilinos)",
  default="")
parser.add_option("--remote-base-dir", dest="remoteBaseDir", type="string",
  help="Name of remote directory (e.g. somemachine:~/base)",
  default="")
parser.add_option("--no-op", dest="noOp", action="store_true",
  help="Don't do anything, just show the command to be run.",
  default=False)
(options, args) = parser.parse_args()

dirName = options.dirName
assert(dirName)
remoteBaseDir = options.remoteBaseDir
assert(remoteBaseDir)

assert(len(args) >= 1)
direction = args[0]

subdir = ""
to_subdir = ""
if len(args) == 2:
  subdir = args[1]
  to_subdir = "/".join(subdir.split("/")[1:])

if  direction == "from":
  if subdir:
    src_dest=remoteBaseDir+"/"+dirName+"/"+subdir+" "+dirName+"/"+to_subdir
  else:
    src_dest=remoteBaseDir+"/"+dirName+" ."
elif direction == "to":
  if subdir:
    src_dest=dirName+"/"+subdir+" "+remoteBaseDir+"/"+dirName+"/"+to_subdir
  else:
    src_dest=dirName+" "+remoteBaseDir
else:
  print "Error, direction (first argument) must be 'to' or 'from'"
  sys.exit(1)

cmnd = "rsync -e ssh -rlptcvz "+src_dest
print "Running: "+cmnd
if not options.noOp:
  os.system(cmnd)
