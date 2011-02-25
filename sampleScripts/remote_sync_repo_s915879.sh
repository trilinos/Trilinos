#!/bin/bash
#
# Usage: ./remote_sync_repo.sh to|from [subdir] [--no-op]
#
# Script that I use to sync between local Trilinos repo on cygwin
# laptop and remote repo on Linux workstation brain.sandia.gov (see
# notes below).
#
# This must be run from the base directory of ./Trilinos
#

EXTRA_ARGS=$@

REPO_NAME="Trilinos"
REMOTE_BASE_DIR="brain:~/PROJECTS/Trilinos.laptop.mirror.base"

./Trilinos/commonTools/python/remote_sync_dir.py \
--dir-name=Trilinos \
--remote-base-dir=brain:~/PROJECTS/Trilinos.laptop.mirror.base \
$EXTRA_ARGS

#
# 2011/02/25
# 
# Here I am using rsync to sync the Trilinos directory tree between
# ./Trilinos and a temp repo on brain.sandia.gov using the
# remote_sync_repo.sh script.  This is to avoid expensive git operations
# that I just can't seem to afford under cygwin.
# 
# The idea is to make changes to the source either only here under
# ./Trilinos on in the remote copy of Trilinos on brain.sandia.gov and
# then use the script remote_sync_repo.sh to copy files back and forth.
# All of the git operations would be performed on brain.sandia.gov where
# they will be super fast.
# 
# Once finished making modifications locally here you would do:
# 
#   $ ./remote_sync_repo.sh to [subdir]
# 
# Then you go to the remote repo on brain.sandia.gov and do git stuff,
# make changes etc.  Once finished making changes on the remote repo you
# sync back from this machine as:
# 
#   $ ./remote_sync_repo.sh from [subdir]
# 
# Note that on the local repo, you can run 'eg diff -- FILE_NAME' and
# that will run pretty quickly.  It is just global git operations that
# take forever.
# 
