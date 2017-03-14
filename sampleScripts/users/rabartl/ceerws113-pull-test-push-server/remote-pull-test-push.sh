#!/bin/bash -e
#
# Invoke a remote pull/test/push from the client machine $localmachine (see
# below).
#
# Call as:
#
#   ./remote-pull-test-push.sh [<localbranch>]
#
# If <localbranch> is not specified, then 'develop' is used.
#
# This script can be called from any directory but the log file
# remote-pull-test-push.log will get generated in the current working
# directory.  Therefore, one should call this script not from inside of the
# Trilinos local git repo.
#

#
# A) Change these for your local and remote machine names
#

# Name of git remote for local machine on remote machine.  (NOTE: This can be
# the name of any git remote really)
localmachine=crf450

# Name of the remote machine w.r.t. SSH.
remotemachine=ceerws1113

# The rest of this will work for any setup, no need to change!

localbranch=$1
if [ "$localbranch" == "" ] ; then
  localbranch=develop
fi

echo "Doing non-blocking remote pull/tets/push on $remotemachine in background.  See log file remote-pull-test-push.log!"

ssh -q $remotemachine \
  "~/remote-pull-test-push-server.sh --extra-pull-from=$localmachine:$localbranch" \
  &> remote-pull-test-push.log &
