#!/bin/bash -e
#
# Invoke a remote pull/test/push from the client machine (see below) using the
# remote server <remote_pull_test_push_server>.
#
# Call as:
#
#  cd <base_dir>/
#  ./Trilinos/sampleScripts/remote-pull-test-push/cee-lan/remote-pull-test-push.sh \
#     <remote_pull_test_push_server>
#
# where <base_dir> is the base directory that the Trilinos git repos exists under:
#
#  <base_dir>/
#    Trilinos/
#

# Get constants and inputs

basedir=$PWD

remote_pull_test_push_sever_base_dir=/scratch/$USER/TRILINOS_PUSH_SERVER

remote_pull_test_push_server=$1
if [ "$remote_pull_test_push_server" == "" ] ; then
  echo "Error: Must pass in remote_pull_test_push_server!"
  exit 1
fi

echo
echo "**************************************************************************"
echo "*** Running remote pull/test/push using server $remote_pull_test_push_server"

cd $basedir/Trilinos/
local_branch_name=`git rev-parse --abbrev-ref HEAD`

echo
echo "***"
echo "*** 1) Force push local branch '$local_branch_name' to 'intermediate-repo'"
echo "***"
echo

cd $basedir/Trilinos/
git push -f intermediate-repo $local_branch_name

echo
echo "***"
echo "*** 2) Hard reset the 'develop' branch in Trilinos repo on $remote_pull_test_push_server"
echo "***"
echo

ssh -q $remote_pull_test_push_server \
  "cd $remote_pull_test_push_sever_base_dir/Trilinos && git checkout develop && git reset --hard @{u}" 

echo
echo "***"
echo "*** 3) Doing non-blocking remote pull/tets/push on $remote_pull_test_push_server ...  See log file checkin-test-remote.log!"
echo "***"


# E) Run remote checkin-test.py script non-blocking 

cd $basedir

ssh -q $remote_pull_test_push_server \
  "cd $remote_pull_test_push_sever_base_dir/Trilinos/CHECKIN && ./checkin-test-sems.sh --extra-pull-from=intermediate-repo:$local_branch_name --do-all --no-rebase --push" &> checkin-test-remote.log &

echo
echo "You may now keep working on your local machine and wait for email notifications!"
