#!/bin/bash -e
#
# Invoke a remote pull/test/push from the client machine (see below) using the
# remote server <remote_pull_test_push_server> with Trilinos soruce tree on
# the remote at <remote_trilinos_base_dir>/Trilinos.
#
# Call as:
#
#  $ cd <local_trilinos_base_dir>/
#  $ ./Trilinos/cmake/std/sems/remote-pull-test-push.sh \
#       <remote_pull_test_push_server>  <remote_trilinos_base_dir>
#
# where <local_trilinos_base_dir> is the base directory for Trilinos on the
# local machine (i.e. <local_trilinos_base_dir>/Trilinos).
#
 
echo
echo "**************************************************************************"
echo "*** Running remote pull/test/push using server $remote_pull_test_push_server"
echo

local_trilinos_base_dir=$PWD

remote_pull_test_push_server=$1
if [ "$remote_pull_test_push_server" == "" ] ; then
  echo "Error: Must pass in remote_pull_test_push_server!"
  exit 1
fi
echo "Remote server: '$remote_pull_test_push_server'"

remote_trilinos_base_dir=$2
if [ "$remote_trilinos_base_dir" == "" ] ; then
  echo "Error: Must pass in remote_trilinos_base_dir!"
  exit 2
fi
echo "Remote base dir: '$remote_trilinos_base_dir'"

blocking_or_nonblocking=$3
if [ "$blocking_or_nonblocking" == "" ] ; then
  # Default is nonblocking
  blocking_or_nonblocking=nonblocking
elif [ "$blocking_or_nonblocking" == "blocking" ] ; then
  # Valid value!
  blocking_or_nonblocking=blocking
elif [ "$blocking_or_nonblocking" == "nonblocking" ] ; then
  # Valid value!
  blocking_or_nonblocking=nonblocking
else
  # Invalid value!
  echo "Error: Third argument value '$blocking_or_nonblocking' is not acceptable!  Must pass in 'blocking' or 'nonblocking'!"
  exit 3
fi
echo "Blocking or nonblocking: '$blocking_or_nonblocking'"

cd $local_trilinos_base_dir/Trilinos/
local_branch_name=`git rev-parse --abbrev-ref HEAD`

echo
echo "***"
echo "*** 1) Force push local branch '$local_branch_name' to 'intermediate-repo'"
echo "***"
echo

cd $local_trilinos_base_dir/Trilinos/
git push -f intermediate-repo $local_branch_name:$local_branch_name

#
# Run the remote commands blocking or nonblocking mode
#

remote_branch_update_cmnds="cd $remote_trilinos_base_dir/Trilinos && git checkout develop && git reset --hard @{u} && git fetch intermediate-repo && git merge --no-ff intermediate-repo/$local_branch_name"
# NOTE: Above, we do the fetch and merge of the branch before running the
# checkin-test.py script on the remote machine.  That is needed for some use
# cases where the checkin-test.py needs to be run given the updated code or it
# will not have the right result.  For example, when changing a package from
# ST to PT and modifiying that package in the same push, the checkin-test.py
# script needs to see the updated PackagesList.cmake file before it runs or it
# will not allow the enable of the package.  (This actually happened with the
# Tempus package.)

remote_checkin_test_cmnds="cd $remote_trilinos_base_dir/Trilinos/CHECKIN && ./checkin-test-sems.sh --do-all --no-rebase --push" 

cd $local_trilinos_base_dir

if [ "$blocking_or_nonblocking" == "blocking" ] ; then

  # blocking

  ssh -q $remote_pull_test_push_server \
    "echo" \
    " && echo \"***\"" \
    " && echo \"*** 2) Hard reset the 'develop' branch and merge in the\"" \
    " \" '$local_branch_name' branch in Trilinos repo on\"" \
    " \" $remote_pull_test_push_server\"" \
    " && echo \"***\"" \
    " && echo" \
    " && $remote_branch_update_cmnds" \
    " && echo" \
    " && echo \"***\"" \
    " && echo \"*** 3) Doing non-blocking remote test/push on\"" \
    " \" '$remote_pull_test_push_server'\"" \
    " && echo \"***\"" \
    " && echo" \
    " && $remote_checkin_test_cmnds" \

else

  # nonblocking
  
  echo
  echo "***"
  echo "*** 2) Hard reset the 'develop' branch and merge in the '$local_branch_name' branch in Trilinos repo on $remote_pull_test_push_server"
  echo "***"
  echo
  
  ssh -q $remote_pull_test_push_server "$remote_branch_update_cmnds"
  
  echo
  echo "***"
  echo "*** 3) Doing non-blocking remote test/push on '$remote_pull_test_push_server' (see log file checkin-test-remote.log)"
  echo "***"
  
  ssh -q $remote_pull_test_push_server "$remote_checkin_test_cmnds" \
    &> checkin-test-remote.log \
    && echo && echo && echo "***" \
    && echo "*** Final test/push results from '$remote_pull_test_push_server':" \
    && echo "***" && echo \
    && tail -n10 checkin-test-remote.log &
  
  echo
  echo "You may now keep working on your local machine and wait for email notifications!  (final result will be printed when complete)"
  echo

fi
