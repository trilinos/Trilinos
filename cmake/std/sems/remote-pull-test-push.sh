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
#       <remote_pull_test_push_server> \
#       <remote_trilinos_base_dir> \
#       <blocking_mode>
#
# where <local_trilinos_base_dir> is the base directory for Trilinos on the
# local machine (i.e. <local_trilinos_base_dir>/Trilinos) and <blocking_mode>
# is either 'blocking' or 'nonblocking'.
#
# If one cannot get public/private SSH access to the remote machine to work
# without typing a password, then one must use 'blocking' mode.  This will
# then only the user to type the password for the remote machine one time.
# Otherwise, for 'nonblocking' mode, the second remote ssh -q invocation that
# requires a password will not work since it is backgrounded.
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

enable_all_packages=$4
if [ "$enable_all_packages" == "" ] ; then
  # Don't enable all packages
  enable_all_packages_arg=""
elif [ "$enable_all_packages" == "all" ] ; then
  enable_all_packages_arg="--enable-all-packages=on"
else
  # Invalid value!
  echo "Error: Forth argument value '$enable_all_packages' is not acceptable!  Must pass in 'all' or empty ''!"
  exit 4
fi
echo "enable_all_packages_arg = '$enable_all_packages_arg'"

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

remote_branch_update_cmnds="cd $remote_trilinos_base_dir/Trilinos && git checkout develop && git fetch && git reset --hard @{u} && git fetch intermediate-repo && git merge --no-ff intermediate-repo/$local_branch_name"
# Some notes about the above commands:
#
# First, the remote tracking 'develop' branch is updated to what is in github
# 'develop'.  That is done with a fetch followed by a reset --hard.  That make
# sure that helps to avoid another trivial merge commit from github 'develop'
# when the checkin-test.py script runs.
#
# Second, we do the fetch and merge of the new branch before running the
# checkin-test.py script on the remote machine.  That is needed for some use
# cases where the checkin-test.py needs to be run given the updated code or it
# will not have the right result.  For example, when changing a package from
# ST to PT and modifiying that package in the same push, the checkin-test.py
# script needs to see the updated PackagesList.cmake file before it runs or it
# will not allow the enable of the package.  (This actually happened with the
# Tempus package.)

remote_checkin_test_cmnds="cd $remote_trilinos_base_dir/Trilinos/CHECKIN && ./checkin-test-sems.sh $enable_all_packages_arg --do-all --no-rebase --push" 

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
