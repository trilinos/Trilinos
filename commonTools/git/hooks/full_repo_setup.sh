#!/bin/bash

# This script sets up the email hooks for a new repo in one shot.
# Only the permissions on the repo must be changed after test.
#
# This script should be run from the GIT_REPO_BASE dir as:
#
#   $ cd GIT_REPO_BASE
#   $ THIS_DIR/full_repo_setup.sh  DEFAULT_EMAIL  EMAIL_SUMMARY_PREFIX  
#

DEFAULT_EMAIL=$1; shift
EMAIL_SUMMARY_PREFIX=$1; shift

_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/full_repo_setup.sh/\1/g"`

GIT_REPO_BASE=$PWD
THIS_DIR=$_SCRIPT_DIR

echo "1) Copy driver files ..."
cd $GIT_REPO_BASE_DIR/hooks
$THIS_DIR/copy_hooks_scripts.sh

echo "2) Copy configure file and customize ..."
cd $GIT_REPO_BASE_DIR
cp $THIS_DIR/config .
eg config --add hooks.mailinglist "$DEFAULT_EMAIL"
eg config --add hooks.emailprefix "[$EMAIL_SUMMARY_PREFIX] "

echo "3) Create initial dirs_to_emails file with default email address ..."

cd $GIT_REPO_BASE_DIR
echo "
   .+ $DEFAULT_EMAIL
   " >> hooks/dirs_to_emails
