#!/bin/sh
# This shell script automates updating of one git branch from another
# one with some amount of sanity checking.

# Update from develop if nothing else is said
if [ $# -eq 0 ]; then
    FROMBRANCH=develop
elif [ $# -eq 1 ]; then
    FROMBRANCH=$1
else
    echo "Syntax: update_branch.sh [<from_branch>]"
    exit 1
fi

CURRBRANCH=`git rev-parse --abbrev-ref HEAD`

# Do we have any modified files
LINECOUNT=`git status -s | wc -l`
if [ $LINECOUNT -ne 0 ]; then
    echo "ERROR: There are modified files in the current branch.  Please check them in";
    exit -1;
fi


echo "Updating $CURRBRANCH from $FROMBRANCH"

git checkout $FROMBRANCH
if [ $? -ne 0 ]; then echo "ERROR in checkout on $FROMBRANCH"; exit -2;fi
git pull --rebase
if [ $? -ne 0 ]; then echo "ERROR in pull on $FROMBRANCH"; exit -3;fi
git checkout $CURRBRANCH
if [ $? -ne 0 ]; then echo "ERROR in checkout on $CURRBRANCH"; exit -4;fi
git merge $FROMBRANCH
if [ $? -ne 0 ]; then echo "ERROR in merge"; exit -5;fi