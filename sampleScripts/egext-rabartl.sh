#!/bin/bash

# Run eg recursively over extra repos

EXTRA_ARGS=$@

BASE_DIR=$PWD

# Must be in the correct order
EXTRA_REPOS_FULL_LIST="Panzer LIMEExt preCopyrightTrilinos"

#
# Run the eg command recursively
#

echo
eg $EXTRA_ARGS

for extra_repo in $EXTRA_REPOS_FULL_LIST; do
  if [ -d $extra_repo ]; then
    echo
    echo "***"
    echo "*** Git Repo: $extra_repo"
    echo "***"
    echo
    cd $extra_repo
    eg $EXTRA_ARGS
    cd $BASE_DIR
  fi
done
