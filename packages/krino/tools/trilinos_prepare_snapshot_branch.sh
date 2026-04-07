#!/usr/bin/env bash

# This is intended to work with a fork of trilinos, ie. https://github.com/drnobleabq/Trilinos

exe() {
  stdbuf -o0 -e0 echo "% $@" ;
  eval "$@" ;
  if [ $? -ne 0 ] ; then
    echo "'$@' failed.";
    exit 1;
  fi
}

verify_clean_repo() {
  if [[ -n $(git status -s) ]] ; then
    echo "Cannot proceed: repository at $1 is not clean"
    exit 2;
  fi
}

export TRILINOS=${TRILINOS:-/fgs/$USER/projects/Trilinos}
export TRILINOS_BRANCH=develop

export SNAPSHOT_BRANCH=krino-snapshot

echo "SIERRA: $SIERRA"
echo "TRILINOS: $TRILINOS"

exe module load sierra-devel

exe cd $TRILINOS
verify_clean_repo $TRILINOS

#Pull request workflow
exe git fetch --all
exe git checkout master
exe git merge upstream/master
exe git push origin master

exe git branch -f $SNAPSHOT_BRANCH
exe git checkout $SNAPSHOT_BRANCH
exe git reset --hard upstream/$TRILINOS_BRANCH
