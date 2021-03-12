#!/usr/bin/env bash

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

verify_no_local_commits() {
  if [[ $(git rev-list --left-only --count $1...origin/$1) -ne 0 ]] ; then
    echo "Cannot proceed: branch $1 has local commits"
    exit 3;
  fi
}

export SIERRA=${SIERRA:-/scratch/$USER/trilinos-snapshot/code}
export SIERRA_BRANCH=master
export TRILINOS=${TRILINOS:-/scratch/$USER/trilinos-snapshot/Trilinos}
export TRILINOS_BRANCH=develop

export COMMIT_MESSAGE="STK: Snapshot $(date +'%m-%d-%y %H:%M')"
export SNAPSHOT_BRANCH=stk-snapshot

echo "SIERRA: $SIERRA"
echo "TRILINOS: $TRILINOS"

exe module load sierra-devel

exe cd $SIERRA
verify_clean_repo $SIERRA
verify_no_local_commits $SIERRA_BRANCH
exe git checkout $SIERRA_BRANCH
exe repo sync

exe cd $TRILINOS
verify_clean_repo $TRILINOS
verify_no_local_commits $TRILINOS_BRANCH
exe git checkout $TRILINOS_BRANCH
exe git pull

exe git checkout $SNAPSHOT_BRANCH
exe git rebase $TRILINOS_BRANCH
exe rm -rf packages/stk
exe cp -r $SIERRA/stk packages/stk
exe git add --all packages/stk
exe git commit -am '"'$COMMIT_MESSAGE'"'
exe git push origin $SNAPSHOT_BRANCH
