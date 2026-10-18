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

verify_no_local_commits() {
  if [[ $(git rev-list --left-only --count $1...origin/$1) -ne 0 ]] ; then
    echo "Cannot proceed: branch $1 has local commits"
    exit 3;
  fi
}

# Copies only the files git tracks in Sierra, so ignored build or editor files (e.g. .cache,
# compile_commands.json, __pycache__) are not copied into Trilinos.
update_package() {
  exe rm -rf packages/$1
  exe "git -C $SIERRA archive --format=tar HEAD $1 | tar -x -C packages"
  exe git add --all packages/$1
}

export SIERRA=${SIERRA:-/fgs/$USER/projects/code_votd}
export SIERRA_BRANCH=master
export TRILINOS=${TRILINOS:-/fgs/$USER/projects/Trilinos}
export TRILINOS_BRANCH=develop

export SNAPSHOT_BRANCH=krino-snapshot

echo "SIERRA: $SIERRA"
echo "TRILINOS: $TRILINOS"

exe module load sierra-devel

exe cd $SIERRA
exe git checkout $SIERRA_BRANCH
verify_clean_repo $SIERRA
exe git fetch origin
verify_no_local_commits $SIERRA_BRANCH
exe git merge --ff-only origin/$SIERRA_BRANCH
STK_VERSION_STRING=$(./stk/stk_util/stk_util/registry/stk_version_gen.sh)

exe cd $TRILINOS

update_package krino
exe git rm -rf --ignore-unmatch packages/krino/krino_sierra packages/krino/Jamfile packages/krino/.clang-format packages/krino/AGENTS.md

export COMMIT_MESSAGE="Krino: Snapshot $(date +'%m-%d-%y %H:%M') from Sierra $STK_VERSION_STRING"
exe git commit -s -m "\"$COMMIT_MESSAGE\""
exe git push origin $SNAPSHOT_BRANCH
