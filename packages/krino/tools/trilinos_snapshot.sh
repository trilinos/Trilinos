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

update_package() {
  exe rm -rf packages/$1
  exe cp -r $SIERRA/$1 packages/$1
  exe git add --all packages/$1
}

export SIERRA=${SIERRA:-/fgs/$USER/trilinos-snapshot/code}
export SIERRA_BRANCH=master
export TRILINOS=${TRILINOS:-/fgs/$USER/trilinos-snapshot/Trilinos}
export TRILINOS_BRANCH=develop

export SNAPSHOT_BRANCH=krino-snapshot

echo "SIERRA: $SIERRA"
echo "TRILINOS: $TRILINOS"

exe module load sierra-devel

exe cd $SIERRA
verify_clean_repo $SIERRA
verify_no_local_commits $SIERRA_BRANCH
exe git checkout $SIERRA_BRANCH
exe repo sync
STK_VERSION_STRING=$(./stk/stk_util/stk_util/registry/stk_version_gen.sh)

exe cd $TRILINOS
verify_clean_repo $TRILINOS
verify_no_local_commits $TRILINOS_BRANCH

#Pull request workflow
exe git fetch --all
exe git checkout master
exe git merge upstream/master
exe git push origin master
exe git checkout develop
exe git merge upstream/develop
exe git push origin develop

exe git checkout $TRILINOS_BRANCH
exe git pull

exe git checkout $SNAPSHOT_BRANCH
exe git reset --hard upstream/$TRILINOS_BRANCH

update_package krino
exe git rm -rf packages/krino/krino_sierra packages/krino/Jamfile packages/krino/.clang-format

export COMMIT_MESSAGE="Krino: Snapshot $(date +'%m-%d-%y %H:%M') from Sierra $STK_VERSION_STRING"
exe git commit -am '"'$COMMIT_MESSAGE'"'
exe git push origin $SNAPSHOT_BRANCH
