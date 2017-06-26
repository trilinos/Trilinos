#!/bin/bash -e

BLOCKING_NONBLOCKING=$1
if [ "$BLOCKING_NONBLOCKING" == "" ] ; then 
  BLOCKING_NONBLOCKING=nonblocking
fi

ENABLE_ALL_PACKAGES=$2

cd $HOME/Trilinos.base/

./Trilinos/cmake/std/sems/remote-pull-test-push.sh \
  ceerws1113 \
  /scratch/$USER/TRILINOS_PUSH_SERVER \
  $BLOCKING_NONBLOCKING \
  $ENABLE_ALL_PACKAGES
