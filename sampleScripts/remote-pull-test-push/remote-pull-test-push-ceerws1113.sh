#!/bin/bash -e

cd $HOME/Trilinos.base/

./Trilinos/cmake/std/sems/remote-pull-test-push.sh \
  ceerws1113 \
  /scratch/$USER/TRILINOS_PUSH_SERVER \
  nonblocking
