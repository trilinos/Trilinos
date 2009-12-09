#!/bin/bash

EXTRA_ARGS=$@

cd ~/
source .bash_profile

cd ~/PROJECTS/Trilinos.base.checkin/BUILDS/CHECKIN/
./checkin-test-godel.sh --do-all \
$EXTRA_ARGS

# ToDo: Put --push back!
#
# --push

echo `date` > checkin-test-godel-remote.finished.out
