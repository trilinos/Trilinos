#!/bin/bash

EXTRA_ARGS=$@

cd ~/
source .bash_profile

cd ~/PROJECTS/Trilinos.base.checkin/BUILDS/CHECKIN/
./checkin-test-godel.sh --do-all --push \
$EXTRA_ARGS

echo `date` > checkin-test-godel-remote.finished.out
