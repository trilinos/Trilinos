#!/bin/bash

EXTRA_ARGS=$@

cd ~/
source .bash_profile

cd ~/Trilinos-remote-testing/Trilinos/CHECKIN/
./checkin-test-trilinos-test2-jw.sh --do-all --push \
$EXTRA_ARGS

echo `date` > checkin-test-trilinos-test2-jw-remote.finished.out
