#!/bin/bash

EXTRA_ARGS=$@

cd ~/
source .bash_profile

cd ~/PROJECTS/Trilinos.base.checkin/BUILDS/CHECKIN/
./checkin-test-godel-remote-test-push.sh "--extra-pull-from=\"gabriel master\"" \
--do-all --push \
$EXTRA_ARGS

echo `date` > checkin-test-godel-remote-test-push-gabriel-remote-driver.finished.out
