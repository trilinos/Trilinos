#!/bin/tcsh

set EXTRA_ARGS="$*"

cd ~/
source .cshrc

cd ~/TBUILDS/CHECKIN/
./checkin-test-orthanc-remote-test-push.sh "--extra-pull-from=\"zan master\"" \
--do-all --push \
$EXTRA_ARGS

echo `date` > checkin-test-orthanc-remote-test-push-zan-remote-driver.finished.out
