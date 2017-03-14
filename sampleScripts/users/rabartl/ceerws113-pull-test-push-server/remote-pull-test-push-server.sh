#!/bin/bash -e
#
# Script that gets called from a client machine using:
#
#   ssh -q <remotemachine> "~/remote-pull-test-push-server.sh ..."

cd /scratch/$USER/Trilinos.base/BUILDS/CHECKIN
./checkin-test-sems.sh --do-all --push "$@"
