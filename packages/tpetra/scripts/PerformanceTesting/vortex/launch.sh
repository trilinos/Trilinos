#!/bin/bash

#Run configure/build on 1 node, but testing on 4 nodes (for 16 total GPUs/ranks)
#Allow 1 hour for building, and 15 minutes for testing
bsub -q normal -Is -nnodes 1 -W 60 $WORKSPACE/VortexPerfScripts/configBuild.sh
bsub -q normal --shared-launch -Is -nnodes 4 -W 25 $WORKSPACE/VortexPerfScripts/test.sh

