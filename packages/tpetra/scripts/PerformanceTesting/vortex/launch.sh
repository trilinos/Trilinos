#!/bin/bash -e

#Run configure/build on 1 node, but testing on 4 nodes (for 16 total GPUs/ranks)
bsub -q normal -Is -nnodes 1 -W 180 $WORKSPACE/PerfScripts/configBuild.sh
bsub -q normal --shared-launch -Is -nnodes 4 -W 240 $WORKSPACE/PerfScripts/test.sh

