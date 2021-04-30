#!/bin/bash

#Allow 90 minutes for all steps
salloc -N1 --time=90 -p short,batch --account=FY150090 $WORKSPACE/EclipsePerfScripts/configBuildTest.sh
