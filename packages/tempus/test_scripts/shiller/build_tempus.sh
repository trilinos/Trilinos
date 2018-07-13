#!/bin/bash


cd $WORKSPACE
if [ ! -d build_drekar ]; then 
    mkdir build_drekar
fi

cd build_drekar
rm -rf C*

$WORKSPACE/Trilinos/tempus/test_scripts/shiller/configure-tempus-shiller.sh 
if [ "$?" -gt 0 ]; then
    echo "Error in configuration of drekar, trying to continue"
fi
make -j16 install
if [ "$?" -gt 0 ]; then
    echo "Error in build of drekar"
    exit 2
fi
