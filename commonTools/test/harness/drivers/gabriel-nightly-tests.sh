#!/bin/bash

# Set up the environment for the sake of cron
cd /home/rabartl
source .bash_profile
./ssh-passkey-setup.sh

# Test cases

TRILINOS_BASE_DIR=/mnt/disk2/rabartl/Trilinos.nightly-tests/Trilinos

echo
echo "Running build for Trilinos serial debug"
echo
date
echo

cd  $TRILINOS_BASE_DIR/commonTools/test/harness
perl runharness --trilinos-dir=$TRILINOS_BASE_DIR --build-name=gabriel-nighly-serial-debug

echo
echo "Running build for Trilinos mpi optimized"
echo
date
echo

cd  $TRILINOS_BASE_DIR/commonTools/test/harness
perl runharness --trilinos-dir=$TRILINOS_BASE_DIR --build-name=gabriel-nighly-mpi

echo
date
echo
