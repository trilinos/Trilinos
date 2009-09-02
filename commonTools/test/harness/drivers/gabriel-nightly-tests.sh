#!/bin/bash

TRILINOS_BUILD=$1; shift

# Set up the environment for the sake of cron
cd /home/rabartl
source .bash_profile
./ssh-passkey-setup.sh

# Test cases

TRILINOS_BASE_DIR=/mnt/disk2/rabartl/Trilinos.nightly-tests
TRILINOS_DIR=$TRILINOS_BASE_DIR/Trilinos


echo
echo "Checking out TrilinosData"
echo
date
echo

cd  $TRILINOS_BASE_DIR
cvs -d :ext:software:/space/CVS co TrilinosData


if [ "$TRILINOS_BUILD" == "serial" ]; then

  echo
  echo "Running build for Trilinos serial debug"
  echo
  date
  echo
  
  cd  $TRILINOS_DIR/commonTools/test/harness
  perl runharness --trilinos-dir=$TRILINOS_DIR --build-name=gabriel-nighly-serial-debug
  
fi

if [ "$TRILINOS_BUILD" == "mpi" ]; then

  echo
  echo "Running build for Trilinos mpi optimized"
  echo
  date
  echo
  
  cd  $TRILINOS_DIR/commonTools/test/harness
  perl runharness --trilinos-dir=$TRILINOS_DIR --build-name=gabriel-nighly-mpi
  
fi

echo
date
echo
