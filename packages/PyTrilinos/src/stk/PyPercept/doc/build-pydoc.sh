#!/bin/sh

# run this from you Trilinos build directory 
TRILINOS_CODE=$1
if [[ "$1" == "" ]] ; then
   echo "must specify path to Trilinos code base"
   exit 1
fi
curwd=`pwd`

cd $TRILINOS_CODE/packages/PyTrilinos/doc/Doxygen/
make Percept_dox.i

cd $curwd
touch $TRILINOS_CODE/packages/PyTrilinos/src/stk/PerceptMesh.i
make -j8

pydoc -w ./packages/PyTrilinos/src/stk/PyPercept/PerceptMesh.py

cp PerceptMesh.html $TRILINOS_CODE/packages/PyTrilinos/src/stk/PyPercept/doc


