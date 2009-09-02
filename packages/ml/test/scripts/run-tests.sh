#! /bin/sh

# Simple script to invoke the correct command line to test ML.
# You might need to modify the following variables, and the
# category.

TRILINOS_HOME=${HOME}/Trilinos
TRILINOS_INST=${TRILINOS_HOME}/LINUX_MPI
COMM=mpi
MPI_GO="mpirun -np "

${TRILINOS_HOME}/commonTools/test/utilities/runtests \
  --trilinos-dir=${TRILINOS_HOME} \
  --comm=$COMM \
  --build-dir=${TRILINOS_INST} \
  --category=MLExamples \
  --mpi-go="${MPI_GO}"
