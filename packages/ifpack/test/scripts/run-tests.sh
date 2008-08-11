#! /bin/sh

# Simple script to invoke the correct command line to test IFPACK.
# You might need to modify the following variables, and the
# category.

TRILINOS_HOME=${HOME}/Trilinos
TRILINOS_INST=${TRILINOS_HOME}/LINUX_MPI
COMM=mpi
MPI_GO="mpirun -n "

${TRILINOS_HOME}/commonTools/test/utilities/runtests \
  --trilinos-dir=${TRILINOS_HOME} \
  --comm=$COMM \
  --build-dir=${TRILINOS_INST} \
  --category=IFPACKTest \
  --mpi-go="${MPI_GO}"
