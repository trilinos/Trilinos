
# CMake code to include in the base-level CMakeLists.txt file for each
# Trilinos package to add the targets 'runtests-serial' and 'runtests-mpi'

FUNCTION(TRILINOS_ADD_CUSTOM_TEST_TARGETS)
ADD_CUSTOM_TARGET(
  runtests-serial-${PROJECT_DIR_NAME}
   ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
  --trilinos-dir=${TRILINOS_HOME_DIR}
  --comm=serial
  --build-dir=${TRILINOS_BUILD_DIR}
  --category=${TRILINOS_TEST_CATEGORY}
  --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
  --verbosity=1
  --packages=${PROJECT_DIR_NAME}
  )


ADD_CUSTOM_TARGET(
  runtests-mpi-${PROJECT_DIR_NAME}
   ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
  --trilinos-dir=${TRILINOS_HOME_DIR}
  --comm=mpi
  --mpi-go="${TRILINOS_MPI_GO}"
  --max-proc=${MPIEXEC_MAX_NUMPROCS}
  --build-dir=${TRILINOS_BUILD_DIR}
  --category=${TRILINOS_TEST_CATEGORY}
  --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
  --verbosity=1
  --packages=${PROJECT_DIR_NAME}
  )
ENDFUNCTION()