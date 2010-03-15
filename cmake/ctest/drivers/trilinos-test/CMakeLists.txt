# trilinos-test driver

set(CTEST_EXE /home/trilinos/cmake/bin/ctest)
set(EG_EXE /home/trilinos/git/bin/eg)
set(BASEDIR /home/bmpersc/nightly/Trilinos.base/release_10)
set(BASEDATADIR /home/bmpersc/nightly)
set(DRIVER_SCRIPT_DIR "${TrilinosDriver_SOURCE_DIR}/trilinos-test")
set(workingdir "${TrilinosDriver_BINARY_DIR}")
set(BRANCH "trilinos-release-10-0-branch")
set(TRILINOS_REPOSITORY_LOCATION "software.sandia.gov:/space/git/Trilinos")

#export CMAKE_LIBRARY_PATH="/home/trilinos/tpl/gcc4.1.2/exodusII_4.84/lib:/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/lib"
#export CMAKE_INCLUDE_PATH="/home/trilinos/tpl/gcc4.1.2/exodusII_4.84/include:/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/include"
#export TRILINOSDATADIRECTORY="${BASEDATADIR}/TrilinosData"

#get the date for use in log files
execute_process(COMMAND date "+%m-%d-%Y" OUTPUT_VARIABLE DATE)

#checkout the trilinos data directory. Needed for a few tests to run
#execute_process(COMMAND cvs -q -d :ext:software:/space/CVS co TrilinosData
#  WORKING_DIRECTORY ${BASEDATADIR})

add_test(mpi_release
  ${CTEST_EXE}
  -S
  "${DRIVER_SCRIPT_DIR}/ctest_linux_nightly_mpi_release_trilinos-test.cmake"
  -VV
  --output-log
  "${workingdir}/MPI_RELEASE_10.0_${DATE}.log"
  )
  
add_test(mpi_release_shared
  ${CTEST_EXE}
  -S
  "${DRIVER_SCRIPT_DIR}/ctest_linux_nightly_mpi_release_shared_trilinos-test.cmake"
  -VV
  --output-log
  "${workingdir}/MPI_RELEASE_10.0_SHARED_${DATE}.log"
  )

add_test(serial_release
  ${CTEST_EXE}
  -S
  "${DRIVER_SCRIPT_DIR}/ctest_linux_nightly_serial_release_trilinos-test.cmake"
  -VV
  --output-log
  "${workingdir}/SERIAL_RELEASE_10.0_${DATE}.log"
  )
