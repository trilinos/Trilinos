# kallikrates driver

# GIT needs to know how to use ssh
set(GIT_SSH "C:\Users\bmpersc\Documents\plink_wrap.bat")

# Set location of CTEST_EXE, and GIT_EXE
set(GIT_EXE "C:\Program Files (x86)\Git\cmd\git")
set(BRANCH trilinos-release-10-0-branch)
set(TRILINOS_REPOSITORY_LOCATION software.sandia.gov:/space/git/nightly/Trilinos)

set(workingdir "${TrilinosDriver_BINARY_DIR}")

# setup the environment for command line cl to run
add_custom_target(vsvars ALL
  COMMAND "C:\Program Files\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat")

# checkout the basics from Trilinos needed to run the dashboard including
# this script.
if(NOT EXISTS "${workingdir}/Trilinos")
  message("Cloning the repository because none exists yet.")
  execute_process(COMMAND "${GIT_EXE}" clone ${TRILINOS_REPOSITORY_LOCATION}
    WORKING_DIRECTORY "${workingdir}")
  execute_process(COMMAND "${GIT_EXE}" checkout --track origin/${BRANCH}
    WORKING_DIRECTORY "${workingdir}/Trilinos")
endif()

# Add targets to update the local repository
add_custom_target(gitcheckout ALL
  COMMAND "%GIT_EXE%" checkout ${BRANCH}
  WORKING_DIRECTORY "${workingdir}/Trilinos"
  COMMENT "Doing update of an existing directory")
add_custom_target(gitpull ALL
  COMMAND "%GIT_EXE%" pull
  WORKING_DIRECTORY "${workingdir}/Trilinos"
  DEPENDS gitcheckout)

# Now run ctest on each of the ctest build scripts for this machine

add_test(serial_release
  "${CTEST_EXE}"
  -S
  "${workingdir}/Trilinos/cmake/ctest/drivers/kallikrates/ctest_windows_nightly_serial_release.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_msvc_nightly_serial_release_kallikrates.out"
)

# disabling temporarily since the mpi build is unstable on windows right now.
#add_test(mpi_release
#  "${CTEST_EXE}"
#  -S
#  "${workingdir}/Trilinos/cmake/ctest/drivers/kallikrates/ctest_windows_nightly_mpi_release.cmake"
#  -VV
#  --output-log
#  "${workingdir}/ctest_msvc_nightly_mpi_release_kallikrates.out"
#)

