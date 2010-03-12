# dash2win64 driver

# CVS needs to know how to use ssh
set(CVS_RSH "C:/Dashboards/Support/Trilinos-ssh/ssh.cmd")
# Set location of CVS_EXE
set(CVS_EXE "C:/Program Files (x86)/CVSNT/cvs.exe")

set(workingdir "${TrilinosDriver_BINARY_DIR}")

# setup the environment for command line cl to run
add_custom_target(vsvars ALL
  COMMAND "C:/Program Files/Microsoft Visual Studio 9.0/Common7/Tools/vsvars32.bat")

# checkout the basics from Trilinos needed to run the dashboard including
# this script.
add_custom_target(cvs_checkout ALL
  COMMAND "${CVS_EXE}" -q -d :ext:software.sandia.gov:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake
  WORKING_DIRECTORY "${workingdir})

# Now add test for each of the ctest build scripts for this machine

add_test(serial_release
  "${DOWNLOADED_CTEST}"
  -S
  "${TrilinosDriver_SOURCE_DIR}/dash2win64/ctest_windows_nightly_serial_release.cmake"
  -VV
  --output-log
  "${TrilinosDriver_BINARY_DIR}/ctest_msvc_nightly_serial_optimized_dash2win64.out"
)
