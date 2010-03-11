# corrin driver

# cygwin CVS needs to know to use ssh
set(CVS_RSH "C:/cygwin/bin/ssh.exe")
# Driver script for dashboards on Corrin
# Set location of CTEST_EXE, and CVS_EXE
set(CVS_EXE "C:\cygwin\bin\cvs.exe")

# Set the base directory which is one above where Trilinos will be 
# checked out.


# setup the environment for command line cl to run
add_custom_target(vsvars ALL
  COMMAND "C:\Program Files\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat")

# checkout the basics from Trilinos needed to run the dashboard including
# this script.
add_custom_target(cvs_checkout ALL
  COMMAND "${CVS_EXE}" -q -d :ext:software.sandia.gov:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake

add_test(corrin_driver ${DOWNLOADED_CTEST} -S ${TrilinosDriver_SOURCE_DIR}/corrin/ctest_windows_nightly_serial_release.cmake -VV --output-log ${TrilinosDriver_BINARY_DIR}/logs/ctest_msvc_nightly_serial_optimized_corrin.out)
