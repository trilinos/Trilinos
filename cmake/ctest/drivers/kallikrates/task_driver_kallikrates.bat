rem Driver script for dashboards on kallikrates

rem Call this script and pass in the base directory for the dashboard
rem save state before changing anything
setlocal

rem CVS needs to know how to use ssh
set CVS_RSH=C:\T10\davcole\ssh.cmd

rem Set location of CTEST_EXE, and CVS_EXE
set CVS_EXE=C:\Program Files (x86)\CVSNT\cvs.exe
set CTEST_EXE=C:\T10\cmake-2.7.20090918-win32-x86\bin\ctest.exe

rem Set the base directory which is one above where Trilinos will be 
rem checked out.

set BASEDIR=%~1

rem setup the environment for command line cl to run
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat"

rem change into the basedir
cd "%BASEDIR%"

rem checkout the basics from Trilinos needed to run the dashboard including
rem this script.
call "%CVS_EXE%" -q -d :ext:software.sandia.gov:/space/CVS co -r trilinos-release-10-0-branch Trilinos/cmake Trilinos/CTestConfig.cmake

rem Now run ctest on each of the ctest build scripts for this machine

call "%CTEST_EXE%" -S "%BASEDIR%\Trilinos\cmake\ctest\drivers\kallikrates\ctest_windows_nightly_serial_release.cmake" -VV >"%BASEDIR%\ctest_msvc_nightly_serial_release_kallikrates.out" 2>&1

rem disabling temporarily since the mpi build is unstable on windows right now.
rem call "%CTEST_EXE%" -S "%BASEDIR%\Trilinos\cmake\ctest\drivers\kallikrates\ctest_windows_nightly_mpi_release.cmake" -VV >"%BASEDIR%\ctest_msvc_nightly_mpi_release_kallikrates.out" 2>&1

endlocal
