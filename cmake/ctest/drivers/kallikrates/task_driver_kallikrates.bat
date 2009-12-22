rem Driver script for dashboards on kallikrates

rem Call this script and pass in the base directory for the dashboard
rem save state before changing anything
setlocal

rem GIT needs to know how to use ssh
set GIT_SSH=C:\Users\bmpersc\Documents\plink_wrap.bat

rem Set location of CTEST_EXE, and GIT_EXE
set GIT_EXE=C:\Program Files (x86)\Git\cmd\git
set CTEST_EXE=C:\Program Files (x86)\CMake 2.8.0 rc5\bin\ctest.exe
set BRANCH=trilinos-release-10-0-branch
set TRILINOS_REPOSITORY_LOCATION=software.sandia.gov:/space/git/nightly/Trilinos

rem Set the base directory which is one above where Trilinos will be 
rem checked out.

set BASEDIR=%~1

rem setup the environment for command line cl to run
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat"

rem change into the basedir
cd "%BASEDIR%"

rem checkout the basics from Trilinos needed to run the dashboard including
rem this script.
if exist Trilinos goto update else goto checkout

:update
  echo "Doing update of an existing directory"
  cd Trilinos
  call "%GIT_EXE%" checkout %BRANCH%
  call "%GIT_EXE%" pull
  cd ..
  goto endif

:checkout
  echo "Cloning the repository because none exists yet."
  call "%GIT_EXE%" clone %TRILINOS_REPOSITORY_LOCATION%
  cd Trilinos
  call "%GIT_EXE%" checkout --track origin/%BRANCH%

:endif

rem Now run ctest on each of the ctest build scripts for this machine

call "%CTEST_EXE%" -S "%BASEDIR%\Trilinos\cmake\ctest\drivers\kallikrates\ctest_windows_nightly_serial_release.cmake" -VV >"%BASEDIR%\ctest_msvc_nightly_serial_release_kallikrates.out" 2>&1

rem disabling temporarily since the mpi build is unstable on windows right now.
rem call "%CTEST_EXE%" -S "%BASEDIR%\Trilinos\cmake\ctest\drivers\kallikrates\ctest_windows_nightly_mpi_release.cmake" -VV >"%BASEDIR%\ctest_msvc_nightly_mpi_release_kallikrates.out" 2>&1

endlocal
