rem Driver script for dashboards on kallikrates

rem Call this script and pass in the base directory for the dashboard
rem save state before changing anything
setlocal

rem Setting the path to have git on it.
set PATH=%PATH%;C:\Program Files (x86)\Git\cmd

rem Set location of CTEST_EXE, and GIT_EXE
set GIT_EXE=C:\Program Files (x86)\Git\cmd\git
set CTEST_EXE=C:\Program Files (x86)\CMake 2.8\bin\ctest.exe

rem Set the base directory which is one above where Trilinos will be 
rem checked out.

rem setup the environment for command line cl to run
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat"

rem Now run ctest on each of the ctest build scripts for this machine

call "%CTEST_EXE%" -S "Trilinos\cmake\ctest\drivers\kallikrates\ctest_windows_nightly_serial_release.cmake" -VV 

endlocal
