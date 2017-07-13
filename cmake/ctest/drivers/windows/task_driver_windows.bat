rem Driver script for dashboards on trilinos-win1

setlocal

rem Set the location of Git, Ninja, etc.
set SEMS_DIR=C:\projects\sems\install\win-x86_64
set NINJA_DIR=%SEMS_DIR%\utility\ninja\1.7.2
set CTEST_EXE=%SEMS_DIR%\utility\cmake\3.8.1\bin\ctest.exe
set GIT_EXE=%SEMS_DIR%\utility\git\2.13.0\cmd\git.exe
set PATH=%NINJA_DIR%;C:\projects\sems\install\win-x86_64\compiler\Microsoft MPI\8.1.12438.1084\Bin\;C:\Windows\system32;C:\Windows;C:\Windows\System32\Wbem;C:\Windows\System32\WindowsPowerShell\v1.0\;C:\Program Files\Microsoft SQL Server\130\Tools\Binn\;C:\Program Files\Microsoft DNX\Dnvm\;C:\projects\sems\install\win-x86_64\utility\git\2.13.0\cmd;C:\Program Files (x86)\Windows Kits\10\Windows Performance Toolkit\
set TRILINOS_REPOSITORY_LOCATION=git@github.com:trilinos/Trilinos.git

set BASEDIR=%~1
cd %BASEDIR%

set SCRIPT_DIR=%BASEDIR%\Trilinos\cmake\ctest\drivers\windows

call %SEMS_DIR%\compiler\VisualStudio\14.0\VC\vcvarsall.bat amd64

rem Run the release build and tests
call %CTEST_EXE% -S "%SCRIPT_DIR%\ctest_windows_mpi_release.cmake" -VV > "%BASEDIR%\ctest_msvc_mpi_windows_release.out" 2>&1

rem Add the /bigobj compiler flag to get the debug build to work properly. Then run the debug build and tests.
set CL=/bigobj
call %CTEST_EXE% -S "%SCRIPT_DIR%\ctest_windows_mpi_debug.cmake" -VV > "%BASEDIR%\ctest_msvc_mpi_windows_debug.out" 2>&1

endlocal
