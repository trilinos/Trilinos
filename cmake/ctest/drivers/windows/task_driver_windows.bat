rem Driver script for dashboards on trilinos-win1

setlocal

rem Set the location of Git, Ninja, etc.
set SEMS_DIR=C:\projects\sems\install\win-x86_64
set NINJA_DIR=%SEMS_DIR%\utility\ninja\1.7.2
set CTEST_EXE=%SEMS_DIR%\utility\cmake\3.12.0\bin\ctest.exe
set GIT_EXE=%SEMS_DIR%\utility\git\2.13.0\cmd\git.exe
set PATH=%NINJA_DIR%;%PATH%
set TRILINOS_REPOSITORY_LOCATION=https://gitlab-ex.sandia.gov/trilinos-project/Trilinos.git

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
