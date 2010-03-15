:# Windows batch file
:#
:# Prerequisites:
:# - eg or git is in the PATH
:# - this machine is set up to pull/clone Trilinos using eg or git
:# - python is in the PATH and is version 2.6 or later
:#     (the Windows version requires python 2.6 for its
:#      tarfile/zipfile extractall method)

setlocal

set MACHINE=%COMPUTERNAME%

:# SCRIPT_DIR is the directory where *this* script is:
:#   ( note: %~dp0 ends in a trailing \ )
:#
set SCRIPT_DIR=%~dp0

:# BASE_DIR is the parent directory of the "Trilinos" source tree,
:# which we compute as relative to SCRIPT_DIR:
:#
set BASE_DIR=%SCRIPT_DIR%..\..\..\..

:# LOGFILE is located in $BASE_DIR/logs. Use current day of week in
:# log file name (most recent 7 saved before cycling around to the
:# same log file name...)
:#
for /f "tokens=1 delims= " %%i in ("%DATE%") do set WEEKDAY=%%i

set LOGFILE=%BASE_DIR%\logs\cron_driver_%MACHINE%_%WEEKDAY%.log

del /q "%LOGFILE%"

:# Make sure logs and tools directories exist:
:#
if not exist "%BASE_DIR%\logs" (
  mkdir "%BASE_DIR%\logs"
  if not exist "%BASE_DIR%\logs" (
    echo error: could not create directory "%BASE_DIR%\logs"
    exit /b 1
  )
)

if not exist "%BASE_DIR%\tools" (
  mkdir "%BASE_DIR%\tools"
  if not exist "%BASE_DIR%\tools" (
    echo error: could not create directory "%BASE_DIR%\tools" >>"%LOGFILE%"
    exit /b 2
  )
)

:# No 'which' on older windows (it's called 'where' on Vista and later...)
:# to get full paths to eg/git and python: use simple 'git' and 'python'
:# instead... -- override by setting EG_EXE and/or PYTHON_EXE in the
:# environment before calling this script.
:#
if "x%EG_EXE%" equ "x" (
  set EG_EXE=git
)

if "x%PYTHON_EXE%" equ "x" (
  set PYTHON_EXE=python
)

if "x%EG_EXE%" equ "x" (
  echo error: no eg or git in the PATH >>"%LOGFILE%"
  exit /b 3
)

if "x%PYTHON_EXE%" equ "x" (
  echo error: no python in the PATH >>"%LOGFILE%"
  exit /b 4
)

set TRILINOS_REPOSITORY=software.sandia.gov:/space/git/Trilinos


echo MACHINE=+%MACHINE%+ >>"%LOGFILE%"
echo SCRIPT_DIR=+%SCRIPT_DIR%+ >>"%LOGFILE%"
echo BASE_DIR=+%BASE_DIR%+ >>"%LOGFILE%"
echo WEEKDAY=+%WEEKDAY%+ >>"%LOGFILE%"
echo LOGFILE=+%LOGFILE%+ >>"%LOGFILE%"
echo EG_EXE=+%EG_EXE%+ >>"%LOGFILE%"
echo PYTHON_EXE=+%PYTHON_EXE%+ >>"%LOGFILE%"
echo TRILINOS_REPOSITORY=+%TRILINOS_REPOSITORY%+ >>"%LOGFILE%"


:# Begin
:#
echo. >>"%LOGFILE%"
echo "Starting nightly TrilinosDriver dashboard on %MACHINE%: %DATE% %TIME%" >>"%LOGFILE%"
echo. >>"%LOGFILE%"

echo >>"%LOGFILE%"
echo "Checking out / updating the scripts:" >>"%LOGFILE%"
echo >>"%LOGFILE%"

:# Checkout / update the Trilinos repository for the latest *.cmake scripts
:#
cd "%BASE_DIR%"
if exist Trilinos (
  echo Doing an update of existing directory >>"%LOGFILE%"
  cd Trilinos
  call "%EG_EXE%" pull >>"%LOGFILE%"
  cd ..
) else (
  echo Cloning the repository because none exists yet >>"%LOGFILE%"
  call "%EG_EXE%" clone %TRILINOS_REPOSITORY% >>"%LOGFILE%"
)

:# Download and install CMake/CTest 'release' build
:#
echo Downloading/installing CMake >>"%LOGFILE%"

rmdir /q /s "%BASE_DIR%\tools\cmake-release"

call "%PYTHON_EXE%" "%BASE_DIR%\Trilinos\cmake\python\download-cmake.py" --skip-detect "--install-dir=%BASE_DIR%\tools\cmake-release" --installer-type=release >>"%LOGFILE%"

for /f "usebackq delims=" %%i in (`dir /s /b "%BASE_DIR%\tools\cmake-release\ctest.exe"`) do set CTEST_EXE=%%i

echo CTEST_EXE=+%CTEST_EXE%+ >>"%LOGFILE%"

if not exist "%CTEST_EXE%" (
  echo error: ctest not found after installation... >>"%LOGFILE%"
  exit /b 5
)

:# Run a single TrilinosDriver dashboard on this machine:
:#
for /f "usebackq delims=" %%i in (`"%CTEST_EXE%" --version`) do set CTEST_VERSION=%%i
echo CTEST_VERSION=+%CTEST_VERSION%+ >>"%LOGFILE%"

:#set TDD_CRON_DRIVER_LOGFILE=%LOGFILE%
set TDD_CRON_DRIVER_SCRIPT=%~f0

echo Running ctest -S TrilinosDriverDashboard.cmake >>"%LOGFILE%"

call "%CTEST_EXE%" -S "%SCRIPT_DIR%TrilinosDriverDashboard.cmake" -VV >>"%LOGFILE%" 2>&1
set CTEST_RESULT=%ERRORLEVEL%

echo CTEST_RESULT=+%CTEST_RESULT%+ >>"%LOGFILE%"

if "x%CTEST_RESULT%" neq "x0" (
  echo error: ctest returned non-zero error value, script will exit with %CTEST_RESULT% >>"%LOGFILE%"
)

:# Done
:#
echo >>"%LOGFILE%"
echo "Ending nightly TrilinosDriver dashboard on %MACHINE%: %DATE% %TIME%" >>"%LOGFILE%"
echo >>"%LOGFILE%"

:# Propagate ctest return value
:#
exit /b %CTEST_RESULT%
