rem Driver script for dashboards on kallikrates

rem Call this script and pass in the base directory for the dashboard
rem save state before changing anything
setlocal

set BASEDIR=%~1

rem Setting the path to have git on it.
set PATH=%PATH%;C:\Program Files (x86)\Git\cmd

rem Set location of CTEST_EXE, and GIT_EXE
set GIT_EXE=C:\Program Files (x86)\Git\cmd\git
set CTEST_EXE=C:\Program Files (x86)\CMake 2.8\bin\ctest.exe

rem Have to set the path so that the tests can find the dlls during runtime. We need a better solution than this long term. Something like -rpath for gnu.
set PATH=%PATH%;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\epetra\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\anasazi\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\epetra\test\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\test\FancyOutputting;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\test\ParameterList


python %BASEDIR%\Trilinos\cmake\ctest\drivers\cron_driver.py
