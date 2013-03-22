rem Driver script for dashboards on kallikrates

rem Call this script and pass in the base directory for the dashboard
rem save state before changing anything
setlocal

set TDD_CTEST_TEST_TYPE=Nightly

set BASEDIR=%~1

rem Setting the path to have git on it.
set PATH=%PATH%;C:\Program Files (x86)\Git\cmd

rem Set location of GIT_EXE
set GIT_EXE=C:\Program Files (x86)\Git\cmd\git

rem Have to set the path so that the tests can find the dlls during runtime. We need a better solution than this long term. Something like -rpath for gnu.
set PATH=%PATH%;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\epetra\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\epetra\test\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\anasazi\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\core\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\parameterlist\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\parameterlist\test\FancyOutputting;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\parameterlist\test\ParameterList;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\comm\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\comm\test\ParameterList;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\numerics\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\remainder\src


python %BASEDIR%\Trilinos\cmake\ctest\drivers\cron_driver.py
