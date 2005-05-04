#!/bin/csh
# see README in this directory for details on how to run this script
#
# Warning, to use this script you must first define TRILINOS_HOME and add the
# base path for the scrpts  'string-replace-list' and 'string-replace-list-r',
# which can be found in $TRILINOS_HOME/commonTools/refactoring, to your path.
#
# The invocation of string-replace-list-r makes the changes specified by 
# the file NewpackageToYourpackagename-string-replace.list
#
# 
#
mv src/New_Package_ConfigDefs.h src/Yourpackagename_ConfigDefs.h
mv src/New_Package_Version.h  src/Yourpackagename_Version.h
mv src/New_Package_config.h.in src/Yourpackagename_config.h.in
mv src/Newp_Hello.cpp src/Yourpackagename_Hello.cpp
mv src/Newp_Hello.h src/Yourpackagename_Hello.h 
mv src/Newp_Jambo.cpp src/Yourpackagename_Jambo.cpp
mv src/Newp_Jambo.h src/Yourpackagename_Jambo.h
mv test/scripts/daily/mpi/new_packageTestAllMpi test/scripts/daily/mpi/YourpackagenameTestAllMpi
mv test/scripts/daily/serial/new_packageTestAllSerial test/scripts/daily/serial/YourpackagenameTestAllSerial
mv python/src/New_Package.i python/src/Yourpackagename.i
mv python/test/testNewPackage.py python/test/testYourpackagename.py
#
#
#
mv ConvertToYourpackagename-string-replace.list ..
string-replace-list-r $TRILINOS_HOME/packages/ConvertToYourpackagename-string-replace.list


