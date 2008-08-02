#!/bin/bash

source /sw/bin/init.sh
export PATH=$PATH:/usr/local/bin
export CVSROOT=:ext:jmwille@software.sandia.gov:/space/CVS
export CVS_RSH=ssh

source ~/.keychain/s858352.srn.sandia.gov-sh > /dev/null

cd /Users/jmwille/TrilinosTestHarness/TrilinosDevelopment/Trilinos/commonTools/test/harness

export DYLD_LIBRARY_PATH=/Users/jmwille/TrilinosTestHarness/TrilinosDevelopment/Trilinos/SERIAL/packages/PyTrilinos/shared:$DYLD_LIBRARY_PATH
perl runharness --trilinos-dir=/Users/jmwille/TrilinosTestHarness/TrilinosDevelopment/Trilinos --build-name=s858352-nightly-serial

