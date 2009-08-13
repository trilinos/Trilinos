#!/usr/bin/env python

from CheckinTest import *


#
# Read in the commandline arguments
#

usageHelp = r"""checkin-test.py [OPTIONS]

This tool does checkin testing for Trilinos with CMake/CTest.

There are two basic configurations that are tested: MPI_DEBUG and
SERIAL_RELEASE.  Several configure options are varied in these two builds to
try to catch as much conditional configration behavior has possible.  If
nothing else, please at least do the MPI_DEBUG build since that will cover the
most code and best supports day-to-day development efforts.  However, if you
are changing code that might break the serial build, please allow the
SERIAL_RELEASE build to be run as well.

The following steps are performed:

1) Do a CVS update of the code (and detect the files that have changed
locally).

2) Select the list of packages to enable forward based on the directories
where there are changed files (or from a list passed in by the user).

3) For each build case (e.g. MPI_DEBUG, SERIAL_RELEASE, etc.)

  3.a) Configure an MPI DEBUG build in a standard way for all of the packages
  that have changed and all of the packages that depend on these packges
  forward.
  
  3.b) Build all configured code with 'make' (with -jN).
  
  3.c) Run all tests for enabled packages
  
  3.d) Analyze the results of the CVS update, configure, build, and tests and
  send email about results.

4) [Optionally] Commit the code given a commit message.

The recommended way to use this tool is to create a new base directory apart
from your standard build directories such as:

  $ cd SOME_OTHER_BASE_DIR
  $ mkdir PRECHECKIN

The most basic way to do the checkin test is:

  $ cd SOME_OTHER_BASE_DIR/PRECHECKIN
  $ $TRILINOS_HOME/cmake/python/checkin-test.py

If your MPI installation, other compilers, and standard TPLs (i.e. BLAS and
LAPACK) can be found automatically, then that is all you will need to do.
However, if the setup can not be determined automatically, then you can append
a set of CMake variables that will get read in the files:

  SOME_OTHER_BASE_DIR/PRECHECKIN/SERIAL_RELEASE.config
  SOME_OTHER_BASE_DIR/PRECHECKIN/MPI_DEBUG.config
  SOME_OTHER_BASE_DIR/PRECHECKIN/COMMON.config

Actually, skeletons of these files will automatically be written out with all
of the CMake cache variables commented out.  Any CMake cache variables in
these files will be read into and passed on the configure line to 'cmake'.

WARNING: Please do not add any CMake cache variables in the *.config files
that will alter what packages are built or what tests are run.  The goal of
these configureation files is to allow you to specifiy the minimum environment
to find MPI and your compilers and basic TPLs.

NOTE: Before running this script, you should first do an CVS update and
examine what files are changed to make sure you want to commit what you have.
Also, please look out for unknown files that you many need to add to the VC
repositoory with 'cvs add'.  Your working directory needs to be 100% ready to
commit before running this script.

NOTE: You don't need to run this script if you have not changed any files that
affect the build or the tests.  For example, if all you have changed is
document files, then you don't need to run this script before committing.

Common use cases for using this script are as follows:

(*) Basic full testing without commit:

   [no options]

(*) Basic full testing with commit:

   --do-commit --commit-msg-header-file=<SOME_FILE_NAME>

(*) The minimum acceptable testing where code has been changed:

  --no-enable-fwd-packages --without-serial-release

  NOTE: This will do only an MPI DEBUG build and will only build and run the
  tests for the packages

(*) Test only the packages modified and not forward dependent packages:

  --no-enable-fwd-packages

  NOTE: This is a safe thing to do when only tests in the given modified
  packages are changed and not library code.  This can speed up the testing
  process and is to be preferred over not running this script at all

(*) Test only a specific set of packages and no others:

  --enable-packages=PACKAGEA,PACKAGEB,PACKAGEC --no-enable-fwd-packages

  NOTE: Using this option is greatly prefered to not running this script at
  all.

(*) MPI DEBUG only (run at least this if nothing else):

  --without-serial-release

(*) Seeing the default option values without doing anything:

   --show-defaults

"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--skip-update", dest="doUpdate", action="store_false",
  help="Skip the update of Trilinos" )
clp.add_option(
  "--do-update", dest="doUpdate", action="store_true",
  help="Do the CVS update of Trilinos", default=True )

clp.add_option(
  "--update-command", dest="updateCommand", type="string",
  default="cvs -q update -dP",
  help="Command used to update the working copy of Trilinos." )

clp.add_option(
  "--skip-configure", dest="doConfigure", action="store_false",
  help="Skip the configure of Trilinos" )
clp.add_option(
  "--do-configure", dest="doConfigure", action="store_true",
  help="Do the configure of Trilinos", default=True )

clp.add_option(
  "--enable-packages", dest="enablePackages", type="string", default="",
  help="List of comma separated Trilinos packages to test changes for" \
  +" (example, 'Teuchos,Epetra').  If this list of packages is empty, then" \
  +" the list of packages to enable will be determined automatically by examining" \
  +" the set of modified files from the version control update log." )

clp.add_option(
  "--no-enable-fwd-packages", dest="enableFwdPackages", action="store_false",
  help="Do not enable forward Trilinos packages" )
clp.add_option(
  "--enable-fwd-packages", dest="enableFwdPackages", action="store_false",
  help="Enable forward Trilinos packages", default=True )

clp.add_option(
  "--skip-build", dest="doBuild", action="store_false",
  help="Skip the build of Trilinos" )
clp.add_option(
  "--do-build", dest="doBuild", action="store_true",
  help="Do the build of Trilinos", default=True )

clp.add_option(
  "--make-options", dest="makeOptions", type="string", default="",
  help="The options to pass to make (e.g. -j4)." )

clp.add_option(
  "--skip-test", dest="doTest", action="store_false",
  help="Skip the running the enabled Trilinos tests" )
clp.add_option(
  "--do-test", dest="doTest", action="store_true",
  help="Do the running of the enabled Trilinos tests", default=True )

clp.add_option(
  "--ctest-options", dest="ctestOptions", type="string", default="",
  help="Extra options to pass to 'ctest'." )

clp.add_option(
  "--no-show-all-tests", dest="showAllTests", action="store_false",
  help="Don't show all of the test results in the summary email" )
clp.add_option(
  "--show-all-tests", dest="showAllTests", action="store_true",
  help="Show all of the tests in the summary email", default=False )

clp.add_option(
  "--skip-commit", dest="doCommit", action="store_false",
  help="Skip the commit at the end if everything works out" )
clp.add_option(
  "--do-commit", dest="doCommit", action="store_true",
  help="Do the commit at the end if everything works out", default=False)

clp.add_option(
  "--commit-msg-header-file", dest="commitMsgHeaderFile", type="string", default="",
  help="Custom commit message file if commiting with --do-commit." \
  + "  If an relative path is given, this is expected to be with respect to the" \
  +" base source directory for Trilinos." )

clp.add_option(
  "--without-mpi-debug", dest="withMpiDebug", action="store_false",
  help="Skip the mpi debug build." )
clp.add_option(
  "--with-mpi-debug", dest="withMpiDebug", action="store_true",
  help="Do the mpi debug build.", default=True )

clp.add_option(
  "--without-serial-release", dest="withSerialRelease", action="store_false",
  help="Skip the serial release build." )
clp.add_option(
  "--with-serial-release", dest="withSerialRelease", action="store_true",
  help="Do the serial release build.", default=True )

clp.add_option(
  "--from-scratch", dest="rebuild", action="store_false",
  help="Blow everything away and rebuild from scratch." )
clp.add_option(
  "--rebuild", dest="rebuild", action="store_true",
  help="Keep and existing build tree and simply rebuild on top of it.", default=True )

clp.add_option(
  "--send-email-to", dest="sendEmailTo", type="string",
  default=os.environ["USER"]+"@sandia.gov",
  help="List of comma-separated email addresses to send notification to." )

clp.add_option(
  "--show-defaults", dest="showDefaults", action="store_true",
  help="Show the default option values and do nothing at all.",
  default=False )

(options, args) = clp.parse_args()

# NOTE: Above, in the pairs of boolean options, the *last* add_option(...) 
# takes effect!  That is whay the commands are ordered the way they are!


#
# Echo the commandline
#

print ""
print "**************************************************************************"
print "Script: checkin-test.py \\"
if options.doUpdate:
  print "  --do-update \\"
else:
  print "  --skip-update \\"
print "  --update-command='"+options.updateCommand+"' \\"
if options.doConfigure:
  print "  --do-configure \\"
else:
  print "  --skip-configure \\"
print "  --enable-packages='"+options.enablePackages+"' \\"
if options.enableFwdPackages:
  print "  --enable-fwd-packages \\"
else:
  print "  --no-enable-fwd-packages \\"
if options.doBuild:
  print "  --do-build \\"
else:
  print "  --skip-build \\"
print "  --make-options='"+options.makeOptions+"' \\"
if options.doTest:
  print "  --do-test \\"
else:
  print "  --skip-test \\"
print "  --ctest-options='"+options.ctestOptions+"' \\"
if options.showAllTests:
  print "  --show-all-tests \\"
else:
  print "  --no-show-all-tests \\"
if options.doCommit:
  print "  --do-commit \\"
else:
  print "  --skip-commit \\"
print "  --commit-msg-header-file='"+options.commitMsgHeaderFile+"' \\"
if options.withMpiDebug:
  print "  --with-mpi-debug \\"
else:
  print "  --without-mpi-debug \\" 
if options.withSerialRelease:
  print "  --with-serial-release \\"
else:
  print "  --without-serial-release \\" 
if options.rebuild:
  print "  --rebuild \\"
else:
  print "  --from-scratch \\"
print "  --send-email-to='"+options.sendEmailTo+"' \\"

#
# Execute
#

import time

if not options.showDefaults:

  print "\nStarting time:", getCmndOutput("date",True)
  
  t1 = time.time()
  success = checkinTest(options)
  t2 = time.time()
  print "\nTotal time for checkin-test.py =", (t2-t1)/60.0, "minutes"
  
  print "\nFinal time:", getCmndOutput("date",True)
  
  if success:
    exit(0)
  else:
    exit(1)
