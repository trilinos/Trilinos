#!/usr/bin/env python

from CheckinTest import *


#
# Read in the commandline arguments
#

usageHelp = r"""checkin-test.py [OPTIONS]

This tool does checkin testing for Trilinos with CMake/CTest and can actually
do the checkin itself in a safe way.


Quickstart:
-----------

In order to do a solid checkin, do the following:

1) Do a 'cvs -nq update -dP' to review the changes that you have made:

  $ cd $TRILINOS_HOME
  $ cvs -nq update -dP > update.out ; cat update.out

  NOTE: If you see any files/directories that you need to add in '?' lines,
  then please do a 'cvs add' on these first.

  NOTE: TRILINOS_HOME is just a mock variable here. You can of course just
  replace $TRILINOS_HOME with the absoute or relative path to the base
  Trilinos source directory.

2) Create a commit log file in the main source directory:

  $ cd $TRILINOS_HOME
  $ xemacs -nw checkin_message

  NOTE: Fill out this checkin message listing what you have changed.  Please
  use the Trilinos template for this file.

3) Do the checkin test (and commit) with:

  $ cd SOME_BASE_DIR
  $ mkdir CHECKIN
  $ cd CHECKIN
  $ $TRILINOS_HOME/cmake/python/checkin-test.py \
      --do-all --commit --commit-msg-header-file=checkin_message

  NOTE: The above will automatically enable the correct packages and then
  build the code, run the tests, send you emails about what happened, and then
  do the commit if everything passed.

  NOTE: Once you start running the checkin-test.py script, you can go off and
  do something else and just check your email to see if all the builds and
  tests passed and if the commit happened or not.

  NOTE: You need to have SSH public/private keys set up to software.sandia.gov
  for the commit to happen automatically.

  NOTE: You can actually finish this commit message file while the
  checkin-test.py script is running.  Just make sure you finish it in time or
  don't pass in --commit and do the commit manually later (or run
  --commit to append the results to the checkin message).

  NOTE: You can do the commit in a second step with a follow-up --commit.

  NOTE: For more details on using this script, see below.


Documentation:
--------------

There are two basic configurations that are tested: MPI_DEBUG and
SERIAL_RELEASE.  Several configure options are varied in these two builds to
try to catch as much conditional configuration behavior has possible.  If
nothing else, please at least do the MPI_DEBUG build since that will cover the
most code and best supports day-to-day development efforts.  However, if you
are changing code that might break the serial build or break non-debug code,
please allow the SERIAL_RELEASE build to be run as well.  Note that the
MPI_DEBUG build actually uses -DCMAKE_BUILD_TYPE=RELEASE with
-DTrilinos_ENABLE_DEBUG=ON to use optimized compiler options but with runtime
debug checking turned on.  This helps to make the tests run faster but still
builds and runs the runtime debug checking code.  Therefore, you should not
use the MPI_DEBUG confgiure options when building a debug version for yourself
to do debugging.

The following steps are performed by this script:

1) Do a CVS update of the code (and detect the files that have changed
locally).  (if --update is set.)

2) Select the list of packages to enable forward based on the directories
where there are changed files (or from a list passed in by the user).  NOTE:
This can be overridden with the options --enable-packages, --disable-packages
and --no-enable-fwd-packages.

3) For each build case (e.g. MPI_DEBUG, SERIAL_RELEASE, etc.)

  3.a) Configure a build directory in a standard way for all of the packages
  that have changed and all of the packages that depend on these packages
  forward. What gets eanbled can be modified (see the enable options above).
  (if --configure is set.)
  
  3.b) Build all configured code with 'make' (e.g. with -jN set through
  --make-options).  (if --build is set.)
  
  3.c) Run all tests for enabled packages.  (if --test is set.)
  
  3.d) Analyze the results of the CVS update, configure, build, and tests and
  send email about results.  (emails only sent out if --send-emails-to is not
  set to ''.)

4) Commit the code given a commit message.  (if --commit is set.)

The recommended way to use this script is to create a new base directory apart
from your standard build directories such as:

  $ cd SOME_BASE_DIR
  $ mkdir CHECKIN

The most basic way to do the checkin test is:

  $ cd SOME_BASE_DIR/CHECKIN
  $ $TRILINOS_HOME/cmake/python/checkin-test.py --do-all

If your MPI installation, other compilers, and standard TPLs (i.e. BLAS and
LAPACK) can be found automatically, then this is all you will need to do.
However, if the setup can not be determined automatically, then you can append
a set of CMake variables that will get read in the files:

  SOME_BASE_DIR/CHECKIN/SERIAL_RELEASE.config
  SOME_BASE_DIR/CHECKIN/MPI_DEBUG.config
  SOME_BASE_DIR/CHECKIN/COMMON.config

Actually, skeletons of these files will automatically be written out with
typical CMake cache variables commented out.  Any CMake cache variables listed
in these files will be read into and passed on the configure line to 'cmake'.

WARNING: Please do not add any CMake cache variables in the *.config files
that will alter what packages are built or what tests are run.  The goal of
these configuration files is to allow you to specify the minimum environment
to find MPI, your compilers, and the basic TPLs.  If you need to fudge what
packages are enabled, please use the script arguments --enable-packages,
--disable-packages, --no-enable-fwd-packages, and/or --enable-all-packages.

NOTE: Before running this script, you should first do an CVS update and
examine what files are changed to make sure you want to commit what you have.
Also, please look out for unknown files that you may need to add to the VC
repository with 'cvs add'.  Your working directory needs to be 100% ready to
commit before running this script.

NOTE: You don't need to run this script if you have not changed any files that
affect the build or the tests.  For example, if all you have changed are
document files, then you don't need to run this script before committing.

NOTE: You can directly call the scripts 'do-configure.base' that get created
by this script in your other build directories in order to match the basic
configuration options that are used for checkin testing.  You just need to set
the other various package and TPL enables for your local work.

Common use cases for using this script are as follows:

(*) Basic full testing without commit:

   --do-all

(*) Basic full testing with commit:

   --do-all --commit --commit-msg-header-file=<SOME_FILE_NAME>

(*) Check commit readiness status:

   [no arguments]

   NOTE: This will examine results for the last testing process and send out
   an email stating if the a commit is ready to perform or not.

(*) Commit after a completed set of tests:

   --commit --commit-msg-header-file=<SOME_FILE_NAME>

   NOTE: This will pick up the results for the last completed test runs and
   append the results of those tests to the checkin-message.

(*) Test only the packages modified and not the forward dependent packages:

  --do-all --no-enable-fwd-packages

  NOTE: This is a safe thing to do when only tests in the modified packages
  are changed and not library code.  This can speed up the testing process and
  is to be preferred over not running this script at all.  It would be very
  hard to make this script automatically determine if only test code has
  changed because every Trilinos package does not follow a set pattern for
  tests and test code.

(*) MPI DEBUG only (run at least this if nothing else):

  --do-all --without-serial-release

(*) The minimum acceptable testing where code has been changed:

  --do-all --no-enable-fwd-packages --without-serial-release

  NOTE: This will do only an MPI DEBUG build and will only build and run the
  tests for the packages that have directly been changed and not forward
  packages.

(*) Test only a specific set of packages and no others:

  --do-all --enable-all-packages=off \
  --enable-packages=<PACKAGEA>,<PACKAGEB>,<PACKAGEC> --no-enable-fwd-packages
  
  NOTE: This will override all logic in the script about which packages will
  be enabled and only the given packages will be enabled.

  NOTE: Using this option is greatly preferred to not running this script at
  all!

(*) See the default option values without doing anything:

  --show-defaults

  NOTE: This is the easiest to figure out what all of the default options are.

"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)


clp.add_option(
  "--update-command", dest="updateCommand", type="string",
  default="cvs -q update -dP",
  help="Command used to update the working copy of Trilinos." )

clp.add_option(
  "--enable-packages", dest="enablePackages", type="string", default="",
  help="List of comma separated Trilinos packages to test changes for" \
  +" (example, 'Teuchos,Epetra').  If this list of packages is empty, then" \
  +" the list of packages to enable will be determined automatically by examining" \
  +" the set of modified files from the version control update log." )

clp.add_option(
  "--disable-packages", dest="disablePackages", type="string", default="",
  help="List of comma separated Trilinos packages to explicitly disable" \
  +" (example, 'Tpetra,NOX').  This list of disables will be appended after" \
  +" all of the listed enables no mater how they are determined (see" \
  +" --enable-packages option).  NOTE: Only use this option to remove packages" \
  +" that will not build for some reason.  You can disable tests that run" \
  +" by using the CTest option -E passed through the --ctest-options argument" \
  +" in this script." )

addOptionParserChoiceOption(
  "--enable-all-packages", "enableAllPackages", ('default', 'on', 'off'), 0,
  "Determine if all Trilinos packages are enabled 'on', or 'off', or let" \
  +" other logic decide 'default'.  Setting to 'off' is appropriate when" \
  +" the logic in this script determines that a global build file has changed" \
  +" but you know that you don't need to rebuild every Trilinos package for" \
  +" a reasonable test.",
  clp )

clp.add_option(
  "--enable-fwd-packages", dest="enableFwdPackages", action="store_true",
  help="Enable forward Trilinos packages." )
clp.add_option(
  "--no-enable-fwd-packages", dest="enableFwdPackages", action="store_false",
  help="Do not enable forward Trilinos packages.", default=True )

clp.add_option(
  "--extra-cmake-options", dest="extraCmakeOptions", type="string", default="",
  help="Extra options to pass to 'cmake' after all other options." \
  +" This should be used only as a last result.  To disable packages, instead use" \
  +" --disable-packages." )

clp.add_option(
  "--make-options", dest="makeOptions", type="string", default="",
  help="The options to pass to make (e.g. -j4)." )

clp.add_option(
  "--ctest-options", dest="ctestOptions", type="string", default="",
  help="Extra options to pass to 'ctest'." )

clp.add_option(
  "--show-all-tests", dest="showAllTests", action="store_true",
  help="Show all of the tests in the summary email." )
clp.add_option(
  "--no-show-all-tests", dest="showAllTests", action="store_false",
  help="Don't show all of the test results in the summary email.",
  default=False )

clp.add_option(
  "--commit", dest="doCommit", action="store_true",
  help="Do the commit at the end if everything works out." \
  + "  Note: You must have SSH public/private keys set up with" \
  + " software.sandia.gov for the commit to happen without having to" \
  + " type your password." )
clp.add_option(
  "--skip-commit", dest="doCommit", action="store_false", default=False,
  help="Skip the commit." )

clp.add_option(
  "--force-commit", dest="forceCommit", action="store_true",
  help="Force the commit even if there are errors.  WARNING: Only do this" \
  +" when you are 100% certain that the errors are not caused by your code" \
  +" changes." )
clp.add_option(
  "--no-force-commit", dest="forceCommit", action="store_false", default=False,
  help="Do not force a commit." )

clp.add_option(
  "--commit-msg-header-file", dest="commitMsgHeaderFile", type="string", default="",
  help="Custom commit message file if commiting with --commit." \
  + "  If an relative path is given, this is expected to be with respect to the" \
  +" base source directory for Trilinos.  The very first line of this file should" \
  +" be the summary line that will be used for the commit." )

clp.add_option(
  "--with-mpi-debug", dest="withMpiDebug", action="store_true",
  help="Do the mpi debug build." )
clp.add_option(
  "--without-mpi-debug", dest="withMpiDebug", action="store_false",
  help="Skip the mpi debug build.", default=True )

clp.add_option(
  "--with-serial-release", dest="withSerialRelease", action="store_true",
  help="Do the serial release build." )
clp.add_option(
  "--without-serial-release", dest="withSerialRelease", action="store_false",
  help="Skip the serial release build.", default=True )

clp.add_option(
  "--rebuild", dest="rebuild", action="store_true",
  help="Keep and existing build tree and simply rebuild on top of it." )
clp.add_option(
  "--from-scratch", dest="rebuild", action="store_false", default=True,
  help="Blow everything away and rebuild from scratch." )

clp.add_option(
  "--send-email-to", dest="sendEmailTo", type="string",
  default=os.environ["USER"]+"@sandia.gov",
  help="List of comma-separated email addresses to send email notification to." \
  +"  By default, this is your Sandia User ID.  In order to turn off email" \
  +" notification, just set --send-email-to='' and no email will be sent." )

clp.add_option(
  "--update", dest="doUpdate", action="store_true",
  help="Do the CVS update of Trilinos.", default=False )

clp.add_option(
  "--configure", dest="doConfigure", action="store_true",
  help="Do the configure of Trilinos.", default=False )

clp.add_option(
  "--build", dest="doBuild", action="store_true",
  help="Do the build of Trilinos.", default=False )

clp.add_option(
  "--test", dest="doTest", action="store_true",
  help="Do the running of the enabled Trilinos tests.", default=False )

clp.add_option(
  "--do-all", dest="doAll", action="store_true",
  help="Do update, configure, build, and test (same as --update --configure" \
  +" --build --test)", default=False )

clp.add_option(
  "--show-defaults", dest="showDefaults", action="store_true",
  help="Show the default option values and do nothing at all.",
  default=False )

(options, args) = clp.parse_args()

# NOTE: Above, in the pairs of boolean options, the *last* add_option(...) 
# takes effect!  That is whay the commands are ordered the way they are!

if options.doCommit and not options.commitMsgHeaderFile:
  print "\nError, if you specify --commit you must also set --commit-msg-header-file!\n"

#
# Echo the commandline
#

print ""
print "**************************************************************************"
print "Script: checkin-test.py \\"
print "  --update-command='"+options.updateCommand+"' \\"
print "  --enable-packages='"+options.enablePackages+"' \\"
print "  --disable-packages='"+options.disablePackages+"' \\"
print "  --enable-all-packages='"+options.enableAllPackages+"'\\"
if options.enableFwdPackages:
  print "  --enable-fwd-packages \\"
else:
  print "  --no-enable-fwd-packages \\"
print "  --extra-cmake-options='"+options.extraCmakeOptions+"' \\"
print "  --make-options='"+options.makeOptions+"' \\"
print "  --ctest-options='"+options.ctestOptions+"' \\"
if options.showAllTests:
  print "  --show-all-tests \\"
else:
  print "  --no-show-all-tests \\"
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
if options.doUpdate:
  print "  --update \\"
if options.doConfigure:
  print "  --configure \\"
if options.doBuild:
  print "  --build \\"
if options.doTest:
  print "  --test \\"
if options.doAll:
  print "  --do-all \\"
if options.doCommit:
  print "  --commit \\"
else:
  print "  --skip-commit \\"
if options.forceCommit:
  print " --force-commit \\"
else:
  print " --no-force-commit \\"


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
    sys.exit(0)
  else:
    sys.exit(1)
