#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# Imports
#

import os
import sys
import traceback
from optparse import OptionParser

if os.environ.get("TRIBITS_CHECKIN_TEST_DEBUG_DUMP", "") == "ON":
  debugDump = True
else:
  debugDump = False

if debugDump:
  print("NOTE: TRIBITS_CHECKIN_TEST_DEBUG_DUMP=ON set in env, doing debug dump ...")

thisFilePath = __file__
if debugDump: print("\nthisFilePath =", thisFilePath)

thisFileRealAbsBasePath = os.path.dirname(os.path.abspath(os.path.realpath(thisFilePath)))
if debugDump: print("\nthisFileRealAbsBasePath = '"+thisFileRealAbsBasePath+"'")

from CheckinTest import *
from GeneralScriptSupport import *
from gitdist import addOptionParserChoiceOption

#
# Utility classes and functions.
#

usageHelp = r"""checkin-test.py [OPTIONS]

This tool does testing of a TriBITS-based project using CTest and this script
can actually do the push itself using git in a safe way.  In fact, it is
recommended that one uses this script to push since it will amend the last
commit message with a (minimal) summary of the builds and tests run with
results and/or send out a summary email about the builds/tests performed.


QUICKSTART
-----------

In order to do a safe push, perform the following recommended workflow
(different variations on this workflow are described in the COMMON USE CASES
section below):

1) Commit changes in the local repo:

  # 1.a) See what files are changed, newly added, etc.
  $ git status

  # 1.b) Stage the files you want to commit
  $ git stage <files you want to commit>

  # 1.c) Create your local commits
  $ git commit -- SOMETHING
  $ git commit -- SOMETHING_ELSE
  ...

  # 1.d) Stash whatever changes are left you don't want to test/push
  $ git stash

  NOTE: You can group your commits any way that you would like (see the basic
  git documentation).

  NOTE: When multiple repos are involved, use the 'gitdist' command instead of
  'git'.  This script is provided at tribits/python_utils/gitdist.  See
  gitdist --help for details.

2) Review the changes that you have made to make sure it is safe to push:

  $ cd $PROJECT_HOME
  $ git local-stat | less               # Look at the full status of local repo
  $ git diff --name-status HEAD ^@{u}   # [Optional] Look at the files that have changed

  NOTE: The command 'local-stat' is a git alias that can be installed with the
  script tribits/python_utils/git-config-alias.sh.  This command is
  recommended over just a raw 'git status' or 'git log' to review commits
  before attempting to test/push commits.  If you have not installed these
  alias, then run the following commands instead:

    $ git status
    $ git log --oneline --name-status HEAD ^@{u}

  NOTE: If you see any files/directories that are listed as 'unknown' returned
  from 'git local-stat', then you will need to do an 'git add' to track them or
  add them to an ignore list *before* you run the checkin-test.py script.
  The git script will not allow you to push if there are new 'unknown' files or
  uncommitted changes to tracked files.

  NOTE: When multiple repos are involved, use 'gitdist-mod-status' to see the
  state of your repos before pushing.  See gitdist --help for details.

3) Set up the checkin base build directory (first time only):

  $ cd $PROJECT_HOME
  $ echo CHECKIN >> .git/info/exclude
  $ mkdir CHECKIN
  $ cd CHECKIN

  NOTE: You may need to set up some configuration files if CMake cannot find
  the right compilers, MPI, and TPLs by default (see detailed documentation
  below).

  NOTE: You might want to set up a simple shell driver script for your common
  use cases.

  NOTE: You can set up a CHECKIN directory of any name in any location
  you want.  If you create one outside of the main source dir, then
  you will not have to add the git exclude shown above.

4) Do the pull, configure, build, test, and push:

  $ cd $PROJECT_HOME
  $ cd CHECKIN
  $ ../checkin-test.py -j4 --do-all --push

  NOTE: The above will: a) pull updates from the tracking branch, b)
  automatically enable the correct packages based on changed files, c)
  configure and build the changed and downstream packages, d) run the tests,
  e) send you emails about what happened, f) do a final pull from the global
  repo, g) optionally amend the last local commit with the test results, and
  h) finally push local commits to the tracking branch if everything passes.

  NOTE: The repo must be on a branch (not a detached head state) with a
  tracking branch '<remoterepo>/<remotebranch>' so that a raw 'git pull' can
  be performed to get updates and to push changes.  Also, the push is done
  explicitly to the tracking branch using:

    git push <remoterepo> <remotebranch>

  NOTE: You must not have any uncommitted changes or the script will stop
  right away.  To run the script, you will may need to first use 'git stash'
  to stash away your unstagged/uncommitted changes *before* running this
  script.

  NOTE: You need to have SSH public/private keys set up to the remote repo
  machines for the git commands invoked in the script to work without you
  having to type a password.  If the 'git pull' fails, the detailed output is
  generally in a pull*.out file (see the checkin-test.out log file for
  details).

  NOTE: You can do the final push in a second invocation of the script with a
  follow-up run with --push and removing --do-all (it will remember the
  results from the build/test cases just ran).  For more details, see detailed
  documentation below.

  NOTE: Once you start running the checkin-test.py script, you can go off and
  do something else and just check your email to see if all the builds and
  tests passed and if the push happened or not.

  NOTE: The commands 'cmake', 'ctest', and 'make' must be in your default path
  before running this script.

  NOTE: Defaults like -j4 can be set using a local-checkin-test-defaults.py
  file (see below).

For more details on using this script, see the detailed documentation below.


DETAILED DOCUMENTATION
-----------------------

The following approximate steps are performed by this script:

----------------------------------------------------------------------------

1) Check to see if the local repo(s) are clean:

  $ git status

  NOTE: If any modified or any unknown files are shown, the process will be
  aborted.  The local repo(s) working directory must be clean and ready to
  push *everything* that is not stashed away.

2) Do a 'git pull' to update the repo (done if --pull or --do-all is set):

  NOTE: If not doing a pull, use --allow-no-pull or --local-do-all.

3) Select the list of packages to enable forward/downstream based on the
package directories where there are changed files (or from a list of packages
passed in by the user).

  NOTE: The automatic enable behavior can be overridden or modified using the
  options --enable-all-packages=on/off, --enable-packages=<p0>,<p1>,... ,
  --disable-packages=<p0>,<p1>,... , and/or --no-enable-fwd-packages.

4) For each build/test case <BUILD_NAME> (e.g. MPI_DEBUG, SERIAL_RELEASE,
extra builds specified with --st-extra-builds and --extra-builds):

  4.a) Configure a build directory <BUILD_NAME> in a standard way for all of
  the packages that have changed and all of the packages that depend on these
  packages forward/downstream. You can manually select which packages get
  enabled (see the enable options above).  (done if --configure, --do-all, or
  --local-do-all is set.)

  4.b) Build all configured code with 'make' (e.g. with -jN set through
  -j or --make-options).  (done if --build, --do-all, or --local-do-all is set.)

  4.c) Run all BASIC tests for enabled packages.  (done if --test, --do-all,
  or --local-do-all is set.)

  4.d) Analyze the results of the pull, configure, build, and tests and send
  email about results.  (emails only sent out if --send-emails-to!="")

5) Do final pull and rebase, append test results to last commit message, and
push (done if --push or --do-all is set)

  5.a) Do a final 'git pull' (done if --pull or --do-all is set)

  5.b) Do 'git rebase <remoterepo>/<remotebranch>' (done if --rebase is set)

  5.c) Amend commit message of the most recent commit with the summary of the
  testing performed.  (done if --append-test-results is set.)

  5.d) Push the local commits to the global repo (done if --push is set)

6) Send out final on actions (i.e. 'DID PUSH' email if a push occurred).
(done if --send-email-to!="" is set and --send-email-only-on-failure is *not*
set)

----------------------------------------------------------------------------

The recommended way to use this script is to create a new base CHECKIN test
directory apart from your standard build directories such as with:

  $ $PROJECT_HOME
  $ mkdir CHECKIN
  $ echo CHECKIN >> .git/info/exclude

The most basic way to do pre-push testing is with:

  $ cd CHECKIN
  $ ../checkin-test.py --do-all [other options]

If your MPI installation, other compilers, and standard TPLs can be found
automatically, then this is all you will need to do.  However, if the setup
cannot be determined automatically, then you can add a set of CMake variables
that will get read in the files:

  COMMON.config
  MPI_DEBUG.config
  SERIAL_RELEASE.config

(or whatever your standard --default-builds are).

Actually, for built-in build/test cases, skeletons of these files will
automatically be written out with typical CMake cache variables (commented
out) that you would need to set out.  Any CMake cache variables listed in
these files will be read into and passed on the configure line to 'cmake'.

WARNING: Please do not add any extra CMake cache variables than what are
needed to get the Primary Tested (PT) --default-builds builds to work.  Adding
other enables/disables will make the builds non-standard and can break these
PT builds.  The goal of these configuration files is to allow you to specify
the minimum environment to find MPI, your compilers, and the required TPLs
(e.g. BLAS, LAPACK, etc.).  If you need to fudge what packages are enabled,
please use the script arguments --enable-packages, --enable-extra-pacakges,
--disable-packages, --no-enable-fwd-packages, and/or --enable-all-packages to
control this, not the *.config files!

WARNING: Please do not add any CMake cache variables in the *.config files
that will alter what packages or TPLs are enabled or what tests are run.
Actually, the script will not allow you to change TPL enables in these
standard *.config files because to do so deviates from a consistent build
configuration for Primary Tested (PT) Code.

NOTE: All tentatively-enabled TPLs (e.g. Pthreads and BinUtils) are hard
disabled in order to avoid different behaviors between machines where they
would be enabled and machines where they would be disabled.

NOTE: If you want to add extra build/test cases that do not conform to
the standard build/test configurations described above, then you need
to create extra builds with the --extra-builds and/or
--st-extra-builds options (see below).

NOTE: Before running this script, you should first do an 'git status' and 'git
diff --name-status HEAD ^@{u}' and examine what files are changed to make sure
you want to push what you have in your local working directory.  Also, please
look out for unknown files that you may need to add to the git repository with
'git add' or add to your ignores list.  There cannot be any uncommitted
changes in the local repo before running this script.

NOTE: You don't need to run this script if you have not changed any files that
affect the build or the tests.  For example, if all you have changed are
documentation files, then you don't need to run this script before pushing
manually.

NOTE: To see detailed debug-level information, set
TRIBITS_CHECKIN_TEST_DEBUG_DUMP=ON in the env before running this script.


COMMON USE CASES (EXAMPLES):
----------------------------

(*) Basic full testing with integrating with global repo(s) without push:

  ../checkin-test.py --do-all

  NOTE: This will result in a set of emails getting sent to your email address
  for the different configurations and an overall push readiness status email.

  NOTE: If everything passed, you can follow this up with a --push (see
  below).

(*) Basic full testing with integrating with local repo and push:

  ../checkin-test.py --do-all --push

  NOTE: By default this will rebase your local commits and amend the last
  commit with a short summary of test results.  This is appropriate for
  pushing commits that only exist in your local repo and are not shared with
  any remote repo.

(*) Push to global repo after a completed set of tests have finished:

  ../checkin-test.py [other options] --push

  NOTE: This will pick up the results for the last completed test runs with
  [other options] and append the results of those tests to the log of the most
  recent commit.

  NOTE: Take the action options for the prior run and replace --do-all with
  --push but keep all of the rest of the options the same.  For example, if
  you did:

    ../checkin-test.py --enable-packages=Blah --default-builds=MPI_DEBUG --do-all

  then follow that up with:

    ../checkin-test.py --enable-packages=Blah --default-builds=MPI_DEBUG --push

  NOTE: This is a common use case when some tests are failing which aborted
  the initial push but you determine it is okay to push anyway and do so with
  --force-push.

(*) Test only the packages modified and not the forward dependent packages:

  ../checkin-test.py --do-all --no-enable-fwd-packages

  NOTE: This is a safe thing to do when only tests in the modified packages
  are changed and not library code.  This can speed up the testing process and
  is to be preferred over not running this script at all.  It would be very
  hard to make this script automatically determine if only test code has
  changed because every package does not follow a set pattern for
  tests and test code.

(*) Run the most important default (e.g. MPI_DEBUG) build/test only:

  ../checkin-test.py --do-all --default-builds=MPI_DEBUG

(*) The minimum acceptable testing when code has been changed:

  ../checkin-test.py \
    --do-all --enable-all-packages=off --no-enable-fwd-packages \
    --default-builds=MPI_DEBUG

  NOTE: This will do only an MPI DEBUG build and will only build and run the
  tests for the packages that have directly been changed and not any forward
  packages.  Replace "MPI_DEBUG" with whatever your most important default
  build is.

(*) Test only a specific set of packages and no others:

  ../checkin-test.py \
    --enable-packages=<P0>,<P1>,<P2> --no-enable-fwd-packages \
    --do-all

  NOTE: This will override all logic in the script about which packages will
  be enabled based on file changes and only the given packages will be
  enabled.  When there are tens of thousands of changed files and hundreds of
  defined packages, this auto-detection algorithm can be very expensive!

  NOTE: You might also want to pass in --enable-all-packages=off in case the
  script wants to enable all the packages (see the output in the
  checkin-test.py log file for details) and you think it is not necessary to
  do so.

  NOTE: Using these options is greatly preferred to not running this script at
  all and should not be any more expensive than the testing you would already
  do manually before a push.

(*) Test changes locally without pulling updates:

  ../checkin-test.py --local-do-all

  NOTE: This will just configure, build, test, and send an email notification
  without updating or changing the status of the local git repo in any way and
  without any communication with the global repo.  Hence, you can have
  uncommitted changes and still run configure, build, test without having to
  commit or having to stash changes.

  NOTE: This will determine what packages to enable and test based on changes
  w.r.t. to the tracking branch.  If not on a tracking branch, or in a
  detached head state, see below.

  NOTE: This is typically not a sufficient level of testing in order to push
  the changes to a shared branch because you have not fully integrated your
  changes yet with other developers.  However, this would be a sufficient
  level of testing in order to do a commit on the local machine and then pull
  to a remote machine for further testing and a push (see below).

(*) Local test of repo version on a detached head or with no tracking branch:

  ../checkin-test.py --enable-all-packages=[on|off] \
    --enable-packages=<P0>,... --local-do-all

  By specifying what packages are enabled and not doing a pull or push, the
  script allows the repo(s) to be in a detached head state or on a branch that
  does not have a tracking branch.  This allows the checkin-test.py script to
  be used, for example, to test versions using 'git bisect'.

(*) Adding extra build/test cases:

  Often you will be working on Secondary Tested (ST) Code or Experimental (EX)
  Code and want to include the testing of this in your pre-push testing
  process along with the standard --default-builds build/test cases which can
  only include Primary Tested (PT) Code.  In this case you can run with:

    ../checkin-test.py --extra-builds=<BUILD1>,<BUILD2>,... [other options]

  For example, if you have a build that enables the TPL CUDA you would do:

    echo "
    -DTPL_ENABLE_MPI:BOOL=ON
    -DTPL_ENABLE_CUDA:BOOL=ON
    " > MPI_DEBUG_CUDA.config

  and then run with:

    ../checkin-test.py --extra-builds=MPI_DEBUG_CUDA --do-all

  This will do the standard --default-builds (e.g. MPI_DEBUG and
  SERIAL_RELEASE) build/test cases along with your non-standard MPI_DEBUG_CUDA
  build/test case.

  NOTE: You can disable the default build/test cases with --default-builds="".
  However, please only do this when you are not going to push because you need
  at least one default build/test case (the most important default PT case,
  e.g. MPI_DEBUG) to do a safe push.

(*) Including extra repos and extra packages:

  You can also use the checkin-test.py script to continuously integrate
  multiple git repos containing add-on packages. To do so, just run:

    ../checkin-test.py --extra-repos=<REPO1>,<REPO2>,... [options]

  NOTE: You have to create local commits in all of the extra repos where there
  are changes or the script will abort.

  NOTE: Extra repos can be specified with more flexibility using the
  --extra-repos-file and --extra-repos-type arguments (also see
  --ignore-missing-extra-repos).

  NOTE: Each of the last local commits in each of the changed repos will get
  amended with the appended summary of what was enabled in the build/test (if
  --append-test-results is set).

(*) Avoid changing any of the local commit SHA1s:

  If you are pushing commits from a shared branch, it is critical that you do
  not change any of the SHA1s of the commits.  Changing the SHA1s for any of
  the commits will mess up various multi-repo, multi-branch workflows.  To
  avoid changing any of the SHA1s of the local commits, one must run with:

    ../checkin-test.py --no-rebase --no-append-test-results [options]

(*) Performing a remote test/push:

  If you develop on a slow machine like your laptop, doing an appropriate
  level of testing can take a long time.  In this case, you can pull the
  changes to another faster remote workstation and do a more complete set of
  tests and push from there.  If you are knowledgeable with git, this will be
  easy and natural to do, without any help from this script.  However, this
  script can still help and automate the steps and can do so in one command
  invocation on the part of the developer.

  On your slow local development machine 'mymachine', do the limited testing
  with:

    ../checkin-test.py --do-all --no-enable-fwd-packages

  On your fast remote test machine, do a full test and push with:

    ../checkin-test.py \
      --extra-pull-from=<remote-repo>:master \
      --do-all --push

  where <remote-name> is a git remote repo name pointing to
  mymachine:/some/dir/to/your/src (see 'git help remote').

  NOTE: You can of course adjust the packages and/or build/test cases that get
  enabled on the different machines.

  NOTE: Once you invoke the checkin-test.py script on the remote test machine
  and it has pulled the commits from mymachine, then you can start changing
  files again on your local development machine and just check your email to
  see what happens on the remote test machine.

  NOTE: If something goes wrong on the remote test machine, you can either
  work on fixing the problem there or you can fix the problem on your local
  development machine and then do the process over again.

  NOTE: If you alter the commits on the remote machine (such as squashing
  commits), you will have trouble merging back on our local machine.
  Therefore, if you have to to fix problems, make new commits and don't alter
  the ones you pulled from your local machine (but rebasing them should be
  okay as long as the local commits on mymachine are not pushed to other
  repos).

  NOTE: Git will resolve the duplicated commits when you pull the commits
  pushed from the remote machine.  Git knows that the commits are the same and
  will do the right thing when rebasing (or just merging).

  NOTE: This would also work for multiple repos if the remote name
  '<remote-repo>' pointed to the right remote repo in all the local repos.

(*) Check push readiness status:

  ../checkin-test.py

  This will examine results for the last testing process and send out an email
  stating if the a push is ready to perform or not.

(*) See the default option values without doing anything:

  ../checkin-test.py --show-defaults

  This is the easiest way to figure out what all of the default options are.

Hopefully the above documentation, the example use cases, the documentation of
the command-line arguments below, and some experimentation will be enough to
get you going using this script for all of your pre-push testing and pushes.
If that is not sufficient, send email to your development support team to ask
for help.


LOCAL DEFAULT COMMAND LINE DEFAULTS
-----------------------------------

If the file local-checkin-test-defaults.py exists in the current directory,
then it will be read in and will change the project defaults for the
command-line arguments.  For example, a valid local-checkin-test-defaults.py
file would look like:

  defaults = [
    "-j10",
    "--no-rebase",
    "--ctest-options=-E '(PackageA_Test1|PackageB_Test2)'"
    ]

Any of the project's checkin-test.py command-line argument defaults can be
changed in this way.  The updated defaults can be observed by running:

  ./checkin-test.py --show-defaults

Any command-line arguments explicitly passed in will override these local
defaults.


HANDLING OF PT, ST, AND EX CODE IN BUILT-IN AND EXTRA BUILDS:
-------------------------------------------------------------

This script will only process PT (Primary Tested) packages in the
--default-builds (e.g. MPI_DEBUG and SERIAL_RELEASE) builds.  This is to avoid
problems of side-effects of turning on ST packages that would impact PT
packages (e.g. an ST package getting enabled that enables an ST TPL which
turns on support for that TPL in a PT package producing different code which
might work but the pure PT build without the extra TPL may actually be broken
and not know it).  Therefore, any non-PT packages that are enabled (either
implicitly through changed files or explicitly by listing in --enable-packages)
will be turned off in the --default-builds builds.  If none of the enabled
packages are PT, then they will all be disabled and the --default-builds
builds will be skipped.

In order to better support the development of ST and EX packages, this script
allows you to define some extra builds that will be invoked and used to
determine overall pass/fail before a potential push.  The option
--st-extra-builds is used to specify extra builds that will test ST packages
(and also PT packages if any are enabled).  If only PT packages are enabled
then the builds specified in --st-extra-builds will still be run.  The
reasoning is that PT packages may contain extra ST features and therefore if
the goal is to test these ST builds it is desirable to also run these builds
because they also my impact downstream ST packages.

Finally, the option --extra-builds will test all enabled packages, including
EX packages, regardless of their test group.  Therefore, when using
--extra-builds, be careful that you watch what packages are enabled.  If you
change an EX package, it will be enabled in --extra-builds builds.

A few use cases might help better demonstrate the behavior.  Consider
the following input arguments specifying extra builds

  --st-extra-builds=MPI_DEBUG_ST --extra-builds=INTEL_DEBUG

with the packages Teuchos, Phalanx, and Meros where Teuchos is PT, Phalanx is
ST, and Meros is EX.

Here is what packages would be enabled in each of the builds:

  --default-builds=MPI_DEBUG,SERIAL_RELEASE \
  --st-extra-builds=MPI_DEBUG_ST \
  --extra-builds=INTEL_DEBUG

and which packages would be excluded:

A) --enable-packages=Teuchos:
   MPI_DEBUG:       [Teuchos]
   SERIAL_RELEASE:  [Teuchos]
   MPI_DEBUG_ST:    [Teuchos]
   INTEL_DEBUG:     [Teuchos]     Always enabled!

B) --enable-packages=Phalanx:
   MPI_DEBUG:       []            Skipped, no PT packages!
   SERIAL_RELEASE:  []            Skipped, no PT packages!
   MPI_DEBUG_ST:    [Phalanx]
   INTEL_DEBUG:     [Phalanx]

C) --enable-packages=Meros:
   MPI_DEBUG:       []            Skipped, no PT packages!
   SERIAL_RELEASE:  []            Skipped, no PT packages!
   MPI_DEBUG_ST:    []            Skipped, no PT or ST packages!
   INTEL_DEBUG:     [Meros]

D) --enable-packages=Teuchos,Phalanx:
   MPI_DEBUG:       [Teuchos]
   SERIAL_RELEASE:  [Teuchos]
   MPI_DEBUG_ST:    [Teuchos,Phalanx]
   INTEL_DEBUG:     [Teuchos,Phalanx]

E) --enable-packages=Teuchos,Phalanx,Meros:
   MPI_DEBUG:       [Teuchos]
   SERIAL_RELEASE:  [Teuchos]
   MPI_DEBUG_ST:    [Teuchos,Phalanx]
   INTEL_DEBUG:     [Teuchos,Phalanx,Meros]

The --extra-builds=INTEL_DEBUG build is always performed with all of the
enabled packages.  This logic given above must be understood in order to
understand the output given in the script.


CONVENTIONS FOR COMMAND-LINE ARGUMENTS:
---------------------------------------

The command-line arguments are segregated into three broad categories: a)
action commands, b) aggregate action commands, and c) others.

a) The action commands are those such as --build, --test, etc. and are shown
with [ACTION] in their documentation.  These action commands have no off
complement.  If the action command appears, then the action will be performed.

b) Aggregate action commands such as --do-all and --local-do-all turn on sets
of other action commands and are shown with [AGGR ACTION] in their
documentation.  The sub-actions that these aggregate action commands turn on
cannot be disabled with other arguments.

c) Other arguments are those that are not marked with [ACTION] or [AGGR
ACTION] tend to either pass in data and turn control flags on or off.


EXIT CODE:
---------

This script returns 0 if the actions requested are successful.  This does not
necessarily imply that it is okay to do a push or that a push was done.  For
example, if only --pull is passed in and is successful, then 0 will be
returned but that does *not* mean that it is okay to do a push.  Therefore, a
return value of is a necessary but not sufficient condition for readiness to
push, it depends on the requested actions.

"""

# ToDo: Break up the above huge documentation block into different "topics"
# and then display those topics with --help-topic=<topic>.  Also provide a
# --help-all that will combine all of the --help-topic documentation with the
# standard documentation to produce where is there now.

def runProjectTestsWithCommandLineArgs(commandLineArgs, configuration = {}):

  clp = ConfigurableOptionParser(configuration.get('defaults', {}), usage=usageHelp)

  clp.add_option(
    "--project-configuration", dest="projectConfiguration", type="string", default="",
    help="Custom file to provide configuration defaults for the project." \
      + "  By default, the file project-checkin-test-config.py is looked for" \
      + " in <checkin-test-path> (in case it is symlinked into <projectDir>/checkin-test.py)" \
      + " if not found there, then it is looked for in <checkin-test-path>/../../.." \
      +"  (assuming default TriBITS snapshot <projectDir>/cmake/tribits/ci_support/)" \
      + " If this file is set to a location that is not in the" \
      + " project's base directory, then --src-dir must be set to point to the" \
      + " project's base directory."
    )

  clp.add_option(
    "--show-defaults", dest="showDefaults", action="store_true",
    help="Show the default option values and do nothing at all.",
    default=False )

  clp.add_option(
    "--project-name", dest="projectName", action="store",
    help="Set the project's name. This is used to locate various files."+\
      "  If not set, then it reads the project name from the PROJECT_NAME"+\
      " variable set in the file SRCDIR/ProjectName.cmake.",
    default=None)

  clp.add_option(
    '--src-dir', dest="srcDir", type="string", default="",
    help="The source base directory for code to be tested.  The default is determined" \
     +" by the location of the found project-checkin-test-config.py file." )

  configuredBuilds = [build for build, unused in
                      configuration.get('cmake', {}).get('default-builds', [])]
  clp.add_option(
    '--default-builds', dest='defaultBuilds', type='string',
    default=','.join(configuredBuilds),
    help="Comma separated list of builds that should always be run by default.")

  clp.add_option(
    "--extra-repos-file", dest="extraReposFile", type="string", default="",
    help="File path to an extra repositories list file.  If set to 'project', then " \
    +"<project_dir>/cmake/ExtraRepositoriesList.cmake is read.  See the argument " \
    +"--extra-repos for details on how this list is used (default empty '')")

  g_extraRepoTypesList = [""]
  g_extraRepoTypesList.extend(g_knownTribitsTestRepoTypes)

  addOptionParserChoiceOption(
    "--extra-repos-type", "extraReposType", g_extraRepoTypesList, 0,
    "The test type of repos to read from <extra_repos_file>.",
    clp )

  clp.add_option(
    "--extra-repos", dest="extraRepos", type="string", default="",
    help="List of comma separated extra repositories " \
    +"containing extra  packages that can be enabled.  The order these repos is "
    +"listed in not important.  This option overrides --extra-repos-file.")

  clp.add_option(
    "--ignore-missing-extra-repos", dest="ignoreMissingExtraRepos", action="store_true",
    help="If set, then extra repos read in from <extra_repos_file> will be ignored " \
    +"and removed from list.  This option is not applicable if <extra_repos_file>=='' " \
    +"or <extra_repos_type>==''." )
  clp.add_option(
    "--require-extra-repos-exist", dest="ignoreMissingExtraRepos", action="store_false",
    default=False,
    help="If set, then all listed extra repos must exist or the script will exit. [default]" )

  clp.add_option(
    "--with-cmake", dest="withCmake", type="string", default="cmake",
    help="CMake executable to use with cmake -P scripts internally (only set" \
    +" by unit testing code).")

  clp.add_option(
    "--skip-deps-update", dest="skipDepsUpdate", action="store_true",
    help="If set, skip the update of the dependency XML file.  If the package structure" \
      " has not changed since the last invocation, then it is safe to use this option.",
    default=False )

  clp.add_option(
    "--enable-packages", dest="enablePackages", type="string", default="",
    help="List of comma separated packages to test changes for" \
    +" (example, 'Teuchos,Epetra').  If this list of packages is empty, then" \
    +" the list of packages to enable will be determined automatically by examining" \
    +" the set of modified files from the version control update log.  Note that"\
    +" this will skip the auto-detection of changed packages based on changed"\
    +" files." )

  clp.add_option(
    "--enable-extra-packages", dest="enableExtraPackages", type="string", default="",
    help="List of comma separated packages to test in addition to the packages" \
    +" that are enabled determined automatically by examining" \
    +" the set of modified files from the version control update log.  This option"\
    +" is mostly just used in ACI sync servers." )

  clp.add_option(
    "--disable-packages", dest="disablePackages", type="string", default="",
    help="List of comma separated packages to explicitly disable" \
    +" (example, 'Tpetra,NOX').  This list of disables will be appended after" \
    +" all of the listed enables no matter how they are determined (see" \
    +" --enable-packages option).  NOTE: Only use this option to remove packages" \
    +" that will not build for some reason.  You can disable tests that run" \
    +" by using the CTest option -E passed through the --ctest-options argument" \
    +" in this script." )

  addOptionParserChoiceOption(
    "--enable-all-packages", "enableAllPackages", ('auto', 'on', 'off'), 0,
    "Determine if all packages are enabled 'on', or 'off', or 'auto'" \
    +" (let other logic decide).  Setting to 'off' is appropriate when" \
    +" the logic in this script determines that a global build file has changed" \
    +" but you know that you don't need to rebuild and test every package for" \
    +" a reasonable test.  Setting --enable-packages effectively disables this" \
    +" option.  Setting this to 'off' does *not* stop the forward enabling" \
    +" of downstream packages for packages that are modified or set by --enable-packages."\
    +" Setting this to 'on' will skip the automatic detection of changed packages"\
    +" based on changed files.  It can be helpful to stop the auto-detection changed"\
    +" packages when there are thousands of changed files and hundreds of defined"\
    +" packages." ,
    clp )

  clp.add_option(
    "--enable-fwd-packages", dest="enableFwdPackages", action="store_true",
    help="Enable forward packages. [default]" )
  clp.add_option(
    "--no-enable-fwd-packages", dest="enableFwdPackages", action="store_false",
    help="Do not enable forward packages.", default=True )

  clp.add_option(
    "--continue-if-no-updates", dest="abortGracefullyIfNoChangesPulled", action="store_false",
    help="If set, then the script will continue if no updates are pulled from any repo. [default]",
    default=False )
  clp.add_option(
    "--abort-gracefully-if-no-changes-pulled", dest="abortGracefullyIfNoChangesPulled", action="store_true",
    help="If set, then the script will abort gracefully if no updates are pulled from any repo.",
    default=False )

  clp.add_option(
    "--continue-if-no-changes-to-push", dest="abortGracefullyIfNoChangesToPush", action="store_false",
    help="If set, then the script will continue if no changes to push from any repo. [default]",
    default=False )
  clp.add_option(
    "--abort-gracefully-if-no-changes-to-push", dest="abortGracefullyIfNoChangesToPush", action="store_true",
    help="If set, then the script will abort gracefully if no changes to push from any repo.",
    default=False )

  clp.add_option(
    "--continue-if-no-enables", dest="abortGracefullyIfNoEnables", action="store_false",
    help="If set, then the script will continue if no packages are enabled. [default]",
    default=False )
  clp.add_option(
    "--abort-gracefully-if-no-enables", dest="abortGracefullyIfNoEnables", action="store_true",
    help="If set, then the script will abort gracefully if no packages are enabled.",
    default=False )

  clp.add_option(
    "--extra-cmake-options", dest="extraCmakeOptions", type="string",
    default=configuration.get('extra-cmake-options', ''),
    help="Extra options to pass to 'cmake' after all other options." \
    +" This should be used only as a last resort.  To disable packages, instead use" \
    +" --disable-packages.  To change test categories, use --test-categories." )

  clp.add_option(
    "--test-categories", dest="testCategories", type="string",
    default="BASIC",
    help="." \
    +" Change the test categories.  Can be 'BASIC', 'CONTINUOUS', " \
      " 'NIGHTLY', or 'HEAVY' (default set by project, see --show-defaults)." )

  clp.add_option(
    "-j", "--parallel", dest="overallNumProcs", type="string", default="",
    help="The options to pass to make and ctest (e.g. -j4)." )

  clp.add_option(
    "--use-makefiles", dest="useNinja", action="store_false",
    help="If set, then -G'Unix Makfiles' used for backend build tool." \
    +" Note: The command 'make' must be in the default path. [default]",
    default=False )
  clp.add_option(
    "--use-ninja", dest="useNinja", action="store_true",
    help="If set, then -GNinja used for backend build tool." \
    +" Note: The command 'ninja' must be in the default path." ,
    default=False )

  clp.add_option(
    "--make-options", dest="makeOptions", type="string", default="",
    help="The options to pass to 'make' (e.g. -j4) or ninja" \
    +" (if --use-ninja given)." )

  clp.add_option(
    "--ctest-options", dest="ctestOptions", type="string", default="",
    help="Extra options to pass to 'ctest' (e.g. -j2)." )

  clp.add_option(
    "--ctest-timeout", dest="ctestTimeOut", type="float", default=300,
    help="timeout (in seconds) for each single 'ctest' test (e.g. 180" \
    +" for three minutes).  This sets the CMake cache var DART_TESTING_TIMEOUT"
    +" which becomes the default timeout for tests, even when running raw"
    +" ctest.  This value can be overridden using the ctest argument --timeout."
    +"  Individual tests may have their own timeouts set which will not be"
    +" impacted by this default global timeout.  See the configure variable"
    +" <Project>_SCALE_TEST_TIMEOUT to scale up timeouts for"
    +" all tests, even those that have individuals timeouts set." )

  clp.add_option(
    "--show-all-tests", dest="showAllTests", action="store_true",
    help="Show all of the tests in the summary email and in the commit message" \
    +" summary (see --append-test-results)." )
  clp.add_option(
    "--no-show-all-tests", dest="showAllTests", action="store_false",
    help="Don't show all of the test results in the summary email. [default]",
    default=False )

  clp.add_option(
    "--without-default-builds", dest="withoutDefaultBuilds", action="store_true",
    default=False,
      help="Skip the default builds (same as --default-builds='')." \
      +"  You would use option along with --extra-builds=BUILD1,BUILD2,... to run your own" \
      +" local custom builds." )

  clp.add_option(
    "--st-extra-builds", dest="stExtraBuilds", type="string", default="",
    help="List of comma-separated ST extra build names.  For each of the build names in" \
    +" --st-extra-builds=<BUILD1>,<BUILD2>,..., there must be a file <BUILDN>.config in" \
    +" the local directory along side the COMMON.config file that defines the special" \
    +" build options for the extra build." )

  clp.add_option(
    "--ss-extra-builds", dest="ssExtraBuilds", type="string", default="",
    help="DEPRECATED!  Use --st-extra-builds instead!. (Default empty "")" )

  clp.add_option(
    "--extra-builds", dest="extraBuilds", type="string", default="",
    help="List of comma-separated extra build names.  For each of the build names in" \
    +" --extra-builds=<BUILD1>,<BUILD2>,..., there must be a file <BUILDN>.config in" \
    +" the local directory along side the COMMON.config file that defines the special" \
    +" build options for the extra build." )

  clp.add_option(
    "--log-file", dest="logFile", type="string", default="checkin-test.out",
    help="File used for detailed log info." )

  clp.add_option(
    "--send-email-to", dest="sendEmailTo", type="string",
    default=getCmndOutput("git config --get user.email", True, False),
    help="List of comma-separated email addresses to send email notification to" \
    +" after every build/test case finishes and at the end for an overall summary" \
    +" and push status." \
    +"  By default, this is the email address you set for git returned by" \
    +" `git config --get user.email`.  In order to turn off email" \
    +" notification, just set --send-email-to='' and no email will be sent." )

  clp.add_option(
    "--skip-case-send-email", dest="skipCaseSendEmail", action="store_true",
    help="If set then if a build/test case is skipped for some reason (i.e." \
    +" because no packages are enabled) then an email will go out for that case." \
    +" [default]" )
  clp.add_option(
    "--skip-case-no-email", dest="skipCaseSendEmail", action="store_false",
    help="If set, then if a build/test case is skipped for some reason (i.e." \
    +" because no packages are enabled) then no email will go out for that case." \
    +" (opposite of --skip-case-send-email) [default]",
    default=True )

  addOptionParserChoiceOption(
    "--send-build-case-email", "sendBuildCaseEmail",
    ('always', 'only-on-failure', 'never'), 0,
    "Determines when email goes out to --send-email-to=<email> for a build" \
    +" case.  But the final status email will still go out if --send-email-to=<email>" \
    +" is not empty. [default = 'always']",
    clp )

  clp.add_option(
    "--send-email-for-all", dest="sendEmailOnlyOnFailure", action="store_false",
    help="If set, then emails will get sent out for all operations. [default]" )
  clp.add_option(
    "--send-email-only-on-failure", dest="sendEmailOnlyOnFailure", action="store_true",
    help="If set, then emails will only get sent out for failures.",
    default=False )

  clp.add_option(
    "--send-email-to-on-push", dest="sendEmailToOnPush", type="string",
    default=configuration.get('SendEmailOnPush', ''),
    help="List of comma-separated email addresses to send email notification to" \
    +" on a successful push.  This is used to log pushes to a central list." \
    +"  In order to turn off this email" \
    +" notification, just set --send-email-to-on-push='' and no email will be sent" \
    +" to these email lists." )

  clp.add_option(
    "--force-push", dest="forcePush", action="store_true",
    help="Force the local push even if there are build/test errors." \
    +" WARNING: Only do this when you are 100% certain that the errors are not" \
    +" caused by your code changes.  This only applies when --push is specified" \
    +" and this script.")
  clp.add_option(
    "--no-force-push", dest="forcePush", action="store_false", default=False,
    help="Do not force a push if there are failures. [default]" )

  clp.add_option(
    "--do-push-readiness-check", dest="doPushReadinessCheck", action="store_true",
    help="Check the push readiness status at the end and send email if not actually" \
    +" pushing. [default]" )
  clp.add_option(
    "--skip-push-readiness-check", dest="doPushReadinessCheck", action="store_false",
    default=True,
    help="Skip push status check." )

  clp.add_option(
    "--rebase", dest="rebase", action="store_true",
    help="Rebase the local commits on top of <remoterepo>/<remotebranch> before amending" \
    +" the last commit and pushing.  Rebasing keeps a nice linear commit" \
    +" history like with CVS or SVN and will work perfectly for the basic" \
    +" workflow of adding commits to the 'master' branch and then syncing" \
    +" up with <remoterepo>/<remotebranch> before the final push. [default]" )
  clp.add_option(
    "--no-rebase", dest="rebase", action="store_false",
    help="Do *not* rebase the local commits on top of <remoterepo>/<remotebranch> before" \
    +" amending the final commit and pushing.  This allows for some more " \
    +" complex workflows involving local branches with multiple merges." \
    +"  However, this will result in non-linear history and will allow for" \
    +" trivial merge commits with <remoterepo>/<remotebranch> to get pushed.  This mode" \
    +" should only be used in cases where the rebase mode will not work or " \
    +" when it is desired to use a merge commit to integrate changes on a" \
    +" branch that you wish be able to easily back out.  For sophisticated" \
    +" users of git, this may in fact be the preferred mode.",
    default=True )

  clp.add_option(
    "--append-test-results", dest="appendTestResults", action="store_true",
    help="Before the final push, amend the most recent local commit by appending a" \
    +" summary of the test results.  This provides a record of what builds" \
    +" and tests were performed in order to test the local changes.  This is only " \
    +" performed if --push is also set.  NOTE: If the same" \
    +" local commit is amended more than once, the prior test summary sections will be" \
    +" overwritten with the most recent test results from the current run. [default]" )
  clp.add_option(
    "--no-append-test-results", dest="appendTestResults", action="store_false",
    help="Do not amend the last local commit with test results.  NOTE: If you have" \
    +" uncommitted local changes that you do not want this script to commit then you" \
    +" must select this option to avoid this last amending commit.  Also, if you are" \
    +" pushing commits from a shared branch and don't want to change any of the SHA1s" \
    +" for the commits, then you must set this option!",
    default=True )

  clp.add_option(
    "--extra-pull-from", dest="extraPullFrom", type="string", default="",
    help="Optional extra git pull(s) to merge in changes from after" \
    +" pulling in changes from the tracking branch.  The format of this argument is:" \
    +" ...,<local-repoi>:<remote-repoi>:<remote-branchi>,... where each pull" \
    +" specification gives the name (not the directory) of the local repo <local-repoi>" \
    +", the remote repo name <remote-repoi>, and the branch in the remote repo to pull" \
    +" <remote-branchi>.  If only two semicolons ':' are given then an pull field takes" \
    +" the form ...,<remote-repo>:<remote-branch>,... where the remote <remote-name>" \
    +" must be defined in all the repos and the branch <remote-branch> must exist" \
    +" in all the remote repos.  If the <remote-repoi> is empty such as with" \
    +" ...,:<remote-repoi>:<remote-branchi>,... then this matches the base git repo." \
    +"  The extra pull(s) are only done if --pull is also specified.  NOTE: when using" \
    +" --extra-repos=<repo0>,<repo1>,... the <local-repoi> must be a named repository" \
    + " that is present in all of the git repos or it will be an error." )

  clp.add_option(
    "--allow-no-pull", dest="allowNoPull", action="store_true", default=False,
    help="Allowing for there to be no pull performed and still doing the other actions." \
    +"  This option is useful for testing against local changes without having to" \
    +" get the updates from the global repo.  However, if you don't pull, you can't" \
    +" push your changes to the global repo.  WARNING: This does *not* stop a pull" \
    +" attempt from being performed by --pull or --do-all!" )

  clp.add_option(
    "--wipe-clean", dest="wipeClean", action="store_true", default=False,
    help="[ACTION] Blow existing build directories and build/test results.  The action can be" \
    +" performed on its own or with other actions in which case the wipe clean will be" \
    +" performed before any other actions. NOTE: This will only wipe clean the builds" \
    +" that are specified and will not touch those being ignored (e.g. SERIAL_RELEASE" \
    +" will not be removed if --default-builds=MPI_DEBUG is specified)." )

  clp.add_option(
    "--pull", dest="doPull", action="store_true", default=False,
    help="[ACTION] Do the pull from the tracking branch and optionally also" \
      +" merge in changes from the repo pointed to by --extra-pull-from.")

  clp.add_option(
    "--configure", dest="doConfigure", action="store_true", default=False,
    help="[ACTION] Do the configure step." )

  clp.add_option(
    "--build", dest="doBuild", action="store_true", default=False,
    help="[ACTION] Do the build step." )

  clp.add_option(
    "--test", dest="doTest", action="store_true", default=False,
    help="[ACTION] Do the running of the enabled tests." )

  clp.add_option(
    "--local-do-all", dest="localDoAll", action="store_true", default=False,
    help="[AGGR ACTION] Do configure, build, and test with no pull (same as setting" \
    +" --allow-no-pull ---configure --build --test)." \
    +"  This is the same as --do-all except it does not do --pull and also allows for no pull." )

  clp.add_option(
    "--do-all", dest="doAll", action="store_true", default=False,
    help="[AGGR ACTION] Do update, configure, build, and test (same as --pull --configure" \
    +" --build --test).  NOTE: This will do a --pull regardless if --allow-no-pull" \
    +" is set or not.  To avoid the pull, use --local-do-all." )

  clp.add_option(
    "--push", dest="doPush", action="store_true", default=False,
    help="[ACTION] Push the committed changes in the local repo into to remote repo" \
      +" pointed to by the tracking branch." )

  clp.add_option(
    "--execute-on-ready-to-push", dest="executeOnReadyToPush", type="string", default="",
    help="[ACTION] A command to execute on successful execution and 'READY TO PUSH'" \
    +" status from this script.  This can be used to do a remote SSH invocation to a" \
    +" remote machine to do a remote pull/test/push after this machine finishes." )

  (options, args) = clp.parse_args(args=commandLineArgs)

  # NOTE: Above, in the pairs of boolean options, the *last* add_option(...)
  # takes effect!  That is why the commands are ordered the way they are!


  #
  # Echo the command-line
  #

  print("")
  print("**************************************************************************")
  print("Script: checkin-test.py \\")

  print("  --src-dir='" + options.srcDir+"' \\")
  print("  --default-builds='" + options.defaultBuilds + "' \\")
  print("  --extra-repos-file='"+options.extraReposFile+"' \\")
  print("  --extra-repos-type='"+options.extraReposType+"' \\")
  print("  --extra-repos='"+options.extraRepos+"' \\")
  if options.ignoreMissingExtraRepos:
    print("  --ignore-missing-extra-repos \\")
  else:
    print("  --require-extra-repos-exist \\")
  if options.skipDepsUpdate:
    print("  --skip-deps-update \\")
  print("  --enable-packages='"+options.enablePackages+"' \\")
  print("  --enable-extra-packages='"+options.enableExtraPackages+"' \\")
  print("  --disable-packages='"+options.disablePackages+"' \\")
  print("  --enable-all-packages='"+options.enableAllPackages+"'\\")
  if options.enableFwdPackages:
    print("  --enable-fwd-packages \\")
  else:
    print("  --no-enable-fwd-packages \\")
  if options.abortGracefullyIfNoChangesPulled:
    print("  --abort-gracefully-if-no-changes-pulled \\")
  else:
    print("  --continue-if-no-updates \\")
  if options.abortGracefullyIfNoChangesToPush:
    print("  --abort-gracefully-if-no-changes-to-push \\")
  else:
    print("  --continue-if-no-changes-to-push \\")
  if options.abortGracefullyIfNoEnables:
    print("  --abort-gracefully-if-no-enables \\")
  else:
    print("  --continue-if-no-enables \\")
  print("  --extra-cmake-options='"+options.extraCmakeOptions+"' \\")
  print("  --test-categories='"+options.testCategories+"' \\")
  if options.overallNumProcs:
    print("  -j"+options.overallNumProcs+" \\")
  if options.useNinja:
    print("  --use-ninja \\")
  else:
    print("  --use-makefiles \\")
  print("  --make-options='"+options.makeOptions+"' \\")
  print("  --ctest-options='"+options.ctestOptions+"' \\")
  print("  --ctest-timeout="+str(options.ctestTimeOut)+" \\")
  if options.showAllTests:
    print("  --show-all-tests \\")
  else:
    print("  --no-show-all-tests \\")
  if options.withoutDefaultBuilds:
    print("  --without-default-builds \\")
  print("  --st-extra-builds='"+options.stExtraBuilds+"' \\")
  print("  --extra-builds='"+options.extraBuilds+"' \\")
  print("  --log-file='"+options.logFile+"' \\")
  print("  --send-email-to='"+options.sendEmailTo+"' \\")
  if options.skipCaseSendEmail:
    print("  --skip-case-send-email \\")
  else:
    print("  --skip-case-no-email \\")
  print("  --send-build-case-email="+str(options.sendBuildCaseEmail)+" \\")
  if not options.sendEmailOnlyOnFailure:
    print("  --send-email-for-all \\")
  else:
    print("  --send-email-only-on-failure \\ ")
  print("  --send-email-to-on-push='"+options.sendEmailToOnPush+"' \\")
  if options.forcePush:
    print("  --force-push \\")
  else:
    print("  --no-force-push \\")
  if options.doPushReadinessCheck:
    print("  --do-push-readiness-check \\")
  else:
    print("  --skip-push-readiness-check \\")
  if options.rebase:
    print("  --rebase \\")
  else:
    print("  --no-rebase \\")
  if options.appendTestResults:
    print("  --append-test-results \\")
  else:
    print("  --no-append-test-results \\")
  if options.extraPullFrom:
    print("  --extra-pull-from='"+options.extraPullFrom+"' \\")
  if options.allowNoPull:
    print("  --allow-no-pull \\")
  if options.wipeClean:
    print("  --wipe-clean \\")
  if options.doPull:
    print("  --pull \\")
  if options.doConfigure:
    print("  --configure \\")
  if options.doBuild:
    print("  --build \\")
  if options.doTest:
    print("  --test \\")
  if options.localDoAll:
    print("  --local-do-all \\")
  if options.doAll:
    print("  --do-all \\")
  if options.doPush:
    print("  --push \\")
  if options.executeOnReadyToPush:
    print("  --execute-on-ready-to-push=("+options.executeOnReadyToPush+") \\")
  if options.ssExtraBuilds:
    print("  --ss-extra-builds='"+options.ssExtraBuilds+"' \\")
    print("\nWARNING: --ss-extra-builds is deprecated!  Use --st-extra-builds instead!")
    if options.stExtraBuilds:
      print("ERROR: Can't set deprecated --ss-extra-builds and --st-extra-builds together!")
      sys.exit(3)
    options.stExtraBuilds = options.ssExtraBuilds


  #
  # Check and adjust the input arguments
  #

  if not os.path.isabs(options.srcDir):
    options.srcDir = os.path.abspath(options.srcDir)

  if options.doAll and options.localDoAll:
    print("\nError, you can not use --do-all and --local-do-all together!  Use on or the other!")
    sys.exit(1)

  if options.doAll and options.allowNoPull:
    print("\nError, you can not use --do-all and --allow-no-pull together! (see the" \
      " documentation for the --do-all, --local-do-all, and --allow-no-pull arguments.)")
    sys.exit(2)

  if options.extraPullFrom:
    for extraPullFromArg in options.extraPullFrom.split(","):
       # Will validate form of argument
      getLocalRepoRemoteRepoAndBranchFromExtraPullArg(extraPullFromArg)


  #
  # Execute the checkin test guts
  #

  import time

  if not options.showDefaults:

    print("\nStarting time:", getCmndOutput("date",True))

    tribitsDir = os.path.abspath(getCompleteFileDirname(__file__)+"/..") # TriBITS dir!

    t1 = time.time()
    success = checkinTest(tribitsDir, options, configuration)
    t2 = time.time()
    print("\nTotal time for checkin-test.py =", formatMinutesStr((t2-t1)/60.0))

    print("\nFinal time:", getCmndOutput("date",True))

    if options.ssExtraBuilds:
      print("\n***")
      print("*** FINAL WARNING: stop using deprecated --ss-extra-builds!  Use --st-extra-builds instead!")
      print("***")

    if success:
      print("\nREQUESTED ACTIONS: PASSED\n")
      return True
    else:
      print("\nREQUESTED ACTIONS: FAILED\n")
      return False
  else:
    return True


def getConfigurationSearchPaths():
  """
  Gets a list of paths to search for the configuration. If this file
  was invoked from a symlink, look in the directory that contains the
  symlink. The returned list will always contain at least one element.
  """
  result = []

  # First, look for the checkin-test.py file's base directory path. It is
  # common practice to symbolically link the checkin-test.py script into the
  # project's base source directory.  NOTE: Don't use realpath here!  We don't
  # want to follow symbolic links!
  result.append(os.path.dirname(os.path.abspath(__file__)))

  # Second, look for the configuration file assuming the checkin-test.py
  # script is run out of the standard snapshotted tribits directory
  # <project-root>/cmake/tribits/ci_support
  result.append(os.path.join(thisFileRealAbsBasePath, '..', '..', '..'))

  return result


def loadConfigurationFile(filepath):
  if os.path.exists(filepath):
    sys_path_old = sys.path
    try:
      modulePath = os.path.dirname(filepath)
      moduleFile = os.path.basename(filepath)
      moduleName, extension = os.path.splitext(moduleFile)
      sys.path = [modulePath] + sys_path_old
      try:
        if debugDump:
          print("\nLoading project configuration from %s..." % filepath)
          print("\nsys.path =", sys.path)
        configuration = __import__(moduleName).configuration
        if debugDump:
          print("\nSetting the default --src-dir='"+modulePath+"'")
        configuration.get("defaults").update({"--src-dir" : modulePath})
        return configuration
      except Exception as e:
        print(e)
        raise e
    finally:
      sys.path = sys_path_old
      if debugDump:
        print("\nsys.path =", sys.path)
  else:
    raise Exception('The file %s does not exist.' % filepath)


def getLocalCmndLineArgs():
  localDefaults = []
  checkinTestDir = os.getcwd()
  localProjectDefaultsBaseName = "local-checkin-test-defaults"
  localProjectDefaultsFile = checkinTestDir+"/"+localProjectDefaultsBaseName+".py"
  if os.path.exists(localProjectDefaultsFile):
    sys_path_old = sys.path
    try:
      sys.path = [checkinTestDir] + sys_path_old
      if debugDump:
        print("\nLoading local default command-line args from "+localProjectDefaultsFile+"...")
        print("\nsys.path =", sys.path)
      localDefaults = __import__(localProjectDefaultsBaseName).defaults
    finally:
      sys.path = sys_path_old
      if debugDump:
        print("\nsys.path =", sys.path)
  return localDefaults


def locateAndLoadConfiguration(path_hints = []):
  """
  Locate and load a module called checkin_test_project_configuration.py. The
  path_hints argument can be used to provide location hints at which to locate
  the file. Returns a configuration dictionary. If the module is not found,
  this dictionary will be empty.
  """
  for path in path_hints:
    candidate = os.path.join(path, "project-checkin-test-config.py")
    if debugDump: print("\nLooking for candidate configuration file '%s'" % candidate)
    if os.path.exists(candidate):
      return loadConfigurationFile(candidate)
  return {}


#
# Main
#

def main(cmndLineArgs):

  # See if the help option is set or not
  helpOpt = len( set(cmndLineArgs) & set(("--help", "-h")) ) > 0

  # See if --show-defaults was set or not
  showDefaultsOpt = len( set(cmndLineArgs) & set(("--show-defaults", "dummy")) ) > 0


  if (not helpOpt) and (not showDefaultsOpt):
    # Get the name of the log file from the input arguments
    logFileName = "checkin-test.out"
    for cmndLineArg in cmndLineArgs:
      cmndLineSplit = cmndLineArg.split("=")
      if cmndLineSplit[0].strip() == "--log-file":
        logFileName = cmndLineSplit[1].strip()
    logFile = open(logFileName, "w", buffering=1) # 1 == line-by-line buffering
    # NOTE: Above, we need to set biffering=1 to make sure that each line that
    # gets written is written out the file.  The is critical so that the user
    # can to a `tail -f <log-file>` to see what is going one.
  else:
    logFile = None

  # There are a lot of print(statements in the implementation. It's
  # easier to reset sys.stdout and sys.stderr to a TeeOutput object
  # than to replace them.
  teeOutput = TeeOutput(logFile)
  originalStdout = sys.stdout
  originalStderr = sys.stderr
  try:
    sys.stdout = teeOutput
    sys.stderr = teeOutput
    try:
      # See if there is a configuration file override.
      configuration = None
      for arg in cmndLineArgs:
        if arg.startswith('--project-configuration='):
          print("Found configuration override %s..." % arg)
          configuration = loadConfigurationFile(arg.split('=')[1])
        elif not configuration and arg.startswith('--src-dir='):
          configuration = locateAndLoadConfiguration([arg.split('=')[1]])
      if not configuration:
        configuration = locateAndLoadConfiguration(getConfigurationSearchPaths())
      localCmndLineArgs = getLocalCmndLineArgs()
      if localCmndLineArgs:
        if debugDump:
          print("\ncmndLineArgs =", cmndLineArgs)
          print("\nlocalCmndLineArgs =", localCmndLineArgs)
        cmndLineArgs = localCmndLineArgs + cmndLineArgs
        if debugDump:
          print("\ncmndLineArgs =", cmndLineArgs)
      if debugDump:
        print("\nConfiguration loaded from configuration file =", configuration)
      success = runProjectTestsWithCommandLineArgs(cmndLineArgs, configuration)
    except SystemExit as e:
      # In Python 2.6, SystemExit inherits Exception, but for proper exit
      # behavior the SystemExit exception must propagate all the way to the top
      # of the call stack. It cannot get handled by the catch Exception below.
      raise e
    except Exception as e:
      success = False
      traceback.print_exc(file=teeOutput)
  finally:
    # Close the log file
    if logFile: logFile.close()
    # Reset stdout and stderr
    sys.stdout = originalStdout
    sys.stderr = originalStderr

  if success:
    return 0
  else:
    return 1

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
