#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


from CheckinTest import *


#
# Read in the command-line arguments
#

usageHelp = r"""checkin-test.py [OPTIONS]

This tool does checkin testing for Trilinos with CMake/CTest and can actually
do the push itself using eg/git in a safe way.  In fact, it is recommended
that you use this script to push since it will amend the last commit message
with a (minimal) summary of the builds and tests run with results.


Quickstart:
-----------

In order to do a solid checkin, perform the following recommended workflow
(different variations on this workflow are described below):

1) Commit changes in the local repo:

  # 1.a) See what files are changed, newly added, etc. that need to be committed
  # or stashed.
  $ eg status

  # 1.b) Stage the files you want to commit (optional)
  $ eg stage <files you want to commit>

  # 1.c) Create your local commits
  $ eg commit -- SOMETHING
  $ eg commit -- SOMETHING_ELSE
  ...

  # 1.d) Stash whatever changes are left you don't want to test/push (optional)
  $ eg stash

  NOTE: You can group your commits any way that you would like (see the basic
  eg/git documentation).

  NOTE: If not installed on your system, the eg script can be found at
  Trilinos/commonTools/git/eg.  Just add it to your path.

2) Review the changes that you have made to make sure it is safe to push:

  $ cd $TRILINOS_HOME
  $ eg local-stat                  # Look at the full status of local repo
  $ eg diff --name-status origin   # [Optional] Look at the files that have changed

  NOTE: The command 'local-stat' is an alias that can be installed with the
  script Trilinos/commonTools/git/git-config-alias.sh.  It is highly
  recommended over just a raw 'eg status' or 'eg log' to review commits before
  attempting to test/push commits.

  NOTE: If you see any files/directories that are listed as 'unknown' returned
  from 'eg local-stat', then you will need to do an 'eg add' to track them or
  add them to an ignore list *before* you run the checkin-test.py script.
  The eg script will not allow you to push if there are new 'unknown' files or
  uncommitted changes to tracked files.

3) Set up the checkin base build directory (first time only):

  $ cd $TRILINOS_HOME
  $ echo CHECKIN >> .git/info/exclude
  $ mkdir CHECKIN
  $ cd CHECKIN

  NOTE: You may need to set up some configuration files if CMake can not find
  the right compilers, MPI, and TPLs by default (see detailed documentation
  below).

  NOTE: You might want to set up a simple shell driver script.  See some
  examples in the files:

    Trilinos/sampmleScripts/checkin-test-*

  NOTE: You can set up a CHECKIN directory of any name in any location you
  want.  If you create one outside of the main Trilinos source dir, then you
  will not have to add the git exclude shown above.

4) Do the checkin build, test, and push:

  $ cd $TRILINOS_HOME
  $ cd CHECKIN
  $ ../checkin-test.py -j4 --do-all --push

  NOTE: The above will: a) pull updates from the global repo, b) automatically
  enable the correct packages, c) build the code, d) run the tests, e) send
  you emails about what happened, f) do a final pull to from the global repo,
  g) amend the last local commit with the test results, and h) finally push
  local commits to the global repo if everything passes.

  NOTE: You must have installed the official versions of eg/git with the
  install-git.py script in order to run this script.  If you don't, the script
  will die right away with an error message telling you what the problem is.

  NOTE: The current branch will be used to pull and push to.  A raw 'eg pull'
  is performed which will get all of the branches from 'origin'.  This means
  that your current branch must be a tracking branch so that it will get
  updated correctly.  The branch 'master' is the most common branch but
  release tracking branches are also common.

  NOTE: You must not have any uncommitted changes or the 'eg pull && eg rebase
  --against origin' command will fail on the final pull/rebase before the push
  and therefore the whole script will fail.  To still run the script, you will
  may need to use 'eg stash' to stash away your unstaged/uncommitted changes
  *before* running this script.

  NOTE: You need to have SSH public/private keys set up to software.sandia.gov
  for the git commands invoked in the script to work without you having to
  type a password.

  NOTE: You can do the final push in a second invocation of the script with a
  follow-up run with --push and removing --do-all (it will remember the
  results from the build/test cases just ran).  For more details, see detailed
  documentation below.

  NOTE: Once you start running the checkin-test.py script, you can go off and
  do something else and just check your email to see if all the builds and
  tests passed and if the push happened or not.

For more details on using this script, see the detailed documentation below.


Detailed Documentation:
-----------------------

There are two basic configurations that are tested by default:
MPI_DEBUG and SERIAL_RELEASE.  Both of these configurations only test
Primary Stable Code (see --extra-builds, --ss-extra-builds and
'Handling of PS and SS Code' for more details)).  Several configure
options are varied in these two builds to try to catch as much
conditional configuration behavior as possible.  If nothing else,
please at least do the MPI_DEBUG build since that will cover the most
code and best supports day-to-day development efforts.  However, if
you are changing code that might break the serial build or break
non-debug code, please allow the SERIAL_RELEASE build to be run as
well.  Note that the MPI_DEBUG build actually uses
-DCMAKE_BUILD_TYPE=RELEASE with -DTrilinos_ENABLE_DEBUG=ON to use
optimized compiler options but with runtime debug checking turned on.
This helps to make the tests run faster but still builds and runs the
runtime debug checking code.  Therefore, you should not use the
MPI_DEBUG configure options when building a debug version for yourself
to do debugging.


The following approximate steps are performed by this script:


----------------------------------------------------------------------------


1) Check to see if the local repo is clean:

  $ eg status

  NOTE: If any modified or any unknown files are shown, the process will be
  aborted.  The local repo working directory must be clean and ready to push
  *everything* that is not stashed away.

2) Do a 'eg pull' to update the code (done if --pull or --do-all is set):

  NOTE: If not doing a pull, use --allow-no-pull or --local-do-all.

3) Select the list of packages to enable forward based on the package
directories where there are changed files (or from a list of packages passed
in by the user).

  NOTE: The automatic enable behavior can be overridden or modified using the
  options --enable-packages, --disable-packages, and/or
  --no-enable-fwd-packages.

4) For each build/test case <BUILD_NAME> (e.g. MPI_DEBUG, SERIAL_RELEASE,
extra builds specified with --extra-builds):

  4.a) Configure a build directory <BUILD_NAME> in a standard way for all of
  the packages that have changed and all of the packages that depend on these
  packages forward. You can manually select which packages get enabled (see
  the enable options above).  (done if --configure, --do-all, or
  --local-do-all is set.)
  
  4.b) Build all configured code with 'make' (e.g. with -jN set through
  -j or --make-options).  (done if --build, --do-all, or --local-do-all is set.)
  
  4.c) Run all BASIC tests for enabled packages.  (done if --test, --do-all,
  or --local-do-all is set.)

    NOTE: By default, only Trilinos_TEST_CATEGORIES=BASIC tests are enabled.
  
  4.d) Analyze the results of the update, configure, build, and tests and send
  email about results.  (emails only sent out if --send-emails-to is not set
  to ''.)

5) Do final pull, append test results to last commit message, and push (done
if --push is set)

  5.a) Do a final 'eg pull && eg rebase --against origin/<current_branch>'
  (done if --pull or --do-all is set)

    NOTE: The final 'eg rebase --against origin/<current_branch>' is required
    to avoid trival merge commits that the Trilinos global get repo will
    reject on the push.
  
  5.b) Amend commit message of the most recent commit with the summary of the
  testing performed.  (done if --append-test-results is set.)
  
  5.c) Push the local commits to the global repo.


----------------------------------------------------------------------------


The recommended way to use this script is to create a new base CHECKIN test
directory apart from your standard build directories such as with:

  $ $TRILINOS_HOME
  $ mkdir CHECKIN
  $ echo CHECKIN >> .git/info/exclude

The most basic way to do the checkin test is:

  $ cd CHECKIN
  $ ../checkin-test.py --do-all [other options]

If your MPI installation, other compilers, and standard TPLs (i.e. BLAS and
LAPACK) can be found automatically, then this is all you will need to do.
However, if the setup cannot be determined automatically, then you can add a
set of CMake variables that will get read in the files:

  COMMON.config
  MPI_DEBUG.config
  SERIAL_RELEASE.config

Actually, for built-in build/test cases, skeletons of these files will
automatically be written out with typical CMake cache variables (commented
out) that you would need to set out.  Any CMake cache variables listed in
these files will be read into and passed on the configure line to 'cmake'.

WARNING: Please do not add any CMake cache variables than what are needed to
get the MPI_DEBUG and SERIAL_RELEASE builds to work.  Adding other
enables/disables will make the builds non-standard and break the Primary
Stable build.  The goal of these configuration files is to allow you to
specify the minimum environment to find MPI, your compilers, and the required
TPLs (e.g. BLAS, LAPACK, etc.).  If you need to fudge what packages are
enabled, please use the script arguments --enable-packages,
--disable-packages, --no-enable-fwd-packages, and/or --enable-all-packages to
control this, not the *.config files!

WARNING: Please do not add any CMake cache variables in the *.config files
that will alter what packages or TPLs are enabled or what tests are run.
Actually, the script will not allow you to change TPL enables in these
standard *.config files because to do so deviates from a consistent build
configuration for Primary Stable Code.

NOTE: All tentatively-enabled TPLs (e.g. Pthreads and BinUtils) are hard
disabled in order to avoid different behaviors between machines where they
would be enabled and machines where they would be disabled.

NOTE: If you want to add extra build/test cases that do not conform to
the standard build/test configurations described above, then you need
to create extra builds with the --extra-builds and/or
--ss-extra-builds options (see below).

NOTE: Before running this script, you should first do an 'eg status' and 'eg
diff --name-status origin..' and examine what files are changed to make sure
you want to push what you have in your local working directory.  Also, please
look out for unknown files that you may need to add to the git repository with
'eg add' or add to your ignores list.  There cannot be any uncommitted changes
in the local repo before running this script.

NOTE: You don't need to run this script if you have not changed any files that
affect the build or the tests.  For example, if all you have changed are
documentation files, then you don't need to run this script before pushing
manually.


Common Use Cases (examples):
----------------------------

(*) Basic full testing with integrating with global repo without push:

  ../checkin-test.py --do-all

  NOTE: This will result in a set of emails getting sent to your email address
  for the different configurations and an overall push readiness status email.

  NOTE: If everything passed, you can follow this up with a --push (see
  below).

(*) Basic full testing with integrating with local repo and push:

  ../checkin-test.py --do-all --push

(*) Push to global repo after a completed set of tests have finished:

  ../checkin-test.py [other options] --push

  NOTE: This will pick up the results for the last completed test runs with
  [other options] and append the results of those tests to the checkin-message
  of the most recent commit.

  NOTE: Take the action options for the prior run and replace --do-all with
  --push but keep all of the rest of the options the same.  For example, if
  you did:

    ../checkin-test.py --enable-packages=Blah --without-serial-release --do-all

  then follow that up with:

    ../checkin-test.py --enable-packages=Blah --without-serial-release --push

  NOTE: This is a common use case when some tests are failing which aborted
  the initial push but you determine it is okay to push anyway and do so with
  --force-push (or just --force for short).

(*) Test only the packages modified and not the forward dependent packages:

  ../checkin-test.py --do-all --no-enable-fwd-packages

  NOTE: This is a safe thing to do when only tests in the modified packages
  are changed and not library code.  This can speed up the testing process and
  is to be preferred over not running this script at all.  It would be very
  hard to make this script automatically determine if only test code has
  changed because every Trilinos package does not follow a set pattern for
  tests and test code.

(*) Run the MPI_DEBUG build/test only:

  ../checkin-test.py --do-all --without-serial-release

(*) The minimum acceptable testing when code has been changed:

  ../checkin-test.py \
    --do-all --enable-all-packages=off --no-enable-fwd-packages \
    --without-serial-release

  NOTE: This will do only an MPI DEBUG build and will only build and run the
  tests for the packages that have directly been changed and not any forward
  packages.

(*) Test only a specific set of packages and no others:

  ../checkin-test.py \
    --enable-packages=<PACKAGEA>,<PACKAGEB>,<PACKAGEC> --no-enable-fwd-packages \
    --do-all
  
  NOTE: This will override all logic in the script about which packages will
  be enabled and only the given packages will be enabled.

  NOTE: You might also want to pass in --enable-all-packages=off in case the
  script wants to enable all the packages (see the output in the
  checkin-test.py log file for details) and you think it is not necessary to
  do so.

  NOTE: Using these options is greatly preferred to not running this script at
  all and should not be any more expensive than the testing you already do.

(*) Test changes locally without pulling updates:

  ../checkin-test.py --local-do-all

  NOTE: This will just configure, build, test, and send an email notification
  without updating or changing the status of the local git repo in any way and
  without any communication with the global repo.  Hence, you can have
  uncommitted changes and still run configure, build, test without having to
  commit or having to stash changes.

  NOTE: This is not a sufficient level of testing in order to push the changes
  to the global repo because you have not fully integrated your changes yet
  with other developers.  However, this would be a sufficient level of testing
  in order to do a commit on the local machine and then pull to a remote
  machine for further testing and a push (see below).

(*) Adding extra build/test cases:

  Often you will be working on Secondary Stable Code or Experimental Code and
  want to include the testing of this in your pre-checkin testing along with
  the standard MPI_DEBUG and SERIAL_RELEASE build/test cases which can only
  include Primary Stable Code.  In this case you can run with:
  
    ../checkin-test.py --extra-builds=<BUILD1>,<BUILD2>,... [other options]
  
  For example, if you have a build that enables the TPL CUDA for Tpetra you
  would do:
  
    echo "
    -DTPL_ENABLE_MPI:BOOL=ON
    -DTPL_ENABLE_CUDA:BOOL=ON
    " > MPI_DEBUG_CUDA.config
  
  and then run with:
  
    ../checkin-test.py \
      --enable-packages=Tpetra --extra-builds=MPI_DEBUG_CUDA --do-all
  
  This will do the standard MPI_DEBUG and SERIAL_RELEASE build/test cases
  along with your non-standard MPI_DEBUG_CUDA build/test case.

  NOTE: You can disable the default build/test cases with
  --without-default-builds.  However, please only do this when you are not
  going to push because we need at least one default build/test case to be
  safe to push.

(*) Including extra repos:

  You can also use the checkin-test.py script to continuously integrate with
  other external extra git repos containing add-on Trilinos packages (such as
  preCopyrightTrilinos).  To do so, just run:
    
    ../checkin-test.py --extra-builds=REPO1,REPO2,... [options]

  where REPOI = preCopyrightTrilinos, LIMEExt, etc.

  NOTE: You have to create local commits in all of the extra repos where there
  are changes or the script will abort.

  NOTE: Each of the last local commits in each of the changed repos will get
  amended with the appended summary of what was enabled in the build/test.

(*) Performing a remote test/push:

  If you develop on a slow machine like your laptop, doing an appropriate
  level of testing can take a long time.  In this case, you can pull the
  changes to another faster remote workstation and do a more complete
  set of tests and push from there.

  On your slow local development machine 'mymachine', do the limited testing
  with:

    ../checkin-test.py --do-all --no-enable-fwd-packages
  
  On your fast remote test machine, do a full test and push with:
  
    ../checkin-test.py \
      --extra-pull-from=mymachine:/some/dir/to/your/trilinos/src:master \
      --do-all --push
  
  NOTE: You can of course adjust the packages and/or build/test cases that get
  enabled on the different machines.
  
  NOTE: Once you invoke the checkin-test.py script on the remote test machine,
  you can start changing files again on your local development machine and
  just check your email to see what happens.
  
  NOTE: If something goes wrong on the remote test machine, you can either
  work on fixing the problem there or you can fix the problem on your local
  development machine and then do the process over again.

  NOTE: If you alter the commits on the remote machine (such as squashing
  commits), you will have trouble merging back on our local machine.
  Therefore, if you have to to fix problems, make new commits and don't alter
  the ones you pulled from your local machine.

  NOTE: Git will resolve the duplicated commits when you pull the commits
  pushed from the remote machine.  Git knows that the commits are the same and
  will do the right thing.
  
(*) Check push readiness status:

  ../checkin-test.py

  NOTE: This will examine results for the last testing process and send out an
  email stating if the a push is ready to perform or not.

(*) See the default option values without doing anything:

  ../checkin-test.py --show-defaults

  NOTE: This is the easiest way to figure out what all of the default options
  are.

Hopefully the above documentation, the example use cases, the documentation of
the command-line arguments below, and some experimentation will be enough to
get you going using this script for all of pre-checkin testing and pushes.  If
that is not sufficient, send email to trilinos-framework@software.sandia.gov
to ask for help.


Handling of PS and SS Code:
---------------------------

This script will only process PS (Primary Stable) packages in the
default MPI_DEBUG and SERIAL_RELEASE builds.  This is to avoid
problems of side-effects of turning on SS packages that would impact
PS packages (e.g. SS Phalanx getting enabled that enables SS Boost
which turns on support for Boost in PS Teuchos producing different
code which might work but the pure PS build without Boost of Teuchos
may actually be broken and not know it).  Therefore, any non PS
packages that are enabled (either implicity through changed files or
explicitly in --enable-packages) will be turned off in the MP_DEBUG
and SERIAL_RELEASE builds.  If none of the enabled packages are PS,
then they will all be disabled and the MPI_DEBUG and SERIAL_RELEASE
builds will be skipped.

In order to better support the development of SS and EX packages, this
script allows you to define some extra builds that will be invoked and
used to determine overall pass/fail before a potential push.  The
option --ss-extra-builds is used to specify extra builds that will
test SS packages (and also PS packages if any are enabled).  However,
if no SS packages are enabled, the builds are skipped.

Finally, the option --extra-builds will test all enabled packages
reguardless of their categorization.

A few use cases might help better demonstrate the behavior.  Consider
the following input arguments specifying extra builds:

  --ss-extra-packages=MPI_DEBUG_SS --extra-packages=INTEL_DEBUG

and the packages Techos, Phalanx, and Meros where Teuchos is PS,
Phalanx is SS, and Meros is EX.

Here is what packages would be enabled in each of the builds
MPI_DEBUG, SERIAL_RELEASE, MPI_DEBUG_SS, and INTEL_DEBUG and which
builds would be skipped:

A) --enable-packages=Teuchos:
   MPI_DEBUG:       [Teuchos]
   SERIAL_RELEASE:  [Teuchos]
   MPI_DEBUG_SS:    [Teuchos]     Skipped, no SS packages!
   INTEL_DEBUG:     [Teuchos]     Always enabled!

B) --enable-packages=Phalanx:
   MPI_DEBUG:       []            Skipped, no PS packages!
   SERIAL_RELEASE:  []            Skipped, no PS packages!
   MPI_DEBUG_SS:    [Phalanx]
   INTEL_DEBUG:     [Phalanx]

C) --enable-packages=Meros:
   MPI_DEBUG:       []            Skipped, no PS packages!
   SERIAL_RELEASE:  []            Skipped, no PS packages!
   MPI_DEBUG_SS:    []            Skipped, no SS packages!
   INTEL_DEBUG:     [Meros]

D) --enable-packages=Teuchos,Phalanx:
   MPI_DEBUG:       [Teuchos]
   SERIAL_RELEASE:  [Teuchos]
   MPI_DEBUG_SS:    [Teuchos,Phalanx]
   INTEL_DEBUG:     [Teuchos,Phalanx]

E) --enable-packages=Teuchos,Phalanx,Meros:
   MPI_DEBUG:       [Teuchos]
   SERIAL_RELEASE:  [Teuchos]
   MPI_DEBUG_SS:    [Teuchos,Phalanx]
   INTEL_DEBUG:     [Teuchos,Phalanx,Meros]

Above, the --ss-enable-builds=MPI_DEBUG_SS build was skipped in use
case 'A' because it did not contain any SS packages beyond what was in
the set of PS packages.  However, the --extra-builds=INTEL_DEBUG build
is always performed with all of the enabled packages.  This logic
given above in order to understand the output given in the script.


Conventions for Command-Line Arguments:
---------------------------------------

The command-line arguments are segregated into three broad categories: a)
action commands, b) aggregate action commands, and c) others.

a) The action commands are those such as --build, --test, etc. and are shown
with [ACTION] in their documentation.  These action commands have no off
complement.  If the action command appears, then the action will be performed.

b) Aggregate action commands such as --do-all and --local-do-all turn on sets
of other action commands and are shown with [AGGR ACTION] in their
documentation.  The sub-actions that these aggregate action commands turn on
and cannot be disabled with other arguments.

c) Other arguments are those that are not [ACTION] or [AGGR ACTION] arguments
and tend to either pass in data and turn control flags on or off.


Exit Code:
---------

This script returns 0 if the actions requested are successful.  This does not
necessarily imply that it is okay to do a push.  For example, if only --pull
is passed in and is successful, then 0 will be returned but that does *not*
mean that it is okay to do a push.  A 0 return value is a necessary but not
sufficient condition for readiness to push.

"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--show-defaults", dest="showDefaults", action="store_true",
  help="Show the default option values and do nothing at all.",
  default=False )

clp.add_option(
  "--eg-git-version-check", dest="enableEgGitVersionCheck", action="store_true",
  help="Enable automatic check for the right versions of eg and git. [default]" )
clp.add_option(
  "--no-eg-git-version-check", dest="enableEgGitVersionCheck", action="store_false",
  help="Do not check the versions of eg and git, just trust they are okay.",
  default=True )

clp.add_option(
  "--trilinos-src-dir", dest="trilinosSrcDir", type="string",
  default='/'.join(getScriptBaseDir().split("/")[0:-2]),
  help="The Trilinos source base directory for code to be tested." )

clp.add_option(
  "--extra-repos", dest="extraRepos", type="string", default="",
  help="List of comma separated extra Trilinos repos " \
  +"(e.g. 'preCopyrightTrilinos,LIMEExt') containing extra " \
  +"trilinos packages that can be enabled.  The order these repos is "
  +"listed in not important.")

clp.add_option(
  "--skip-deps-update", dest="skipDepsUpdate", action="store_true",
  help="If set, skip the update of the dependency XML file (debug only).",
  default=False )

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
  "--enable-all-packages", "enableAllPackages", ('auto', 'on', 'off'), 0,
  "Determine if all Trilinos packages are enabled 'on', or 'off', or let" \
  +" other logic decide 'auto'.  Setting to 'off' is appropriate when" \
  +" the logic in this script determines that a global build file has changed" \
  +" but you know that you don't need to rebuild every Trilinos package for" \
  +" a reasonable test.  Setting --enable-packages effectively disables this" \
  +" option.  NOTE: Setting this to 'off' does *not* stop the forward enabling" \
  +" of downstream packages for packages that are modified or set by --enable-packages." \
  +"  This only stops the setting of -DTrilinos_ENABLE_ALL_PAKCAGES:BOOL=ON.",
  clp )

clp.add_option(
  "--enable-fwd-packages", dest="enableFwdPackages", action="store_true",
  help="Enable forward Trilinos packages. [default]" )
clp.add_option(
  "--no-enable-fwd-packages", dest="enableFwdPackages", action="store_false",
  help="Do not enable forward Trilinos packages.", default=True )

clp.add_option(
  "--continue-if-no-updates", dest="abortGracefullyIfNoUpdates", action="store_false",
  help="If set, then the script will continue if no updates are pulled from any repo. [default]",
  default=False )
clp.add_option(
  "--abort-gracefully-if-no-updates", dest="abortGracefullyIfNoUpdates", action="store_true",
  help="If set, then the script will abort gracefully if no updates are pulled from any repo.",
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
  "--extra-cmake-options", dest="extraCmakeOptions", type="string", default="",
  help="Extra options to pass to 'cmake' after all other options." \
  +" This should be used only as a last resort.  To disable packages, instead use" \
  +" --disable-packages." )

clp.add_option(
  "-j", dest="overallNumProcs", type="string", default="",
  help="The options to pass to make and ctest (e.g. -j4)." )

clp.add_option(
  "--make-options", dest="makeOptions", type="string", default="",
  help="The options to pass to make (e.g. -j4)." )

clp.add_option(
  "--ctest-options", dest="ctestOptions", type="string", default="",
  help="Extra options to pass to 'ctest' (e.g. -j2)." )

clp.add_option(
  "--ctest-timeout", dest="ctestTimeOut", type="float", default=300,
  help="timeout (in seconds) for each single 'ctest' test (e.g. 180" \
  +" for three minutes)." )

clp.add_option(
  "--show-all-tests", dest="showAllTests", action="store_true",
  help="Show all of the tests in the summary email and in the commit message" \
  +" summary (see --append-test-results)." )
clp.add_option(
  "--no-show-all-tests", dest="showAllTests", action="store_false",
  help="Don't show all of the test results in the summary email. [default]",
  default=False )

clp.add_option(
  "--without-mpi-debug", dest="withMpiDebug", action="store_false",
  help="Skip the mpi debug build.", default=True )

clp.add_option(
  "--without-serial-release", dest="withSerialRelease", action="store_false",
  help="Skip the serial release build.", default=True )

clp.add_option(
  "--without-default-builds", dest="withoutDefaultBuilds", action="store_true",
  default=False,
    help="Skip the default builds (same as --without-mpi-debug --without-serial-release)." \
    +"  You would use option along with --extra-builds=BUILD1,BUILD2,... to run your own" \
    +" local custom builds." )

clp.add_option(
  "--ss-extra-builds", dest="ssExtraBuilds", type="string", default="",
  help="List of comma-separated SS extra build names.  For each of the build names in" \
  +" --ss-extra-builds=<BUILD1>,<BUILD2>,..., there must be a file <BUILDN>.config in" \
  +" the local directory along side the COMMON.config file that defines the special" \
  +" build options for the extra build." )

clp.add_option(
  "--extra-builds", dest="extraBuilds", type="string", default="",
  help="List of comma-separated extra build names.  For each of the build names in" \
  +" --extra-builds=<BUILD1>,<BUILD2>,..., there must be a file <BUILDN>.config in" \
  +" the local directory along side the COMMON.config file that defines the special" \
  +" build options for the extra build." )

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
  help="If set then if a build/test case is skipped for some reason (i.e." \
  +" because no packages are enabled) then no email will go out for that case." \
  +" [default]",
  default=True )

clp.add_option(
  "--send-email-for-all", dest="sendEmailOnlyOnFailure", action="store_false",
  help="If set, then emails will get sent out for all operations. [default]" )
clp.add_option(
  "--send-email-only-on-failure", dest="sendEmailOnlyOnFailure", action="store_true",
  help="If set, then emails will only get sent out for failures.",
  default=False )

clp.add_option(
  "--send-email-to-on-push", dest="sendEmailToOnPush", type="string",
  default="trilinos-checkin-tests@software.sandia.gov",
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
  help="Rebase the local commits on top of origin/master before amending" \
  +" the last commit and pushing.  Rebasing keeps a nice linear commit" \
  +" history like with CVS or SVN and will work perfectly for the basic" \
  +" workflow of adding commits to the 'master' branch and then syncing" \
  +" up with origin/master before the final push. [default]" )
clp.add_option(
  "--no-rebase", dest="rebase", action="store_false",
  help="Do not rebase the local commits on top of origin/master before" \
  +" amending the final commit and pushing.  This allows for some more " \
  +" complex workflows involving local branches with multiple merges." \
  +"  However, this will result in non-linear history and will allow for" \
  +" trivial merge commits with origin/master to get pushed.  This mode" \
  +" should only be used in cases where the rebase mode will not work or " \
  +" when it is desired to use a merge commit to integrate changes on a" \
  +" branch that you wish be able to easily back out.  For sophisticated" \
  +" users of git, this may in fact be the prefered mode.",
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
  +" must select this option to avoid this last amending commit.",
  default=True )

clp.add_option(
  "--extra-pull-from", dest="extraPullFrom", type="string", default="",
  help="Optional extra git pull '<repository>:<branch>' to merge in changes from after" \
  +" pulling in changes from 'origin'.  This option uses a colon with no spaces in between" \
  +" <repository>:<branch>' to avoid issues with passing arguments with spaces." \
  +"  For example --extra-pull-from=machine:/base/dir/repo:master." \
  +"  This extra pull is only done if --pull is also specified.  NOTE: when using" \
  +" --extra-repo=REPO1,REPO2,... the <repository> must be a named repository that is" \
  +" present in all of the git repos or it will be an error." )

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
  +" will not be removed if --without-serial-release is specified)." )

clp.add_option(
  "--pull", dest="doPull", action="store_true", default=False,
  help="[ACTION] Do the pull from the default (origin) repository and optionally also" \
    +" merge in changes from the repo pointed to by --extra-pull-from." )

clp.add_option(
  "--configure", dest="doConfigure", action="store_true", default=False,
  help="[ACTION] Do the configure of Trilinos." )

clp.add_option(
  "--build", dest="doBuild", action="store_true", default=False,
  help="[ACTION] Do the build of Trilinos." )

clp.add_option(
  "--test", dest="doTest", action="store_true", default=False,
  help="[ACTION] Do the running of the enabled Trilinos tests." )

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
  help="[ACTION] Push the committed changes in the local repo into to global repo" \
    +" 'origin' for the current branch.  Note: If you have uncommitted changes this" \
    +" command will fail.  Note: You must have SSH public/private keys set up with" \
    +" the origin machine (e.g. software.sandia.gov) for the push to happen without" \
    +" having to type your password." )

clp.add_option(
  "--execute-on-ready-to-push", dest="executeOnReadyToPush", type="string", default="",
  help="[ACTION] A command to execute on successful execution and 'READY TO PUSH'" \
  +" status from this script.  This can be used to do a remote SSH invocation to a" \
  +" remote machine to do a remote pull/test/push after this machine finishes." )

(options, args) = clp.parse_args()

# NOTE: Above, in the pairs of boolean options, the *last* add_option(...) 
# takes effect!  That is why the commands are ordered the way they are!


#
# Echo the command-line
#

print ""
print "**************************************************************************"
print "Script: checkin-test.py \\"

if options.enableEgGitVersionCheck:
  print "  --eg-git-version-check \\"
else:
  print "  --no-eg-git-version-check \\"
print "  --trilinos-src-dir='"+options.trilinosSrcDir+"' \\"
print "  --extra-repos='"+options.extraRepos+"' \\"
if options.skipDepsUpdate:
  print "  --skip-deps-update \\"
print "  --enable-packages='"+options.enablePackages+"' \\"
print "  --disable-packages='"+options.disablePackages+"' \\"
print "  --enable-all-packages='"+options.enableAllPackages+"'\\"
if options.enableFwdPackages:
  print "  --enable-fwd-packages \\"
else:
  print "  --no-enable-fwd-packages \\"
if options.abortGracefullyIfNoUpdates:
  print "  --abort-gracefully-if-no-updates \\"
else:
  print "  --continue-if-no-updates \\"
if options.abortGracefullyIfNoEnables:
  print "  --abort-gracefully-if-no-enables \\"
else:
  print "  --continue-if-no-enables \\"
print "  --extra-cmake-options='"+options.extraCmakeOptions+"' \\"
if options.overallNumProcs:
  print "  -j"+options.overallNumProcs+" \\"
print "  --make-options='"+options.makeOptions+"' \\"
print "  --ctest-options='"+options.ctestOptions+"' \\"
print "  --ctest-timeout="+str(options.ctestTimeOut)+" \\"
if options.showAllTests:
  print "  --show-all-tests \\"
else:
  print "  --no-show-all-tests \\"
if not options.withMpiDebug:
  print "  --without-mpi-debug \\" 
if not options.withSerialRelease:
  print "  --without-serial-release \\" 
if options.withoutDefaultBuilds:
  print "  --without-default-builds \\" 
print "  --ss-extra-builds='"+options.ssExtraBuilds+"' \\"
print "  --extra-builds='"+options.extraBuilds+"' \\"
print "  --send-email-to='"+options.sendEmailTo+"' \\"
if options.skipCaseSendEmail:
  print "  --skip-case-send-email \\"
else:
  print "  --skip-case-no-email \\"
if not options.sendEmailOnlyOnFailure:
  print "  --send-email-for-all \\"
else:
  print "  --send-email-only-on-failure \\ "
print "  --send-email-to-on-push='"+options.sendEmailToOnPush+"' \\"
if options.forcePush:
  print "  --force-push \\"
else:
  print "  --no-force-push \\"
if options.doPushReadinessCheck:
  print "  --do-push-readiness-check \\"
else:
  print "  --skip-push-readiness-check \\"
if options.rebase:
  print "  --rebase \\"
else:
  print "  --no-rebase \\"
if options.appendTestResults:
  print "  --append-test-results \\"
else:
  print "  --no-append-test-results \\"
if options.extraPullFrom:
  print "  --extra-pull-from='"+options.extraPullFrom+"' \\"
if options.allowNoPull:
  print "  --allow-no-pull \\"
if options.wipeClean:
  print "  --wipe-clean \\"
if options.doPull:
  print "  --pull \\"
if options.doConfigure:
  print "  --configure \\"
if options.doBuild:
  print "  --build \\"
if options.doTest:
  print "  --test \\"
if options.localDoAll:
  print "  --local-do-all \\"
if options.doAll:
  print "  --do-all \\"
if options.doPush:
  print "  --push \\"
if options.executeOnReadyToPush:
  print "  --execute-on-ready-to-push=("+options.executeOnReadyToPush+") \\"


#
# Check the input arguments
#

if options.doAll and options.localDoAll:
  print "\nError, you can not use --do-all and --local-do-all together!  Use on or the other!"
  sys.exit(1)

if options.doAll and options.allowNoPull:
  print "\nError, you can not use --do-all and --allow-no-pull together! (see the" \
    " documentation for the --do-all, --local-do-all, and --allow-no-pull arguments.)"
  sys.exit(2)

if options.extraPullFrom:
  getRepoSpaceBranchFromOptionStr(options.extraPullFrom) # Will validate form


#
# Execute the checkin test guts
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
    print "\nREQUESTED ACTIONS: PASSED\n"
    sys.exit(0)
  else:
    print "\nREQUESTED ACTIONS: FAILED\n"
    sys.exit(1)
