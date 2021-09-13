# CDash and Jenkins drivers for ATDM Trilinos builds

The files and scripts in this directory drive builds of Trilinos on the
various ATDM test platforms using Jenkins jobs and submit results to the
Trilinos CDash site.

**Outline:**
* <a href="#base-ctest-cdash-config">Base CTest/CDash configuration</a>
* <a href="#system-specific-drivers">System-specific driver files</a>
* <a href="#split-ctest-s-drivers">Split ctest -S drivers and submits to CDash</a>
* <a href="#run-locally-and-debug">Running locally and debugging</a>
* <a href="#installing">Installing as a byproduct of of running ctest -S drivers</a>
* <a href="#setup-jenkins-jobs">Setting up Jenkins jobs</a>
* <a href="#specific-system-directories">Specific <system_name> directories</a>
* <a href="#howto-add-new-system">How add a new system</a>


<a name="base-ctest-cdash-config"/>

## Base CTest/CDash configuration

The base directory:

```
  Trilinos/cmake/ctest/drivers/atdm/
```

contains files that are used for driving ATDM builds on various machines.
These files define common behavior and reduces duplication to ease
maintenance.

This directory contains the file:

```
  Trilinos/cmake/ctest/drivers/atdm/TrilinosCTestDriverCore.atdm.cmake
```

which sets up the CTest -S script driver options that are used for all
automated ATDM builds of Trilinos.  It uses the configuration of Trilinos
given in the files defined in the directory:

```
  Trilinos/cmake/std/atdm/
```

That file reads `JOB_NAME` from the environment and uses it to set the CDash build
name.  (Therefore, the Jenkins `JOB_NAME` is the same as the CDash build name
for all of these Trilinos ATDM builds.)  It also other CMake and CTest options
that are pulled out of the environment set by the
`cmake/std/atdm/<system_name>/environment.sh` script.  (See `$ENV{<varName>}`
to see what variables are pulled out of the environment.)

This directory contains a CTest -S driver script:

```
  Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake
```

which is called using:

```
  ctest -S <base-dir>/Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake
```

which can be run using any way desired and it will clone a new Trilinos git
repo (if not cloned already).  (But this file gets directly run by the
universal driver script `ctest-s-driver.sh` described above.)

This directory contains the file:

```
  Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
```

which sets up and runs `ctest -S .../atdm/ctest-driver.cmake`.  (Note, there
are also `ctest-s-driver-config-build.sh` and `ctest-s-driver-test.sh` drivers
as explained <a href="#split-ctest-s-drivers">below</a>.)

This base directly also contains the script:

```
  Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh
```

which just runs:

```
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/$JOB_NAME.sh
```

for that current system for the given job name.


<a name="system-specific-drivers"/>

## System-specific driver files

Each system `<system_name>` is given a sub-directory under this directory:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/
```

Each `atdm/<system_name>/` directory minimally contains a local driver script:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/local-driver.sh
```

which launches the `ctest-s-driver.sh` script in the appropriate way for that
system (such as using a batch running system like Slum).  This script assumes
a directory structure as set up by Jenkins but does not really require Jenkins
to run it.

The sub-directory:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
```

contains specific drivers with the file names of the Jenkins build names:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/$JOB_NAME.sh
```

This file sets some tweaks for that particular build `$JOB_NAME` link such as
which CDash Track/Group results are sent to and other tweaks like this.
Having this file and using `smart-jenkins-driver.sh` allows customizing almost
anything about a particular ATDM build of Trilinos without having to touch the
Jenkins job configuration (which is not under any type of version control).


<a name="split-ctest-s-drivers"/>

## Split ctest -S drivers and submits to CDash

On some machines (e.g. the SNL HPC machines), the update, configure and build
must be done on a login/compile node and running the tests must be done on a
compute node.  To accommodate this, the update, configure, build and submit of
those results to CDash can be done in one `ctest -S` driver invocation and
then running the tests and submitting those results to CDash can be done in a
later `ctest -S` driver invocation.

To use this approach, the `<system-name>/local-driver.sh` script should be set
up to run as follows:

```
source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-config-build.sh
<comamnd-to-run-on-compute-node> \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh
```

The first script `ctest-s-driver-config-build.sh` must be sourced instead of
just run because it sets up the environment and changes to a subdir.  The environment must be
loaded on the login/compile node and not the compute node and therefore, it
must be done this way.

With newer versions of CMake (3.10+), this will result in all of the results
going to the same build (i.e. row) on CDash.  For older versions of CMake, the
test results show up as a separate build (same `site` and `build name` but
different `build stamp`).


<a name="run-locally-and-debug"/>

## Running locally and debugging

To test out locally, first set up a local directory and symlink as:

```
$ cd <some_base_build_dir>/
$ ln -s <some_base_dir>/Trilinos/cmake/std/atdm/ctest-s-local-test-driver.sh .
````

Once that directory and the symlinked script `ctest-s-local-test-driver.sh`
are set up, then one can drive and test out builds (for a subset of packages)
as:

```
$ env \
    Trilinos_PACKAGES=Kokkos,Teuchos,Tpetra \
    CTEST_DO_SUBMIT=OFF \
  ./ctest-s-local-test-driver.sh <build-base-name-0> <build-base-name-1> ...
```

where the build names `<build-base-name-i>` (e.g. `gnu-opt-debug`) must match
the entries given in variable `ATDM_CONFIG_ALL_SUPPORTED_BUILDS` in the file
`cmake/std/atdm/<system_name>/all_supported_builds.sh` for the local system.
This will not submit to CDash due to `CTEST_DO_SUBMIT=OFF` so to see the
status of the build and tests, see the generated files:

```
  <some_base_build_dir>/<full_build_name>/smart-jenkins-driver.out
```

(e.g. `<full_build_name>` = `Trilinos-atdm-<system_name>-gnu-opt-debug`) and
also examine the generated `*.xml` configure, build, and test files created
under:

```
  <some_base_build_dir>/<full_build_name>/SRC_AND_BUILD/BUILD/Testing/
```

to see if it is doing the right thing.

To test the submit to CDash (Experimental Track/Group) without a complete
rebuild, run again with:

```
$ env \
    Trilinos_PACKAGES=Kokkos,Teuchos,Tpetra \
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE \
    CTEST_DO_SUBMIT=ON \
  ./ctest-s-local-test-driver.sh <build-base-name-0> <build-base-name-1> ...
```

To test that all of the builds specified in the file
`cmake/std/atdm/<system_name>/all_supported_builds.sh` have matching driver
scripts in the directory `cmake/ctest/drivers/atdm/<system_name>/drivers/`,
run with `all` like:

```
$ env \
    Trilinos_PACKAGES=Kokkos \
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE \
    CTEST_DO_SUBMIT=OFF \
  ./ctest-s-local-test-driver.sh all
```

and examine the generated `<full_build_name>/smart-jenkins-driver.out` files
to ensure that they were all found correctly.

If things look good after all of that testing, then the builds are ready to be
set up as Jenkins (or GitLab CI or cron, etc.) jobs.

NOTE: When one is running on a loaded/shared machine and therefore needs to
use less processes to build and test, one can use the environment variables
`ATDM_CONFIG_BUILD_COUNT_OVERRIDE` and
`ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERIDE` and use them as, for example:

```
$ env \
    ATDM_CONFIG_BUILD_COUNT_OVERRIDE=8 \
    ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERIDE=12 \
    Trilinos_PACKAGES=Kokkos,Teuchos,Tpetra \
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE \
    CTEST_DO_SUBMIT=OFF \
  ./ctest-s-local-test-driver.sh <build-base-name>
```

That can also be handy for specializing the automated builds on specific SEMS
and CEE RHEL7 machines, for example, that may have more or less hardware
cores.


<a name="installing"/>

## Installing as a byproduct of running ctest -S drivers

These scripts support installing Trilinos as a byproduct of running the `ctest
-S` driver scripts.  These installations are often done as the `jenkins`
entity account from the `jenkins-srn.sandia.gov` site or as the
'atdm-devops-admin' entity account (e.g. from a cron job).  In order to
properly set up installations of Trilinos on all of the supported systems such
that the `atdm-devops-admin` entity account can edit and remote the installs,
some features of TriBITS are used to run `chgrp` and `chmod` on the installed
directories and files.  In addition, automatic `<date>/` subdirectories are
created which for each testing day.

The following (bash) environment variables determine the behavior of the ATDM
`ctest -S` scripts for building and installing Trilinos using this scheme:

* `ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE=<install-prefix-base>`:
  Defines the base directory installs of Trilinos under
  `<install-prefix-base>/<date>/<system-build-name>`.  This directory
  `<install-prefix-base>` must be owned by the user 'atdm-devops-admin', the
  group 'wg-run-as-atdm-devops' and it should be world readable and group
  read/writable.  (If `ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE==""`
  and `ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT!=""` and
  `ATDM_CONFIG_USE_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=="1"`, then
  `ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE` is set to
  `${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT}`.)  (The var
  `ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE` is given a default value
  in the file `cmake/std/atdm/atdm_devops_install_defaults.sh` if it is not
  already set in `cmake/std/atdm/<system_name>/environment.sh`.)

* `ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR=<base-dir>`:
  Defines the base directory for setting the group and permissions.  If not
  set in the env already and if
  `${ATDM_CONFIG_USE_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR_DEFAULT}==1`,
  then this will bet set to `<install-prefix-base>/<date>`.

* `ATDM_CONFIG_MAKE_INSTALL_GROUP`: Defines the group that will get set on all
  files and dirs under
  `${ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}`.  If not
  already set in the env, then this will be set to
  `${ATDM_CONFIG_MAKE_INSTALL_GROUP_DEFAULT}` if
  `${ATDM_CONFIG_USE_MAKE_INSTALL_GROUP_DEFAULT}==1`.  (The var
  `ATDM_CONFIG_MAKE_INSTALL_GROUP_DEFAULT` is given a default value in the
  file `cmake/std/atdm/atdm_devops_install_defaults.sh` if it is not already
  set in `cmake/std/atdm/<system_name>/environment.sh`.)

* `ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS=[0|1]`: Set to '1' to use the
  defaults for the variables that have defaults (i.e. this sets the
  environment variables
  `ATDM_CONFIG_USE_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=1`,
  `ATDM_CONFIG_USE_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR_DEFAULT=1`,
  and `ATDM_CONFIG_USE_MAKE_INSTALL_GROUP_DEFAULT=1`.)

The defaults for some of these can be set for all systems in the file
`cmake/std/atdm/atdm_devops_install_defaults.sh`.  These defaults can then be
overridden in the `cmake/std/atdm/<system_name>/environment.sh` for each
system.

Then the cron or jenkins driver jobs can activate the usage of these defaults
and perform standard installs as a bi-product of the testing process as
follows:

```
export ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS=1
export CTEST_DO_INSTALL=ON
${WORKSPACE}/Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh
```

That will result in the install of Trilinos under:

```
<install-prefix-base>/<date>/<system-build-name>/
```

where all of the files and directories `<install-prefix-base>/<date>` on down
will be owned by the group `${ATDM_CONFIG_MAKE_INSTALL_GROUP}` and will be
given group read/write and "other" read access.

NOTE:

* The `<date>` in the format `YYYY-MM-DD` is automatically determined to
  correspond to the CDash `date=<date>` field for the given build of Trilinos
  (assuming that `ctest_start()` is called almost immediately which it should
  be within a second or less).

* The build name `<system-build-name>` is taken from the full build name
  stored in the environment variable `${JOB_NAME}` (with `Trilinos-atdm-`
  removed from the beginning of the Jenkins job name).

Internally, for each build, the environment variable
`ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX` is set to this full install path
(which then gets picked up in the `ATDMDevEnvSettings.cmake` file during the
CMake configure step).

**WARNING:** Do **NOT** directly set the environment variable
`ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX`.  That would result in every Trilinos
build getting installed on top of each other in the same installation
directory!


## Installing as the 'atdm-devops-admin' account using the 'jenkins' entity account

In order to protect the installation Trilinos from other 'jenkins' jobs, a
strategy has been implemented that allows performs the final install using the
`atdm-devops-admin` account using a setuid program called
`run-as-atdm-devops-admin` that in installed on each supported system.  The
setup of that program under the `atdm-devops-admin` user account is described
in:

* https://gitlab.sandia.gov/atdm-devops-admin/run-as-atdm-devops-admin/blob/master/README.md

This documentation below assumes that the program 'run-as-atdm-devops-admin'
is correctly installed on each given system.

The following additional (bash) environment variables determine the behavior
of the ATDM `ctest -S` scripts for building and installing Trilinos using this
scheme:

* `ATDM_CONFIG_WORKSPACE_BASE=<workspace-base>`: Defines a different base
  workspace directory under which the subdir `SRC_AND_BUILD` is created and
  used (and the scripts 'cd' into that workspace).  This directory
  `<workspace-base>` must be owned and be writable by the
  `wg-run-as-atdm-devops` group and must be given the sticky group bit
  `chmod g+s <workspace-base>` so that the 'jenkins' account can create files
  and directories under this directory.  If not set, then `WORKSPACE` (set by
  the Jenkins job) is used as the base working directory. (If
  `ATDM_CONFIG_WORKSPACE_BASE==""` and
  `ATDM_CONFIG_WORKSPACE_BASE_DEFAULT!=""` and
  `ATDM_CONFIG_USE_WORKSPACE_BASE_DEFAULT=="1"`, then
  `ATDM_CONFIG_WORKSPACE_BASE` is set to
  `${ATDM_CONFIG_WORKSPACE_BASE_DEFAULT}`.)

* `ATDM_CONFIG_INSTALL_PBP_RUNNER=<base-dir>/run-as-atdm-devops-admin`:
  Defines an executable that is used to run the install command in the target
  `install_package_by_package`.  This allows inserting the
  `run-as-atdm-devops-admin` setuid program to run the install command as the
  'atdm-devops-admin' user.  (If `ATDM_CONFIG_INSTALL_PBP_RUNNER==""` and
  `ATDM_CONFIG_INSTALL_PBP_RUNNER_DEFAULT!=""` and
  `ATDM_CONFIG_USE_INSTALL_PBP_RUNNER_DEFAULT=="1"`, then
  `ATDM_CONFIG_INSTALL_PBP_RUNNER` is set to
  `${ATDM_CONFIG_INSTALL_PBP_RUNNER_DEFAULT}`)

The variables `ATDM_CONFIG_WORKSPACE_BASE_DEFAULT` and
`ATDM_CONFIG_INSTALL_PBP_RUNNER_DEFAULT` are meant to be set in the
`atdm/<system_name>/environment.sh` file as, for example:

```
export ATDM_CONFIG_WORKSPACE_BASE_DEFAULT=/home/atdm-devops-admin/jenkins
export ATDM_CONFIG_INSTALL_PBP_RUNNER_DEFAULT=/home/atdm-devops-admin/tools/run-as-atdm-devops-admin
```

Running with:

```
export ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS=1
export CTEST_DO_INSTALL=ON
${WORKSPACE}/Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh
```

will result in the alternate workspace directory being create as:

```
export WORKSPACE=${ATDM_CONFIG_WORKSPACE_BASE}/${ATDM_CONFIG_SYSTEM_NAME}/${JOB_NAME}
```

The inner clone of Trilinos and the build of Trilinos will be performed under
that subdir.

That will result in the install of Trilinos as the 'atdm-devops-admin' user
under:

```
<install-prefix-base>/<date>/<system-build-name>/
```


<a name="setup-jenkins-jobs"/>

## Setting up Jenkins jobs

To set up a Jenkins build, you must set the following in the Jenkins build
configuration GUI:

* "Source Code Management"
  * "Git"
    * "Repository URL": `https://github.com/trilinos/Trilinos.git
    * "Branch": `develop`
* "Build Triggers"
  * "Build Periodically" (**checked**)
    * "Schedule": `H 0 ***`
* "Build Environment"
  * "Abort of the build if it's stuck" (**checked**)
    * "Time-out strategy": `Deadline`
      * "Deadline time": `20:45`
      * "Deadline tolerance in minutes": `1`
* "Build"
  * "Execute shell"
    * "Command": `Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh`

But be careful not to put any additional settings that what is absolutely
necessary because Jenkins configurations are not under version control and
there is no tractability for changes in these settings!


<a name="specific-system-directories"/>

## Specific <system_name> directories

The following `<system_name>` sub-directories exist (in alphabetical order):

* `cee-rhel7/`: Contains files to drive builds on CEE LAN RHEL7 machines with
  the 'sparc-dev' modules.

* `tlcc2/`: Contains files to drive builds on the SRN HPC TLCC-2 machines
  (e.g. 'chama', 'skybridge', etc.).

* `waterman/`: Contains files to drive builds on the SRN Test Bed machine
  `waterman`.


<a name="howto-add-new-system"/>

## How to add a new system

To add a new system, first add a new failing unit test for that system in the
file:


```
  Trilinos/cmake/std/atdm/test/unit_tests/get_system_info_unit_tests.sh
```

That unit test should fail!

Second, update the file:

```
  Trilinos/cmake/std/atdm/utils/get_known_system_info.sh
```

by adding the new system name to the list variable:

```
ATDM_KNOWN_SYSTEM_NAMES_LIST=(
  ...
  )
```

Then, if the system selection is done by matching to the `hostname`, add a new
`elif` statement for that set of machines.  Note that more than one `hostname`
machine may map to the same `<new_system_name>`.

However, if adding a new system type that will run on many machines and not
looking at the `hostname` on the machine, then add a new `if` block to the
section for the logic.  For an example, see how the system types `tlcc2`,
`sems-rhel7`, and `cee-rhel7` are handled.

The variable `ATDM_HOSTNAME` (set to exported variable
`ATDM_CONFIG_CDASH_HOSTNAME`) is used for the CDash site name.  This makes it
so that any node `white05`, `white12`, etc. just says `white` on CDash.  This
is important for the CDash 'next' and 'previous' relationships to work.  (But
for `CTEST_TEST_TYPE=Experimental` builds, the real `hostname` is used which
is stored in the exported environment variable `ATDM_CONFIG_REAL_HOSTNAME`.  This
ensures that queries with `cdash/queryTests.php` don't accidentally pick up
tests from "Experimental" builds.)

The variable `ATDM_SYSTEM_NAME` (set to the exported variable
`ATDM_CONFIG_SYSTEM_NAME`) must be set to `<new_system_name>` which is
selected for this new system type.

Make sure the new unit test(s) in the file `get_system_info_unit_tests.sh`
pass and add more unit tests to cover full behavior for the new system if
needed.

Then, create a new directory for the new system called `<new_system_name>`:

```
  Trilinos/cmake/std/atdm/<new_system_name>/
```

and fill in the file:

```
  Trilinos/cmake/std/atdm/<new_system_name>/environment.sh
```

The file `<new_system_name>/environment.sh` contains all of the
system-specific settings to get an ATDM Trilinos build to work on this new
platform.  For tips, see examples of these files in other
`cmake/std/atdm/<system_name>/` directories.  And see how these environment variables are
read and interpreted in the file:

```
  Trilinos/cmake/std/atdm/ATDMDevEnvSettings.cmake
```

to understand their impact on the configuration, build and testing.

In addition, a custom set of build configuration options can be set in the
file:

```
  Trilinos/cmake/std/atdm/<new_system_name>/custom_builds.sh
```

This file can be used to put in special logic for special compilers and
compiler versions and other types of logic.  This file gets sourced before the
standard logic is executed in the file `atdm/utils/set_build_options.sh`.  (To
see an example of the usage of this file, see
`atdm/cee-rhel7/custom_builds.sh` and `atdm/cee-rhel7/environment.sh`.)

A few of the environment variables that need to be set in the file
`<new_system_name>/environment.sh` worth specifically discussing are:

* `ATDM_CONFIG_USE_NINJA`: If Ninja is available on the system and will be
  loaded into the environment in this script, set this to `TRUE`.  Otherwise, if
  makefiles are to be used, set to `FALSE`.  When possible, it is best to use
  Ninja but this needs to be a Fortran-enabled Ninja (see [Installing Ninja
  from Source](
  https://tribits.org/doc/TribitsBuildReference.html#installing-ninja-from-source)).

* `ATDM_CONFIG_BUILD_COUNT`: Set to a positive integer to control how many
  processes will be used when building.  When using
  `ATDM_CONFIG_USE_NINJA=TRUE`, this can be set to `0` or `-1` to allow
  `ninja` to be run without a `-j<N>` option and therefore use all of the
  cores (see [Building in parallel with Ninja](
  https://tribits.org/doc/TribitsBuildReference.html#building-in-parallel-with-ninja)).
  But when using Makefiles, this must be set to a positive integer.  Don't set
  too high so that it will overload the machine.

* `ATDM_CONFIG_CTEST_PARALLEL_LEVEL`: This determines the parallel level when
  running `ctest -j<N>` (see [CTEST_PARALLEL_LEVEL](
  https://tribits.org/doc/TribitsDevelopersGuide.html#determining-what-testing-related-actions-are-performed-tribits-ctest-driver)).
  When using threading with OpenMP, note that this will need to be reduced
  based on the value set for `OMP_NUM_THREADS`.

* `ATDM_CONFIG_KOKKOS_ARCH`: This must be set correctly for each compiler
  supported.  Knowing the correct value to set is beyond the scope of this
  document.

After creating the file:

```
  Trilinos/cmake/std/atdm/<new_system_name>/environment.sh
```

one can do some local configures and builds of Trilinos packages as documented
in the file

```
  Trilinos/cmake/std/atdm/README.md
```

Then, add the file:

```
  Trilinos/cmake/std/atdm/<new_system_name>/all_supported_builds.sh
```

and fill in the list of builds that are going to be supported on this machine.  For example, this looks like:

```
  export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-<system_name>-

  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
    gnu-debug-openmp
    gnu-opt-openmp
    ...
    )
```

The names
`${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[i]}.sh`
must match the ctest -S driver files under:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
```

The variable `ATDM_CONFIG_ALL_SUPPORTED_BUILDS` is used by
`./checkin-test-sems.sh all ...` to test all configurations on a given machine
using the `checkin-test.py` script.

NOTE: The array variable `ATDM_CONFIG_ALL_SUPPORTED_BUILDS` can also be
specified as:

```
  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS="gnu-debug-openmp gnu-opt-openmp ..."
```

which works as well.

Once local configurations are complete and some local testing if finished,
then create the ctest -S / Jenkins driver directory:

```
  Trilinos/cmake/ctest/drivers/atdm/<new_system_name>/
```

and create the basic driver file:

```
  Trilinos/cmake/ctest/drivers/atdm/<new_system_name>/local-driver.sh
```

and one or more smart Jenkins driver files:

```
  Trilinos/cmake/ctest/drivers/atdm/<new_system_name>/drivers/$JOB_NAME.sh
```

Use examples from other `Trilinos/cmake/ctest/drivers/atdm/<system_name>/`
directories for inspiration.  Then test the configurations one at a time as
described <a href="#run-locally-and-debug">above</a>.

Once the basic configurations seem to be working, then commit your changes on
a topic branch and [submit a pull request (PR) to
Trilinos](https://github.com/trilinos/Trilinos/wiki/Submitting-a-Trilinos-Pull-Request).
Make sure and mention `@fryeguy52` and add the label `ATDM`.

Once the PR has been has been merged, then set up the Jenkins jobs to run the
builds as described <a href="#setup-jenkins-jobs">above</a>. Please note
that the variable `JOB_NAME` is set by Jenkins and is the name of currently
running job.  Your driver files in
`Trilinos/cmake/ctest/drivers/atdm/<new_system_name>/drivers/` must be named
to exactly match the Jenkins `JOB_NAME` variable.

ToDo: Fill in more detail, add an FAQ, etc.
