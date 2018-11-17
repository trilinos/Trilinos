# CDash and Jenkins drivers for ATDM Trilinos builds

The files and scripts in this directory drive builds of Trilinos on the
various ATDM test platforms using Jenkins jobs and submit results to the
Trilinos CDash site.

**Outline:**
* <a href="#base-ctestcdash-configuration">Base CTest/CDash configuration</a>
* <a href="#system-specific-driver-files">System-specific driver files</a>
* <a href="#split-ctest--s-drivers-and-submits-to-cdash">Split ctest -S drivers and submits to CDash</a>
* <a href="#running-locally-and-debugging">Running locally and debugging</a>
* <a href="#setting-up-jenkins-jobs">Setting up Jenkins jobs</a>
* <a href="#specific-system_name-directories">Specific <system_name> directories</a>
* <a href="#how-add-a-new-system">How add a new system</a>


## Base CTest/CDash configuration

The base directory:

```
  Trilinos/cmake/ctest/drivers/atdm/
```

contains files that are used for driving builds on all machines.  These files
define common behavior and reduces duplication to ease maintenance.

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
to see what variables are pulled out of the env.)

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
as explained <a href="#split-ctest--s-drivers-and-submits-to-cdash">below</a>.)

This base directly also contains the script:

```
  Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh
```

which just runs:

```
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/$JOB_NAME.sh
```

for that current system for the given job name.


## System-specific driver files

Each system `<system_name>` is given a subdirectory under this directory:

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
use less processes to build and test, one can use the env vars
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
and CEE RHEL6 machines, for example, that may have more or less hardware
cores.


## Setting up Jenkins jobs

To set up a Jenkins build, you must set the following in the Jenkins build
configuration GUI:

* "Source Code Management"
  * "Git"
    * "Repository URL": `https://github.com/trilinos/Trilinos.git` (until we
         can figure out how to clone from `sofware.sandia.gov:/git/nightly/Trilinos` with Jenkins)
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


## Specific <system_name> directories

The following `<system_name>` sub-directories exist (in alphabetical order):
 
* `cee-rhel6/`: Contains files to drive builds on CEE LAnL RHEL6 machines with
  a SEMS environment.

* `chama/`: Contains files to drive builds on the SRN HPC machine `chama`.

* `mutrino/`: Contains files to drive builds on SNL machine mutrino.

* `ride/`: Contains the files to drive builds on the SRN test bed machine
  `ride` which also can be run on the SON machine `white`.
 
* `sems-rhel6/`: Contains files to drive builds on rhel6 machines with the SEMS
  environment.

* `sems_gcc-7.2.0/`: Contains driver scripts for an on-off GCC 7.2.0 build
  based on the SEMS system.  This build really does not fit into the system
  described above but it put in this directory since it is targeted to support
  ATDM.  It also shows that a given system can have its own driver files if it
  needs to.

* `serrano/`: Contains files to drive builds on the SRN HPC machine `serrano`.

* `shiller/`: Contains the files to drive builds on the SRN test bed machine
  `shiller` which also can be run on the SON machine `hansen`.

* `waterman/`: Contains files to drive builds on the SRN Test Bed machine
  `waterman`.


## How add a new system

To add a new system, first add a new `elseif` statement for the new system in
the file:

```
  Trilinos/cmake/std/atdm/utils/get_known_system_name.sh
```

Note that more than one `hostname` machine may map to the same
`<new_system_name>` (e.g. both `white` and `ride` machines map to the system
`ride`).

The variable `ATDM_HOSTNAME` (set to exported variable
`ATDM_CONFIG_KNOWN_HOSTNAME`) is used for the CDash site name.  This makes it
so that any node `white05`, `white12`, etc. just says `white` on CDash.  This
is important for the CDash 'next' and 'previous' relationships to work.  (But
for `CTEST_TEST_TYPE=Experimental` builds, the real `hostname` is used which
is stored in the exported env variable `ATDM_CONFIG_REAL_HOSTNAME`.  This
ensures that queries with `cdash/queryTests.php` don't accidentally pick up
tests from "Experimental" builds.)

The variable `ATDM_SYSTEM_NAME` (set to the exported variable
`ATDM_CONFIG_KNOWN_SYSTEM_NAME`) must be set to `<new_system_name>` which is
selected for this new system type.

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
standard logic is executed in in `atdm/utils/set_build_options.sh`.  (To see
an example of the usage of this file, see `atdm/cee-rhel6/custom_builds.sh`
and `atdm/cee-rhel6/environment.sh`.)

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
`"${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[i]}.sh"
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
described <a href="#running-locally-and-debugging">above</a>.

Once the basic configurations seem to be working, then commit your changes on
a topic branch and [submit a pull request (PR) to
Trilinos](https://github.com/trilinos/Trilinos/wiki/Submitting-a-Trilinos-Pull-Request).
Make sure and mention `@fryeguy52` and add the label `ATDM`.

Once the PR has been has been merged, then set up the Jenkins jobs to run the
builds as described <a href="#setting-up-jenkins-jobs">above</a>. Please note that the 
variable `JOB_NAME` is set by Jenkins and is the name of currently running job.  Your 
driver files in `Trilinos/cmake/ctest/drivers/atdm/<new_system_name>/drivers/` must 
be named to exactly match the Jenkins `JOB_NAME` variable.

ToDo: Fill in more detail, add an FAQ, etc.
