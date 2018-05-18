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

That file reads `JOB_NAME` from the env and uses it to set the CDash build
name.  (Therefore, the Jenkikns `JOB_NAME` is the same as the CDash build name
for all of these Trilinos ATDM builds.)  It also other CMake and CTest options
that are pulled out of the env set by the
`cmake/std/atdm/<system_name>/environment.sh` script.  (See `$ENV{<varName>}`
to see what vars are pulled out of the env.)

This directory contains a CTest -S driver script:

```
  Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake
```

which is called using:

```
  ctest -S <base-dir>/Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake
```

which can be run using any way desired and it will clone a new Trilinos git
repo (if not cloned already).  (But this file get directly run by the
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
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh
```

The first script `ctest-s-driver-config-build.sh` must be sourced instead of
just run because it sets up the env and changes to a subdir.  The env must be
loaded on the login/compile node and not the compute node and therefore, it
must be done this way.

With newer versions of CMake (3.10+), this will result in all of the results
going to the same build (i.e. row) on CDash.  For older versions of CMake, the
test results show up as a separate build (same `site` and `build name` but
different `build stamp`).

## Running locally and debugging

To run locally, first set up a mock Jenkins workspace directory structure and
set up symlinks to the local Trilinos git repo as:

```
$ cd <some_base_build_dir>/
$ mkdir MOCK_jenkins_driver
$ cd MOCK_jenkins_driver/
$ ln -s <some_base_dir>/Trilinos .
$ mkdir SRC_AND_BUILD
$ cd SRC_AND_BUILD/
$ ln -s <some_base_dir>/Trilinos .
$ cd ..
```

Then any of these builds can be tested locally without submitting to CDash with:

```
$ cd <some_base_build_dir>/MOCK_jenkins_driver/
$ time env \
    JOB_NAME=<some-build-name> \
    WORKSPACE=$PWD \
    Trilinos_PACKAGES=Kokkos,Teuchos,Tpetra \
    CTEST_TEST_TYPE=Experimental \
    CTEST_DO_SUBMIT=OFF \
    CTEST_DO_UPDATES=OFF \
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE \
  <some_base_dir>/Trilinos/cmake/ctest/drivers/atdm/<system_name>/local-driver.sh \
    &> console.out
```

(Where it is **CRITICAL** that you set `CTEST_DO_UPDATES=OFF` or it will hard
reset your local Trilinos git repo!)

Or if a `<system_name>/drivers/<some-build-name>.sh` file
(e.g. `Trilinos-atdm-hansel-shiller-gnu-debug-openmp.sh`) already exists,
instead use:

```
$ cd <some_base_build_dir>/MOCK_jenkins_driver/
$ time env \
    JOB_NAME=<some-build-name> \
    WORKSPACE=$PWD \
    Trilinos_PACKAGES=Kokkos,Teuchos,Tpetra \
    CTEST_TEST_TYPE=Experimental \
    CTEST_DO_SUBMIT=OFF \
    CTEST_DO_UPDATES=OFF \
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE \
  <some_base_dir>/Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh \
    &> console.out
```

Then you can loot at the `console.out` file and examine the `*.xml` configure,
build, and test files created under:

```
  <some_base_build_dir>/MOCK_jenkins_driver/SRC_AND_BUILD/BUILD/Testing/
```

to see if it is doing the right thing.

If that looks good, then you can do an experimental submit (avoiding a rebulid)
with:

```
$ cd <some_base_build_dir>/MOCK_jenkins_driver/
$ time env \
    JOB_NAME=<some-build-name> \
    WORKSPACE=$PWD \
    Trilinos_PACKAGES=Kokkos,Teuchos,Tpetra \
    CTEST_TEST_TYPE=Experimental \
    CTEST_DO_SUBMIT=ON \
    CTEST_DO_UPDATES=OFF \
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE \
  <some_base_dir>/Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh \
    &> console.out
```

If that submit looks good, then the job is ready to set up as a Jenkins job.

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
      * "Deadline tolerence in minues": `1`
* "Build"
  * "Execute shell"
    * "Command": `Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh`

But be careful not to put any additional settings that what is absolutely
necessary because Jenkins configurations are not under version control and
there is no tractability for changes in these settings!

## Specific <system_name> directories

The following `<system_name>` sub-directories exist (in alphabetical order):
* `mutrino/`: Contains files to drive builds on SNL machine mutrino.
 
* `rhel6/`: Contains files to drive builds on rhel6 machines with the SEMS
  environment.

* `ride/`: Contains the files to drive builds on the SRN test bed machine
  `ride` which also can be run on the SON machine `white`.

* `sems_gcc-7.2.0/`: Contains driver scripts for a on-off GCC 7.2.0 build
  based on the SEMS system.  This build really does not fit into the system
  described above but it put in this directory since it is targeted to support
  ATDM.  It also shows that a given system can have its own driver files if it
  needs to.

* `shiller/`: Contains the files to drive builds on the SRN test bed machine
  `shiller` which also can be run on the SON machine `hansen`.
  
* `toss3/`: Contains files to drive builds on the SRN HPC machines `serrano`
  and `chama`.

## How add a new system

ToDo: Fill in!
