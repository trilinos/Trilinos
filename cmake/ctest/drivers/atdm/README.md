# CDash and Jenkins drivers for ATDM Trilinos builds

The files and scripts in this directory drive builds of Trilinos on the
various ATDM test platforms using Jenkins jobs and submit results to the
Trilinos CDash site.

**Outline:**
* <a href="#base-ctestcdash-configuration">Base CTest/CDash configuration</a>
* <a href="#system-specific-driver-files">System-specific driver files</a>
* <a href="#running-locally-and-debugging">Running locally and debugging</a>
* <a href="#setting-up-jenkins-jobs">Setting up Jenkins jobs</a>
* <a href="#specific-system_name-directories">Specific <system_name> directories</a>
* <a href="#how-add-a-new-system">How add a new system</a>

## Base CTest/CDash configuration

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
for all of these Trilinos ATDM builds.)

This directory also contains the file:

```
  Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
```

which sets up and runs the system-specific ctest -S script.

## System-specific driver files

Each system `<system_name>` is given a subdirectory under this directory:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/
```

That directory contains a CTest -S driver script:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/ctest-driver.cmake
```

which is directly run with `ctest -S`.

A good example of this is:

```
  Trilinos/cmake/ctest/drivers/atdm/shiller/ctest-driver.cmake
```

The `<system_name>/ctest-driver.cmake` files sets any CMake or CTest options
that are speicfic to that system like:

* `CTEST_NOTES_FILES`: Append any extra notes files to CDash.

* `CTEST_BUILD_FLAGS`: The flags to pass to the build (e.g. `"-j32 -k 999999"`
  for Ninja)

* `CTEST_PARALLEL_LEVEL`: Set the parallel level for runnign tests (e.g. `16`
  if using 2 MPI threads per core).

* `CTEST_CMAKE_GENERATOR`: Set to `Ninja` if using Ninja on that system.

* `CTEST_SITE`: Set to the name of the site to show up on CDash (i.e. to avoid
  different site names for different nodes like `shiller01`, `shiller02`, or
  `shiller03` messing up previous/next on CDash.

Given `<system_name>/ctest-driver.cmake`, then:

```
  ctest -S <base-dir>/Trilinos/cmake/ctest/drivers/atdm/<system_name>/ctest-driver.cmake
```

can be run using any way desired and it will clone a new Trilinos git repo (if
not cloned already).

In addition, each `atdm/<system_name>/` directory contains a local driver
script:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/local-driver.sh
```

which runs `ctest -S` in the appropriate way for that system.  This script
assumes a directory structure as set up by Jenkins but does not really require
Jenkins to run it.

The sub-directory:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
```

contains specific drivers with the names of the Jenkins build names:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/$JOB_NAME.sh
```

This file sets some tweaks for that particular build `$JOB_NAME` link which
CDash Track/Group results are sent to and other tweaks like this.

To run these files automatically from the Jenkins build configuration, the
file:

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/smart-jenkins-driver.sh
```

is provided which just runs:

```
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/drivers/$JOB_NAME.sh
```

## Running locally and debugging

To run locally, first set up a mock Jenkins workspace directory structure and
set up symlinks to local

```
$ cd <some_base_build_dir>
$ mkdir MOCK_jenkins_driver
$ cd MOCK_jenkins_driver/
$ ln -s <some_base_dir>/Trilinos .
$ mkdir SRC_AND_BUILD
$ cd SRC_AND_BUILD/
$ ln -s <some_base_dir>/Trilinos .
$ cd ..
```

Then any of these builds can be first tested out with:

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

or if the `<system_name>/drivers/$JOB_NAME.sh` file already exists, instead use:

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
  <some_base_dir>/Trilinos/cmake/ctest/drivers/atdm/<system_name>/smart-jenkins-driver.sh \
    &> console.out
```

Then you can examine the `*.xml` configure, build, and test files created
under:

  <some_base_build_dir>/MOCK_jenkins_driver/SRC_AND_BUILD/BUILD/Testing/

to see if it is doing the right thing.

If that looks good, you can do an experimental submit (avoiding a rebulid)
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
  <some_base_dir>/Trilinos/cmake/ctest/drivers/atdm/<system_name>/smart-jenkins-driver.sh \
    &> console.out
```

## Setting up Jenkins jobs

To set up a Jenkins build, you must set the following in the Jenkins build
configuation GUI:

* "Source Code Management"
  * "Git"
    * "Repository URL": `https://github.com/trilinos/Trilinos.git` (until we
         can figure out how to clone from `sofware.sandia.gov:/git/Trilinos`)
* "Build Triggers"
  * "Build Periodically" (checked)
    * "Schedule": `H 0 ***`
* "Build Environment"
  * "Abort of the build if it's stuck" (checked)
    * "Time-out strategy": `Deadline`
      * "Deadline time": `20:45`
      * "Deadline tolerence in minues": `1`
* "Build"
  * "Execute shell"
    * "Command": `Trilinos/cmake/ctest/drivers/atdm/<system_name>/smart-jenkins-driver.sh`

## Specific <system_name> directories

The following `<system_name>` sub-directories exist:

* `shiller/`: The SRN test bed machine `shiller` which also can be run on the
  SON machine `hansen`.

* ???

## How add a new system

ToDo: Fill in!







