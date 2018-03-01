# Local ATDM builds of Trilinos

This directory `cmake/std/atdm/` contains a set of `*.sh` shell scripts and
`*.cmake` files that define a standard ATDM-focused build of Trilinos on a
number of different ATDM platforms.

This is used to define a set of automated builds of Trilinos that are run with
Jenkins and post to the main Trilinos CDash site using build names that start
with `Trilinos-atdm`.  The CTest driver scripts and Jenkins driver scripts
that use this ATDM configuration are defined in the Trilinos directory
`cmake/ctest/drivers/atdm/` (see the `README.md` file in that directory for
details) but these details are not necessary in order to just reproduce a
build locally as described below.

**Outline:**
* <a href="#quick-start">Quick-start</a>
* <a href="#specific-instructions-for-each-system">Specific instructions for each system</a>
* <a href="#directory-structure-and-contents">Directory structure and contents</a>
* <a href="#specific-systems-supported">Specific systems supported</a>

## Quick-start

After [cloning the Trilinos git
repo](https://github.com/trilinos/Trilinos/wiki/VC-%7C-Initial-Git-Setup) on
one of the supported ATDM machines, a local configure of Trilinos enabling a
few packages is performed as:

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh <job-name>

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_<Package1>=ON \
  $TRILINOS_DIR

$ make NP=16  # Uses ninja -j16

$ ctest -j16  # Might need to be run with srun or some other command, see below
```

The command:

```
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh <job-name>
```

determines what machine you are on (using `hostname`) and then loads the
correct environment automatically for that machine and for the build options
passed through in `<job-name>` (or errors out if the current machine is not
one of the supported machines).

The `<job-name>` argument is a single string of the form
`XXX-<keyword0>-<keyword1>-...`.

The following `<job-name>` keywords specify the `<COMPILER>`:

* `gnu`: Use the GCC compilers (`<COMPILER>=GNU`)
* `intel`: Use the Intel compilers (`<COMPILER>=INTEL`)
* `clang`: Use the LLVM Clang compilers (`<COMPILER>=CLANG`)
* `cuda`: Do a CUDA build for that system (`<COMPILER>=CUDA`, `NODE_TYPE=CUDA`)

The following `<job-name>` keywords specify debug or optimized `<BUILD_TYPE>
`(used for the CMake cache var `CMAKE_BUILD_TYPE`with default
`<BUILD_TYPE>=DEBUG`):

* `opt`: Use `<BUILD_TYPE>=RELEASE`

The following `<job-name>` keywords determine the Kokkos threading model
`<NODE_TYPE>` (default is `<NODE_TYPE>=SERIAL` unless `<COMPILER>=CUDA`):

* `openmp`: Use OpenMP for host threading (`NODE_TYPE=OPENMP`)
* `pthread`: Use Pthreads for host threading (`NODE_TYPE=THREAD`)
* `serial`: Use no host threading (`NODE_TYPE=SERIAL`, DEFAULT)

If `cuda` is given, then `<NODE_TYPE>` is automatically set to `CUDA`.

All other strings in `<job-name>` are ignored but are allowed for
informational purposes.  The reason that a `<job-name>` string is defined in
this form is that this can be used as the Jenkins job name and the Trilinos
build name that shows up on CDash.  This makes it very easy to define the
configuration options and maintain the Jenkins build jobs.  The combination
`<COMPILER>-<BUILD_TYPE>-<NODE_TYPE>` is used to define the CMake variable
`ATDM_JOB_NAME_KEYS_STR` that is used to uniquely define a build on a
particular system (see below).

Some examples of `<job-name>` keyword sets used on various platforms include:
* `gnu-debug-openmp`
* `gnu-opt-openmp`
* `intel-debug-openmp`
* `intel-opt-openmp`
* `cuda-debug` (`<NODE_TYPE>` is implicitly `CUDA`)
* `cuda-opt` (`<NODE_TYPE>` is implicitly `CUDA`)

The script `cmake/std/atdm/load-env.sh` when sourced sets a set of bash
environment variables that are prefixed with `ATDM_CONFIG_` and other standard
variables.

The file `ATDMDevEnv.cmake` pulls bash environment variables set by the
sourced `atdm/load-env.sh` script and sets up a number of CMake cache
variables such as the compilers, compiler options, and sets up a number of
Trilinos configuration variables (many of which are echoed to the STDOUT when
running `cmake`).  This also default enables all of the standard TPLs used by
ATDM Application codes through Trilinos and sets up locations for the include
directories and libraries for each of these TPLs.

When included, the file `ATDMDevEnv.cmake` also disables many packages and
subpackages not used for the ATDM configuration of Trilinos.  This uses a
so-called back-listing approach which allows one to only directly enable the
packages you want to use with `Triinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON` and
then disable the ones you don't want.  This is a much more flexible way to
define a standard configuration of Trilinos that allows different sets of
actual packages to be enabled based on the needs of the different ATDM
application customers and Trilinos developers just needing to enable a subset
of packages.

When `ATDMDevEnv.cmake` is being processed, if there is a "tweaks" file
defined for a build, then it will be picked up in the CMake cache var <a
href="#ATDM_TWEAKS_FILES">ATDM_TWEAKS_FILES</a> and that file will be read in
using `INCLUDE()` to process the extra options contained within it.

## Specific instructions for each system

* <a href="#hansenshiller">hansen/shiller</a>
* ???

### hansen/shiller

Once logging on to `hansen` (on the SON) or `shiller` (on the SRN), one can
directly configure and build on the login node (being careful not to overload
the node).  But to run the tests, one must run on the compute nodes using the
`srun` command.  For example, to configure, build and run the tests for say
`MueuLu` on `hansen`, (after cloning Trilinos on the `develop` branch) one
would do:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-opt-openmp

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ srun ctest -j16
```

## Directory structure and contents

This base directory:

```
  cmake/std/atdm/
```

contains the following files:

* **ATDMDevEnv.cmake**: Reads vars out of the env (loaded with `load-env.sh`)
  to set compilers and other options for the Trilinos build.

* **ATDMDisables.cmake**: Disables a bunch of Trilinos packages and
  subpackages not used by ATDM application customers.  This file gets included
  automatically in `ATDMDevEnv.cmake` (so you don't need to list it in local
  configures of Trilinos).  But this file is also included in the outer `ctest
  -S` driver script code for ATDM builds of Trilinos which is needed for
  correct package-by-package testing of Trilinos.

* **ATDMDevEnvUtils.cmake**: Defines some simple macros and functions used in
  the above `*.cmake` files.

Each supported ATDM system `<system-name>` has its own sub-directory with the
contents:

```
  <system-name>/
    environment.sh  # Load env for the given system based on $JOB_NAME keys
    tweaks/
       <COMPILER0>-<BUILD_TYPE0>-<NODE_TYPE0>.cmake
       <COMPILER1>-<BUILD_TYPE1>-<NODE_TYPE1>.cmake
       ...
```

<a name="ATDM_TWEAKS_FILES"/>

The files in the `cmake/std/atdm/<system-name>/tweaks/` directory contain
special settings for specific builds for a specific system.  Typically, this
file contains (temporary) disables for tests for that given build.  When a
configure is performed, the variable `ATDM_JOB_NAME_KEYS_STR` set to
`<COMPILER>-<BUILD_TYPE>-<NODE_TYPE>` (printed to STDOUT) is used to define a
default file name:

```
  Trilinos/cmake/std/atdm/<system-name>/tweaks/$ATDM_JOB_NAME_KEYS_STR.cmake
```

If that file exists, then it is set as the default for the cmake cache var
`ATDM_TWEAKS_FILES` (printes to STDOUT) and that file is included and its
options are read.

## Specific systems supported

The specific `cmake/std/atdm/<system-name>/` sub-directories and the systems
they support are:

* `shiller/`: Supports GNU, Intel, and CUDA builds on both the SRN machine
  `shiller` and the mirror SON machine `hansen`.

* ???
