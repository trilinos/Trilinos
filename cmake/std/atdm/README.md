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
* <a href="#checkin-test-atdmsh">checkin-test-atdm.sh</a>
* <a href="#specific-instructions-for-each-system">Specific instructions for each system</a>
* <a href="#troubleshooting-configuration-problems">Troubleshooting configuration problems</a>
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
* `cuda`: Do a CUDA build (`<COMPILER>=CUDA`, `NODE_TYPE=CUDA`)
  - `cuda-8.0`: Use CUDA 8.0
  - `cuda-9.0`: Use CUDA 9.0

When using `gnu`, `intel`, `clang`, and `cuda` without specifying a version
(e.g. `cuda-9.0`), then a default version of the compilers for that system
will be chosen (see the loaded env for the default chosen version).  Each
system may only support a subset of these compilers; see the
`cmake/std/atdm/<system-name>/environment.sh` file for details on which
compilers and which versions are supported.  If you choose a compiler that is
not supported, an error message will be provided.  If `default` is used, then
the default compiler for the system will be selected.

The following `<job-name>` keywords specify debug or optimized `<BUILD_TYPE>
`(used for the CMake cache var `CMAKE_BUILD_TYPE`with default
`<BUILD_TYPE>=DEBUG`):

* `debug`: Use `<BUILD_TYPE>=DEBUG`
* `opt`: Use `<BUILD_TYPE>=RELEASE`

The following `<job-name>` keywords determine the Kokkos threading model
`<NODE_TYPE>` (default is `<NODE_TYPE>=SERIAL` unless `<COMPILER>=CUDA`):

* `openmp`: Use OpenMP for host threading (`NODE_TYPE=OPENMP`)
* `pthread`: Use Pthreads for host threading (`NODE_TYPE=THREAD`)
* `serial`: Use no host threading (`NODE_TYPE=SERIAL`, DEFAULT)

If `cuda` (or `cuda-8.0`, `cuda-9.0`, etc.) is given, then `<NODE_TYPE>` is
automatically set to `CUDA`.

All other strings in `<job-name>` are ignored but are allowed for
informational purposes.  The reason that a `<job-name>` string is defined in
this form is that this can be used as the Jenkins job name and the Trilinos
build name that shows up on CDash.  This makes it very easy to define the
configuration options and maintain the Jenkins build jobs.  The combination
`<COMPILER>-<BUILD_TYPE>-<NODE_TYPE>` is used to define the CMake variable
`ATDM_JOB_NAME_KEYS_STR` that is used to uniquely define a build on a
particular system and to manage a system of tweaks for each of the supported
builds (see below).

Some examples of `<job-name>` keyword sets used on various platforms include:
* `gnu-debug-openmp`
* `gnu-opt-openmp`
* `intel-debug-openmp`
* `intel-opt-openmp`
* `cuda-debug` (`<NODE_TYPE>` is implicitly `CUDA`)
* `cuda-opt` (`<NODE_TYPE>` is implicitly `CUDA`)
* `cuda-8.0-debug` (`<NODE_TYPE>` is implicitly `CUDA`)
* `cuda-9.0-opt` (`<NODE_TYPE>` is implicitly `CUDA`)

The script `cmake/std/atdm/load-env.sh` when sourced sets some bash
environment variables that are prefixed with `ATDM_CONFIG_` and other standard
variables.

The file `ATDMDevEnv.cmake` pulls bash environment variables set by the
sourced `atdm/load-env.sh` script and sets up a number of CMake cache
variables such as the compilers, compiler options, and sets up a number of
Trilinos configuration variables (many of which are echoed to the STDOUT when
running `cmake`).  This also enables by default all of the standard TPLs used
by ATDM Application codes through Trilinos and sets up locations for the
include directories and libraries for each of these TPLs.

When included, the file `ATDMDevEnv.cmake` also disables many packages and
subpackages not used for the ATDM configuration of Trilinos.  This uses a
so-called back-listing approach which allows one to only directly enable the
packages you want to use with `Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON` and
then disable the ones you don't want.  This is a much more flexible way to
define a standard configuration of Trilinos that allows different sets of
actual packages to be enabled based on the needs of the different ATDM
application customers and Trilinos developers just needing to enable a subset
of packages.  But if package X does get enabled, then it will always have the
same configuration options independent of any other packages that are enabled.

When `ATDMDevEnv.cmake` is being processed, if there is a "tweaks" file
defined for a build, then it will be picked up in the CMake cache var <a
href="#ATDM_TWEAKS_FILES">ATDM_TWEAKS_FILES</a> and that file will be read in
using `INCLUDE()` to process the extra options contained within it.


## checkin-test-atdm.sh

For those Trilinos developers comfortable with using the
`Trilinos/checkin-test.py` script, multiple local builds and testing on a
system can also be driven with the provided `checkin-test-atdm.sh` wrapper
script.  This can be used to drive a number of builds on system as:

```
$ cd <some_build_dir>/

$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-sems.sh .

$ ./checkin-test-sems.sh <job-name-0> <job-name-1> ... \
  --enable-packages=<Package> --local-do-all
```

That will configure, build, and run tests for each specified build
`<job-name-0>` and send a summary email when complete.  See comments at the
top of the script `checkin-test-atdm.sh` for more details.  The parallel level
for building and running tests are determined by the env vars
`ATDM_CONFIG_BUILD_COUNT` and `ATDM_CONFIG_CTEST_PARALLEL_LEVEL`,
respectfully, as set by default for the given system.  These can be overridden
by setting the env vars `ATDM_CONFIG_BUILD_COUNT_OVERRIDE` and
`ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERIDE`, respectfully as, for example:

```
$ env \
  ATDM_CONFIG_BUILD_COUNT_OVERRIDE=8 \
  ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERIDE=12 \
  ./checkin-test-sems.sh ...
```

A value of `ATDM_CONFIG_BUILD_COUNT_OVERRIDE=0` or less than `0` is allowed
when using Ninja (i.e. `ATDM_CONFIG_USE_NINJA=ON`) in which case `ninja` will
be run with non `-j<N>` argument, and therefore all of the non-loaded cores
will be used.

Note that to run tests for a CUDA build or to run tests on platforms that must
run on a compute node one will need to run these on a compute node on the
system that has a GPU.  On such a system one would run:

```
$ ./checkin-test-sems.sh <job-name-0> <job-name-1> ... \
  --enable-packages=<Package> --configure --build \
  && \
  <command-to-run-on-compute-node> \
  ./checkin-test-sems.sh <job-name-0> <job-name-1> ... \
  --enable-packages=<Package> --test
```

See <a href="#specific-instructions-for-each-system">Specific instructions for
each system</a> for details.

Note that one can create a `local-checkin-test-defaults.py` file to set
defaults like:

```
defaults = [
  "--no-enable-fwd-packages",
  "--enable-all-packages=off",
  ]
```

and then run:

```
$ ./checkin-test-sems.sh <job-name-0> <job-name-1> ... \
  --enable-packages=<Package> --local-do-all
```

However, a default `local-checkin-test-defaults.py` is created the first time
the `checkin-test-atdm.sh` script is run.

## Specific instructions for each system

* <a href="#ridewhite">ride/white</a>
* <a href="#shillerhansen">shiller/hansen</a>
* <a href="#chamaserrano">chama/serrano</a>
* <a href="#mutrino">mutrino</a>
* <a href="#sems-rhel6-environment">SEMS rhel6 environment</a>


### ride/white

Once logged on to `white` (on the SON) or `ride` (on the SRN), one can
directly configure and build on the login node (being careful not to overload
the node).  But to run the tests, one must run on the compute nodes using the
`bsub` command to run if using a CUDA build.  For example, to configure, build
and run the tests for the `cuda-debug` build for say `MueLu` on `white`,
(after cloning Trilinos on the `develop` branch) one would do:

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cuda-debug

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ bsub -x -Is -q rhel7F -n 16 ctest -j16
```

The ATDM configuration of Trilinos is set up to run on the Firestone nodes
(Dual-Socket POWER8, 8 cores per socket, K80 GPUs).  This confiugration will
not work on the other GPU nodes currently.

Note that one can also run the same build a tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-sems.sh .
$ bsub -x -I -q rhel7F -n 16 \
  ./checkin-test-sems.sh cuda-debug \
  --enable-all-packages=off --no-enable-fwd-packages \
  --enable-packages=MueLu \
  --local-do-all
```


### shiller/hansen

Once logged on to `hansen` (on the SON) or `shiller` (on the SRN), one can
directly configure and build on the login node (being careful not to overload
the node).  But to run the tests, one must run on the compute nodes using the
`srun` command.  For example, to configure, build and run the tests for say
`MueLu` on `hansen`, (after cloning Trilinos on the `develop` branch) one
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

Note that one can also run the same build a tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-sems.sh .
$ ./checkin-test-sems.sh intel-opt-openmp \
  --enable-all-packages=off --no-enable-fwd-packages \
  --enable-packages=MueLu \
  --configure --build
$ srun ./checkin-test-sems.sh intel-opt-openmp \
  --enable-all-packages=off --no-enable-fwd-packages \
  --enable-packages=MueLu \
  --test
```

When using `srun`, for some reason, one must configure and build on the login
node but then use `srun` to run tests on the compute node.  For some reason,
the `ninja` program in the `PATH` can't be found on the compute node when
using `srun` with the checkin-test-atdm.sh script.

WARNING: One must use `srun` and **not** `salloc` to allocate and run on a
compute node.  The way that these machines are currently set up, running with
`salloc` the command `squeue` will seem to show that the job is allocated and
running on a compute node but the job is actually running on the login node!
(Running `ps -AF | grep <script-name>` and `top` on the login node will
clearly show that that the job is actually running on the login node and
presumably the actual compute node is idle.)


### chama/serrano

Once logged on to `chama` or `serrano`, one can directly configure and build
on the login node (being careful not to overload the node).  But to run the
tests, one must run on the compute nodes using the `srun` command.  For
example, to configure, build and run the tests for say `MueLu` on `serrano`
or `chama`, (after cloning Trilinos on the `develop` branch) one would do:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-opt-openmp

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ salloc -N1 --time=0:20:00 --account=<YOUR_WCID> ctest -j16
```

To get information on <YOUR_WCID> used above, there is a WC tool tab on
computing.sandia.gov

To use the checkin-test-atdm.sh script, you must split running the tests from
the configure and build as with:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-sems.sh .
$ ./checkin-test-sems.sh intel-opt-openmp \
  --enable-packages=MueLu --allow-no-pull --configure --build
$ salloc -N1 --time=0:20:00 --account=<YOUR_WCID> \
  ./checkin-test-sems.sh intel-opt-openmp \
  --enable-packages=MueLu --test
```

### mutrino

Once logged on to `mutrino`, one can directly configure and build
on the login node (being careful not to overload the node).  But to run the
tests, one must run on the compute nodes using the `salloc` command.  For
example, to configure, build and run the tests for say `MueLu` on `mutrino`, 
(after cloning Trilinos on the `develop` branch) one would:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-opt-openmp

$ cmake \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make -j16

$ salloc -N 1 -p standard -J $JOB_NAME ctest -j16
```

### SEMS rhel6 environment

Once logged on to a rhel6 machine with the sems NFS, one can directly
configure, build, and run tests.  For example, to configure, build and run the
tests for `MueLu` one would clone Trilinos on the `develop` branch and then
do the following:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-opt-openmp

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ ctest -j16 \
```

## Troubleshooting configuration problems

There are situations were a particular developer on a particular system will
not be able to cleanly reproduce a Trilinos configuration as described above.
Some of the issues that can cause that are described here.

The first set of problems that might occur are that the version of `git` that
is loaded by sourcing `Trilinos/cmake/std/atdm/load-env.sh` for a given system
may conflict with the developer's global `~/.gitconfig` file on that system.
It may be that an older or newer version of git is loaded by the script
(i.e. `<atdm-loaded-git>`) than the git version that the developer has been
using on the system (i.e. `<developer-git>`) and therefore some of the
settings in the developer's `~/.gitconfig` file may not be compatible with the
different version of git.  One approach to address this is run `module list`
to see the version of git that is loaded by this script and then to swap out
that version of git for the developer's regular version of git using `module
swap <atdm-loaded-git> <developer-git>`.  Or, the developer can just call
`export PATH=<path-to-developers--git>/bin:$PATH` after sourcing
`Trilinos/cmake/std/atdm/load-env.sh` to select the developer's regular
version of git.  The git commands that are used by TriBITS, Trilinos, and the
ATDM Trilinos configuration and testing scripts should work with just about
any version of git 2.0+.  Another approach would be for the user to
(temporarily) edit their `~/.gitconfig` file to address the problems.

## Directory structure and contents

This base directory:

```
  cmake/std/atdm/
```

contains the following files:

* **ATDMDevEnvSettings.cmake**: Reads vars out of the env (loaded with
  `load-env.sh`) to set compilers and other options for the Trilinos build.
  It also default enables all of the TPLs used in ATDM configurations of
  Trilinos.  However, this file does **not** include `ATDMDisables.cmake`.
  Therefore, this file can be used to set up builds of Trilinos that enable
  arbitrary sets of packages.

* **ATDMDisables.cmake**: Disables a bunch of Trilinos packages and
  subpackages not used by ATDM application customers.  This file gets included
  automatically in `ATDMDevEnv.cmake` (so you don't need to list it in local
  configures of Trilinos).  But this file is also included in the outer `ctest
  -S` driver script code for ATDM builds of Trilinos which is needed for
  correct package-by-package testing of Trilinos.

* **ATDMDevEnv.cmake**: Includes `ATDMDevEnvSettings.cmake` and
  `ATDMDisables.cmake`.

* **checkin-test-atdm.sh**: Uses the `Trilinos/checkin-test.py` script to
  drive builds and tests on the given platform.  (See comments in the top of
  the script for instructions.)

Each supported ATDM system `<system-name>` has its own sub-directory with the
contents:

```
  <system-name>/
    environment.sh  # Load env for the given system based on $JOB_NAME keys
    all_supported_builds.sh  # [Optional] List of all supported builds
    tweaks/
       <COMPILER0>-<BUILD_TYPE0>-<NODE_TYPE0>.cmake  # [Optional]
       <COMPILER1>-<BUILD_TYPE1>-<NODE_TYPE1>.cmake  # [Optional]
       ...
```

The optional file `<system-name>/all_supported_builds.sh` contains a list of
all of the supported builds on the system.  This sets the environment variable `ATDM_CONFIG_ALL_SUPPORTED_BUILDS` as:

```
  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS="gnu-debug-openmp gnu-opt-openmp ..."
```

This is used in the `checkin-test-atdm.sh` script to run all of the builds for
a system with `checkin-test-atdm.sh all [other options]`.

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
* `chama/`: Supports SNL HPC machine `chama`.

* `mutrino/`: Supports SNL HPC machine `mutrino`.

* `rhel6/`: RHEL6 systems with the SEMS NFS environment

* `ride/`: Supports GNU and CUDA builds on both the SRN machine `ride` and the
  mirror SON machine `white`.

* `serrano/`: Supports SNL HPC machine `serrano`.

* `shiller/`: Supports GNU, Intel, and CUDA builds on both the SRN machine
  `shiller` and the mirror SON machine `hansen`.
  

