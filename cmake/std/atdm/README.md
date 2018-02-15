# Local ATDM builds of Trilinos

This directory `cmake/std/atdm/` contains a set of `*.sh` shell scripts and
`*.cmake` files that define a standard ATDM-focused build of Trilinos on a
number of different ATDM platforms.

This is used to define a set of automated builds of Trilinos that are run with
Jenkins and post to the main Trilinos CDash site using build names that start
with `Trilinos-atdm`.  The CTest driver scripts and Jenkins driver scripts
that use this ATDM configuration are defined in the Trilinos directory
`cmake/ctest/drivers/atdm/` (see the `READM.md` file in that directory for
details) but these details are not necessary in order to just reproduce a
build locally as described below.

**Outline:**
* <a href="#quick-start">Quick-start</a>
* <a href="#specific-instructions-for-each-system">Specific instructions for each system</a>
* <a href="#directory-structure-and-contents">Directory structure and contents</a>

## Quick-start

After cloning the Trilinos git repo one of the supported ATDM machines, a
local configure of Trilinos enabling a few packages is performed as:

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
`XXX-<keyword0>-<keyword1>-...` where the following keywords are recognized:

* `gnu`: Use the GCC compilers
* `intel`: Use the Intel compilers
* `cuda`: Do a CUDA build for that system
* `openmp`: Use OpenMP for threading

All other keywords are ignored but are allowed for informational purposes.
The reason that a `<job-name>` string is used in this form is that this can be
used as the Jenkins job name and the Trilinos build name that shows up on
CDash.  This makes it very easy to define the configuration options and
maintain the Jenkins build jobs.

The file `ATDMDevEnv.cmake` pulls vars out of the env set by the
`atdm/load-env.sh` script and sets up a number of CMake cache variables such
as the compilers, compiler options, and sets up a number of Trilinos
configuration variables (many of which are echoed to the STDOUT when running
`cmake`).

Note that the file `ATDMDevEnv.cmake` also disables many packages and
subpackages not used for the ATDM configuration of Trilinos.  This uses a
so-called back-listing approach which allows one to only directly enable the
packages you want to use with `Triinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON` and
then disable the ones you don't want.  This is a much more flexible way to
define a standard configuration of Trilinos that allows different sets of
actual packages to be enabled based on the needs of the different ATDM
application customers and Trilinos developers just needing to enable a subset
of packages.

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

* **ATDMDevEnv.cmake**: Reads vars out of the env (loaded with `load-env.sh`) to
  set compilers and other options for the Trilinos build.

* **ATDMDisables.cmake**: Disables a bunch of Trilinos packages and subpackages
  not used by ATDM application customers.  This file gets included
  automatically in `ATDMDevEnv.cmake` (so you don't need to list it in local
  configures of Trilinos).  But this file is also included in the outer `ctest
  -S` driver script code for ATDM builds of Trilinos which is needed for
  correct package-by-package testing of Trilinos.

* **ATDMDevEnvUtils.cmake**: Defines some simple macros and functions used in
  the above `*.cmake` files.

Each supported ATDM system has its own sub-directory that contains the
`environment.sh` script for that system.  The sub-directories and the systems
they support are:

* `shiller/`: Supports GNU, Intel, and CUDA builds on both the SRN machine
  `shiller` and the mirror SON machine `hansen`.

* ???
