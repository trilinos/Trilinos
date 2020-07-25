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
* <a href="#parallel-build-and-test-processes">Parallel build and test processes</a>
* <a href="#installation-and-usage">Installation and usage</a>
* <a href="#checkin-test-atdmsh">checkin-test-atdm.sh</a>
* <a href="#ctest-s-local-test-driversh">ctest-s-local-test-driver.sh</a>
* <a href="#specific-instructions-for-each-system">Specific instructions for each system</a>
* <a href="#building-and-installing-trilinos-for-atdm-applications">Building and installing Trilinos for ATDM Applications</a>
* <a href="#troubleshooting-configuration-problems">Troubleshooting configuration problems</a>
* <a href="#directory-structure-and-contents">Directory structure and contents</a>
* <a href="#disabling-failing-tests">Disabling failing tests</a>
* <a href="#specific-systems-supported">Specific systems supported</a>
* <a href="#custom-systems-and-configurations">Custom systems and configurations</a>


## Quick-start

After [cloning the Trilinos git
repo](https://github.com/trilinos/Trilinos/wiki/VC-%7C-Initial-Git-Setup) on
one of the supported ATDM machines, a local configure of Trilinos enabling a
few packages is performed using the bash shell (or opening a new bash shell
using `bash -l` if bash is not the user's default shell) as:

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh <build-name>

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_<Package1>=ON \
  $TRILINOS_DIR

$ make NP=16  # Uses ninja -j16

$ ctest -j4  # Might need to be run with srun or some other command, see below
```

The command:

```
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh <build-name>
```

determines what machine you are on (using `hostname`) and then loads the
correct environment automatically for that machine and for the build options
passed through in `<build-name>` (or errors out if the current machine is not
one of the supported machines) (see
[build-name-keywords](#build-name-keywords)).  An example of the output for
this command (on 'white') is:

```
$ source cmake/std/atdm/load-env.sh gnu-openmp-debug
Hostname 'white11' matches known ATDM host 'white' and system 'ride'
Setting compiler and build options for build name 'gnu-openmp-debug'
Using white/ride compiler stack GNU to build DEBUG code with Kokkos node type OPENMP and KOKKOS_ARCH=Power8
```

<a name="build-name-keywords"/>

**[build-name-keywords]** The `<build-name>` argument is a single string of
keywords of the form `XXX-<keyword0>-<keyword1>-...-YYY` (or
`XXX_<keyword0>_<keyword1>_..._YYY`, either separator is supported).  The
typical order and format of this string is:

    <system_name>-<use_mpi>-<kokkos_arch>-<compiler>-<kokkos_thread>-<rdc>-<complex>-<shared_static>-<release_debug>-<pt>

but any order of these keywords is supported.  Also, the keywords are
case-insensitive All of these keywords, except for `<compiler>` (which can be
`default`), are optional.  All of the other keywords have default values.  Any
strings not recognized as keywords are ignored.  (Therefore, misspelling the
name of a keyword is ignored.) See some examples of build name strings
[below](#build-name-examples).

Each of these keywords [`<system_name>`](#system_name),
[`<use_mpi>`](#use_mpi), [`<kokkos_arch>`](#kokkos_arch),
[`<compiler>`](#compiler), [`<kokkos_thread>`](#kokkos_thread),
[`<rdc>`](#rdc), [`<fpic>`](#fpic), [`<complex>`](#complex),
[`<shared_static>`](#shared_static), [`<release_debug>`](#release_debug), and
[`<pt>`](#pt), are described below.

<a name="system_name"/>

**`<system_name>`**: Typically, the system name is determined automatically by
examining the `hostname`, standard system env vars, or other files on the
machine and matching to a known supported system.  Therefore, it is typically
not necessary to specify `<system_name>` in the `<build-name>` keys string.
But there are some cases where more then one `<system_name>` env are supported
on the same machine.  For example, on CEE LAN RHEL6 machines, both the <a
href="#sems-rhel6-environment">sems-rhel6</a> and <a
href="#cee-rhel6-and-rhel7-environment">cee-rhel6</a> environments are
supported.  On these CEE LAN RHEL6 machines, when `cee-rhel6` is included in
`<build-name>`, then the `cee-rhel6` env will be selected.  But if
`sems-rhel6` is included in the build name (or no system name is listed in the
build name), then the `sems-rhel6` env will be selected by default on such
machines.  The same is true for CEE LAN RHEL7 machines with the <a
href="#sems-rhel6-environment">sems-rhel7</a> and <a
href="#cee-rhel6-and-rhel7-environment">cee-rhel6</a> environments.  And if
`spack-rhel` is included in `<build-name>`, then the <a
href="#spack-rhel-environment">spack-rhel</a> will attempted to be loaded.
(In that case, one must ensure that the ATDM Spack modules have been defined.)

<a name="kokkos_arch"/>

**`<kokkos_arch>`:** The `<build-name>` string can also contain keywords to
determine the `KOKKOS_ARCH` option of the build.  This is the case-insensitive
architecture name that is recognized by the CMake
[KOKKOS_ARCH](https://trilinos.org/docs/files/TrilinosBuildReference.html#configuring-with-kokkos-and-advanced-back-ends)
configure option for Trilinos and Kokkos.  Some common supported Kokkos
architectures for the host node include `BDW`, `HSW`, `Power8`, `Power9`, and
`KNL`.  When a GPU is present, some common Kokkos architecture options include
`Kepler37`, `Pascal60`, and `Volta70`.  (Note that the `KOKKOS_ARCH` keywords
are case-insensitive and can be `hsw`, 'volta70`, etc.)  If one selects a
`KOKKOS_ARCH` value that is not supported by the current system or selected
compiler, then the `load-env.sh` script will return an error message listing
the valid choices for `KOKKOS_ARCH` for each supported compiler.

This setup does not currently support specifying multiple `KOKKOS_ARCH` values
(since there is no example yet where that would be needed or useful) but such
functionality could be supported in the future if needed.

<a name="compiler"/>

**`<compiler>`:** The following **lower case** `<build-name>` keywords specify
the `<COMPILER>` variable include:

* `default`: Auto-select the default compiler and version for the given system
* `gnu`: Use the GCC compilers (`<COMPILER>=GNU`)
  - `gnu-4.9.3`: Use GNU GCC 4.9.3 compilers
  - `gnu-7.2.0`: Use GNU GCC 7.2.0 compilers
  - See what other GNU versions that may be supported on the given system
* `intel`: Use the Intel compilers (`<COMPILER>=INTEL`)
  - See what Intel versions are supported on the given system
* `clang`: Use the LLVM Clang compilers (`<COMPILER>=CLANG`)
  - See what Clang versions are supported on the given system
* `cuda`: Do a CUDA build (`<COMPILER>=CUDA`, `NODE_TYPE=CUDA`)
  - `cuda-8.0`: Use CUDA 8.0
  - `cuda-9.0`: Use CUDA 9.0
  - `cuda-9.2`: Use CUDA 9.2
  - See what CUDA versions are supported on the given system

When using `gnu`, `intel`, `clang`, and `cuda` without specifying a version
(e.g. `cuda-9.2`), then a default version of the compilers for that system
will be chosen (see the loaded env for the default chosen version).  To see
what compilers and compiler versions are supported for a given system, run
`source cmake/std/atdm/load-env.sh nothing` which will then print out an error
message listing them.  Each system may only support a subset of these
compilers and or may support multiple versions of a specific compiler.  (See
the optional file `cmake/std/atdm/<system-name>/custom_builds.sh` and the
required file `cmake/std/atdm/<system-name>/environment.sh` for details on
which compilers and which versions are supported for a given system.)  If
`default` is used, then the default compiler for the system will be selected.
Carefully examine STDOUT after running `source cmake/std/atdm/load-env
<build-name>` to see what compiler gets selected.

<a name="use_mpi"/>

**`<use_mpi>`:** The following `<build-name>` keywords determine if MPI is
  enabled or not in Trilinos (i.e. the value of `TPL_ENABLE_MPI`):

* `mpi`: Enable MPI and use MPI compler wrappers (`TPL_ENABLE_MPI=ON`, DEFAULT)
* `no-mpi`: Don't enable MPI and use raw compilers (except for CUDA builds that use `nvcc_wrapper` for the C++ compiler) (`TPL_ENABLE_MPI=OFF`)

NOTE: Setting `no-mpi` also switches to some non-MPI TPL builds and disables other TPLs like SuperLUDist that require MPI.

<a name="kokkos_thread"/>

**`<kokkos_thread>`:** The following `<build-name>` keywords determine the
Kokkos threading / back-end model variable `<NODE_TYPE>` (default is
`<NODE_TYPE>=SERIAL` unless `<COMPILER>=CUDA`):

* `serial`: Use no host threading (`NODE_TYPE=SERIAL`, DEFAULT)
* `openmp`: Use OpenMP for host threading (`NODE_TYPE=OPENMP`)

If the compiler is set as `cuda` (or `cuda-8.0`, `cuda-9.2`, etc.), then
`<NODE_TYPE>` is automatically set to `CUDA`.

<a name="rdc"/>

**`<rdc>`:** The following `<build-name>` keywords determine the value for the
Trilinos CMake cache var `Kokkos_ENABLE_Cuda_Relocatable_Device_Code` in CUDA
builds (does nothing in non-CUDA builds):

* `rdc`: Set `Kokkos_ENABLE_Cuda_Relocatable_Device_Code=ON`
* `no-rdc`: Set `Kokkos_ENABLE_Cuda_Relocatable_Device_Code=OFF`

NOTE: Setting `rdc` also currently adds the `nvcc_wrapper` option
`--remove-duplicate-link-files` as well.

<a name="fpic"/>

**`<fpic>`:** The following `<build-name>` keyword will result in `-fPIC`
being added to `CMAKE_CXX_FLAGS`:

* `fpic`: Add `-fPIC` to `CMAKE_CXX_FLAGS`

<a name="complex"/>

**`<complex>`:** The following `<build-name>` keywords determine if support
for the `complex<double>` scalar type is built into the code and is tested or
not:

* `complex`: Enable support for `complex<double>` scalar type (set
  `Trilinos_ENABLE_COMPLEX=ON`)
* `no-complex`: Do not enable support for `complex<double>` scalar type (set
  `Trilinos_ENABLE_COMPLEX=ON`) (DEFAULT)

NOTE: Setting `Trilinos_ENABLE_COMPLEX=ON` only enables `complex<double>` not
`complex<float>` by default.

<a name="shared_static"/>

**`<shared_static>`:** The following `<build-name>` keywords specify debug if
a shared or static library build of Trilinos is to be created (which also
impacts if shared or stack TPL libs are linked to on some system):

* `static`: `BUILD_SHARED_LIBS=OFF` (DEFAULT)
* `shared`: `BUILD_SHARED_LIBS=ON`

<a name="release_debug"/>

**`<release_debug>`:** The following `<build-name>` keywords specify debug or
optimized build and the `<BUILD_TYPE>` variable (used to set the CMake cache
var `CMAKE_BUILD_TYPE=[DEBUG|RELEASE]` and turn on or off runtime debug
checking `Trilinos_ENABLE_DEBUG=ON`, e.g. array bounds checking, pointer
checking etc.):

* `release-debug` or `opt-dbg` (or using `_`): (`<BUILD_TYPE>=RELEASE-DEBUG`)
  * Set `CMAKE_BULD_TYPE=RELEASE` (i.e. `-O3` compiler options)
  * Turn **ON** runtime debug checking
  * NOTE: This build runs runtime checks to catch developer and user mistakes
    but still runs fairly fast.
* `release` or `opt`: (`<BUILD_TYPE>=RELEASE`)
  * Set `CMAKE_BULD_TYPE=RELEASE` (i.e. `-O3` compiler options)
  * Turn **OFF** runtime debug checking
  * NOTE: This build runs fast with minimal checks (i.e. production).
* `debug` or `dbg`: (`<BUILD_TYPE>=DEBUG`, DEFAULT)
  * Set `CMAKE_BULD_TYPE=DEBUG` (i.e. `-O0 -g` compiler options)
  * Turn **ON** runtime debug checking
  * NOTE: This build supports running in a debugger.

<a name="pt"/>

**`<pt>`:** The following `<build-name>` keywords specify all Primary Tested
(pt) packages are allowed to be enabled.  The default is to allow enable of
Secondary Tested packages with disables for packages that the ATDM APPs do not
use (as specified in the file `ATDMDisables.cmake`):

* `pt`: Allow enable if all PT packages, don't disable any PT packages by
  default and enable Fortran.

All other strings in `<build-name>` are ignored but are allowed for
informational purposes.  The reason that a `<build-name>` string is defined in
this form is that this can be used as the Jenkins job name and the Trilinos
build name that shows up on CDash.  This makes it very easy to define the
configuration options and maintain the Jenkins build jobs.  The combination
`<COMPILER>-<BUILD_TYPE>-<NODE_TYPE>-<KOKKOS_ARCH>` is used to define the
CMake variable `ATDM_BUILD_NAME_KEYS_STR` that is used to uniquely define a
build on a particular system and to manage a system of tweaks for each of the
supported builds (see below).

<a name="build-name-examples"/>

**[build-name-examples]** Some examples of `<build-name>` keyword sets used on
various platforms include:
* `gnu-openmp-debug`
* `gnu-openmp-opt`
* `intel-openmp-debug`
* `intel-openmp-opt`
* `sems-rhel6-gnu-openmp-debug`
* `cee-rhel6-gnu-openmp-debug`
* `sems-rhel6-intel-openmp-opt`
* `cee-rhel6-intel-openmp-opt`
* `cee-rhel6-gnu-7.2.0-openmpi-1.10.2-debug-openmp`
* `intel-debug-openmp-KNL`
* `intel-opt-openmp-HSW`
* `cuda-debug` (`<NODE_TYPE>` is implicitly `CUDA`)
* `cuda-opt` (`<NODE_TYPE>` is implicitly `CUDA`)
* `cuda-8.0-debug` (`<NODE_TYPE>` is implicitly `CUDA`)
* `cuda-9.0-opt-Kepler37` (`<NODE_TYPE>` is implicitly `CUDA`)

The script `cmake/std/atdm/load-env.sh` when sourced sets some bash
environment variables that are prefixed with `ATDM_CONFIG_` and other standard
variables.  This includes setting the var `ATDM_CONFIG_BUILD_NAME` which stores
the input `<build-name>` which is used in other parts of the system.

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

When `ATDMDevEnv.cmake` (or `ATDMDevEnvSettings.cmake`) is being processed, if
there are "tweaks" files defined for a build, then they will be picked up in
the CMake cache var <a href="#ATDM_TWEAKS_FILES">ATDM_TWEAKS_FILES</a> and
those files will be read in using `INCLUDE()` to process the extra options
contained within it.


## Installation and usage

When including the file `ATDMDevEnv.cmake` (or `ATDMDevEnvSettings.cmake`) in
the CMake configuration as described above, the cmake configure automatically
sets up to install a script:

```
  <install-prefix>/load_matching_env.sh
```

which after installation can then be sourced by clients using:

```
$ source <install-prefix>/load_matching_env.sh
```

Sourcing this file loads the compilers, MPI, and TPLs and sets up the various
`ATDM_CONG_` environment variables described above.  It also sets the
environment variable `export Trilinos_ROOT=<install-prefix>` that clients can
use to point back to the Trilinos installation directory.  Also, this variable
will allow `find_package(Trilinos)` to automatically find the right Trilinos
installation with no further action.

The install location `<install-prefix>` can be set using the CMake cache
variable:

```
  -D CMAKE_INSTALL_PREFIX=<install-prefix> \
```

or by setting the environment variable:

```
$ export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX=<install-prefix>
```

when configuring Trilinos.

If the environment variable `ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX` is set,
then it will be used to set `CMAKE_INSTALL_PREFIX` internally and override any
value that might be passed in or set otherwise.  (This is a `FORCE`D cache
variable set on `CMAKE_INSTALL_PREFIX` so this value will appear in the
`CMakeCache.txt` file.)

If permissions need to be set on a base dir of `CMAKE_INSTALL_PREFIX`, then one
can set a base dir of this through:

```
$ export ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR=<install-base-dir>
```

where `<install-base-dir>` must be a base directory of `<install-prefix>`.  If
not explicitly set in this way, then it is assumed to be whatever the set
value is for `CMAKE_INSTALL_PREFIX`.

By default, every file and directory created during the install under
`<install-base-dir>` will be made explicitly group read/write and "other"
readable.  (That can be changed by setting CMake Cache vars starting with
`Trilinos_MAKE_INSTALL_`.)

The owning group for everything under `<install-base-dir>` can be set using:

```
$ export ATDM_CONFIG_MAKE_INSTALL_GROUP=<owning-group>
```

Otherwise, the owning group will be set by the group sticky bit or by the
default user's group (on systems that don't support the group sticky bit).

The name of the installed script `load_matching_env.sh` can be changed at
configure-time using the CMake cache variable:

```
  -D ATDM_INSTALLED_ENV_LOAD_SCRIPT_NAME=<load-matching-env-sh> \
```

where if the CMake cache variable `ATDM_INSTALLED_ENV_LOAD_SCRIPT_NAME` is not
specified, then it is given the name `load_matching_env.sh` by default.


## Parallel build and test processes

By default, each system's `<system_name>/environment.sh` script automatically
selects the parallel build and test jobs by setting the env vars:

* `ATDM_CONFIG_BUILD_COUNT`: Number of default parallel build jobs passed to
  `ninja -j ${ATDM_CONFIG_BUILD_COUNT}`

* `ATDM_CONFIG_PARALLEL_COMPILE_JOBS_LIMIT`: Max number of parallel processes
  allowed to build object files with ninja (default is usually empty).

* `ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT`: Max number of parallel processes
  allowed to link libraries and executables (default is usually empty).

* `ATDM_CONFIG_CTEST_PARALLEL_LEVEL`: Number passed to `ctest -j ${ATDM_CONFIG_CTEST_PARALLEL_LEVEL}`

These values can be overridden by setting the following env vars before running
`source cmake/std/atdm/load-env.sh <build-name>`:

* `ATDM_CONFIG_BUILD_COUNT_OVERRIDE`
* `ATDM_CONFIG_PARALLEL_COMPILE_JOBS_LIMIT_OVERRIDE`
* `ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT_OVERRIDE`
* `ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERRIDE`


## checkin-test-atdm.sh

For those Trilinos developers comfortable with using the
`Trilinos/checkin-test.py` script, multiple local builds and testing on a
system can also be driven with the provided `checkin-test-atdm.sh` wrapper
script.  This can be used to drive a number of builds on system as:

```
$ cd <some_build_dir>/

$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .

$ ./checkin-test-atdm.sh <build-name-0> <build-name-1> ... \
  --enable-packages=<Package> --local-do-all
```

That will configure, build, and run tests for each specified build
`<build-name-0>` and send a summary email when complete.  All of the supported
builds on the local system can be run by using `all` instead of
`<build-name-0> <build-name-1> ...`.  See comments at the top of the script
`checkin-test-atdm.sh` for more details.

The parallel level for building and running tests are determined by the env
vars `ATDM_CONFIG_BUILD_COUNT` and `ATDM_CONFIG_CTEST_PARALLEL_LEVEL`,
respectfully, as set by default for the given system by the `atdm/load-env.sh`
script for the given machine.  (On most machines, these are fixed but on
generic systems like <a href="#sems-rhel6-environment">sems-rhel6</a>, they
are computed from the number of cores on that machine).  These values can be
overridden by setting the env vars `ATDM_CONFIG_BUILD_COUNT_OVERRIDE` and
`ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERIDE`, respectfully as, for example:

```
$ env \
  ATDM_CONFIG_BUILD_COUNT_OVERRIDE=8 \
  ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERIDE=12 \
  ./checkin-test-atdm.sh [options]
```

A value of `ATDM_CONFIG_BUILD_COUNT_OVERRIDE=0` or less than `0` is allowed
when using Ninja (i.e. `ATDM_CONFIG_USE_NINJA=ON`) in which case `ninja` will
be run with non `-j<N>` argument, and therefore all of the non-loaded cores
will be used.

Alternatively, one can override the parallel build and test running levels and
set other make/ninja and ctest options using the checkin-test arguments
`--make-options` and `--ctest-options`.  For example, to use 20 processes to
build with Nina, have Ninja keep going even if there are build errors, and run
ctest with 10 processes, one can use:

```
$ ./checkin-test-atdm.sh \
  --make-options="-j20 -k 99999999" \
  --ctest-options=-j10 \
  [other options]
```

To run tests for a CUDA build or to run tests on platforms that must run on a
compute node one will need to run these on a compute node on the system that
has a GPU.  On such a system one would run:

```
$ ./checkin-test-atdm.sh <build-name-0> <build-name-1> ... \
  --enable-packages=<Package> --configure --build \
  && \
  <command-to-run-on-compute-node> \
  ./checkin-test-atdm.sh <build-name-0> <build-name-1> ... \
  --enable-packages=<Package> --test
```

See <a href="#specific-instructions-for-each-system">Specific instructions for
each system</a> for details for how to run on the compute node for that
system.

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
$ ./checkin-test-atdm.sh <build-name-0> <build-name-1> ... \
  --enable-packages=<Package> --local-do-all
```

However, a default `local-checkin-test-defaults.py` is created the first time
the `checkin-test-atdm.sh` script is run and will set these as the defaults
(after which can be modified).


## ctest-s-local-test-driver.sh

When one wants to run local builds to test a branch and submit results to
CDash (so that they are archived and for others to see), then a simple way to
that is to use the provided
[`ctest-s-local-test-driver.sh`](https://github.com/trilinos/Trilinos/blob/develop/cmake/std/atdm/ctest-s-local-test-driver.sh)
script.  This script uses the CTest -S Jenkins driver system in the directory
`Trilinos/cmake/ctest/drivers/atdm/` and the specific Jenkins driver files in
the directory

```
  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
```

to run builds and submit to the Experimental CDash Track/Group.

To use this, first set up a local directory and symlink as:

```
$ cd <some_base_build_dir>/
$ ln -s <some_base_dir>/Trilinos/cmake/std/atdm/ctest-s-local-test-driver.sh .
````

Then one can run any of the builds with defined driver files listed under:

```
    Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
      <full-build-name-1>.sh
      <full-build-name-2>.sh
      ...
```

using:

```
$ env \
    Trilinos_PACKAGES=<pkg0>,<pkg1>,... \
    ATDM_CTEST_S_USE_FULL_BUILD_NAME=1 \
  ./ctest-s-local-test-driver.sh <build-base-name-0> <build-base-name-1> ...
```

(Or leave out `Trilinos_PACKAGES` to test all of the ATDM packages.)  That
will submit results to the Trilinos CDash project to the "Experimental" Group
(the CDash group cannot be changed).  This will automatically allocate nodes
and run just like it was running as a Jenkins job so the details of how this
is done are completely taken care of by the existing setup for the current
system.

To run all of the supported builds listed in the variable
`ATDM_CONFIG_ALL_SUPPORTED_BUILDS` in the file
`cmake/std/atdm/<system_name>/all_supported_builds.sh`, use:

```
$ env \
    Trilinos_PACKAGES=<pkg0>,<pkg1>,... \
  ./ctest-s-local-test-driver.sh all
```

One can examine the progress of the builds and tests locally by looking at the
generated files:

```
  <some_base_build_dir>/<full_build_name>/smart-jenkins-driver.out
```

and also examine the generated `*.xml` configure, build, and test files
created under:

```
  <some_base_build_dir>/<full_build_name>/SRC_AND_BUILD/BUILD/Testing/
```

To avoid submitting results to CDash and rebuilding (instead of blowing away
the build dir each time), run:

```
$ env \
    Trilinos_PACKAGES=<pkg0>,<pkg1>,... \
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE \
    CTEST_DO_SUBMIT=OFF \
  ./ctest-s-local-test-driver.sh <build-base-name-0> <build-base-name-1> ...
```

One can also do run local installs using:

```
$ env \
    Trilinos_PACKAGES=<pkg0>,<pkg1>,... \
    ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX=install \
    CTEST_DO_INSTALL=ON \
  ./ctest-s-local-test-driver.sh <build-base-name-0> <build-base-name-1> ...
```

That will install the enabled Trilinos packages under:

```
  <full-build-name>/SRC_AND_BUILD/BUILD/install/
```

For more details, see the help documentation in the script itself
[`ctest-s-local-test-driver.sh`](https://github.com/trilinos/Trilinos/blob/develop/cmake/std/atdm/ctest-s-local-test-driver.sh). Also,
see
[TRIBITS_CTEST_DRIVER()](https://tribits.org/doc/TribitsDevelopersGuide.html#determining-what-testing-related-actions-are-performed-tribits-ctest-driver)
for a description of all of the options that can be set as env vars to, for
example, skip the configure, skip the build, skip running tests, etc.


## Specific instructions for each system

* <a href="#ridewhite">ride/white</a>
* <a href="#shillerhansen">shiller/hansen</a>
* <a href="#tlcc-2-and-cts-1">TLCC-2 and CTS-1</a>
* <a href="#mutrino">mutrino</a>
* <a href="#sems-rhel6-environment">SEMS RHEL6 Environment</a>
* <a href="#sems-rhel7-environment">SEMS RHEL7 Environment</a>
* <a href="#spack-rhel-environment">Spack RHEL Environment</a>
* <a href="#cee-rhel6-and-rhel7-environment">CEE RHEL6 and RHEL7 Environment</a>
* <a href="#waterman">waterman</a>
* <a href="#ats-2">ATS-2</a>
* <a href="#astra-vanguard-arm-system">ASTRA (Vanguard ARM System)</a>
* <a href="#ats-1">ATS-1</a>


### ride/white

Once logged on to 'white' (on the SON) or 'ride' (on the SRN), one can
directly configure and build on the login node (being careful not to overload
the node) using the `ride` env.  But to run the tests, one must run on the
compute nodes using the `bsub` command to run if using a CUDA build.  For
example, to configure, build and run the tests for the `cuda-debug` build for
say `MueLu` on 'white', (after cloning Trilinos on the `develop` branch) one
would do:

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cuda-debug

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ bsub -x -Is -q rhel7F -n 16 ctest -j4
```

The ATDM configuration of Trilinos is set up to run on the Firestone nodes
(Dual-Socket POWER8, 8 cores per socket, K80 GPUs).  This configuration will
not work on the other GPU nodes currently.

**NOTE:** While the above example shows loading the environment, configuring
and building on the login node, one can also do these on the compute nodes as
well.  In fact, that is what the CTest -S drivers do in automated testing on
'white' and 'ride'.

Note that one can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ bsub -x -I -q rhel7F -n 16 \
  ./checkin-test-atdm.sh cuda-debug \
  --enable-packages=MueLu \
  --local-do-all
```


### shiller/hansen

Once logged on to 'hansen' (on the SON) or 'shiller' (on the SRN), one can
directly configure and build on the login node (being careful not to overload
the node) using the `shiller` env.  But to run the tests, one must run on the
compute nodes using the `srun` command.  For example, to configure, build and
run the tests for say `MueLu` on 'hansen', (after cloning Trilinos on the
`develop` branch) one would do:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-opt-openmp

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ srun ctest -j4
```

**NOTE:** While the above example shows loading the environment, configuring
and building on the login node, one can also do these on the compute nodes as
well.  In fact, that is what the CTest -S drivers do in automated testing on
'hansen' and 'shiller'.

Note that one can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ srun ./checkin-test-atdm.sh intel-opt-openmp \
  --enable-packages=MueLu \
  --local-do-all
```


### TLCC-2 and CTS-1

Once logged on to any TLCC2 machine (e.g. 'chama', 'skybridge') or the CTS-1
machine 'serrano', 'eclipse' or 'ghost' machines, one can directly configure
and build on the login node (being careful not to overload the node) using the
`chama` and `serrano` envs, respectively.  But to run the tests, one must run
on the compute nodes using the `srun` command.  For example, to configure,
build and run the tests for say `MueLu` on 'serrano', (after cloning Trilinos
on the `develop` branch) one would do:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-opt-openmp

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ salloc -N1 --time=0:20:00 --account=<YOUR_WCID> ctest -j4
```

To get information on `<YOUR_WCID>` used above, there is a WC tool tab on
computing.sandia.gov

**NOTE:** Unlike some of the other machines, one must load the environment,
configure and build on the login node and then run the test suite on a compute
node on this system.  This is what the CTest -S driver does on TLCC-2 and
CTS-1 systems like 'chama' and 'serrano' in order to drive jobs and submit to
CDash.

To use the checkin-test-atdm.sh script, you must split running the tests from
the configure and build as with:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ ./checkin-test-atdm.sh intel-opt-openmp \
  --enable-packages=MueLu --allow-no-pull --configure --build
$ salloc -N1 --time=0:20:00 --account=<YOUR_WCID> \
  ./checkin-test-atdm.sh intel-opt-openmp \
  --enable-packages=MueLu --test
```


### mutrino

The default environment on mutrino is now ats1. Please see 
<a href="#ats-1">ATS-1</a>.
Once logged on to 'mutrino', one can directly configure and build on the login
node (being careful not to overload the node) using the `mutrino` env.  But to
run the tests, one must run on the compute nodes using the `salloc` command.
For example, to configure, build and run the tests for say `MueLu` on
'mutrino', (after cloning Trilinos on the `develop` branch) one would:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh mutrino-intel-opt-openmp

$ cmake \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make -j16

# to run on the Haswell partition
$ salloc -N 1 -p standard -J $ATDM_CONFIG_BUILD_NAME ctest -j16

# to run on the KNL partition
$ salloc -N 1 -p knl -J $ATDM_CONFIG_BUILD_NAME ctest -j16
```

**NOTE:** Unlike some of the other machines, one must load the environment,
configure and build on the login node and then run the test suite on a compute
node on this system.  This is what the CTest -S driver on 'mutrino' does in
order to drive jobs and submit to CDash.


### SEMS RHEL6 Environment

Once logged on to a SNL COE RHEL6 machine with the sems NFS env, one can
directly configure, build, and run tests using the `sems-rhel6` env.  For
example, to configure, build and run the tests for `MueLu` one would clone
Trilinos on the `develop` branch and then do the following:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh sems-rhel6-intel-opt-openmp

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ ctest -j8
```

NOTE:     Above    including     `sems-rhel6`     in     the    build     name
`sems-rhel6-intel-opt-openmp` is  not necessary but  is recommended when  on a
CEE LAN RHEL6 machine  to be explicit that the SEMS env is  being used and not
the <a href="#cee-rhel6-and-rhel7-environment">CEE RHEL6 env</a>.

One can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ ./checkin-test-atdm.sh sems-rhel6-intel-opt-openmp \
  --enable-packages=MueLu \
  --local-do-all
```

NOTE: The number of parallel build and test processes in this case are
determine automatically from the number of cores on the current machine.  But
this can be overridden by setting the env var
`ATDM_CONFIG_NUM_CORES_ON_MACHINE_OVERRIDE` before running `source
cmake/std/atdm/load-env.sh <build_name>`.


### SEMS RHEL7 Environment

Once logged on to a SNL COE RHEL7 machine with the SEMS NFS env, one can
directly configure, build, and run tests using the `sems-rhel7` env.  For
example, to configure, build and run the tests for `MueLu` one would clone
Trilinos on the `develop` branch and then do the following:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cuda-9.2-Pascal60-release-debug

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ ctest -j8
```

One can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ ./checkin-test-atdm.sh cuda-9.2-Pascal60-release-debug \
  --enable-packages=MueLu \
  --local-do-all
```

NOTE: The number of parallel build and test processes in this case are
determine automatically from the number of cores on the current machine.  But
this can be overridden by setting the env var
`ATDM_CONFIG_NUM_CORES_ON_MACHINE_OVERRIDE` running `source
cmake/std/atdm/load-env.sh <build_name>`.

NOTE: The default Intel compiler license server can be overridden by setting
the env var:

```
$ export ATDM_CONFIG_LM_LICENSE_FILE_OVERRIDE=<some-url>
```

NOTE: If the SEMS modules are not already defined in the current shell
(i.e. the environment variable `SEMS_MODULEFILES_ROOT` is empty), then the
`atdm/load-env.sh` script for this 'sems-rhel7' system will attempt to define
the SEMS modules by running:

```
$ module use /projects/sems/modulefiles/projects
$ export SEMS_MODULEFILES_ROOT=/projects/sems/modulefiles
```

if that directory exists.


### Spack RHEL Environment

The env 'spack-rhel' should work on any Red Hat Enterprise Linux (RHEL) (and
perhaps many other Linux systems) that have the SNL ATDM Spack modules
installed on them.  See the [installation
documentation](https://gitlab.sandia.gov/atdm/atdm-spack-scripts/blob/master/README.md).
**WARNING:** This Spack env is still under development and may change in the
future.

Once logged onto a Linux machine with the SNL ATDM Spack modules installed,
one can directly configure, build, and run tests using the `spack-rhel` env.
For example, to configure, build and run the tests for `MueLu` one would clone
Trilinos on the `develop` branch and then do the following:


```
$ cd <some_build_dir>/

$ source <spack-install-base-dir>/setup-env.sh

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh spack-rhel-gnu-openmp-opt

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ ctest -j8
```

One can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ env ATDM_CHT_DEFAULT_ENV=spack-rhel-default \
  ./checkin-test-atdm.sh spack-rhel-gnu-openmp-opt \
  --enable-packages=MueLu \
  --local-do-all
```

NOTE: Above one must set `ATDM_CHT_DEFAULT_ENV=spack-rhel-default` in the env
when passing in `all` in order for it to select the correct set of supported
builds for the `spack-rhel` env and also to load the correct env to find
Python, etc.


### CEE RHEL6 and RHEL7 Environment

Once logged into any CEE LAN RHEL6 or RHEL7 SRN machine, one can configure,
build, and run tests for any ATDM Trilinos package using the `cee-rhel6` env.
For example, to configure, build and run the tests for the
`cee-rhel6-clang-opt-openmp` build for say `MueLu` on a CEE LAN machine,
(after cloning Trilinos on the `develop` branch) one would do:

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cee-rhel6-clang-opt-openmp

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ ctest -j16
```

NOTE: Above one must include `cee-rhel6` in the build name
`cee-rhel6-clang-opt-openmp` in order to select the `cee-rhel6` env on a CEE
LAN RHEL6 machine or the <a href="#sems-rhel6-environment">sems-rhel6</a> env
will be used by default.

One can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ env ATDM_CHT_DEFAULT_ENV=cee-rhel6-default \
  ./checkin-test-atdm.sh cee-rhel6-clang-opt-openmp \
  --enable-packages=MueLu \
  --local-do-all
```

NOTE: Above one must set `ATDM_CHT_DEFAULT_ENV=cee-rhel6-default` in the env
when passing in `all` in order for it to select the correct set of supported
builds for the `cee-rhel6` env.

NOTE: The number of parallel build and test processes in this case are
determine automatically from the number of cores on the current machine.  But
this can be overridden by setting the env var
`ATDM_CONFIG_NUM_CORES_ON_MACHINE_OVERRIDE` before running `source
cmake/std/atdm/load-env.sh <build_name>`.


### waterman

Once logged on to 'waterman' (SRN), one can directly configure and build on
the login node (being careful not to overload the node) using the `waterman`
env.  But to run the tests, one must run on the compute nodes using the `bsub`
command to run if using a CUDA build.  For example, to configure, build and
run the tests for the default `cuda-debug` build for say `MueLu` (after
cloning Trilinos on the `develop` branch) one would do:

```bash
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cuda-debug

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=20

$ bsub -x -Is -n 20 ctest -j2
```

**NOTE:** While the above example shows loading the environment, configuring
and building on the login node, one can also do these on the compute nodes as
well.  In fact, that is what the CTest -S drivers do in automated testing on
'waterman'.  To get an interactive compute node, do:

```
$ bsub -x -Is -n 20 bash
```

Then one can configure, build, and run tests interactively on that compute
node.

Note that one can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ bsub -x -Is -n 20 \
  ./checkin-test-atdm.sh cuda-debug --enable-packages=MueLu --local-do-all
```


### ATS-2

Once logged on a supported ATS-2 system (called system 'ats2') like 'vortex'
(SRN), one can either build and configure on a login node or a compute node.
But one must always run the tests from the launch node (allocated using
'bsub').  Make sure to setup SSH keys as described in `/opt/VORTEX_INTRO`
before trying to do anything.

For example, to configure, build and run the tests for the default
`cuda-debug` build for `Tpetra` (after cloning Trilinos on the 'develop'
branch), run the following from the login node on 'vortex':

```bash
$ cd <some_build_dir>/

# Load env and configure on the login node
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cuda-debug
$ cmake -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON \
  -DTrilinos_ENABLE_Kokkos=ON \
  $TRILINOS_DIR

# Get a node allocation and log on to launch node
$ bsub -J <job-name> -W 6:00 -Is bash

# Build on the allocated compute node
$ lrun -n 1 make NP=32

# Run the test suite from the launch node
$ ctest -j4
```

CTest runs everything using the `jsrun` command.  You must run jsrun from the
launch node allocated using `bsub`.  That is why the raw `ctest` command is
run on the launch node.

One can also build directly from the login node on a compute node using:

```
$ lalloc 1 -W 4:00 make NP=32
```

Or, to SSH from the launch node (acquired with `bsub -Is bash` as shown above)
to the compute node and run make from there, use:

```bash
$ ssh $(atdm_ats2_get_allocated_compute_node_name)
$ make NP=32
```

One can also directly run the test suite from the login node using:

```bash
$ bsub -Is -W 2:00 ctest -j4
```

The advantage of creating a node allocation first using `bsub -Is bash` and
logging into the launch node is that your builds and test runs from there are
completely interactive since you don't need to wait for a node allocation.
(Also, you can set an LSF job name with `bsub -J <job-name>` which you can't
do with `lalloc`.)

The MPI test tests in Trilinos are actually run by a wrapper script
`trilinos_jsrun` which calls the `jsrun` command and modifies the input
arguments to accommodate the MPI test suite in Trilinos (see the
implementation of the script `trilinos_jsrun` for details on what it does).
By default, the script `trilinos_jsrun` will set `export
TPETRA_ASSUME_CUDA_AWARE_MPI=0` if `TPETRA_ASSUME_CUDA_AWARE_MPI` is unset in
the environment.  Therefore, by default, the tests are run without CUDA-aware
MPI on this system.

To explicitly **disable CUDA-aware MPI** when running the test suite from the
launch node, set the environment variable and run as:

```bash
$ export TPETRA_ASSUME_CUDA_AWARE_MPI=0
$ lrun -n 1 ctest -j4
```

To explicitly **enable CUDA-aware MPI** when running the test suite from the
launch node, set the environment variable and run as:

```bash
$ export TPETRA_ASSUME_CUDA_AWARE_MPI=1
$ lrun -n 1 ctest -j4
```

Note that one must run the function <a
href="#ctest-s-local-test-driversh">ctest-s-local-test-driver.sh</a> from the
login node, not the launch node or the compute node!  It is designed to run
from the login node and it will allocate resources to run on a compute node
when it runs.

One can also use the <a href="#checkin-test-atdmsh">checkin-test-atdm.sh</a>
tool to drive local builds.  To do so, run from the launch node as:

```bash
# Allocate a compute node and log onto the launch node
$ bsub -J <job-name> -W 6:00 -Is bash

# Do the configure and build from the allocated compute node
$ lrun -n 1 ./checkin-test-atdm.sh <buildname0> <buildname1> ... \
    --enable-packages=<pkg0>,<pkg1>,... --configure --build

# Run tests from the launch node which will call jsrun
$ ./checkin-test-atdm.sh <buildname0> <buildname1> ... \
    --enable-packages=<pkg0>,<pkg1>,... --test
```

**NOTES:**
- Do **NOT** do `module purge` before loading the environment. Simply start
  off with a clean default environment (through a fresh login) on 'vortex'.
- Please do not build on the launch node after logging onto it using `bsub -Is
  bash`.  That takes up CPU resources on the launch node that need to be used
  by all users of the cluster to run individual MPI jobs.
- The `lrun` command is just a simpler wrapper around `jsrun` that is meant to
  be more similar to the SLURM `srun` command.  That is, `lrun -n <N> <cmnd>`
  should behave similarly to `jsrun --np <N> <cmnd>`.
- One can directly log into any compute node (independent of if that node is
  allocated to you or not).  That can be useful to examine a running job but
  use caution when doing so as not to disturb the job running.


### ASTRA (Vanguard ARM System)

Once logged onto a supported Vanguard ARM system (called system 'van1-tx2')
system like 'stria', one can build and configure on a login node.

To configure, build and run the tests for the default `arm-20.0` build for
`Kokkos` (after cloning Trilinos on the 'develop' branch), run the following
from a login node on 'stria':

```bash
$ cd <some_build_dir>/

# List available environments
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh help

# Load arm env and configure on the login node
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh arm-20.0
$ cmake -G Ninja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON \
  -DTrilinos_ENABLE_Tpetra=ON \
  $TRILINOS_DIR

$ make NP=8  # This is a shared node!

# Get a node allocation and run ctest
$ salloc -N 1 --time=2:00:00  -p short,batch --account=<YOUR_WCID> ctest -j4
```

One can also get an allocation first and then configure, build on a compute
node, and then run the test suite using:

```bash
$ salloc -N 1 --time=4:00:00 -p short,batch --account=<YOUR_WCID> bash
# NOTE: After the above runs, hostname=stria-login<n> but now a compute node
# has been allocated for immediately usage.

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh arm-20.0

$ cmake -G Ninja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON \
  -DTrilinos_ENABLE_Tpetra=ON \
  $TRILINOS_DIR

$ srun -N 1 make NP=20  # We have the entire compute node to ourselves!

$ ctest -j4
```

The advantage of the latter approach is that one just waits once for a node
allocation and then one can immediately run fast parallel builds on the compute
node (taking up the entire node).  Then one can run the test suite multiple
times without waiting for a new allocation.

One can also directly build on a compute node from the login node in one
command:

```bash
$ srun -N 1 --time=2:00:00 -p short,batch --account=<YOUR_WCID> make NP=20
```

To use the `ctest-s-local-test-driver.sh` script, one must set one's WCID
account using:

```
$ export ATDM_CONFIG_WCID_ACCOUNT=<YOUR_WCID>
```

If `ATDM_CONFIG_WCID_ACCOUNT` is not set, then a default account will be used.
(But if the user is not approved for that account, then the allocation will
fail.)

**NOTES:**
- To get information on <YOUR_WCID> used above, there is a WC tool tab on
  computing.sandia.gov.
- CTest runs everything using the `mpirun` command and this must be run from
  inside a `salloc` or `sbatch` allocation and can **not** be directly
  launched from a compute node.  For example, one cannot get an interactive
  shell directly on a compute node using `srun ... bash' an` then run `mpirun`
  from there.

### ATS-1

Once logged onto a supported ATS-1 system (called system 'ats1')
system like 'mutrino', one can build and configure on a login node.

To configure, build and run the tests for the default `intel-19` build for
`Kokkos` (after cloning Trilinos on the 'develop' branch), run the following
from a login node on 'mutrino':

```bash
$ cd <some_build_dir>/

# List available environments
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh ats1-help

# Load HSW env and configure on the login node
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-19
$ cmake -G Ninja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON \
  -DTrilinos_ENABLE_Tpetra=ON \
  $TRILINOS_DIR

$ make NP=8  # This is a shared node!

# Get a node allocation and run ctest
$ salloc -N 1 --time=2:00:00  -p short,batch --account=<YOUR_WCID> ctest -j1
```

One can also run tests in the background on a compute node using sbatch:
```
echo '#!/bin/bash' > do-test.sh
echo 'ctest -j1' >> do-test.sh
$ sbatch --output=do-test.out -N1 --time=2:00:00 -J ats1-do-test --account=<YOUR_WCID> do-test.sh
```

To drop to a compute node for an interactive session:
```
$ salloc -N1 --account=<YOUR_WCID>
```

To use the `ctest-s-local-test-driver.sh` script, one must set one's WCID
account using:

```
$ export ATDM_CONFIG_WCID_ACCOUNT=<YOUR_WCID>
```

If `ATDM_CONFIG_WCID_ACCOUNT` is not set, then a default account will be used.
(But if the user is not approved for that account, then the allocation will
fail.)

**NOTES:**
- **IMPORTANT** Do not use anything but `ctest -j1` for running tests.
- One cannot run more than one instance of srun at a time on a given compute node.
- To get information on <YOUR_WCID> used above, there is a WC tool tab on
  computing.sandia.gov.
- CTest runs everything using the `srun` command and this must be run from
  inside a `salloc` or `sbatch` allocation.

## Building and installing Trilinos for ATDM Applications

See the following internal SNL wiki page for instructions on building and
testing that ATDM APPs (e.g. SPARC and EMPIRE) against ATDM Trilinos builds:

* [Building ATDM APPs against Trilinos](https://snl-wiki.sandia.gov/display/CoodinatedDevOpsATDM/Building+ATDM+APPs+against+Trilinos)


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

Another problem occurs for developers who don't use bash but instead need to
switch into a bash shell.  In this case, one must use `bash -l` so that the
user's `.bash_profile` file will get sourced.  Without this, one can't load
the modules when sourcing `Trilinos/cmake/std/atdm/load-env.sh`.


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
  subpackages not used by ATDM application customers, disables of test suites
  and individual tests across various builds, etc..  This file gets included
  automatically in `ATDMDevEnv.cmake` (so you don't need to list it in local
  configures of Trilinos).  But this file is also included in the outer `ctest
  -S` driver script code for ATDM builds of Trilinos which is needed for
  correct package-by-package testing of Trilinos.

* **ATDMDevEnv.cmake**: Includes `ATDMDevEnvSettings.cmake` and
  `ATDMDisables.cmake`.

* **checkin-test-atdm.sh**: Uses the `Trilinos/checkin-test.py` script to
  drive builds and tests on the given platform.  (See comments in the top of
  the script for instructions.)

* **ctest-s-local-test-driver.sh**: Uses the script
  `Trilinos/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh` script to drive
  builds and tests on the given platform and submit results to CDash.  (See
  comments in the top of the script for instructions.)

Each supported ATDM system `<system-name>` has its own sub-directory with the
contents:

```
  <system-name>/
    environment.sh  # Load env for the given system based on $ATDM_CONFIG_BUILD_NAME keys
    all_supported_builds.sh  # [Optional] List of all supported builds
    custom_builds.sh  # [Optional] Special logic for compiler keywords, etc.
    tweaks/
       <COMPILER0>_<BUILD_TYPE0>_<NODE_TYPE0>_<KOKKOS_ARCH0>.cmake  # [Optional]
       <COMPILER1>_<BUILD_TYPE1>_<NODE_TYPE1>_<KOKKOS_ARCH0>.cmake  # [Optional]
       ...
       Tweaks.cmake # [Optional]
```

The optional file `<system-name>/all_supported_builds.sh` contains a list of
all of the supported builds on the system.  This sets the bash environment
array variable `ATDM_CONFIG_ALL_SUPPORTED_BUILDS` and the bash var
`ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX` as, for example:

```
  export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-<system_name>-

  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
    gnu-debug-openmp
    gnu-opt-openmp
    ...
    )
```

The variable `ATDM_CONFIG_ALL_SUPPORTED_BUILDS` is used in the
`checkin-test-atdm.sh` script for the `all` argument to run all of the builds
for a system with `checkin-test-atdm.sh all [other options]`.  Both the
variables `ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX` and
`ATDM_CONFIG_ALL_SUPPORTED_BUILDS` are used in the
`ctest-s-local-test-driver.sh` script in order to drive ctest -S Experimental
builds that submit to CDash.

The optional file `<system-name>/custom_builds.sh` contains specialized logic
for compiler versions and other specialized keywords and versions.  (For an
example, see `atdm/cee-rhel6/custom-builds.sh` and
`atdm/cee-rhel6/environment.sh`.)

<a name="ATDM_TWEAKS_FILES"/>

The **ATDM TWEAKS FILES** in the `cmake/std/atdm/<system-name>/tweaks/`
directory contain special settings for specific builds for a specific system.
Typically, these files contains (temporary) disables for tests and test
exectuables for that given build.  When a configure is performed, the internal
CMake variable `ATDM_BUILD_NAME_KEYS_STR` set to
`<COMPILER>_<BUILD_TYPE>_<NODE_TYPE>_<KOKKOS_ARCH>` (printed to STDOUT) is
used to define a default file name:

```
  Trilinos/cmake/std/atdm/<system-name>/tweaks/${ATDM_BUILD_NAME_KEYS_STR}.cmake
```

If that file exists, then it is set as the default for the cmake cache
variable `ATDM_TWEAKS_FILES` (prints to STDOUT) and that file is included and
its options are set as CMake cache variables.  For example, this is what the
output looks like for a build on 'waterman':

```
-- Reading in configuration options from cmake/std/atdm/ATDMDevEnv.cmake ...
-- ATDM_BUILD_NAME_KEYS_STR='CUDA-9.2_RELEASE-DEBUG_CUDA_POWER9_VOLTA70'
-- ATDM_TWEAKS_FILES='<...>/cmake/std/atdm/waterman/tweaks/CUDA-9.2_RELEASE-DEBUG_CUDA_POWER9_VOLTA70.cmake'
-- Including ATDM build tweaks file <...>/cmake/std/atdm/waterman/tweaks/CUDA-9.2_RELEASE-DEBUG_CUDA_POWER9_VOLTA70.cmake ...
```

In addition, if the file:

```
  Trilinos/cmake/std/atdm/<system-name>/tweaks/Tweaks.cmake
```

exists, then it will be included after the above
`tweaks/${ATDM_BUILD_NAME_KEYS_STR}.cmake` file for the matching build.
Disables for all builds on a system or for many related builds on a system can
go into the `Tweaks.cmake` file to avoid having to duplicate disables across
multiple `${ATDM_BUILD_NAME_KEYS_STR}.cmake` files.  Details are in the next
section.


## Disabling failing tests

There are situations where specific tests must be disabled on certain
platforms for certain builds or based on other criteria (see sub-process
[Temporarily disable the failing code or
test](https://snl-wiki.sandia.gov/display/CoodinatedDevOpsATDM/Triaging+and+addressing+ATDM+Trilinos+Failures#TriagingandaddressingATDMTrilinosFailures-5.Makesuretheissueisaddressedinatimelyway:)).
There are various ways to disable tests with the Trilinos TriBITS/CMake-based
build and test system.  First, tests can be disabled in the `CMakeLists.txt`
files that define the tests themselves using various logic.  But the way to
selectively disable tests for the ATDM Trilinos builds that will be described
here will be to only modify files under the `Trilinos/cmake/std/atdm/`
directory.  This will be done by setting the CMake cache variable
`<full-test-name>_DISABLE=ON` for each test to be disabled where
`<full-test-name>` is the full test name.  For example, to disable the test
`MueLu_UnitTestsTpetra_MPI_4`, one would set the CMake cache variable
`MueLu_UnitTestsTpetra_MPI_4_DISABLE=ON`.  The TriBITS/CMake configure system
will then disable the test and print the following line during configuration
(when `Trilinos_TRACE_ADD_TEST=ON` is also set):

```
-- MueLu_UnitTestsTpetra_MPI_4: NOT added test because MueLu_UnitTestsTpetra_MPI_4_DISABLE='ON'!
```

The CMake function `ATDM_SET_ENABLE()` is provided in order to set the disable
such as with:

```
ATDM_SET_ENABLE(MueLu_UnitTestsTpetra_MPI_4_DISABLE ON)
```

Using this function will result in the CMake cache variable
`<full-test-name>_DISABLE` to get set to `ON` by default and will print the
disable to the CMake STDOUT when the `ATDM_SET_ENABLE()` is run, such as:

```
-- Setting default MueLu_UnitTestsTpetra_MPI_4_DISABLE=ON
```

This approach also allows a user to enable the test by adding
`-D<full-test-name>_DISABLE=OFF` to the configure line (which overrides the
default `ON`).

However, before any tests are disabled, there must be a Trilinos GitHub issue
addressing the failing test(s) and the GitHub issue number (e.g. `#2839`) must
put into a comment above the `ATDM_SET_ENABLE()` statement and the issue
number must be added to the git commit message (see [Create a Trilinos GitHub
issue for the
failure](https://snl-wiki.sandia.gov/display/CoodinatedDevOpsATDM/Triaging+and+addressing+ATDM+Trilinos+Failures#TriagingandaddressingATDMTrilinosFailures-3.CreateaTrilinosGitHubissueforthefailure:)).

The best file to add the `ATDM_SET_ENABLE(<full-test-name>_DISABLE ON)`
statement to based on different situations is described in the following
sub-sections:

* <a href="#disable-a-test-for-a-specific-build-on-a-specific-platform">Disable a test for a specific build on a specific platform</a>
* <a href="#disable-a-test-for-several-or-all-builds-on-a-specific-platform">Disable a test for several or all builds on a specific platform</a>
* <a href="#disable-a-test-for-builds-on-all-platforms">Disable a test for builds on all platforms</a>


### Disable a test for a specific build on a specific platform

In order to disable a test for a single build on single platform, one will
want to put the `ATDM_SET_ENABLE()` statement into the [tweaks
file](#ATDM_TWEAKS_FILES) for that build and platform:

```
  Trilinos/cmake/std/atdm/<system-name>/tweaks/<ATDM_BUILD_NAME_KEYS_STR>.cmake
```

The tweak file being looked for is printed out in the CMake configure output
as the line:

```
-- ATDM_TWEAKS_FILES='.../Trilinos/cmake/std/atdm/<system-name>/tweaks/<ATDM_BUILD_NAME_KEYS_STR>.cmake'
```

For example, for the `intel-debug-openmp-KNL` build on 'mutrino', the printout
would be:

```
-- ATDM_TWEAKS_FILES='.../Trilinos/cmake/std/atdm/mutrino/tweaks/INTEL_DEBUG_OPENMP_KNL.cmake'
```

For example, Trilinos commit
[73ae19cf0c](https://github.com/trilinos/Trilinos/commit/73ae19cf0cc7295a7f36de342ea51718226825b7)
shows the disable of the test:

```
# Disable test that times out for some unknown reason (#2925)
ATDM_SET_ENABLE(Stratimikos_test_aztecoo_thyra_driver_MPI_1_DISABLE ON)
```

in both the files `cmake/std/atdm/shiller/tweaks/GNU_DEBUG_SERIAL.cmake` and
`cmake/std/atdm/shiller/tweaks/GNU_RELEASE_SERIAL.cmake` (before `-HSW` was
added to the names).

NOTE: Adding a comment with the Trilinos GitHub Issue ID (`#2925` in this
example) is critical for tractability to remind why the test was disabled and
where to find the disable to remove it later.

To avoid duplicate `ATDM_SET_ENABLE()` statements, one can use the approach of
creating a single `*.cmake` that is included in the different tweak files as
described below.


### Disable a test for several or all builds on a specific platform

It is often the case that a test needs to be disabled for several (or all)
builds for a given platform.  The best way to do this is to disable these in
the file:

```
cmake/std/atdm/<system_name>/tweaks/Tweaks.cmake
```

If that file exists, it will get included after the
`<ATDM_BUILD_NAME_KEYS_STR>.cmake` file as described above.

Typical logic in a `Tweaks.cmake` file may look like:

```
# Disable tests for all builds for this system
ATDM_SET_ENABLE(<full_test_name_1>_DISABLE ON)
ATDM_SET_ENABLE(<full_test_name_2>_DISABLE ON)
...

IF (Trilinos_ENABLE_DEBUG)
  # Disable tests for all debug builds on this system
  ...
ENDIF()

IF (ATDM_COMPILER STREQUAL "GNU-7.2.0")  # See <system_name>/environment.sh
  # Disables for all non-CUDA GNU 7.2.0 builds
  ...
ENDIF()

IF (ATDM_NODE_TYPE STREQUAL "SERIAL")
  # Disable tests for all SERIAL builds for this system
  ...
ELSEIF (ATDM_NODE_TYPE STREQUAL "OPENMP")
  # Disable tests for all OpenMP builds for this system
  ...
ELEIF (ATDM_NODE_TYPE STREQUAL "CUDA")
  # Disable tests for all CUDA builds for this system
  ...
  IF (Trilinos_ENABLE_DEBUG)
    # Disable tests for all CUDA debug builds for this system
    ...
  ENDIF()
  IF (ATDM_CUDA_RDC and Trilinos_ENABLE_DEBUG)
    # Disable tests for all CUDA, RDC, debug builds for this system
    ...
  ENDIF()
ENDIF()
```

Any CMake variable that has been set in the `ATDMDevEnvSettings.cmake` file
before these tweak files are included can be used in if-logic but the
recommended variables are:

* `ATDM_COMPILER` (uppercase)
* `ATDM_KOKKOS_ARCH_JOB_NAME_KEYS` (uppercase and separated by `_`)
* `ATDM_NODE_TYPE` (values `CUDA`, `OPENMP`, `SERIAL`)
* `ATDM_KOKKOS_ARCH` (uppercase and separated by `,`)
* `ATDM_CUDA_RDC` (`ON`/`OFF`)
* `ATDM_FPIC` (`ON`/`OFF`)
* `ATDM_COMPLEX` (`ON`/`OFF`),
* `ATDM_SHARED_LIBS` (`ON`/`OFF`)
* `ATDM_CMAKE_BUILD_TYPE` (values `DEBUG`, `RELEASE` and `RELEASE-DEBUG`)
* `Trilinos_ENABLE_DEBUG` (`ON`/`OFF`)
* `ATDM_PT_PACKAGES` (`ON`/`OFF`)

No other variables should be used in if-logic in these files as other
variables may change in the future.


### Disable a test for builds on all platforms

There are rare cases where a test (or a set of tests) will need to be disabled
for all or a subset of builds across all systems.  To do that using the
above-described approaches would require touching files in every
`cmake/std/atdm/<system-name>/tweaks/` directory.  In rare cases like this,
one can accomplish this by adding a (conditional) `ATDM_SET_ENABLE()`
statement for each test disable directly to the file:

```
  Trilinos/cmake/std/atdm/ATDMDisables.cmake
```

For example, Trilinos commit [5e52db03ff](https://github.com/trilinos/Trilinos/commit/5e52db03ff33acb5b9a0be7ba7507a8bb0de6e30) added the CMake code:

```
# Disable test that fails for all openmp builds (#3035)
IF (ATDM_NODE_TYPE STREQUAL "OPENMP")
  ATDM_SET_ENABLE(MueLu_UnitTestsTpetra_MPI_4_DISABLE ON)
ENDIF()
```

to the file `ATDMDisables.cmake` to disable the test
`MueLu_UnitTestsTpetra_MPI_4` for all OpenMP builds across all platforms.
(Note that that disable was later removed in Trilinos commit
[62fa6663a6](https://github.com/trilinos/Trilinos/commit/62fa6663a6d5a757d786ac87752c3e2074d28414)
after the test was fixed.)


## Specific systems supported

The specific `cmake/std/atdm/<system-name>/` sub-directories and the systems
they support are:

* `cee-rhel6/`: CEE LAN RHEL6 systems with a CEE environment

* `mutrino/`: Supports SNL HPC machine 'mutrino'.

* `ride/`: Supports GNU and CUDA builds on both the SRN machine 'ride' and the
  mirror SON machine 'white'.

* `sems-rhel6/`: SNL COE RHEL6 systems with the SEMS NFS environment

* `sems-rhel7/`: SNL COE RHEL7 systems with the SEMS NFS environment

* `spack-rhel/`: RHEL (and likely other Linux) systems with the SNL ATDM Spack modules installed.

* `cts1/`: Supports SNL HPC CTS-1 machines 'serrano', 'eclipse', and
  'ghost'.

* `shiller/`: Supports GNU, Intel, and CUDA builds on both the SRN machine
  'shiller' and the mirror SON machine 'hansen'.

* `tlcc2/`: Supports SNL HPC TLCC-2 machines 'chama', 'skybridge', etc..

* `waterman/`: Supports GNU and CUDA builds on the SRN machine 'waterman'.

* `ats2/`: Supports GNU and CUDA builds on the SRN machine 'vortex'.


## Custom systems and configurations

In addition to officially defined system configurations described
[above](#specific-instructions-for-each-system), one can also define a custom
system configuration and use that.  To do so, create a new directory with the
contents:

```
<some-base-dir>/<custom-system-name>/
  environment.sh           # Required
  custom_builds.sh         # Optional
  all_supported_builds.sh  # Optional
```

(only the `environment.sh` script is required).

Then one can explicitly use that custom system setup using:

```
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh <build-name> \
  <some-base-dir>/<custom-system-name>
```

The name of the new custom system name is taken from the last directory name
`<custom-system-name>` and the files in that directory are read and treated
like any of the offically defined system configurations.  Also, if the name of
an officially supported system is present in `<build-name>`, it will be
ignored.

A second approach is to register a custom system configuration as:

```
export ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR=<some-base-dir>/<custom-system-name>
```

and then that configuration can be selected by adding the keyword
`<custom-system-name>` to the `<build-name>` passed in as:

```
$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh \
  <custom-system-name>-<other-keywords>
```

In this case, if `<custom-system-name>` is found in the input `<build-name>`
argument, then the custom configuration will be used even if one of the
officially defined configurations would otherwise match.  But if
`<custom-system-name>` is not included in the build name, then the registered
custom system configuration will not be selected.

When a custom configuration is selected, when Trilinos is installed, the
directory `<some-base-dir>/<custom-system-name>` is copied to the install tree
and the installed [`load-matching-env.sh`](#installation-and-usage) script
will automatically load that custom env according to the custom build
configuration script.

A simple example can be seen in:

* [cmake/std/atdm/examples/my-sems-rhel6-config/environment.sh](examples/my-sems-rhel6-config/environment.sh)

which works on any SEMS RHEL6 machine similarly to the offically defined <a
href="#sems-rhel6-environment">SEMS RHEL6 Environment</a>.

To see how things can be specified look at examples of other `environment.sh`
scripts in the offically defined configurations under:

* [cmake/std/atdm/<system_name>/](.)

where `<system_name>` is `ride`, `waterman`, `tlcc2`, etc.
