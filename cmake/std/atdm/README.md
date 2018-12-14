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
* <a href="#installation-and-usage">Installation and usage</a>
* <a href="#checkin-test-atdmsh">checkin-test-atdm.sh</a>
* <a href="#ctest-s-local-test-driversh">ctest-s-local-test-driver.sh</a>
* <a href="#specific-instructions-for-each-system">Specific instructions for each system</a>
* <a href="#building-and-installing-trilinos-for-atdm-applications">Building and installing Trilinos for ATDM Applications</a>
* <a href="#troubleshooting-configuration-problems">Troubleshooting configuration problems</a>
* <a href="#directory-structure-and-contents">Directory structure and contents</a>
* <a href="#disabling-failing-tests">Disabling failing tests</a>
* <a href="#specific-systems-supported">Specific systems supported</a>


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

$ ctest -j16  # Might need to be run with srun or some other command, see below
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
the form `XXX-<keyword0>-<keyword1>-...-YYY` (or
`XXX_<keyword0>_<keyword1>_..._YYY`, either seprator is supported) .  The
typical order and format of this string is:

    <system_name>-<kokkos_arch>-<compiler>-<kokkos_thread>-<shared_static>-<release_debug>

(but almost any order is supported).  All of these keywords, except for
`<compiler>` (which can be `default`), are optional.  All of the other
keywords have reasonable defaults for a given system.  See some examples of
build name strings [below](#build-name-examples).

Each of these keywords [`<system_name>`](#system_name),
[`<kokkos_arch>`](#kokkos_arch), [`<compiler>`](#compiler),
[`<kokkos_thread>`](#kokkos_thread), [`<shared_static>`](#shared_static)
and [`<release_debug>`](#release_debug), is described below.

<a name="system_name"/>

**`<system_name>`**: Typically, the system name is determined automatically by
examining the `hostname` or other files on the system and matching to known
hosts.  Therefore, it is typically not necessary to specify `<system_name>` in
the `<build-name>` keys string.  But there are some cases where more then one
`<system_name>` env are supported on the same machine.  For example, on CEE
LAN RHEL6 machines, both the <a href="#sems-rhel6-environment">SEMS RHEL6
env</a> and <a href="#cee-rhel6-environment">CEE RHEL6 env</a> are supported.
On these CEE LAN RHEL6 machines, when `cee-rhel6` is included in
`<build-name>`, then the `cee-rhel6` env will be selected.  But if
`sems-rhel6` is included in the build name or no system name is given, then
the `sems-rhel6` env will be selected by default on such machines.

<a name="kokkos_arch"/>

**`<kokkos_arch>`:** The `<build-name>` string can also contain keywords to
determine the `KOKKOS_ARCH` option of the build.  This is the case-sensitive
architecture name that is recognized by the CMake
[KOKKOS_ARCH](https://trilinos.org/docs/files/TrilinosBuildReference.html#configuring-with-kokkos-and-advanced-back-ends)
configure option for Trilinos and Kokkos.  Some common supported Kokkos
architectures for the host node include `BDW`, `HSW`, `Power8`, `Power9`, and
`KNL`.  When a GPU is present, some common Kokkos architecture options include
`Kepler37` and `Pascal60`.  If one selects a `KOKKOS_ARCH` value that is not
supported by the current system or selected compiler, then the `load-env.sh`
script will return an error message listing the value choices for
`KOKKOS_ARCH` for each supported compiler.

Note that currently only a single `KOKKOS_ARCH` value is recognized in the
`<build-name>` string and it must be proceeded a dash '-' such as with
`intel-KNL` or `cuda-Kepler37`.  This setup does not currently support
specifying multiple `KOKKOS_ARCH` values (since there is no example yet where
that would be needed or useful) but such functionality could be supported in
the future if needed.

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

<a name="kokkos_thread"/>

**`<kokkos_thread>`:** The following `<build-name>` keywords determine the
Kokkos threading / backend model variable `<NODE_TYPE>` (default is
`<NODE_TYPE>=SERIAL` unless `<COMPILER>=CUDA`):

* `serial`: Use no host threading (`NODE_TYPE=SERIAL`, DEFAULT)
* `pthread`: Use Pthreads for host threading (`NODE_TYPE=THREAD`)
* `openmp`: Use OpenMP for host threading (`NODE_TYPE=OPENMP`)

If `cuda` (or `cuda-8.0`, `cuda-9.2`, etc.) is given, then `<NODE_TYPE>` is
automatically set to `CUDA`.

<a name="shared_static"/>

**`<shared_static>`:** The following `<build-name>` keywords specify debug if a shared or static library build of Trilinos is to be created (which also impacts if shared or stack TPL libs are linked to on some system):

* `static`: `BUILD_SHARED_LIBS=OFF`, DEFAULT
* `shared`: `BUILD_SHARED_LIBS=ON`

<a name="release_debug"/>

**`<release_debug>`:** The following `<build-name>` keywords specify debug or
optimized build and the `<BUILD_TYPE> variable `(used to set the CMake cache
var `CMAKE_BUILD_TYPE=[DEBUG|RELEASE]` and turn on or off runtime debug
checking (e.g. array bounds checking, pointer checking etc.)):

* `release-debug` or `opt-dbg` (or using `_`): (`<BUILD_TYPE>=RELEASE_RELEASE`)
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

All other strings in `<build-name>` are ignored but are allowed for
informational purposes.  The reason that a `<build-name>` string is defined in
this form is that this can be used as the Jenkins job name and the Trilinos
build name that shows up on CDash.  This makes it very easy to define the
configuration options and maintain the Jenkins build jobs.  The combination
`<COMPILER>-<BUILD_TYPE>-<NODE_TYPE>-<KOKKOS_ACRCH>` is used to define the
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

When `ATDMDevEnv.cmake` is being processed, if there is a "tweaks" file
defined for a build, then it will be picked up in the CMake cache var <a
href="#ATDM_TWEAKS_FILES">ATDM_TWEAKS_FILES</a> and that file will be read in
using `INCLUDE()` to process the extra options contained within it.


## Installation and usage

When including the `ATDMDevEnv.cmake` file (or `ATDMDevEnvSettings.cmake`) at
configure time as described above, the cmake configure automatically sets up
to install an environment script:

```
  <install-prefix>/<load-matching-env-sh>
```

where `<install-prefix>` and `<load-matching-env-sh>` are set at
configure-time using:

```
  -D CMAKE_INSTALL_PREFIX=<install-prefix> \
  -D ATDM_INSTALLED_ENV_LOAD_SCRIPT_NAME=<load-matching-env-sh> \
  -D ATDM_TRILINOS_INSTALL_PREFIX_ENV_VAR_NAME=<trilinos-install-prefix-var-name> \
```

* If `ATDM_INSTALLED_ENV_LOAD_SCRIPT_NAME` is not specified then it is given the
name `load_matching_env.sh` by default.

* If `ATDM_TRILINOS_INSTALL_PREFIX_ENV_VAR_NAME` is not specified then it is
given the name `ATDM_TRILINOS_INSTALL_PREFIX` by default.

After installation with `make install`, a client can load the environment to
use this ATDM configuration of Trilinos by running:

```
$ source <install-prefix>/<load-matching-env-sh>
```

Sourcing this file sets all of the various `ATDM_CONG_` environment variables
described above and also sets the environment variable
`<trilinos-install-prefix-var-name>` to `<install-prefix>`.


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
CDash so that they are archived and for others to see, then a simple way to
that is to use the provided `ctest-s-local-test-driver.sh` script.  This
script uses the CTest -S Jenkins driver system in the directory
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

Then run any of the build names (e.g. `gnu-opt-debug`) listed in the variable
`ATDM_CONFIG_ALL_SUPPORTED_BUILDS` in the file
`cmake/std/atdm/<system_name>/all_supported_builds.sh` (or `all` for all of
the defined builds) for the system as:

```
$ env \
    Trilinos_PACKAGES=<pkg0>,<pkg1>,... \
  ./ctest-s-local-test-driver.sh <build-base-name-0> <build-base-name-1> ...
```

That will submit results to the Trilinos CDash project to the "Experimental"
CDash Group (the CDash group can not be changed).  This will automatically
allocate nodes and run just like it was running as a Jenkins job so the
details of how this is done are completely taken care of by the existing setup
for the current system.

One can examine the progress of the builds and tests locally by looking at the
generated files:

```
  <some_base_build_dir>/<full_build_name>/smart-jenkins-driver.out
```

(e.g. `<full_build_name>` = `Trilinos-atdm-<system_name>-gnu-opt-debug`) and
also examine the generated `*.xml` configure, build, and test files created
under:

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

See
[TRIBITS_CTEST_DRIVER()](https://tribits.org/doc/TribitsDevelopersGuide.html#determining-what-testing-related-actions-are-performed-tribits-ctest-driver)
for a description of all of the options that can be set as env vars to, for
example, skip configure, skip the build, skip running tests, etc.


## Specific instructions for each system

* <a href="#ridewhite">ride/white</a>
* <a href="#shillerhansen">shiller/hansen</a>
* <a href="#chamaserrano">chama/serrano</a>
* <a href="#mutrino">mutrino</a>
* <a href="#sems-rhel6-environment">SEMS RHEL6 Environment</a>
* <a href="#sems-rhel7-environment">SEMS RHEL7 Environment</a>
* <a href="#cee-rhel6-environment">CEE RHEL6 Environment</a>
* <a href="#waterman">waterman</a>


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

$ bsub -x -Is -q rhel7F -n 16 ctest -j16
```

The ATDM configuration of Trilinos is set up to run on the Firestone nodes
(Dual-Socket POWER8, 8 cores per socket, K80 GPUs).  This configuration will
not work on the other GPU nodes currently.

**NOTE:** While the above example shows loading the environment, configuring
and building on the login node, one can also do these on the compute nodes as
well.  In fact, that is what the CTest -S drivers do in automated testing on
'white' and 'ride'.

Note that one can also run the same build a tests using the <a
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

$ srun ctest -j16
```

**NOTE:** While the above example shows loading the environment, configuring
and building on the login node, one can also do these on the compute nodes as
well.  In fact, that is what the CTest -S drivers do in automated testing on
'hansen' and 'shiller'.

Note that one can also run the same build a tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ srun ./checkin-test-atdm.sh intel-opt-openmp \
  --enable-packages=MueLu \
  --local-do-all
```


### chama/serrano

Once logged on to 'chama' or 'serrano', one can directly configure and build
on the login node (being careful not to overload the node) using the `chama`
and `serrano` envs, respectively.  But to run the tests, one must run on the
compute nodes using the `srun` command.  For example, to configure, build and
run the tests for say `MueLu` on 'serrano' or 'chama', (after cloning Trilinos
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

$ salloc -N1 --time=0:20:00 --account=<YOUR_WCID> ctest -j16
```

To get information on <YOUR_WCID> used above, there is a WC tool tab on
computing.sandia.gov

**NOTE:** Unlike some of the other machines, one must load the environment,
configure and build on the login node and then run the test suite on a compute
node on this system.  This is what the CTest -S driver on 'chama' and
'serrano' does in order to drive jobs and submit to CDash.

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

Once logged on to 'mutrino', one can directly configure and build on the login
node (being careful not to overload the node) using the `mutrino` env.  But to
run the tests, one must run on the compute nodes using the `salloc` command.
For example, to configure, build and run the tests for say `MueLu` on
'mutrino', (after cloning Trilinos on the `develop` branch) one would:


```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh intel-opt-openmp

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

NOTE: Above including `sems-rhel6` in the build name
`sems-rhel6-intel-opt-openmp` is not necessary but is recommended when on a
CEE LAN RHEL6 machine to be explicit that the SEMS env is being used and not
the <a href="#cee-rhel6-environment">CEE RHEL6 env</a>.

One can also run the same build a tests using the <a
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

One can also run the same build a tests using the <a
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


### CEE RHEL6 Environment

Once logged into any CEE LAN RHEL6 SRN machine, one can configure, build, and
run tests for any ATDM Trilinos package using the `cee-rhel6` env.  For
example, to configure, build and run the tests for the
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

One can also run the same build a tests using the <a
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

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cuda-debug

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=20

$ bsub -x -Is -n 20 ctest -j20
```

**NOTE:** While the above example shows loading the environment, configuring
and building on the login node, one can also do these on the compute nodes as
well.  In fact, that is what the CTest -S drivers do in automated testing on
'waterman'.

Note that one can also run the same build a tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ bsub -x -Is -n 20 \
  ./checkin-test-atdm.sh cuda-debug --enable-packages=MueLu --local-do-all
```


## Building and installing Trilinos for ATDM Applications

The sections below describe how to configure, build, and install Trilinos for
usage by the ATDM applications:

* <a href="#building-and-installing-trilinos-for-empire">Building and installing Trilinos for EMPIRE</a>
* <a href="#building-and-installing-trilinos-for-sparc">Building and installing Trilinos for SPARC</a>


### Building and installing Trilinos for EMPIRE

Configuring, building and installing Trilinos for EMPIRE and then building and testing EMPIRE against that Trilinos installation is an easy process.  To build the ATDM Trilinos configuration on any supported system just do:

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh <build-name>

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake,cmake/std/atdm/apps/empire/EMPIRETrilinosPackagesEnables.cmake \
  -DCMAKE_INSTALL_PREFIX=<trilinos-install-dir> \
  $TRILINOS_DIR

$ make NP=16  # Uses ninja -j16

$ make NP=16 install
```

Once the Trilinos installation is complete, one can configure, build, and test EMPIRE with:

```
$ cd <empire-build-dir>/

$ source <trilinos-install-dir>/load_matching_env.sh

$ rm -r CMake*

$ cmake \
  -GNinja \
  -DTrilinos_INSTALL_DIR=${ATDM_TRILINOS_INSTALL_PREFIX} \
  [other options] \
  <empire-src-dir>

$ ninja -j20

$ ctest -j8
```

The EMPIRE configuration gets everything it needs from the sourced and loaded environment and from the Trilinos installation directory.  Easy as pie!

Details on where to get the EMPIRE source repos and the exact CMake options to use when configuring EMPIRE are given at:

* [Building and Testing EMPIRE against Local Trilinos Installation](https://snl-wiki.sandia.gov/display/CoodinatedDevOpsATDM/Building+ATDM+APPs+Against+Local+Installs+of+Trilinos#BuildingATDMAPPsAgainstLocalInstallsofTrilinos-BuildingandTestingEMPIREagainstLocalTrilinosInstallation)


### Building and installing Trilinos for SPARC


To configure, build, and install Trilinos for usage by SPARC, one must use a very specific name for the installation directory or the CMake configure of SPARC will not pick it up.  Specific instructions and some helper scripts for building and installing Trilinos and then building and testing SPARC against Trilinos installation are given at:

* [Building and Testing SPARC against Local Trilinos Installation](https://snl-wiki.sandia.gov/display/CoodinatedDevOpsATDM/Building+ATDM+APPs+Against+Local+Installs+of+Trilinos#BuildingATDMAPPsAgainstLocalInstallsofTrilinos-BuildingandTestingSPARCagainstLocalTrilinosInstallation)

But for a specific build example on a CEE RHEL6 machine, one wouild configure, build, and install Triilnos for SPARC using:

```
$ cd <some_build_dir>/

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh cee-rhel6-gnu-7.2.0-openmpi-1.10.2-serial-static-dbg

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake,cmake/std/atdm/apps/sparc/SPARCTrilinosPackagesEnables.cmake \
  -DCMAKE_INSTALL_PREFIX=<base-dir>/cee-cpu_gcc-7.2.0_openmpi-1.10.2_serial_static_dbg
  $TRILINOS_DIR

$ make NP=16  # Uses ninja -j16

$ make NP=16 install
```

In addition, configuring and building SPARC requires special CMake configure scripts be used.  Again, consult the above web page and helper scripts.


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
    custom_bulds.sh  # [Optional] Special logic for compiler keywords, etc.
    tweaks/
       <COMPILER0>-<BUILD_TYPE0>-<NODE_TYPE0>-<KOKKOS_ARCH0>.cmake  # [Optional]
       <COMPILER1>-<BUILD_TYPE1>-<NODE_TYPE1>-<KOKKOS_ARCH0>.cmake  # [Optional]
       ...
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
example, see `atdm/cee-rhel6/cutome-builds.sh` and
`atdm/cee-rhel6/environment.sh`.)

<a name="ATDM_TWEAKS_FILES"/>

The **ATDM TWEAKS FILES** in the `cmake/std/atdm/<system-name>/tweaks/`
directory contain special settings for specific builds for a specific system.
Typically, this file contains (temporary) disables for tests for that given
build.  When a configure is performed, the internal CMake variable
`ATDM_BUILD_NAME_KEYS_STR` set to
`<COMPILER>-<BUILD_TYPE>-<NODE_TYPE>-<KOKKOS_ARCH>` (printed to STDOUT) is
used to define a default file name:

```
  Trilinos/cmake/std/atdm/<system-name>/tweaks/${ATDM_BUILD_NAME_KEYS_STR}.cmake
```

If that file exists, then it is set as the default for the cmake cache
variable `ATDM_TWEAKS_FILES` (prints to STDOUT) and that file is included and
its options are read.  For example, this is what the output looks like on
'waterman':

```
-- Reading in configuration options from cmake/std/atdm/ATDMDevEnv.cmake ...
-- ATDM_BUILD_NAME_KEYS_STR='GNU-RELEASE-OPENMP-POWER9'
-- ATDM_TWEAKS_FILES='<...>/Trilinos/cmake/std/atdm/waterman/tweaks/GNU-RELEASE-OPENMP-POWER9.cmake'
-- Including ATDM build tweaks file <...>//Trilinos/cmake/std/atdm/waterman/tweaks/GNU-RELEASE-OPENMP-POWER9.cmake ...
```


## Disabling failing tests

There are situations where specific tests must be disabled on certain
platforms for certain builds or based on other criteria (see sub-process
[Temporarily disable the failing code or
test](https://snl-wiki.sandia.gov/display/CoodinatedDevOpsATDM/Triaging+and+addressing+ATDM+Trilinos+Failures#TriagingandaddressingATDMTrilinosFailures-5.Makesuretheissueisaddressedinatimelyway:)).
There are various ways to disable tests with the Trilinos TriBITS/CMake-based
build and test system.  Tests can be disabled in the `CMakeLists.txt` files
that define the tests themselves using various logic.  But the way to
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
-- ATDM_TWEAKS_FILES='.../Trilinos/cmake/std/atdm/mutrino/tweaks/INTEL-DEBUG-OPENMP-KNL.cmake'
```

For example, Trilinos commit
[73ae19cf0c](https://github.com/trilinos/Trilinos/commit/73ae19cf0cc7295a7f36de342ea51718226825b7)
shows the disable of the test:

```
# Disable test that times out for some unkown reason (#2925)
ATDM_SET_ENABLE(Stratimikos_test_aztecoo_thyra_driver_MPI_1_DISABLE ON)
```

in both the files `cmake/std/atdm/shiller/tweaks/GNU-DEBUG-SERIAL.cmake` and
`cmake/std/atdm/shiller/tweaks/GNU-RELEASE-SERIAL.cmake` (before `-HSW` was
added to the names).

NOTE: Adding a comment with the Trilinos GitHub Issue ID (`#2925` in this
example) is critical for tractability to remind why the test was disabled and
where to find the disable to remove it later.

To avoid duplicate `ATDM_SET_ENABLE()` statements, one can use the approach of
creating a single `*.cmake` that is included in the different tweak files as
described below.


### Disable a test for several or all builds on a specific platform

It is often the case that a test needs to be disabled for several (or all)
builds for a given platform.  An efficient way to do this is to create a new
`*.cmake` file that contains the `ATDM_SET_ENABLE()` statements and then
include that new file in all of the tweaks files on that system where the
tests should be disabled.

For example, the Trilinos commit
[3450efd421](https://github.com/trilinos/Trilinos/commit/3450efd421f1ce2b47700853aa4c5801f667202a)
shows how a set of tests were disabled for all of the CUDA builds on the
system `ride` through the creation of the file:

```
  Trilinos/cmake/std/atdm/ride/tweaks/CUDA_COMMON_TWEAKS.cmake
```

and then the inclusion of that file in the specific tweak files for each CUDA
build:

```
  Trilinos/cmake/std/atdm/ride/tweaks/CUDA-DEBUG-CUDA.cmake
  Trilinos/cmake/std/atdm/ride/tweaks/CUDA-RELEASE-CUDA.cmake
```

(before `-POWER8-KEPLER37` was added to the names) using the inserted CMake
statement:

```
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/CUDA_COMMON_TWEAKS.cmake")
```

An example of using a `*.cmake` file to disable the same set of tests in all
of the builds for a given system is shown in Trilinos commit
[33a933b004](https://github.com/trilinos/Trilinos/commit/33a933b004f88710274906fad612380049e1e82e).
This example shows the creation of the file:

```
  Trilinos/cmake/std/atdm/ride/tweaks/ALL_COMMON_TWEAKS.cmake
```

and then the inclusion of that file in all of the specific tweaks files on
'ride' with the statement:

```
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ALL_COMMON_TWEAKS.cmake")
```

in each of those files.


### Disable a test for builds on all platforms

There are rare cases where a test (or a set of tests) will need to be disabled
for all or a subset of builds across all systems.  To do that using the
above-described approaches would require touching files in every
`cmake/std/atdm/<system-name>/tweaks/` directory.  In rare cases like this,
one can accomplish this by adding a (conditional) `ATDM_SET_ENABLE()`
statement for each test disable directly to the file:

```
  Trilinos/cmake/std/atdm/ATDMDevEnvSettings.cmake
```

For example, Trilinos commit [5e52db03ff](https://github.com/trilinos/Trilinos/commit/5e52db03ff33acb5b9a0be7ba7507a8bb0de6e30) added the CMake code:

```
# Disable test that fails for all openmp builds (#3035)
IF (ATDM_USE_OPENMP)
  ATDM_SET_ENABLE(MueLu_UnitTestsTpetra_MPI_4_DISABLE ON)
ENDIF()
```

to the file `ATDMDevEnvSettings.cmake` to disable the test
`MueLu_UnitTestsTpetra_MPI_4` for all OpenMP builds across all platforms.
(Note that that disable was later removed in Trilinos commit
[62fa6663a6](https://github.com/trilinos/Trilinos/commit/62fa6663a6d5a757d786ac87752c3e2074d28414)
after the test was fixed.)


## Specific systems supported

The specific `cmake/std/atdm/<system-name>/` sub-directories and the systems
they support are:

* `cee-rhel6/`: CEE LANL RHEL6 systems with a CEE environment

* `chama/`: Supports SNL HPC machine `chama`.

* `mutrino/`: Supports SNL HPC machine `mutrino`.

* `ride/`: Supports GNU and CUDA builds on both the SRN machine `ride` and the
  mirror SON machine `white`.

* `sems-rhel6/`: SNL COE RHEL6 systems with the SEMS NFS environment

* `sems-rhel7/`: SNL COE RHEL7 systems with the SEMS NFS environment

* `serrano/`: Supports SNL HPC machine `serrano`.

* `shiller/`: Supports GNU, Intel, and CUDA builds on both the SRN machine
  `shiller` and the mirror SON machine `hansen`.

* `waterman/`: Supports GNU and CUDA builds on the SRN machine `waterman`.
