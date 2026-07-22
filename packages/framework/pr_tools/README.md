# Trilinos Pull Request Testing Tools

## `PullRequestLinuxDriverTest.py`
This script is intended to drive the "core" configure/build/test/report part of the process.  As such, this is the tool that is called as part of most of the GitHub Actions (seen in `.github/workflows/AT2.yml`).  Its primary purpose is to use the abstraction provided by [the `GenConfig` tool](https://github.com/sandialabs/GenConfig) (the configuration name, e.g. `rhel8_gcc-openmpi-openmp_release-debug_static_no-kokkos-arch_no-asan_no-complex_no-fpic_mpi_no-pt_no-rdc_no-uvm_deprecated-on_no-package-enables`) to decode out the individual CMake options and environment settings and set them up / pass them onwards to the CTest driver scripts at `cmake/SimpleTesting` (documentation [here](https://github.com/trilinos/Trilinos/tree/develop/cmake/SimpleTesting/README.md)).  It assumes that the `GenConfig` tool (and said tools dependencies, e.g. `LoadEnv`) are already set up e.g. via a call to the `get_dependencies.sh` script at `packages/framework/get_dependencies.sh`.

### Behavior
1. Call `LoadEnv` to determine environment changes from `genconfig-configuration`
2. Apply environment changes from previous step
3. Determine which packages to enable in `packageEnables.cmake` and the subproject labels for `package_subproject_labels.cmake` via calling `commonTools/framework/get-changed-trilinos-packages.sh`
4. Move to build directory
5. Call `GenConfig` to write CMake config file (`generatedPRFragment.cmake`) for ingestion by CMake later in process
6. Call into `cmake/SimpleTesting/ctest-driver.cmake` with appropriate arguments

## `PullRequestLinuxDriverMerge.py`
This script does the merge-related parts of the CI testing process.  It merges a passed source ref (branch or commit) into a passed target branch.  It is not needed if something else (e.g. GitHub Actions) performs the merge.

### Assumptions
* There is a local clone of Trilinos already existing inside of a passed location (which must be named `Trilinos`)
* `origin` exists as a remote and is set to the target repository remote

## `PullRequestLinuxDriver.sh`
This script is a "driver of drivers".  It is designed to do some shell-specific pre-behavior that is difficult to handle in Python prior to calling the underlying "worker" utilities (`PullRequestLinuxDriver*.py`)

### Behavior
1. Set up environment
2. If testing `develop` version of Kokkos/KokkosKernels, call `SetKokkosDevelop.sh` and go straight to step 5
3. Call merge utility
4. If this script or the merge script (pre-merge things that could influence behavior of the merge) have changed, re-execute using new versions of pre-merge tools
5. Call configure/build/test/report utility

### Assumptions
* The script used will be contained within the repository being tested (i.e. it cannot test a different clone of the Trilinos repository)

## `SetKokkosDevelop.sh`
This script sets up clones of Kokkos and KokkosKernels from their home repositories within Trilinos, with the intention of allowing a given clone of Trilinos to be tested against the current development version of Kokkos/KokkosKernels.

## `LaunchDriver.py`
Highest-level entry point.  Use information about the system (e.g. hostname) to decide to run specific workflows.  In the past, this meant running a different version of `PullRequestLinuxDriver.sh` depending on e.g. whether or not a specific HPC was being used (where that script would run a queue submission command to submit the "real" `PullRequestLinuxDriver.sh` inside of an HPC allocation).  This is not currently necessary because current CI testing all happens on similar machines (i.e. no queue, Linux, x86_64).


# Future plans

`PullRequestLinuxDriverTest.py` will continue to be the primary testing driver.  Much of the work it performs is still important and relevant to the CI testing process.  However, it is fair to question its necessity if one classifies it as merely a Python wrapper for the CTest driver script doing the bulk of the "work" (`cmake/SimpleTesting/ctest-driver.cmake`).  Working with Python is considerably more straightforwards than the CMake language, so keeping this script for now is certainly sensible.

`PullRequestLinuxDriverMerge.py` can be removed as soon as there is no need for a local tool to perform merge operations.  This capability is fully replaced by GitHub Actions for AutoTester2-based builds.  It should be carefully considered if any other testing approaches need to persist (e.g. running AT1-based builds on any platforms), and if so, this should probably remain in-place.

`PullRequestLinuxDriver.sh` can be removed as soon as `PullRequestLinuxDriverMerge.py` is no longer needed, as the primary purpose is to handle the complex re-launch behavior that may be induced by changes to the tooling itself.

`SetKokkosDevelop.sh` can either persist or be re-authored, but the capability is required by doing upstream integration builds of Kokkos/KokkosKernels.  In the long run it probably makes the most sense to move this in-line into a scheduled GitHub Action, further eliminating the need for `PullRequestLinuxDriver.sh`.

`LaunchDriver.py` can be removed as long as different launch types are not needed on different systems.  The primary reason for its existance was to enable the paradigm of launching a tool within a queue allocation, enabling a queue wrapper to be used on specific systems.  There are other ways to handle this particular concern, so it may be possible to eliminate it even if that use case re-occurs (it is not currently used).
