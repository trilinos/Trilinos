CMake File Structure
====================
This CMake structure attempts to capture the structure of the original
`TFW_testing_single_configure_prototype` that is used by the current / old Trilinos
Pull-Request (PR) framework.

**Note**: _Trilinos uses the **`pull_request_changes`** branch from this repository_.

Specifically, this new configuration will replace the existing `simple_testing.cmake`
script.

This diagram describes the current structure of the CMake files and how they
interact with each other.
<center>
     <img src="img/CMake-structure.png" alt="CMake Structure Diagram" height="600"/>
</center>

CMake Files
-----------

| Filename                      | Guarded? | Purpose                                                                  |
|-------------------------------|:--------:|--------------------------------------------------------------------------|
| `ctest-driver.cmake`          |          | Main driver script / entry-point                                         |
| `ctest-common.cmake`          |    X     | Common activities such as option handling, etc.                          |
| `ctest-cdash-setup.cmake`     |    X     | Handle CDash configurations                                              |
| `ctest-stage-configure.cmake` |          | Handle the _configuration_ stage of the test                             |
| `ctest-stage-build.cmake`     |          | Handle the _build_ stage of the test                                     |
| `ctest-stage-test.cmake`      |          | Handle the _test_ stage of the test.                                     |
| `ctest-functions.cmake`       |    X     | Functions and macros                                                     |
| `ctest-CTestConfig.cmake`     |          | Required file that is copied into the build directory at configure time. |

The _guarded_ files use the CMake command [`include_guard()`][1] which should prevent that file
from being include more than once in an include chain.

Options and Variables (`simple_testing.cmake`)
----------------------------------------------
The `simple_testing.cmake` file has a number of optional parameters that can be sent into the
CMake system via `-D<VARNAME>:BOOL=<ON|OFF>` parameters:

| Option                     | Required? | Default                                     | PR Override | Purpose                          |
|----------------------------|:---------:|---------------------------------------------|-------------|----------------------------------|
| `build_name`               |    YES    | N/A                                         | YES         | Sets `CTEST_BUILD_NAME`          |
| `ctest_submit_retry_count` |    NO     | 5                                           |             |                                  |
| `ctest_submit_retry_delay` |    NO     | 3                                           |             |                                  |
| `dashboard_model`          |    NO     | `Experimental`                              | YES but NO  | Set to the same value as default |
| `dashboard_track`          |    NO     | `Experimental`                              | YES but NO  | Set to the same value as default |
| `skip_by_parts_submit`     |    NO     | `ON`                                        | YES         |                                  |
| `skip_clean_build_dir`     |    NO     | `ON`                                        |             |                                  |
| `skip_single_submit`       |    NO     | `ON`                                        |             |                                  |
| `skip_update_step`         |    NO     | `OFF`                                       | YES         |                                  |
| `skip_upload_config_files` |    NO     | `OFF`                                       |             |                                  |
| `PARALLEL_LEVEL`           |    NO     | _num cores_                                 | YES         |                                  |
| `TEST_PARALLEL_LEVEL`      |    NO     | `${PARALLEL_LEVEL}`                         | YES         |                                  |
| `build_root`               |    NO     | `${CTEST_SOURCE_DIRECTORY}/nightly_testing` |             |                                  |
| `build_dir`                |    NO     | `${build_root}/${CTEST_BUILD_NAME}`         | YES         | The CMake build dir              |
| `configure_script`         |    YES    | N/A                                         | YES         | See note below table             |
| `package_enables`          |    YES    | N/A                                         | YES         | `packageEnables.cmake`           |
| `subprojects_file`         |    YES    | N/A                                         | YES         | `package_subproject_list.cmake`  |

1. `configure_script` points to the `cmake/std/PullRequestLinux<COMPILER><COMPILER_VERSION>TestingSettings.cmake` file.
    - Example: `${WORKSPACE}/Trilinos/cmake/std/PullRequestLinuxGCC8.3.0TestingSettings.cmake`

See `TrilinosPRConfigurationStandard.py`[2] for information on what options are set to something
other than the default during normal Trilinos PR operations.

Expected Operation
------------------
The expected operation of this set of files to replace the old `simple_testing.cmake` is to load
the [`ctest-driver.cmake`](ctest-driver.cmake) file in its place.


`ctest-driver.cmake` Options
----------------------------

| Option                     | Type   | Required? | Default                                       | Description                                                        |
|----------------------------|--------|:---------:|-----------------------------------------------|--------------------------------------------------------------------|
| `build_name`               | STRING |    YES    | N/A                                           | The build name, see Jenkins' `${JOB_NAME}` envvar.                 |
| `subprojects_file`         | STRING |    YES    | N/A                                           | This is the `package_subproject_list.cmake` file.                  |
| `source_dir`               | PATH   |    YES    | N/A                                           | Path to the source directory.                                      |
| `configure_script`         | STRING |    YES    | N/A                                           | Test settings CMake script.                                        |
| `package_enables`          | STRING |    YES    | N/A                                           | This is the `packageEnables.cmake` file.                           |
| `ctest_submit_retry_count` | STRING |    NO     | 5                                             | Number of times to retry a ctest submssion.                        |
| `ctest_submit_retry_delay` | STRING |    NO     | 3                                             | Delay (seconds) between attempts to submit to cdash.               |
| `dashboard_model`          | STRING |    NO     | Experimental                                  | CDash model                                                        |
| `dashboard_track`          | STRING |    NO     | Experimental                                  | CDash track                                                        |
| `skip_clean_build_dir`     | BOOL   |    NO     | ON                                            | Skip cleaning the build directory (`ctest_empty_binary_directory`) |
| `skip_update_step`         | BOOL   |    NO     | OFF                                           | Skip the update step (`ctest_update()`) of the repository.         |
| `skip_by_parts_submit`     | BOOL   |    NO     | ON                                            | Skip submission to CDash after each phase.                         |
| `skip_single_submit`       | BOOL   |    NO     | ON                                            | Skip single submission                                             |
| `skip_upload_config_files` | BOOL   |    NO     | OFF                                           | Skip upload config files (???)                                     |
| `build_root`               | STRING |    NO     | `${source_dir}/nightly_testing`               | Used to generate `build_dir` if `build_dir` is not defined.        |
| `build_dir`                | STRING |    NO     | `${build_root}/${CTEST_BUILD_NAME}`           | Path to the build directory.                                       |
| `PARALLEL_LEVEL`           | STRING |    NO     | `<num cores>`                                 |                                                                    |
| `TEST_PARALLEL_LEVEL`      | STRING |    NO     | `${PARALLEL_LEVEL}`                           |                                                                    |

1. It might worthwhile to remove `build_root` since it's only used to create `build_dir` IF `build_dir` is not passed in
   via a `-D` option.
2. Related to (1), we might also change `build_dir` to be `BUILD_DIR` and pass that in.



Example CTest call from a Trilinos PR
-------------------------------------
This is an example, for reference, of how the `ctest` command is invoked in the current/old Trilinos
PR test driver.
```bash
ctest \
   -S simple_testing.cmake \
   -Dsource_dir=${WORKSPACE}/Trilinos \
   -Dbuild_name=PR-9495-test-Trilinos_pullrequest_gcc_8.3.0-5164 \
   -Dskip_by_parts_submit=OFF \
   -Dskip_update_step=ON \
   -Ddashboard_model=Experimental \
   -Ddashboard_track=Pull Request \
   -DPARALLEL_LEVEL=20 \
   -DTEST_PARALLEL_LEVEL=4 \
   -Dbuild_dir=${WORKSPACE}/pull_request_test \
   -Dconfigure_script=${WORKSPACE}/Trilinos/cmake/std/PullRequestLinuxGCC8.3.0TestingSettings.cmake \
   -Dpackage_enables=../packageEnables.cmake \
   -Dsubprojects_file=../package_subproject_list.cmake
```

Additional Notes and Pitfalls
=============================

`ctest_test()` and `CAPTURE_CMAKE_ERROR`
----------------------------------------
For Trilinos testing we should avoid checking the value returned by `ctest_test()`
for the `CAPTURE_CMAKE_ERROR` parameter. The reason for this is that CTest will
flag this as an error (i.e., -1 is returned) if there were _no tests run_.

The _no tests were run_ though is a valid 'success' result for Trilinos PR tests since
this project enables packages dynamically based on what packages have modifications.
This allows some PR's to go through without building Trilinos, which is advantageous
when only documentation or perhaps the testing framework itself is modified and we do
not need to spend O(5 hours) for the test suite to run.

[1]: https://cmake.org/cmake/help/latest/command/include_guard.html
[2]: https://github.com/trilinos/Trilinos/blob/master/cmake/std/trilinosprhelpers/TrilinosPRConfigurationStandard.py
