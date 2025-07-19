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
| `skip_upload_config_files` | BOOL   |    NO     | OFF                                           | Skip upload config files (???)                                     |
| `build_root`               | STRING |    NO     | `${source_dir}/nightly_testing`               | Used to generate `build_dir` if `build_dir` is not defined.        |
| `build_dir`                | STRING |    NO     | `${build_root}/${CTEST_BUILD_NAME}`           | Path to the build directory.                                       |
| `PARALLEL_LEVEL`           | STRING |    NO     | `<num cores>`                                 |                                                                    |
| `TEST_PARALLEL_LEVEL`      | STRING |    NO     | `${PARALLEL_LEVEL}`                           |                                                                    |
| `SKIP_RUN_TESTS`           | BOOL   |    NO     | OFF                                           | Skip running any tests (any tests enabled will still compile)      |


1. It might worthwhile to remove `build_root` since it's only used to create `build_dir` IF `build_dir` is not passed in
   via a `-D` option.
2. Related to (1), we might also change `build_dir` to be `BUILD_DIR` and pass that in.


Example CTest call from a Trilinos PR
-------------------------------------
This is an example, for reference, of how the `ctest` command is invoked in the current Trilinos
PR test driver.

Assumptions:
- Existing generated CMake fragment file (usually generated by GenConfig with `--cmake-fragment=GenConfigFragment.cmake` option)
- Existing `packageEnables.cmake` and `package_subprojects_list.cmake` file in build directory
- `$TRILINOS_DIR` is set to the root source directory for your copy of Trilinos
- `$BUILD_DIR` is set to your created out-of-source directory for building Trilinos

NOTE: `generatedPRFragment.cmake` file can be generated by following Step 8. of [our PR reproduction steps](https://github.com/trilinos/Trilinos/wiki/Reproducing-Pull-Request-Testing-Errors-%E2%80%90%E2%80%90-AutoTester1-(Legacy))

NOTE: `packageEnables.cmake` and `package_subprojects_list.cmake` file can be generated by running
```shell
$TRILINOS_DIR/commonTools/framework/get-changed-trilinos-packages.sh \
    <upstream_remote_name>/develop \
    HEAD \
    packageEnables.cmake \
    package_subproject_list.cmake
```

From your `$BUILD_DIR` with a GenConfig environment loaded, run the following command to launch the Trilinos CTest driver
that will configure, build, and test.
```shell
ctest \
   -S ${TRILINOS_DIR}/cmake/SimpleTesting/cmake/ctest-driver.cmake \
   -Dsource_dir=${TRILINOS_DIR} \
   -Dbuild_name=PR-Manual-Trilinos-CTest-Driver-Run-01 \
   -Dskip_by_parts_submit=OFF \
   -Dskip_update_step=ON \
   -Ddashboard_model=Experimental \
   -Ddashboard_track=Experimental \
   -DPARALLEL_LEVEL=20 \
   -DTEST_PARALLEL_LEVEL=20 \
   -Dbuild_dir=${BUILD_DIR} \
   -Dconfigure_script=${BUILD_DIR}/GenConfigFragment.cmake
   -Dpackage_enables=${BUILD_DIR}/packageEnables.cmake
   -Dsubprojects_file=${BUILD_DIR}/package_subproject_list.cmake
```

See `TrilinosPRConfigurationStandard.py`[2] for information on what options are set to something
other than the default during normal Trilinos PR operations.

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
[2]: https://github.com/trilinos/Trilinos/blob/master/packages/framework/pr_tools/trilinosprhelpers/TrilinosPRConfigurationStandard.py
