----------------------------------------
ChangeLog for TriBITS
----------------------------------------

## 2023-5-03:

* **Added:** Added support for non-fully TriBITS-compatible external packages.
  Now, a `<Package>Config.cmake` file need not define
  `<UpstreamPkg>::all_libs` targets for all of its upstream dependencies.  The
  updated macro `tribits_process_enabled_tpls()` will find any missing
  upstream external packages/TPLs as needed (see updated documentation in the
  section "TriBITS-Compliant External Packages" in the "TriBITS Users Guide"
  and the section "Processing of external packages/TPLs and TriBITS-compliant
  external packages" in the "TriBITS Maintainers Guide").

## 2023-02-24:

* **Changed:** Upgraded minimum required CMake version from 3.17 to 3.23.
  Existing TriBITS projects that have already upgraded to require CMake 3.23+
  should not notice any major changes due to this change.

## 2023-02-21:

* **Added:** Added support for pre-installed internal packages treated as
  external packages.  Now, any set of internally defined TriBITS packages for
  a TriBITS project can be pre-built and pre-installed and the remaining
  packages in the TriBITS project can be configured to point to those by
  setting `-D TPL_ENABLE_<Package>=ON`.  This allows great flexibility in how
  a TriBITS project's packages can be and built, installed, and deployed.
  This technically implements "Use Case 3: Configure/build pointing to a
  subset of already installed TriBITS packages in same repo" in [TriBITS
  #63](https://github.com/TriBITSPub/TriBITS/issues/63).  See the section
  "Building against pre-installed packages" in the updated build reference
  documentation for details.

* **Fixed:** Setting `-D<Project>_ENABLE_<TplName>=ON` for an external
   package/TPL `<TplName>` will not correctly enable and process the TPL.

## 2022-12-07:

* **Changed:** Setting `-D<Project>_ENABLE_<TplName>=ON` now triggers the
  enable an external package/TPL `<TplName>` similar to the way that
  `-DTPL_ENABLE_<TplName>` has always done.  This is technically a change in
  backward compatibility because setting `<Project>_ENABLE_<TplName>=ON` for
  an external package/TPL used to be ignored.  This change was done as part of
  a general refactoring to unify the handling of internal and external
  packages and is a side effect of that refactoring (see [TriBITS
  #63](https://github.com/TriBITSPub/TriBITS/issues/63)).  (Initially, setting
  `-D<Project>_ENABLE_<TplName>=ON` resulted in a failed configure because the
  refactoring was not complete for the handling of external packages/TPL.  But
  this was fixed in a latter update.)

## 2023-10-25:

* **Added:** New option `<Project>_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES`
  skips the install of the project-level `<Project>Config.cmake` file.  The
  default value is ``FALSE`` so as to maintain backward compatibility.  (The
  project can change the default.)

* **Changed:** External packages/TPLs are now processed at the base project
  scope level.  This allows simple `set()` statements in package module files
  or config files included by `find_package()` to have project-level scope for
  the entire TriBITS project.  This is more similar to how a raw CMake project
  would usually behave that calls `find_package()` in the base
  `CMakeLists.txt` file.  Before, calls to `find_package()` were wrapped in a
  CMake `function()` called from the base project directory scope.  So while
  IMPORTED targets created from a `find_package()` command where visible at
  the base directory project-level scope, local variables were not.  With this
  change, now they are.

## 2023-01-10:

* **Added:** Added back support for deprecated variable
  `<Project>_ASSERT_MISSING_PACKAGES` that was removed
  [2022-10-11](#2022-10-11).  When `<Project>_ASSERT_MISSING_PACKAGES` is set
  to a non-null value, it overrides the default value for
  `<Project>_ASSERT_DEFINED_DEPENDENCIES` (but setting
  `<Project>_ASSERT_DEFINED_DEPENDENCIES` in the cache takes precedence).

## 2023-01-06:

* **Changed:** Changed all TPL dependencies back to 'Optional' so that
  disabling an external package/TPL will **not** disable any downstream
  external packages/TPLs that list a dependency on that external package/TPL.
  This undoes the change on [2022-10-20](#2022-10-20) and restores backward
  compatibility to the behavior before that change.

## 2022-12-20:

* **Deprecated:** The macro `set_and_inc_dirs()` is deprecated and replaced by
  `tribits_set_and_inc_dirs()`.  Use the script
  `TriBITS/refactoring/replace_set_and_inc_dirs_r.sh` to update
  `CMakeLists.txt` files.

## 2022-11-03:

* **Deprecated:** The long-deprecated TriBITS function override
  `include_directories()` now emits a deprecated warning.  To replace all
  usages of `include_directories()` that should be
  `tribits_include_directories()`, use the script
  `TriBITS/refactoring/replace_include_directories_r.sh` (see documentation in
  that script).

* **Deprecated:** Many previously deprecated TriBITS features now will trigger
  a CMake DEPRECATION warning message by default (by calling
  `message(DEPRECATION ...)`).  The message printed to the CMake output will
  typically describe how to remove the usage of the deprecated feature.  To
  remove deprecation warnings, change to use the non-deprecated features
  mentioned in the deprecation warning message.  To temporarily disable
  deprecation warnings, configure with `-D
  TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE=IGNORE` (see build reference entry
  for `TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE` for more details).

## 2022-10-20:

* **Changed:** Disabling an external package/TPL will now disable any
  downstream external packages/TPLs that list a dependency on that external
  package/TPL through its
  [`FindTPL<tplName>Dependencies.cmake`](https://tribitspub.github.io/TriBITS/users_guide/index.html#findtpl-tplname-dependencies-cmake)
  file.  Prior to this, disabling an external package/TPL would not disable
  dependent downstream external packages/TPLs (it would only disable
  downstream dependent required internal packages).  To avoid this, simply
  leave the enable status of the upstream external package/TPL empty "" and no
  downstream propagation of disables will take place.

## 2022-10-16:

* **Removed:** Removed the variables `<Project>_LIBRARY_DIRS`,
  `<Project>_TPL_LIST` and `<Project>_TPL_LIBRARIES` from the installed
  `<Project>Config.cmake` file.  These are not needed after the change to
  modern CMake targets `<Package>::all_libs` (see `<Package>::all_libs`
  below).  To determine if a TPL is enabled, check `if (TARGET
  <tplName>::all_libs)`.  To get the libraries and include dirs for a TPL,
  link against the IMPORTED target `<tplName>::all_libs` (see the updated
  TriBITS example APP projects for details).

* **Removed:** Removed the variables `<Package>_PACKAGE_LIST`,
  `<Package>_TPL_LIST`, `<Package>_INCLUDE_DIR`, `<Package>_LIBRARY_DIRS`,
  `<Package>_TPL_INCLUDE_DIRS`, `<Package>_TPL_LIBRARIES` and
  `<Package>_TPL_LIBRARY_DIRS` from the generated `<Package>Config.cmake`
  files.  These are not needed with the move to modern CMake targets (see
  `<Package>::all_libs` below).

* **Changed:** Changed `<Package>_LIBRARIES` in generated
  `<Package>Config.cmake` files from the full list of the package's library
  targets to just `<Package>::all_libs`.  (There is no need to list the
  individual libraries after the move to modern CMake targets.)

## 2022-10-11:

* **Changed:** Added option `<Project>_ASSERT_DEFINED_DEPENDENCIES` to
  determine if listed external package/TPL and internal package dependencies
  are defined within the project or not.  The initial default is `FATAL_ERROR`
  for development mode and `IGNORE` for release mode.  (Previously, undefined
  external package/TPL dependencies where ignore.)  To set a different
  default, set `<Project>_ASSERT_DEFINED_DEPENDENCIES_DEFAULT` to `WARNING`,
  for example, in the project's `ProjectName.cmake` file.

* **Removed:** `<Project>_ASSERT_MISSING_PACKAGES` has been removed and setting
  it will result in a `FATAL_ERROR`.  Instead, use
  `<Project>_ASSERT_DEFINED_DEPENDENCIES` (and make sure all of your project's
  listed TPL dependencies are all defined within the project).

## 2022-10-02:

* **Changed:** The TriBITS FindTPLCUDA.cmake module changed
  `find_package(CUDA)` to `find_package(CUDAToolkit)` (the former is
  deprecated as of CMake 3.17).  This avoids imported target namespace
  conflicts with downstream CMake projects that call
  `find_package(CUDAToolkit)` (see [Trilinos
  #10954](https://github.com/trilinos/Trilinos/issues/10954)).

## 2022-09-16:

* **Changed:** Changed nomenclature for packages and TPLs (see updated
  "Maintainers Guide" section "TriBITS System Data Structures"): "TPLs" =>
  "External Packages/TPLs"; "Packages" => "Internal Top-Level Packages"; "SE
  Packages" => "Internal Packages". This impacted many internal variables as
  well as printed qualities.  Behavior should otherwise be identical
  w.r.t. input state.  The only observable change that users should see is the
  text used to describe the different sets of packages and TPLs.  (This is
  working towards a uniform handling of packages and TPLs (see [TriBITS
  #63](https://github.com/TriBITSPub/TriBITS/issues/63)).

* **Deprecated:** The rarely used input var
  `<Project>_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES` is deprecated
  and the new var name
  `<Project>_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES` should be used
  instead.

## 2022-08-22:

* **Added:** Added support for exporting cache variables for packages in their
    `<Package>Config.cmake` files using the new function
    `tribits_pkg_export_cache_var()`.

## 2022-08-18:

* **Changed:** Made setting parent package tests/examples enable/disable
  correctly propagate down to subpackages in a more intuitive way (see
  [TriBITSPub/TriBITS#268](https://github.com/TriBITSPub/TriBITS/issues/268)).
  This also results in not enabling tests for subpackages that are not
  explicitly enabled or enabled as part of the forward sweep of packages
  enables due to `<Project>_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON`.

## 2022-08-11:

* **Changed:** Fix a bug where the test dependencies for an enabled
  subpackage that resulted in a build error (see
  [trilinos/Trilinos#10842](https://github.com/trilinos/Trilinos/issues/10842)
  and
  [TriBITSPub/TriBITS#510](https://github.com/TriBITSPub/TriBITS/issues/510)).

## 2022-07-20:

* **Changed:** Fixed TriBITS generated and installed `<tplName>Config.cmake`
  files to not point into the build dir but instead point into relative dir
  for upstream TPL's when calling find_dependency() (see
  [TribitsPub/TriBITS#500](https://github.com/TriBITSPub/TriBITS/issues/500)).
  This also makes these files relocatable.

## 2022-07-14:

* **Added:** Added support for `FindTPL<tplName>Dependencies.cmake` with new
  macro `tribits_extpkg_define_dependencies()` that sets
  `<tplName>_LIB_ALL_DEPENDENCIES`.  Now `<tplName>_LIB_ENABLED_DEPENDENCIES`
  is automatically set from `<tplName>_LIB_ALL_DEPENDENCIES` based on what
  TPLs are actually enabled.  This avoids the problem described below from
  directly setting `<tplName>_LIB_ENABLED_DEPENDENCIES` without reguard to
  what TPLs are actually enabled.  This maintains backward compatibility for
  existing configure scripts where an upstream TPL may not be enabled in some
  strange configure scripts (see
  [TriBITSPub/TriBITS#494](https://github.com/TriBITSPub/TriBITS/issues/494)).

## 2022-05-25:

* **Changed:** Dependencies between external packages (TPLs) must now be
  specified in order for correct linkage.  No longer will listing the external
  packages (TPLs) in the correct order in the `<repoDir>/TPLsList.cmake` file
  and listing all upstream external packages (TPLs) in
  `<packageDir>/cmake/Dependencies.cmake` be sufficient.  For now,
  dependencies between external packages (TPLs) can be set in the
  `<packageDir>/TPLsList.cmake` file by setting the cache var
  `<tplName>_LIB_ENABLED_DEPENDENCIES` for each downstream external package
  (TPL).  (See
  [`TribitsExampleProject2/TPLsList.cmake`](https://github.com/TriBITSPub/TriBITS/blob/master/tribits/examples/TribitsExampleProject2/TPLsList.cmake)
  for an example.)  Later, a more scalable approach for setting these vars
  will be devised.  However, this means it is no longer necessary for a
  package to list all of its upstream external packages/TPLs, only its direct
  dependencies.

* **Changed:** All compliant external packages (TPLs) must now set the
  imported target `<tplName>::all_libs` in their `FindTPL<tplName>.cmake`
  files instead of the variables `TPL_<tplName>_INCLUDE` and
  `TPL_<tplName>_LIBRARIES`.

* **Changed:** The project-level variables `<Project>_INCLUDE_DIRS`,
  `<Project>_TPL_INCLUDE_DIRS`, `<Project>_LIBRARY_DIRS` and
  `<Project>_TPL_LIBRARY_DIRS` from the file `<Project>Config.cmake` and
  package-level variables `<Package>_INCLUDE_DIRS`,
  `<Package>_TPL_INCLUDE_DIRS`, `<Package>_LIBRARY_DIRS` and
  `<Package>_TPL_LIBRARY_DIRS` from the installed files
  `<Package>Config.cmake` are now set to empty.  Downstream CMake projects
  will need to link against the exported project-level target
  `<Project>::all_libs` or `<Project>::all_selected_libs` or the package-level
  targets `<Package>::all_libs` in order to get the list of include
  directories that are propagated through those targets by CMake.

* **Changed:** The function `tribits_add_library()` no longer sets the
  directory property `INCLUDE_DIRECTORIES` and instead passes include
  directory information between targets directly using the
  `INTERFACE_INCLUDE_DIRECTORIES` target property.  This may result in
  situations where include directories were getting set through the directory
  property `INCLUDE_DIRECTORIES` and getting picked up in the builds of
  targets that lacked the proper library dependencies.  Therefore, this is
  technically a break in backward compatibility, but not for well-formed
  TriBITS projects.

* **Changed:** Variables `${PACKAGE_NAME}_ENABLE_<depPkg>` are now set to
  `TRUE` for required upstream internal and external packages/TPLs `<depPkg>`
  (in order to simplify internal TriBITS logic).  However, a side-effect of
  this change is that CMake code that was ifed out with an `if
  (${PACKAGE_NAME}_ENABLE_<depPkg>)` statement (because that variable was not
  defined and therefore defaults to `FLASE`) for a required upstream
  dependency `<depPkg>` will now be enabled.  (This mistake can happen when an
  optional dependency `<depPkg>` is changed to a required dependency but the
  `if()` statements based on `${PACKAGE_NAME}_ENABLE_<depPkg>` are not
  removed.  Well-formed TriBITS projects and packages should not experience
  this problem or notice any change.)

## 2022-05-16:

* **Added:** The function `tribits_add_advanced_test(`) now correctly accepts
  a list of regexes for `PASS_REGULAR_EXPRESSION` and
  `FAIL_REGULAR_EXPRESSION` and behave the same as the raw CTest properties of
  the same names.

## 2022-03-10:

* **Changed:** The `tribits_add_advanced_test()` command now correctly reports
  unparsed/unrecognized arguments.  This will break some sloppy usages. (One
  TriBITS test had to be fixed.)

* **Changed:** The `tribits_add_test()` and `tribits_add_advanced_test()`
  behave differently with data passed in with explicit colon arguments.  For
  example, before `PASS_REGULAR_EXPRESSION "<regex0>;<regex1>"` could be used
  to pass in a list of regular expressions but with the new handling of
  function arguments, this now gets set as a single regex
  `"<regex0>\\;<regex1>"` which is not the same.  The fix (that also works
  with older versions of TriBITS) is to pass in multiple regexes as separate
  arguments as `PASS_REGULAR_EXPRESSION "<regex0>" "<regex1>"`.

* **Added:** The `tribits_add_test()`, `tribits_add_advanced_test()`, and
  `tribits_add_executable_and_test()` functions now allow handling semi-colons
  ';' in the quoted arguments to CMND/EXEC `ARGS` and `ENVIRONMENT` variable
  values by adding a `LIST_SEPARATOR <sep>` argument (same as for
  `ExternalProject_Add()`).

* **Changed:** The `tribits_add_test()` and `tribits_add_advanced_test()`
  functions switched over from using `cmake_parse_arguments(... ${ARGN})` to
  using `cmake_parse_arguments(PARSE_ARGV ...)` and, therefore, now will no
  longer ignore empty arguments.  This will break backward compatibility for
  cases where empty quoted arguments like `"${SomeVar}"` are passed in where
  the variable `SomeVar` is empty.

## 2022-03-02:

* **Added:** The project-level cache variable `<Project>_IMPORTED_NO_SYSTEM`
  was added to set the `IMPORTED_NO_SYSTEM` property (CMake versions 3.23+
  only) on the IMPORTED library targets in the installed
  `<Package>Config.cmake` files (see updated TriBITS users guide and build
  reference documentation for `<Project>_IMPORTED_NO_SYSTEM`).  Setting this
  to `ON` results in the include directories for this project's IMPORTED
  library targets to be listed on the compile lines in downstream CMake
  projects using `-I` instead of the default `-isystem` for IMPORTED library
  targets.  Setting this option to `ON` returns backward compatibility for the
  move to modern CMake targets which involved setting the include directories
  on the IMPORTED library targets using `target_include_directories()`
  described below (which changed the include directories from being listed as
  `-I` to `-isystem` by default).<br> **Workaround:** As a workaround for
  CMake versions less than 3.23, downstream CMake projects can set
  `CMAKE_NO_SYSTEM_FROM_IMPORTED=TRUE` in their CMake configure as described
  below.<br> For more details, see
  [TriBITSPub/TriBITS#443](https://github.com/TriBITSPub/TriBITS/issues/443).

## 2021-11-18:

* **Changed:** The default `<Project>_GENERATE_REPO_VERSION_FILE_DEFAULT` will
  be overridden to `OFF` if the 'git' executable cannot be found at configure
  time.  See updated TriBITS Developer's Guide documentation.

* **Changed:** The default value for `<Project>_ENABLE_Fortran` is set to
  `OFF` on WIN32 systems.  (Getting a Fortran compiler for native Windows is
  not typically very easy.)

## 2021-10-11

* **Changed:** The `<Package>Config.cmake` for each enabled package generated
  in the build directory tree have been moved from
  `<buildDir>/packages/<packageDir>/` to
  `<buildDir>/cmake_packages/<PackageName>/`.  (This makes it easy for
  `find_package(<PackageName>)` to find these files by simply adding the
  directory `<buildDir>/cmake_packages` to `CMAKE_PREFIX_PATH` and then
  `<Package>Config.cmake` for any enabled package will be found automatically
  found by CMake.)

* **Added/Changed:** Added the include directories for each library target
  with `target_include_directories()`.  This makes the usage of the variables
  `<Package>_INCLUDE_DIRS`, `<Package>_TPL_INCLUDE_DIRS`,
  `<Project>_INCLUDE_DIRS`, and `<Project>_TPL_INCLUDE_DIRS` unnecessary by
  downstream CMake projects. (See the changes to the
  `TribitsExampleApp/CmakeLists.txt` file that removed calls to
  `include_directories()` involving these variables.)  However, this change
  will also cause downstream CMake projects to pull in include directories as
  `SYSTEM` includes (e.g. using `-isystem` instead of `-I`) from IMPORTED
  library targets.  This changes how these include directories are searched
  and could break some fragile build environments that have the same header
  file names in multiple include directories searched by the compiler.
  Changing to `-isystem` will also silence any regular compiler warnings from
  headers found under these include directories.<br> ***Workarounds:*** One
  workaround for this is for the downstream CMake project to set the cache
  variable `CMAKE_NO_SYSTEM_FROM_IMPORTED=TRUE` which will restore the include
  directories for the IMPORTED library targets for the TriBITS project as
  non-SYSTEM include directories (i.e. `-I`) but it will also cause all
  include directories for all IMPORTED library targets to be non-SYSTEM
  (i.e. `-I`) even if they were being handled as SYSTEM include directories
  using `-isystem` before.  Therefore, that could still break the downstream
  project as it might change what header files are found for these other
  IMPORTED library targets and may expose many new warnings (which may have
  been silenced by their include directories being pulled in using
  `-isystem`).  The other workaround would be to clean up the list of include
  directories or delete some header files in those include directories so that
  only the correct header files can be found (regardless of the include
  directory search order).  For more details, see
  [TriBITSPub/TriBITS#443](https://github.com/TriBITSPub/TriBITS/issues/443).

## 2021-09-13

* **Removed:** Support for generation and installation of `Makefile.export.*`
  files has been removed along with the cache variable
  `<Project>_ENABLE_EXPORT_MAKEFILES`.  This is to allow the refactoring of
  TriBITS to use modern CMake targets that propagate all information and
  removing complex dependency tracking information from TriBITS (see
  TriBITSPub/TriBITS#63 and TriBITSPub/TriBITS#299).

## 2021-06-17

* **Added:** Added tool `tribits/python_utils/lower_case_cmake.py` and driver
  `tribits/refactoring/lower-case-cmake-tree.sh` that can make CMake command
  calls lower-case and make macro and function definition names lower case.
  This was applied to the entire TriBITS repository (and then minor
  modifications to fix a few issues).  This should not impact users of TriBITS
  but it does change the case of commands in the TriBITS error messages.  (See
  TriBITSPub/TriBITS#274)

## 2021-05-15

* **Changed/Fixed:** Changed so that all arguments passed to `ctest -S`
  command in 'dashboard' target are passed if their value is non-empty instead
  of `TRUE`.  Also, `CTEST_BUILD_NAME`, `TRIBITS_2ND_CTEST_DROP_SITE` and
  `TRIBITS_2ND_CTEST_DROP_LOCATION` can now be set and overridden in the CMake
  configure and will get passed to the `ctest -S` command by the 'dashboard'
  target.  This now allows setting these latter vars to `OFF` or `FALSE` which
  allows skipping a submit to a secondary CDash site if the underlying project
  is set up to submit to a secondary CDash site by default (see
  trilinos/Trilinos#9079).  These changes technically break backward
  compatibility because now setting some vars in the CMake configure will now
  get passed on to `ctest -S` command in the 'dashboard' target and some vars
  set in the env when calling 'make dashboard' will now be ignored if they
  were set in the CMake cache.  But many of these vars are only forwarded to
  the ctest -S command in the env if they are set to non-empty in the
  configure.  So existing processes that set `CTEST_BUILD_NAME` in the env
  when running `make dashboard` but not in the CMake configure will still
  behave the same (which is why no existing TriBITS tests had to change).  For
  these reasons, the chance of these changes causing a problem for most users
  should be very small and in fact it restores what most people would consider
  to be logical and useful behavior.

## 2021-03-12

* **Changed:** Upgrade minimum required CMake version from 3.10 to 3.17.  Existing
  TriBITS projects that have already upgraded to require CMake 3.17+ should not
  notice any major changes due to this change.

## 2020-11-12

* **Changed:** The default for `<Project>_ENABLE_EXPLICIT_INSTANTIATION` (ETI)
  was changed from `OFF` to `ON`.  This was turned on in practice for almost
  all users of Trilinos so while this change technically breaks backward
  compatibility of TriBITS, in practice, this will likely not impact any
  important customers.  (And new customers of Trilinos and related TriBITS
  projects will enjoy the benefits of ETI by default.  See
  trilinos/Trilinos#8130.)

## 2020-10-08

* **Changed:** Tests defined with `TRIBITS_ADD_TEST()` and
  `TRIBITS_ADD_ADVANCED_TEST()` with `NUM_MPI_PROCS > 1` will now not be
  added for non-MPI (`TPL_ENABLE_MPI=OFF`) configurations.  Before, the test
  would get added and basically ignore the value of `NUM_MPI_PROCS`.  This
  change together with the change to move to using `add_test(NAME <name>
  COMMAND <command>)` mostly restores backward compatibility for projects
  using TriBITS.  For more details, see trilinos/Trilinos#8110.

## 2020-09-28

* **Changed/Fixed:** Tests defined with `TRIBITS_ADD_TEST()` and
  `TRIBITS_ADD_ADVANCED_TEST()` have been updated from using
  `add_test(<name> <command>)` to using `add_test(NAME <name> COMMAND
  <command>)`.  This now causes CMake to error-out if there is more than one
  test with the same name in the same directory.  Before, CMake would allow
  it (but the behavior in that case was undocumented and undefined).  For
  more details, see trilinos/Trilinos#8110.

## 2020-06-16

* **Added/Deprecated:** The variables `<Package>_FORTRAN_COMPILER` and
  `<Package>_FORTRAN_FLAGS` in the `<Package>Config.cmake` files (in the build
  dir and the install dir) are deprecated.  Please use the new variables
  `<Package>_Fortran_COMPILER` and `<Package>_Fortran_FLAGS`.  (Justification:
  In raw CMake, the compiler is called 'Fortran', not 'FORTRAN'.  Also, the
  name 'Fortran' is already used in the `<Project>Config.cmake` file.)

## 2020-04-16

* **Changed/Removed:** TriBITS Core: `CMAKE_CXX_STANDARD` is now always set to
  at least `11` to enable C++11 support.  The `${PROJECT_NAME}_ENABLE_CXX11`
  option has been removed, but a variable of the same name is hard-coded to ON
  for compatibility.

* **Changed/Removed:** TriBITS Core: `CMAKE_CXX_FLAGS` will no longer contain
  an explicit C++ standard option like `-std=c++11`.  TriBITS clients'
  dependents using CMake will still compile with support for C++11 language
  constructs because libraries come with `INTERFACE_COMPILE_FEATURES`
  including `cxx_std_11`.  However, those using `Makefile.export.<Package>`
  files will no longer get a `-std=c++11` flag.

## 2018-10-10

* **Changed:** TriBITS Core: Changed minimum CMake version from 2.8.11 to
  3.10.0.

## 2017-09-25

* **Added:** TriBITS CTest Driver: Added cache and env vars
  `CTEST_SUBMIT_RETRY_COUNT` and `CTEST_SUBMIT_RETRY_DELAY` to allow the
  number of `ctest_submit()` submit attempts to retry and how long to pause
  between retries.  Before, these were hard-coded to 25 and 120 respectively,
  which means that something like a MySQL insertion error could consume as
  much as 50 minutes before moving on!  The new defaults are set at 5 retries
  with a 3 sec delay (which appear to be the CTest defaults).

## 2017-09-30

* **Added:** TriBITS Core: Added `TEST_<IDX> COPY_FILES_TO_TEST_DIR` block for
    `TRIBITS_ADD_ADVANCED_TEST()`.  This was added in such a way so to avoid
    clashing with existing usages of the script (since the new arguments
    `SOURCE_DIR` and `DEST_DIR` are only parsed if `COPY_FILES_TO_TEST_DIR` is
    listed in the `TEST_<IDX> block`.

## 2017-09-05

* **Changed:** TriBITS Core: Un-parsed and otherwise ignored arguments to many
  TriBITS functions and macros are now flagged (see developers guide
  documentation for `${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`).  The
  default value is `WARNING` which results in simply printing a warning but
  allow configure to complete.  This allows one to see the warnings but for
  the project to continue to work as before.  But this can be changed to
  `SEND_ERROR` or `FATAL_ERROR` that will fail the configure.

## 2017-06-24

* **Added:** TriBITS CTest Driver: Add new all-at-once mode for `ctest -S`
  driver scripts using `TRIBITS_CTEST_DRIVER()` by setting the variable
  `${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE=TRUE`.  This works with older versions
  of CMake/CTest and CDash but, by default, will just return a single glob of
  results, not breaking-out results on a package-by-package basis.  Therefore,
  this is disabled by default and package-by-package mode is used by default.
  But if `${PROJECT_NAME}_CTEST_USE_NEW_AAO_FEATURES=TRUE` is set, then
  TriBITS will take advantage of new CMake, CTest, and CDash features
  (currently on a branch) to display the results on CDash broken down
  package-by-package.  Once these changes are merged to the CMake/CTest and
  CDash 'master' branches, then the default for
  `${PROJECT_NAME}_CTEST_USE_NEW_AAO_FEATURES` will be set to `TRUE`
  automatically when it detects an updated version of CMake/CTest is present.
  In the future, at some point, the TriBITS default for
  `${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE` will change from `FALSE` to `TRUE`
  since that is a much more efficient way to drive automated testing.

## 2017-05-25

* **Added/Deprecated:** TriBITS Core: The usage of `PARSE_ARGUMENTS()` has
  been deprecated and replaced by `CMAKE_PARSE_ARGUMENTS()` everywhere in
  TriBITS. Any call to `PARSE_ARGUMENTS()` will warn users and tell them to
  use `CMAKE_PARSE_ARGUMENTS()` instead.

## 2017-05-17

* **Changed:** TriBITS Core: TriBITS now unconditionally sets
  `${PROJECT_NAME}_ENABLE_Fortran_DEFAULT` to `ON`.  Projects will now need to
  put in special logic to set to `OFF` or `ON` for certain platforms.

## 2017-01-11

* **Changed/Fixed:** TriBITS Core: TriBITS now correctly sets the default
  value for DART_TESTING_TIMEOUT to 1500 seconds and will scale it by
  `<Project>_SCALE_TEST_TIMEOUT` even if `DART_TESTING_TIMEOUT` is not
  explicitly set.

## 2016-12-07

* **Removed:** TriBITS Core: The long deprecated variable
  `${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE` has been removed.  Upgrading
  existing TriBITS projects just requires a simple string replacement of
  `_ENABLE_SECONDARY_STABLE_CODE` with `_ENABLE_SECONDARY_TESTED_CODE` in all
  files.  Since Trilinos has turned on ST code by default and many other
  TriBITS projects don't differentiate between PT and ST code, this change
  should not even break those projects, even if they don't update.

## 2016-11-02

* **Added/Changed/Removed:** TriBITS Python Utils: gitdist now accepts
  `--dist-repos` and `--dist-not-repos` arguments and requires that the base
  repo '.' be explicitly listed in the `.gitdist[.default]` files and in
  `--dist-repos`.  The arguments `--dist-extra-repos`,
  `--dist-not-extra-repos` and `--dist-not-base-repo` are not longer
  supported.  See `gitdist --help` for more details.

* **Changed:** TriBITS projects now install with full RPATH set by default
  (see "Setting install RPATH" in build reference guide).

## 2016-10-22

* **Changed:** TriBITS Core: `TRIBITS_ADD_TEST()` argument for
  `FAIL_REGULAR_EXPRESSION` now works when circular RCP detection is enabled.
  This is technically a break in backward compatibility since now that
  argument will not be ignored and any tests that specified this may change
  behavior.

* **Changed:** TriBITS Core: `TRIBITS_ADD_ADVANCED_TEST()` block `TEST_<IDX>`
  argument `FAIL_REGULAR_EXPRESSION` now works.  Before, it was just being
  ignored.  This is technically a break in backward compatibility since now
  that argument will not be ignored and any tests that specified this may
  change behavior.

* **Added/Changed:** TriBITS Core: Added `TRIBITS_ADD_ADVANCED_TEST()` block
  `TEST_<IDX>` option `WILL_FAIL` that has the same behavior as the built-in
  CTest option `WILL_FAIL`.  Note that this technically can break backward
  compatibility since `WILL_FAIL` may have been interpreted to be a value from
  another argument and now will not.

## 2016-01-22

* **Added/Deprecated:** TriBITS Core: Change test category `WEEKLY` to `HEAVY`
  and depreciate `WEEKLY`.  You can still use `WEEKLY` but it will result in a
  lot of warnings.

## 2015-12-03

* **Added:** TriBITS CI Support: `checkin-test.py`: Added support for tracking
  branches for each repo independently and not assume 'origin' and not assume
  that all of the repos are on the same branch or will be pulling and pushing
  to the same remote branch.  This will make it easier to use the
  `checkin-test.py` script to set up various integration scenarios.  See
  TriBITSPub/TriBITS#15 for details.

## 2015-04-14

* MAJOR: TriBITS Core: When configuring with
    ${PROJECT_NAME}_ENABLE_CXX11=ON, if C++11 support cannot be verified, then
    the configure will fail hard right away.  Before, TriBITS would disable
    C++11 support and continue.

## 2014-11-22

* **Added:** TriBITS Core: Added `${PROJECT_NAME}_TRACE_ADD_TEST`: Now you can
  print a single line to STDOUT if a test got added (and its important
  properties) or not and if not then why the test did not get added.

## 2014-09-22

* **Changed:** TriBITS Core: Changed minimum version of CMake from 2.7 to
  2.8.11.

## 2014-09-21

* **Added:** TriBITS Dashboard Driver: Added support for the env var
  `TRIBITS_TDD_USE_SYSTEM_CTEST` so that if equal to `1`, then the TriBITS
  Dashboard Driver (TDD) system will use the CTest (and CMake) in the env will
  be used instead of being downloaded using `download-cmake.py`.  This not
  only speeds up the automated builds, but it also ensures that the automated
  testing uses exactly the install of CMake/CTest that is used by the
  developers on the system.  Also, it has been found that `download-cmake.py`
  will download and install a 32bit version even on 64bit machines.

<!--

LocalWords: Fortran cmake CMake CMAKE ctest CTest CTEST CDash
LocalWords:  MAKEFILES refactoring tribits TriBITS TriBITSPub trilinos Trilinos
LocalWords: INSTANTIATION STDOUT gitdist

-->
