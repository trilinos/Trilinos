# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include("${CMAKE_CURRENT_LIST_DIR}/../test_support/TribitsAddExecutableTestHelpers.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../test_support/TribitsAddTestHelpers.cmake")
include(TribitsCommonArgsHelpers)
include(TribitsGeneralMacros)
include(TribitsLibIsTestOnly)
include(TribitsReportInvalidTribitsUsage)
include(TribitsDeprecatedHelpers)

include(PrintVar)
include(AppendSet)
include(CMakeParseArguments)


# @FUNCTION: tribits_add_executable()
#
# Function used to create an executable (typically for a test or example),
# using the built-in CMake command ``add_executable()``.
#
# Usage::
#
#   tribits_add_executable(
#     <exeRootName>  [NOEXEPREFIX]  [NOEXESUFFIX]  [ADD_DIR_TO_NAME]
#     SOURCES <src0> <src1> ...
#     [CATEGORIES <category0>  <category1> ...]
#     [HOST <host0> <host1> ...]
#     [XHOST <host0> <host1> ...]
#     [HOSTTYPE <hosttype0> <hosttype1> ...]
#     [XHOSTTYPE <hosttype0> <hosttype1> ...]
#     [EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...]
#     [DIRECTORY <dir>]
#     [TESTONLYLIBS <lib0> <lib1> ...]
#     [IMPORTEDLIBS <lib0> <lib1> ...]
#     [COMM [serial] [mpi]]
#     [LINKER_LANGUAGE (C|CXX|Fortran)]
#     [TARGET_DEFINES -D<define0> -D<define1> ...]
#     [INSTALLABLE]
#     [ADDED_EXE_TARGET_NAME_OUT <exeTargetName>]
#     )
#
# *Sections:*
#
# * `Formal Arguments (tribits_add_executable())`_
# * `Executable and Target Name (tribits_add_executable())`_
# * `Additional Executable and Source File Properties (tribits_add_executable())`_
# * `Install Target (tribits_add_executable())`_
#
# .. _Formal Arguments (tribits_add_executable()):
#
# **Formal Arguments (tribits_add_executable())**
#
#   ``<exeRootName>``
#
#     The root name of the executable (and CMake target) (see `Executable and
#     Target Name (tribits_add_executable())`_).  This must be the first
#     argument.
#
#   ``NOEXEPREFIX``
#
#     If passed in, then ``${PACKAGE_NAME}_`` is not added the beginning of
#     the executable name (see `Executable and Target Name
#     (tribits_add_executable())`_).
#
#   ``NOEXESUFFIX``
#
#     If passed in, then ``${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}`` and
#     not added to the end of the executable name (except for native Windows
#     builds, see `Executable and Target Name (tribits_add_executable())`_).
#
#   ``ADD_DIR_TO_NAME``
#
#     If passed in, the directory path relative to the package's base
#     directory (with "/" replaced by "_") is added to the executable name
#     (see `Executable and Target Name (tribits_add_executable())`_).  This
#     provides a simple way to create unique test executable names inside of a
#     given TriBITS package.  Only test executables in the same directory
#     would need to have unique ``<execRootName>`` passed in.
#
#   ``SOURCES <src0> <src1> ...``
#
#     Gives the source files that will be compiled into the built executable.
#     By default, these sources are assumed to be in the current working
#     directory (or can contain the relative path or absolute path).  If
#     ``<srci>`` is an absolute path, then that full file path is used.  This
#     list of sources (with adjusted directory path) are passed into
#     ``add_executable(<exeTargetName> ... )``.  After calling this function,
#     the properties of the source files can be altered using the built-in
#     CMake command ``set_source_file_properties()``.
#
#   ``DIRECTORY <dir>``
#
#     If specified, then the generated executable ``<exeTargetName>`` is
#     placed in the relative or absolute directory ``<dir>``.  If ``<dir>`` is
#     not an absolute path, then the generated executable is placed in the
#     directory ``${CMAKE_CURRENT_BINARY_DIR}/<dir>/``.  Also, the sources for
#     the executable listed in ``SOURCES <src0> <src1> ...`` are assumed to be
#     in the relative or absolute directory ``<dir>`` instead of the current
#     source directory.  This directory path is prepended to each source file
#     name ``<srci>`` unless ``<srci>`` is an absolute path.  If ``<dir>`` is
#     not an absolute path, then source files listed in ``SOURCES`` are
#     assumed to be in the directory ``${CMAKE_CURRENT_SOURCE_DIR}/<dir>/``.
#
#   ``CATEGORIES <category0> <category1> ...``
#
#     Gives the `Test Test Categories`_ for which this test will be added.
#     See `tribits_add_test()`_ for more details.
#
#   ``HOST <host0> <host1> ...``
#
#     The list of hosts for which to enable the test (see `tribits_add_test()`_).
#
#   ``XHOST <host0> <host1> ...``
#
#     The list of hosts for which **not** to enable the test (see
#     `tribits_add_test()`_).
#
#   ``HOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which to enable the test (see
#     `tribits_add_test()`_).
#
#   ``XHOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which **not** to enable the test (see
#     `tribits_add_test()`_).
#
#   ``EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...``
#
#     If specified, gives the names of CMake variables that must evaluate to
#     true, or the test will not be added (see `tribits_add_test()`_).
#
#   ``TESTONLYLIBS <lib0> <lib1> ...``
#
#     Specifies extra test-only libraries defined in this CMake project that
#     will be linked to the executable using ``target_link_libraries()``.
#     Note that regular libraries (i.e. not ``TESTONLY``) defined in the
#     current package or any upstream packages can *NOT* be listed!  TriBITS
#     automatically links non ``TESTONLY`` libraries in this package and
#     upstream packages to the executable.  The only libraries that should be
#     listed in this argument are either ``TESTONLY`` libraries.
#
#   ``IMPORTEDLIBS <lib0> <lib1> ...``
#
#     Specifies extra external libraries that will be linked to the executable
#     using ``target_link_libraries()``.  This can only be used for libraries
#     that are built external from this CMake project and are not provided
#     through a proper `TriBITS TPL`_.  The latter usage of passing in
#     external libraries is not recommended.  External libraries should be
#     handled as declared `TriBITS TPLs`_.  So far, the only case where
#     ``IMPORTEDLIBS`` has been shown to be necessary is to pass in the
#     standard C math library ``m``.  In every other case, a TriBITS TPL
#     should be used instead.
#
#   ``COMM [serial] [mpi]``
#
#     If specified, selects if the test will be added in serial and/or MPI
#     mode.  See the ``COMM`` argument in the script
#     `tribits_add_test()`_ for more details.
#
#   ``LINKER_LANGUAGE (C|CXX|Fortran)``
#
#     If specified, overrides the linker language used by setting the built-in
#     CMake target property ``LINKER_LANGUAGE``.  TriBITS sets the default
#     linker language as follows::
#
#       if (${PROJECT_NAME}_ENABLE_CXX)
#         set(LINKER_LANGUAGE CXX)
#       elseif (${PROJECT_NAME}_ENABLE_C)
#         set(LINKER_LANGUAGE C)
#       else()
#         # Let CMake set the default linker language it wants based
#         # on source file extensions passed into ``add_executable()``.
#       endif()
#
#     The reason for this logic is that on some platform if you have a Fortran
#     or C main that links to C++ libraries, then you need the C++ compiler to
#     do the final linking.  CMake does not seem to automatically know that it
#     is pulling in C++ libraries and therefore needs to be told use C++ for
#     linking.  This is the correct default behavior for mixed-language
#     projects. However, this argument allows the developer to override this
#     logic and use any linker language desired based on other considerations.
#
#   ``TARGET_DEFINES -D<define0> -D<define1> ...``
#
#     Add the listed defines using
#     ``target_compile_definitions(<exeTargetName> ...)``.  These should only
#     affect the listed sources for the built executable and not other
#     targets.
#
#   ``INSTALLABLE``
#
#     If passed in, then an install target will be added to install the built
#     executable into the ``${CMAKE_INSTALL_PREFIX}/bin/`` directory (see
#     `Install Target (tribits_add_executable())`_).
#
#   ``ADDED_EXE_TARGET_NAME_OUT <exeTargetName>``
#
#     If specified, then on output the variable ``<exeTargetName>`` will be
#     set with the name of the executable target passed to
#     ``add_executable(<exeTargetName> ... )``.  Having this name allows the
#     calling ``CMakeLists.txt`` file access and set additional target
#     properties (see `Additional Executable and Source File Properties
#     (tribits_add_executable())`_).
#
# .. _Executable and Target Name (tribits_add_executable()):
#
# **Executable and Target Name (tribits_add_executable())**
#
# By default, the full name of the executable and target name
# is::
#
#   <exeTargetName> = ${PACKAGE_NAME}_<exeRootName>
#
# If ``ADD_DIR_TO_NAME`` is set, then the directory path relative to the
# package base directory (with "/" replaced with "_"), or ``<relDirName>``, is
# added to the executable name to form::
#
#   <exeTargetName> = ${PACKAGE_NAME}_<relDirName>_<exeRootName>
#
# If the option ``NOEXEPREFIX`` is passed in, then the prefix
# ``${PACKAGE_NAME}_`` is removed.
#
# The executable suffix ``${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}`` will be
# added to the actual executable file name if the option ``NOEXESUFFIX`` is
# *not* passed in but this suffix is never added to the target name.
# (However, note that on Windows platforms, the default ``*.exe`` extension is
# always added because windows will not run an executable in many contexts
# unless it has the ``*.exe`` extension.)
#
# The reason that a default prefix is prepended to the executable and target
# name is because the primary reason to create an executable is typically to
# create a test or an example that is private to the package.  This prefix
# helps to namespace the executable and its target so as to avoid name clashes
# with targets in other packages.  It also helps to avoid clashes if the
# executable gets installed into the install directory (if ``INSTALLABLE`` is
# specified).  For general utility executables on Linux/Unix systems,
# ``NOEXEPREFIX`` and ``NOEXESUFFIX`` should be passed in.  In this case, one
# must be careful to pick ``<exeRootName>`` that will be sufficiently globally
# unique.  Please use common sense when picking non-namespaced names.
#
# .. _Additional Executable and Source File Properties (tribits_add_executable()):
#
# **Additional Executable and Source File Properties (tribits_add_executable())**
#
# Once ``add_executable(<exeTargetName> ... )`` is called and this function
# exists, one can set and change properties on the ``<exeTargetName>``
# executable target using the built-in ``set_target_properties()`` command as
# well as properties on any of the source files listed in ``SOURCES`` using
# the built-in ``set_source_file_properties()`` command just like in any CMake
# project.  IF the executable is added, its name will be returned by the
# argument ``ADDED_EXE_TARGET_NAME_OUT <exeTargetName>``.  For example::
#
#   tribits_add_executable( someExe ...
#     ADDED_EXE_TARGET_NAME_OUT  someExe_TARGET_NAME )
#
#   if (someExe_TARGET_NAME)
#     set_target_properties( ${someExe_TARGET_NAME}
#       PROPERTIES  LINKER_LANGUAGE  CXX )
#   endif()
#
# The ``if(someExe_TARGET_NAME)`` is needed in case the executable does not
# get added for some reason (see `Formal Arguments
# (tribits_add_executable())`_ for logic that can result in the executable
# target not getting added).
#
# .. _Install Target (tribits_add_executable()):
#
# **Install Target (tribits_add_executable())**
#
# If ``INSTALLABLE`` is passed in, then an install target using the built-in
# CMake command ``install(TARGETS <exeTargetName> ...)`` is added to install
# the built executable into the ``${CMAKE_INSTALL_PREFIX}/bin/`` directory
# (actual install directory path is determined by
# ``${PROJECT_NAME}_INSTALL_RUNTIME_DIR``, see `Setting the install prefix`_).
#
function(tribits_add_executable EXE_NAME)

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("")
    message("TRIBITS_ADD_EXECUTABLE: ${EXE_NAME} ${ARGN}")
  endif()

  tribits_add_executable_assert_correct_call_context()

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    "NOEXEPREFIX;NOEXESUFFIX;ADD_DIR_TO_NAME;INSTALLABLE"
    #one_value_keywords
    ""
    #multi_value_keywords
    "SOURCES;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;EXCLUDE_IF_NOT_TRUE;DIRECTORY;TESTONLYLIBS;IMPORTEDLIBS;DEPLIBS;COMM;LINKER_LANGUAGE;TARGET_DEFINES;DEFINES;ADDED_EXE_TARGET_NAME_OUT"
    ${ARGN}
    )

  tribits_check_for_unparsed_arguments()

  # Executable not added by default!
  if(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    set(${PARSE_ADDED_EXE_TARGET_NAME_OUT} "" PARENT_SCOPE)
  endif()

  set(EXE_BINARY_NAME ${EXE_NAME})
  tribits_add_executable_adjust_exe_name(EXE_BINARY_NAME)

  tribits_add_executable_is_skipped(skipAddExecutable)
  if (skipAddExecutable)
    return()
  endif()

  tribits_add_executable_get_adjusted_sources_list(EXE_SOURCES)

  tribits_add_executable_assert_testonlylibs()

  tribits_add_executable_assert_importedlibs()

  tribits_add_executable_convert_from_deplibs()

  #
  # B) Add the executable and set its properties
  #

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("TRIBITS_ADD_EXECUTABLE: add_executable(${EXE_BINARY_NAME} ${EXE_SOURCES})")
  endif()
  add_executable(${EXE_BINARY_NAME} ${EXE_SOURCES})
  append_global_set(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${EXE_BINARY_NAME})

  if(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    set(${PARSE_ADDED_EXE_TARGET_NAME_OUT} ${EXE_BINARY_NAME} PARENT_SCOPE)
  endif()

  if (PARSE_DEFINES)
    tribits_deprecated("Passing extra defines through 'DEFINES' ${PARSE_DEFINES}"
      " is deprecated.  Instead, pass them through 'TARGET_DEFINES'.  The 'DEFINES'"
      " argument was incorrectly implemented by calling add_definitions() which has"
      " directory scope and not function scope as was documented.  This resulted in"
      " confusing behavior.  If one wishes to set defines at the directly level,"
      " just call add_definitions() directly.")
    add_definitions(${PARSE_DEFINES})
  endif()

  if (PARSE_TARGET_DEFINES)
    target_compile_definitions(${EXE_BINARY_NAME} PUBLIC ${PARSE_TARGET_DEFINES})
  endif()

  if(PARSE_NOEXESUFFIX AND NOT WIN32)
    set_target_properties(${EXE_BINARY_NAME} PROPERTIES SUFFIX "")
  else()
    set_target_properties(${EXE_BINARY_NAME} PROPERTIES SUFFIX
      ${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX})
  endif()

  tribits_set_linker_language_from_arg( ${EXE_BINARY_NAME}
    "${PARSE_LINKER_LANGUAGE}" )

  assert_defined(${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
  if (${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
    #message("${EXE_BINARY_NAME}: Adding property LINK_SEARCH_START_STATIC")
    set_property(TARGET ${EXE_BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
  endif()

  if(PARSE_DIRECTORY)
    set_target_properties( ${EXE_BINARY_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${PARSE_DIRECTORY} )
  endif()

  set_property(TARGET ${EXE_BINARY_NAME} APPEND PROPERTY
    LABELS ${PACKAGE_NAME}Exes ${PARENT_PACKAGE_NAME}Exes)

  #
  # C) Link ${EXE_BINARY_NAME} to direct upstream libraries
  #

  target_link_libraries(${EXE_BINARY_NAME} PUBLIC ${${PACKAGE_NAME}_LIBRARIES})
  foreach(depPkg IN LISTS ${PACKAGE_NAME}_LIB_ENABLED_DEPENDENCIES
      ${PACKAGE_NAME}_TEST_ENABLED_DEPENDENCIES
    )
    target_link_libraries(${EXE_BINARY_NAME} PUBLIC ${depPkg}::all_libs)
  endforeach()
  foreach(testOnlyLib ${PARSE_TESTONLYLIBS})
    target_link_libraries(${EXE_BINARY_NAME} PUBLIC
      "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${testOnlyLib}")
  endforeach()

  #
  # D) Install if asked
  #

  if(${PROJECT_NAME}_INSTALL_EXECUTABLES AND PARSE_INSTALLABLE)
    install(
      TARGETS ${EXE_BINARY_NAME}
      EXPORT ${PROJECT_NAME}
        DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}
      COMPONENT ${PACKAGE_NAME}
    )
  endif()

endfunction()


# Assert tribits_add_executable() is called in the correct context
#
# NOTE: This read the variables from the enclosing tribits_add_executable()
# function call scope.
#
function(tribits_add_executable_assert_correct_call_context)

  if (CURRENTLY_PROCESSING_SUBPACKAGE)

    # This is a subpackage being processed

    if(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_subpackage() before tribits_add_executable()"
        " in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

    if(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_executable() before "
        " tribits_subpackage_postprocess() in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

  else()

    # This is a package being processed

    if(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_package() or tribits_package_decl() before"
        " tribits_add_executable() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

    if(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_executable() before "
        " tribits_package_postprocess() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

  endif()

endfunction()


# Modify EXE_BINARY_NAME based in passed-in args
#
# NOTE: This read the variables from the enclosing tribits_add_executable()
# function call scope.
#
function(tribits_add_executable_adjust_exe_name  exeBinaryNameInOut)

  set(exeBinaryName ${${exeBinaryNameInOut}})
  if (PARSE_ADD_DIR_TO_NAME)
    set(dirName "")
    tribits_create_name_from_current_source_directory(dirName)
    set(exeBinaryName ${dirName}_${exeBinaryName})
  endif()

  if (DEFINED PACKAGE_NAME AND NOT PARSE_NOEXEPREFIX)
    set(exeBinaryName ${PACKAGE_NAME}_${exeBinaryName})
  endif()
  set(${exeBinaryNameInOut} ${exeBinaryName} PARENT_SCOPE)

endfunction()


# Check if to skip adding the executable based on different criteria
#
# NOTE: This read the variables from the enclosing tribits_add_executable()
# function call scope.
#
function(tribits_add_executable_is_skipped   skipAddExecutableOut)

  set(skipAddExecutable FALSE)

  set(ADD_THE_TEST FALSE)
  set(TEST_NAME ${EXE_NAME})  # For error message
  tribits_add_test_process_categories(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    set(skipAddExecutable TRUE)
  endif()
  set(TEST_NAME "")

  set(ADD_THE_TEST FALSE)
  tribits_add_test_process_host_hosttype(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    set(skipAddExecutable TRUE)
  endif()

  tribits_process_comm_args(ADD_SERIAL_EXE  ADD_MPI_EXE  ${PARSE_COMM})
  if (NOT ADD_SERIAL_EXE AND NOT ADD_MPI_EXE)
    set(skipAddExecutable TRUE)
  endif()

  if (${EXE_BINARY_NAME}_EXE_DISABLE)
    message("-- "
      "${EXE_BINARY_NAME} EXE NOT being built due to ${EXE_BINARY_NAME}_EXE_DISABLE="
      "'${${EXE_BINARY_NAME}_EXE_DISABLE}'")
    set(skipAddExecutable TRUE)
  endif()

  set(${skipAddExecutableOut} ${skipAddExecutable} PARENT_SCOPE)

endfunction()


# Get adjusted list of source files
#
# NOTE: This read the variables from the enclosing tribits_add_executable()
# function call scope.
#
function(tribits_add_executable_get_adjusted_sources_list  exeSourcesOut)

  set(exeSources "")
  if(PARSE_DIRECTORY )
    foreach( srcFile ${PARSE_SOURCES} )
      if(IS_ABSOLUTE ${srcFile})
        list(APPEND exeSources "${srcFile}")
      else()
        list(APPEND exeSources "${PARSE_DIRECTORY}/${srcFile}")
      endif()
    endforeach( )
  else()
    foreach( srcFile ${PARSE_SOURCES} )
      list(APPEND exeSources "${srcFile}")
    endforeach( )
  endif()

  set(${exeSourcesOut} ${exeSources} PARENT_SCOPE)

endfunction()


# Assert tribits_add_executable() TESTONLYLIBS
#
# NOTE: This read the variables from the enclosing tribits_add_executable()
# function call scope.
#
function(tribits_add_executable_assert_testonlylibs)
  # Assert that TESTONLYLIBS only contains TESTONLY libs!
  foreach(testOnlyLib ${PARSE_TESTONLYLIBS})
    set(prefixedTestOnlyLib "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${testOnlyLib}")
    tribits_lib_is_testonly(${prefixedTestOnlyLib} libIsTestOnlyLib)
    if (NOT libIsTestOnlyLib)
      message(FATAL_ERROR "ERROR: '${testOnlyLib}' in TESTONLYLIBS not a TESTONLY lib!"
        "  If this a regular library in this package or in an dependent upstream"
        " package then TriBITS will link automatically to it.  If you remove this and it"
        " does not link, then you need to add a new package dependency to"
        " this package's dependencies file"
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    elseif(PARSE_INSTALLABLE)
      message(FATAL_ERROR "ERROR: TESTONLY lib '${testOnlyLib}' not allowed with"
        " INSTALLABLE executable!  An INSTALLABLE executable can only depend on"
        " non-TESTONLY libraries!  Otherwise, when shared libs are used, and"
        " TESTONLY library would not be installed and the installed executable"
        " would be unusable!" )
    endif()
  endforeach()
endfunction()


# Assert tribits_add_executable() IMPORTEDLIBS
#
# NOTE: This read the variables from the enclosing tribits_add_executable()
# function call scope.
#
function(tribits_add_executable_assert_importedlibs)

  # Assert that IMPORTEDLIBS are not TESTONLY libs are not regular package
  # libs!
  foreach(importedLib ${PARSE_IMPORTEDLIBS})
    set(prefixedImportedLib "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${importedLib}")
    tribits_lib_is_testonly(${prefixedImportedLib} importedLibIsTestOnlyLib)
    if (libIsTestOnly)
      message(FATAL_ERROR
        "ERROR: Lib '${importedLib}' being passed through"
        " IMPORTEDLIBS is not allowed to be a TESTONLY lib!"
        "  Use TESTONLYLIBS instead!" )
    endif()
    list(FIND ${PACKAGE_NAME}_LIBRARIES "${PACKAGE_NAME}::${prefixedImportedLib}"
      foundPrefixedImportedLibInPkgLibs_idx)
    if (NOT foundPrefixedImportedLibInPkgLibs_idx EQUAL -1)
      message(FATAL_ERROR
        "ERROR: Lib '${importedLib}' in IMPORTEDLIBS is in"
        " this package and is *not* an external lib!"
        "  TriBITS takes care of linking against libs the current"
        " package automatically.  Please remove '${importedLib}' from IMPORTEDLIBS!")
    elseif (TARGET ${prefixedImportedLib})
      message(FATAL_ERROR
        "ERROR: Lib '${importedLib}' being passed through"
        " IMPORTEDLIBS is *not* an external library but instead is a library"
        " defined in this CMake project!"
        "  TriBITS takes care of linking against libraries in dependent upstream"
        " packages.  If you want to link to a library in an upstream"
        " package then add the package name for that library to the appropriate"
        " list in this package's dependencies file"
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    endif()
  endforeach()

endfunction()


# Convert from tribits_add_executable() DEPLIBS to IMPORTEDLIBS and TESTONLYLIBS
#
# NOTE: This is a macro as it updates local variables in the
# tribits_add_executable() scope!
#
macro(tribits_add_executable_convert_from_deplibs)

  # Convert from old DEPLIBS to TESTONLYLIBS and IMPORTEDLIBS
  foreach(depLib ${PARSE_DEPLIBS})
    set(prefixedDepLib "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${depLib}")
    tribits_lib_is_testonly(${prefixedDepLib} depLibIsTestOnlyLib)
    if (depLibIsTestOnlyLib)
      tribits_deprecated("Passing TESTONLY lib '${depLib}' through DEPLIBS"
        " is deprecated!  Instead, please pass through TESTONLYLIBS instead!"
        "  DEPLIBS is deprecated!")
      list(APPEND PARSE_TESTONLYLIBS ${depLib})
    elseif (TARGET ${prefixedDepLib})
      tribits_deprecated("Passing non-TESTONLY lib '${depLib}' through DEPLIBS"
      " is deprecated!  The library '${depLib}' appears to be a"
      " library defined in this CMake project."
      "  TriBITS takes care of linking against libraries in dependent upstream"
      " packages.  Therefore, please remove '${depLib}' from this list."
      "   If you want to link to a library from an upstream"
      " package, then add the package name to the appropriate category"
      " in this package's dependencies file: "
      " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
      # ToDo: Convert the above 'WARNING' to 'SEND_ERROR'
    else()
      message(WARNING "WARNING: Passing external lib '${depLib}' through"
        " DEPLIBS is deprecated!  Instead, pass through IMPORTEDLIBS!"
        "  DEPLIBS is deprecated!"
        "  Please note that only external libs are allowed to be passed through"
        " IMPORTEDLIBS.")
      list(APPEND PARSE_IMPORTEDLIBS ${depLib})
    endif()
  endforeach()

endmacro()
