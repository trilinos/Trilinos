# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER


include(TribitsAddExecutableTestHelpers)
include(TribitsCommonArgsHelpers)
include(TribitsAddTestHelpers)
include(TribitsGeneralMacros)
include(TribitsReportInvalidTribitsUsage)

include(PrintVar)
include(AppendSet)
include(CMakeParseArguments)

###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of the file
### TribitsPackageMacros.cmake before trying to change anything in this file!
###


#
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
#     will be linked to the executable using ``target_link_libraries()``.  Note
#     that regular libraries (i.e. not ``TESTONLY``) defined in the current SE
#     package or any upstream SE packages can *NOT* be listed!  TriBITS
#     automatically links non ``TESTONLY`` libraries in this package and
#     upstream packages to the executable.  The only libraries that should be
#     listed in this argument are either ``TESTONLY`` libraries.  The include
#     directories for each test-only library will automatically be added
#     using::
#
#       include_directories(${<libi>_INCLUDE_DIRS})
#
#     where ``<libi>_INCLUDE_DIRS`` was set by::
#
#       tribits_add_library(<libi> ... TESTONLY ...)
#
#     Therefore, to link to a defined ``TESTONLY`` library in any upstream
#     enabled package, one just needs to pass in the library name through
#     ``TESTONLYLIBS ... <libi> ...`` and that is it!
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
  
  #
  # Confirm that package and subpackage macros/functions have been called inteh correct order
  #
  
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

  if(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    set(${PARSE_ADDED_EXE_TARGET_NAME_OUT} PARENT_SCOPE)
  endif()

  #
  # B) Exclude building the test executable based on some criteria
  #

  set(ADD_THE_TEST FALSE)
  set(TEST_NAME ${EXE_NAME})  # For error message
  tribits_add_test_process_categories(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    return()
  endif()
  set(TEST_NAME)

  set(ADD_THE_TEST FALSE)
  tribits_add_test_process_host_hosttype(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    return()
  endif()

  tribits_process_comm_args(ADD_SERIAL_EXE  ADD_MPI_EXE  ${PARSE_COMM})
  if (NOT ADD_SERIAL_EXE AND NOT ADD_MPI_EXE)
    return()
  endif()

  #
  # C) Add the executable
  #

  set(LIBRARY_NAME_PREFIX "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}")

  if (NOT TRIBITS_ADD_EXECUTABLE_UNIT_TESTING)
    tribits_include_directories(REQUIRED_DURING_INSTALLATION_TESTING
      ${${PACKAGE_NAME}_INCLUDE_DIRS})
#    set_property(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
#      ${${PACKAGE_NAME}_LIBRARY_DIRS})
  endif()

  set (EXE_SOURCES)
  set(EXE_BINARY_NAME ${EXE_NAME})

  # If requested create a modifier for the name that will be inserted between
  # the package name and the given name or exe_name for the test
  if(PARSE_ADD_DIR_TO_NAME)
    set(DIRECTORY_NAME "")
    tribits_create_name_from_current_source_directory(DIRECTORY_NAME)
    set(EXE_BINARY_NAME ${DIRECTORY_NAME}_${EXE_BINARY_NAME})
  endif()

  if(DEFINED PACKAGE_NAME AND NOT PARSE_NOEXEPREFIX)
    set(EXE_BINARY_NAME ${PACKAGE_NAME}_${EXE_BINARY_NAME})
  endif()

  # Exclude the build if requested
  if (${EXE_BINARY_NAME}_EXE_DISABLE)
    message("-- "
      "${EXE_BINARY_NAME} EXE NOT being built due to ${EXE_BINARY_NAME}_EXE_DISABLE="
      "'${${EXE_BINARY_NAME}_EXE_DISABLE}'")
    return()
  endif()

  # If exe is in subdirectory prepend that dir name to the source files
  if(PARSE_DIRECTORY )
    foreach( SOURCE_FILE ${PARSE_SOURCES} )
      if(IS_ABSOLUTE ${SOURCE_FILE})
        set (EXE_SOURCES ${EXE_SOURCES} ${SOURCE_FILE})
      else()
        set (EXE_SOURCES ${EXE_SOURCES} ${PARSE_DIRECTORY}/${SOURCE_FILE})
      endif()
    endforeach( )
  else()
    foreach( SOURCE_FILE ${PARSE_SOURCES} )
      set (EXE_SOURCES ${EXE_SOURCES} ${SOURCE_FILE})
    endforeach( )
  endif()

  # Assert that TESTONLYLIBS only contains TESTONLY libs!
  foreach(TESTONLYLIB ${PARSE_TESTONLYLIBS})
    set(PREFIXED_LIB "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${TESTONLYLIB}")
    if (NOT ${PREFIXED_LIB}_INCLUDE_DIRS)
      message(FATAL_ERROR "ERROR: '${TESTONLYLIB}' in TESTONLYLIBS not a TESTONLY lib!"
        "  If this a regular library in this SE package or in an dependent upstream SE"
        " package then TriBITS will link automatically to it.  If you remove this and it"
        " does not link, then you need to add a new SE package dependency to"
        " this SE package's dependencies file"
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    elseif(PARSE_INSTALLABLE)
      message(FATAL_ERROR "ERROR: TESTONLY lib '${TESTONLYLIB}' not allowed with"
        " INSTALLABLE executable!  An INSTALLABLE executable can only depend on"
        " non-TESTONLY libraries!  Otherwise, when shared libs are used, and"
        " TESTONLY library would not be installed and the installed executable"
        " would be unusable!" )
    endif()
  endforeach()

  # Assert that IMPORTEDLIBS are not TESTONLY libs are not regular package
  # libs!
  foreach(IMPORTEDLIB ${PARSE_IMPORTEDLIBS})
    set(PREFIXED_LIB "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${IMPORTEDLIB}")
    if (${PREFIXED_LIB}_INCLUDE_DIRS)
      message(FATAL_ERROR
        "ERROR: Lib '${IMPORTEDLIB}' being passed through"
        " IMPORTEDLIBS is not allowed to be a TESTONLY lib!"
        "  Use TESTONLYLIBS instead!" )
    endif()
    list(FIND ${PACKAGE_NAME}_LIBRARIES ${PREFIXED_LIB} FOUND_IDX)
    if (NOT FOUND_IDX EQUAL -1)
      message(FATAL_ERROR
        "ERROR: Lib '${IMPORTEDLIB}' in IMPORTEDLIBS is in"
        " this SE package and is *not* an external lib!"
        "  TriBITS takes care of linking against libs the current"
        " SE package automatically.  Please remove '${IMPORTEDLIB}' from IMPORTEDLIBS!")
    elseif (TARGET ${PREFIXED_LIB})
      message(FATAL_ERROR
        "ERROR: Lib '${IMPORTEDLIB}' being passed through"
        " IMPORTEDLIBS is *not* an external library but instead is a library"
        " defined in this CMake project!"
        "  TriBITS takes care of linking against libraries in dependent upstream"
        " SE packages.  If you want to link to a library in an upstream SE"
        " package then add the SE package name for that library to the appropriate"
        " list in this SE package's dependencies file"
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    endif()
  endforeach()

  # Convert from old DEPLIBS to TESTONLYLIBS and IMPORTEDLIBS
  foreach(DEPLIB ${PARSE_DEPLIBS})
    set(PREFIXED_LIB "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${DEPLIB}")
    if (${PREFIXED_LIB}_INCLUDE_DIRS)
      message(WARNING "WARNING: Passing TESTONLY lib '${DEPLIB}' through DEPLIBS"
        " is deprecated!  Instead, please pass through TESTONLYLIBS instead!"
        "  DEPLIBS is deprecated!")
      list(APPEND PARSE_TESTONLYLIBS ${DEPLIB})
    elseif (TARGET ${PREFIXED_LIB})
      message(WARNING "WARNING: Passing non-TESTONLY lib '${DEPLIB}' through DEPLIBS"
      " is deprecated!  The library '${DEPLIB}' appears to be a"
      " library defined in this CMake project."
      "  TriBITS takes care of linking against libraries in dependent upstream"
      " SE packages.  Therefore, please remove '${DEPLIB}' from this list."
      "   If you want to link to a library from an upstream SE"
      " package, then add the SE package name to the appropriate category"
      " in this SE package's dependencies file: "
      " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    else()
      message(WARNING "WARNING: Passing external lib '${DEPLIB}' through"
        " DEPLIBS is deprecated!  Instead, pass through IMPORTEDLIBS!"
        "  DEPLIBS is deprecated!"
        "  Please note that only external libs are allowed to be passed through"
        " IMPORTEDLIBS.")
      list(APPEND PARSE_IMPORTEDLIBS ${DEPLIB})
    endif()
  endforeach()

  foreach(TESTONLYLIB_IN ${PARSE_TESTONLYLIBS})
    set(TESTONLYLIB "${LIBRARY_NAME_PREFIX}${TESTONLYLIB_IN}")
    if (${TESTONLYLIB}_INCLUDE_DIRS)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message(STATUS "Adding include directories ${TESTONLYLIB}_INCLUDE_DIRS ...")
      endif()
      include_directories(${${TESTONLYLIB}_INCLUDE_DIRS})
    endif()
  endforeach()

  if (PARSE_DEFINES)
    message(WARNING "WARNING: Passing extra defines through 'DEFINES' ${PARSE_DEFINES}"
      " is deprecated.  Instead, pass them through 'TARGET_DEFINES'.  The 'DEFINES'"
      " argument was incorrectly implemented by calling add_definitions() which has"
      " directory scope and not function scope as was documented.  This resulted in"
      " confusing behavior.  If one wishes to set defines at the directly level,"
      " just call add_definitions() directly.")
    add_definitions(${PARSE_DEFINES})
  endif()

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("TRIBITS_ADD_EXECUTABLE: add_executable(${EXE_BINARY_NAME} ${EXE_SOURCES})")
  endif()
  add_executable(${EXE_BINARY_NAME} ${EXE_SOURCES})
  append_global_set(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${EXE_BINARY_NAME})

  if(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    set(${PARSE_ADDED_EXE_TARGET_NAME_OUT} ${EXE_BINARY_NAME} PARENT_SCOPE)
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

  set(LINK_LIBS)

  # First, add in the passed in TESTONLY dependent libraries
  if (PARSE_TESTONLYLIBS)
    foreach(LIB ${PARSE_TESTONLYLIBS})
      list(APPEND LINK_LIBS "${LIBRARY_NAME_PREFIX}${LIB}")
    endforeach()
  endif()

  # Second, add the package's own regular libraries
  if(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    list(APPEND LINK_LIBS ${${PACKAGE_NAME}_LIBRARIES})
  else()
    list(APPEND LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})
  endif()

  # Third, add the IMPORTEDLIBS
  if (PARSE_IMPORTEDLIBS)
    list(APPEND LINK_LIBS ${PARSE_IMPORTEDLIBS})
  endif()

  # Call include_directories() and link_directories(...) for upstream
  # dependent Packages and TPLs and accumulate the list of libraries that will
  # need to be linked to.

  if(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    AND NOT ${PACKAGE_NAME}_INCLUDE_DIRS
    )
    # No libraries have been added for this package so
    # add the upstream package and TPL includes and libraries
    tribits_sort_and_append_package_include_and_link_dirs_and_libs(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)
    tribits_sort_and_append_tpl_include_and_link_dirs_and_libs(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)
  endif()

  tribits_sort_and_append_package_include_and_link_dirs_and_libs(
    ${PACKAGE_NAME}  TEST  LINK_LIBS)

  if(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    tribits_sort_and_append_tpl_include_and_link_dirs_and_libs(
      ${PACKAGE_NAME}  TEST  LINK_LIBS)
  else()
    list(APPEND LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_TPL_LIBRARIES})
  endif()

  # Last, add last_lib to get extra link options on the link line
  if (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    list(APPEND LINK_LIBS last_lib)
  endif()

  if (${PROJECT_NAME}_DUMP_LINK_LIBS)
      message("-- ${EXE_NAME}:LINK_LIBS='${LINK_LIBS}'")
  endif()

  target_link_libraries(${EXE_BINARY_NAME} PUBLIC ${LINK_LIBS})

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

  if(${PROJECT_NAME}_INSTALL_EXECUTABLES AND PARSE_INSTALLABLE)
    install(
      TARGETS ${EXE_BINARY_NAME}
      EXPORT ${PROJECT_NAME}
        DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}
      COMPONENT ${PACKAGE_NAME}
    )
  endif()
endfunction()


#
# Setup include directories and library dependencies
#

#if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
#  message("TribitsAddExecutable.cmake")
#  print_var(${PACKAGE_NAME}_INCLUDE_DIRS)
#  print_var(${PACKAGE_NAME}_LIBRARY_DIRS)
#endif()
#
#if (NOT TRIBITS_ADD_EXECUTABLE_UNIT_TESTING)
#  include_directories(REQUIRED_DURING_INSTALLATION_TESTING
#    ${${PACKAGE_NAME}_INCLUDE_DIRS})
#  set_property(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
#    ${${PACKAGE_NAME}_LIBRARY_DIRS})
#endif()
