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


INCLUDE(TribitsAddExecutableTestHelpers)
INCLUDE(TribitsCommonArgsHelpers)
INCLUDE(TribitsAddTestHelpers)
INCLUDE(TribitsGeneralMacros)

INCLUDE(PrintVar)
INCLUDE(AppendSet)
INCLUDE(CMakeParseArguments)

###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of the file
### TribitsPackageMacros.cmake before trying to change anything in this file!
###


#
# @FUNCTION: TRIBITS_ADD_EXECUTABLE()
#
# Function used to create an executable (typically for a test or example),
# using the built-in CMake command ``ADD_EXECUTABLE()``.
#
# Usage::
#
#   TRIBITS_ADD_EXECUTABLE(
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
# * `Formal Arguments (TRIBITS_ADD_EXECUTABLE())`_
# * `Executable and Target Name (TRIBITS_ADD_EXECUTABLE())`_
# * `Additional Executable and Source File Properties (TRIBITS_ADD_EXECUTABLE())`_
# * `Install Target (TRIBITS_ADD_EXECUTABLE())`_
#
# .. _Formal Arguments (TRIBITS_ADD_EXECUTABLE()):
#
# **Formal Arguments (TRIBITS_ADD_EXECUTABLE())**
#
#   ``<exeRootName>``
#
#     The root name of the executable (and CMake target) (see `Executable and
#     Target Name (TRIBITS_ADD_EXECUTABLE())`_).  This must be the first
#     argument.
#
#   ``NOEXEPREFIX``
#
#     If passed in, then ``${PACKAGE_NAME}_`` is not added the beginning of
#     the executable name (see `Executable and Target Name
#     (TRIBITS_ADD_EXECUTABLE())`_).
#
#   ``NOEXESUFFIX``
#
#     If passed in, then ``${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}`` and
#     not added to the end of the executable name (except for native Windows
#     builds, see `Executable and Target Name (TRIBITS_ADD_EXECUTABLE())`_).
#
#   ``ADD_DIR_TO_NAME``
#
#     If passed in, the directory path relative to the package's base
#     directory (with "/" replaced by "_") is added to the executable name
#     (see `Executable and Target Name (TRIBITS_ADD_EXECUTABLE())`_).  This
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
#     ``ADD_EXECUTABLE(<exeTargetName> ... )``.  After calling this function,
#     the properties of the source files can be altered using the built-in
#     CMake command ``SET_SOURCE_FILE_PROPERTIES()``.
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
#     See `TRIBITS_ADD_TEST()`_ for more details.
#
#   ``HOST <host0> <host1> ...``
#
#     The list of hosts for which to enable the test (see `TRIBITS_ADD_TEST()`_).
#
#   ``XHOST <host0> <host1> ...``
#
#     The list of hosts for which **not** to enable the test (see
#     `TRIBITS_ADD_TEST()`_).
#
#   ``HOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which to enable the test (see
#     `TRIBITS_ADD_TEST()`_).
#
#   ``XHOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which **not** to enable the test (see
#     `TRIBITS_ADD_TEST()`_).
#
#   ``EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...``
#
#     If specified, gives the names of CMake variables that must evaluate to
#     true, or the test will not be added (see `TRIBITS_ADD_TEST()`_).
#
#   ``TESTONLYLIBS <lib0> <lib1> ...``
#
#     Specifies extra test-only libraries defined in this CMake project that
#     will be linked to the executable using ``TARGET_LINK_LIBRARY()``.  Note
#     that regular libraries (i.e. not ``TESTONLY``) defined in the current SE
#     package or any upstream SE packages can *NOT* be listed!  TriBITS
#     automatically links non ``TESTONLY`` libraries in this package and
#     upstream packages to the executable.  The only libraries that should be
#     listed in this argument are either ``TESTONLY`` libraries.  The include
#     directories for each test-only library will automatically be added
#     using::
#
#       INCLUDE_DIRECTORIES(${<libi>_INCLUDE_DIRS})
#
#     where ``<libi>_INCLUDE_DIRS`` was set by::
#
#       TRIBITS_ADD_LIBRARY(<libi> ... TESTONLY ...)
#
#     Therefore, to link to a defined ``TESTONLY`` library in any upstream
#     enabled package, one just needs to pass in the library name through
#     ``TESTONLYLIBS ... <libi> ...`` and that is it!
#
#   ``IMPORTEDLIBS <lib0> <lib1> ...``
#
#     Specifies extra external libraries that will be linked to the executable
#     using ``TARGET_LINK_LIBRARY()``.  This can only be used for libraries
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
#     `TRIBITS_ADD_TEST()`_ for more details.
#
#   ``LINKER_LANGUAGE (C|CXX|Fortran)``
#
#     If specified, overrides the linker language used by setting the built-in
#     CMake target property ``LINKER_LANGUAGE``.  TriBITS sets the default
#     linker language as follows::
#
#       IF (${PROJECT_NAME}_ENABLE_CXX)
#         SET(LINKER_LANGUAGE CXX)
#       ELSEIF (${PROJECT_NAME}_ENABLE_C)
#         SET(LINKER_LANGUAGE C)
#       ELSE()
#         # Let CMake set the default linker language it wants based
#         # on source file extensions passed into ``ADD_EXECUTABLE()``.
#       ENDIF()
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
#     ``TARGET_COMPILE_DEFINITIONS(<exeTargetName> ...)``.  These should only
#     affect the listed sources for the built executable and not other
#     targets.
#
#   ``INSTALLABLE``
#
#     If passed in, then an install target will be added to install the built
#     executable into the ``${CMAKE_INSTALL_PREFIX}/bin/`` directory (see
#     `Install Target (TRIBITS_ADD_EXECUTABLE())`_).
#
#   ``ADDED_EXE_TARGET_NAME_OUT <exeTargetName>``
#
#     If specified, then on output the variable ``<exeTargetName>`` will be
#     set with the name of the executable target passed to
#     ``ADD_EXECUTABLE(<exeTargetName> ... )``.  Having this name allows the
#     calling ``CMakeLists.txt`` file access and set additional target
#     properties (see `Additional Executable and Source File Properties
#     (TRIBITS_ADD_EXECUTABLE())`_).
#
# .. _Executable and Target Name (TRIBITS_ADD_EXECUTABLE()):
#
# **Executable and Target Name (TRIBITS_ADD_EXECUTABLE())**
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
# .. _Additional Executable and Source File Properties (TRIBITS_ADD_EXECUTABLE()):
#
# **Additional Executable and Source File Properties (TRIBITS_ADD_EXECUTABLE())**
#
# Once ``ADD_EXECUTABLE(<exeTargetName> ... )`` is called and this function
# exists, one can set and change properties on the ``<exeTargetName>``
# executable target using the built-in ``SET_TARGET_PROPERTIES()`` command as
# well as properties on any of the source files listed in ``SOURCES`` using
# the built-in ``SET_SOURCE_FILE_PROPERTIES()`` command just like in any CMake
# project.  IF the executable is added, its name will be returned by the
# argument ``ADDED_EXE_TARGET_NAME_OUT <exeTargetName>``.  For example::
#
#   TRIBITS_ADD_EXECUTABLE( someExe ...
#     ADDED_EXE_TARGET_NAME_OUT  someExe_TARGET_NAME )
#
#   IF (someExe_TARGET_NAME)
#     SET_TARGET_PROPERTIES( ${someExe_TARGET_NAME}
#       PROPERTIES  LINKER_LANGUAGE  CXX )
#   ENDIF()
#
# The ``IF(someExe_TARGET_NAME)`` is needed in case the executable does not
# get added for some reason (see `Formal Arguments
# (TRIBITS_ADD_EXECUTABLE())`_ for logic that can result in the executable
# target not getting added).
#
# .. _Install Target (TRIBITS_ADD_EXECUTABLE()):
#
# **Install Target (TRIBITS_ADD_EXECUTABLE())**
#
# If ``INSTALLABLE`` is passed in, then an install target using the built-in
# CMake command ``INSTALL(TARGETS <exeTargetName> ...)`` is added to install
# the built executable into the ``${CMAKE_INSTALL_PREFIX}/bin/`` directory
# (actual install directory path is determined by
# ``${PROJECT_NAME}_INSTALL_RUNTIME_DIR``, see `Setting the install prefix at
# configure time`_) .
#
FUNCTION(TRIBITS_ADD_EXECUTABLE EXE_NAME)

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_ADD_EXECUTABLE: ${EXE_NAME} ${ARGN}")
  ENDIF()
  
  #
  # Confirm that package and subpackage macros/functions have been called inteh correct order
  #
  
  IF (CURRENTLY_PROCESSING_SUBPACKAGE)

    # This is a subpackage being processed

    IF(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      MESSAGE(FATAL_ERROR "Must call TRIBITS_SUBPACKAGE() before TRIBITS_ADD_EXECUTABLE()"
        " in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    ENDIF()

    IF(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      MESSAGE(FATAL_ERROR "Must call TRIBITS_ADD_EXECUTABLE() before "
        " TRIBITS_SUBPACKAGE_POSTPROCESS() in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    ENDIF()

  ELSE()

    # This is a package being processed

    IF(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_CALLED)
      MESSAGE(FATAL_ERROR "Must call TRIBITS_PACKAGE() before TRIBITS_ADD_EXECUTABLE()"
        " in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    ENDIF()

    IF(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
      MESSAGE(FATAL_ERROR "Must call TRIBITS_ADD_EXECUTABLE() before "
        " TRIBITS_PACKAGE_POSTPROCESS() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    ENDIF()

  ENDIF()


  #
  # A) Parse the input arguments
  #

  CMAKE_PARSE_ARGUMENTS(
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

  TRIBITS_CHECK_FOR_UNPARSED_ARGUMENTS()

  IF(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    SET(${PARSE_ADDED_EXE_TARGET_NAME_OUT} PARENT_SCOPE)
  ENDIF()

  #
  # B) Exclude building the test executable based on some criteria
  #

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_CATEGORIES(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_HOST_HOSTTYPE(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  TRIBITS_PROCESS_COMM_ARGS(ADD_SERIAL_EXE  ADD_MPI_EXE  ${PARSE_COMM})
  IF (NOT ADD_SERIAL_EXE AND NOT ADD_MPI_EXE)
    RETURN()
  ENDIF()

  #
  # C) Add the executable
  #

  SET(LIBRARY_NAME_PREFIX "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}")

  IF (NOT TRIBITS_ADD_EXECUTABLE_UNIT_TESTING)
    TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
      ${${PACKAGE_NAME}_INCLUDE_DIRS})
#    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
#      ${${PACKAGE_NAME}_LIBRARY_DIRS})
  ENDIF()

  SET (EXE_SOURCES)
  SET(EXE_BINARY_NAME ${EXE_NAME})

  # If requested create a modifier for the name that will be inserted between
  # the package name and the given name or exe_name for the test
  IF(PARSE_ADD_DIR_TO_NAME)
    SET(DIRECTORY_NAME "")
    TRIBITS_CREATE_NAME_FROM_CURRENT_SOURCE_DIRECTORY(DIRECTORY_NAME)
    SET(EXE_BINARY_NAME ${DIRECTORY_NAME}_${EXE_BINARY_NAME})
  ENDIF()

  IF(DEFINED PACKAGE_NAME AND NOT PARSE_NOEXEPREFIX)
    SET(EXE_BINARY_NAME ${PACKAGE_NAME}_${EXE_BINARY_NAME})
  ENDIF()

  # Exclude the build if requested
  IF (${EXE_BINARY_NAME}_EXE_DISABLE)
    MESSAGE("-- "
      "${EXE_BINARY_NAME} EXE NOT being built due to ${EXE_BINARY_NAME}_EXE_DISABLE="
      "'${${EXE_BINARY_NAME}_EXE_DISABLE}'")
    RETURN()
  ENDIF()

  # If exe is in subdirectory prepend that dir name to the source files
  IF(PARSE_DIRECTORY )
    FOREACH( SOURCE_FILE ${PARSE_SOURCES} )
      IF(IS_ABSOLUTE ${SOURCE_FILE})
        SET (EXE_SOURCES ${EXE_SOURCES} ${SOURCE_FILE})
      ELSE()
        SET (EXE_SOURCES ${EXE_SOURCES} ${PARSE_DIRECTORY}/${SOURCE_FILE})
      ENDIF()
    ENDFOREACH( )
  ELSE()
    FOREACH( SOURCE_FILE ${PARSE_SOURCES} )
      SET (EXE_SOURCES ${EXE_SOURCES} ${SOURCE_FILE})
    ENDFOREACH( )
  ENDIF()

  # Assert that TESTONLYLIBS only contains TESTONLY libs!
  FOREACH(TESTONLYLIB ${PARSE_TESTONLYLIBS})
    SET(PREFIXED_LIB "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${TESTONLYLIB}")
    IF (NOT ${PREFIXED_LIB}_INCLUDE_DIRS)
      MESSAGE(FATAL_ERROR "ERROR: '${TESTONLYLIB}' in TESTONLYLIBS not a TESTONLY lib!"
        "  If this a regular library in this SE package or in an dependent upstream SE"
        " package then TriBITS will link automatically to it.  If you remove this and it"
        " does not link, then you need to add a new SE package dependency to"
        " this SE package's dependencies file"
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    ELSEIF(PARSE_INSTALLABLE)
      MESSAGE(FATAL_ERROR "ERROR: TESTONLY lib '${TESTONLYLIB}' not allowed with"
        " INSTALLABLE executable!  An INSTALLABLE executable can only depend on"
        " non-TESTONLY libraries!  Otherwise, when shared libs are used, and"
        " TESTONLY library would not be installed and the installed executable"
        " would be unusable!" )
    ENDIF()
  ENDFOREACH()

  # Assert that IMPORTEDLIBS are not TESTONLY libs are not regular package
  # libs!
  FOREACH(IMPORTEDLIB ${PARSE_IMPORTEDLIBS})
    SET(PREFIXED_LIB "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${IMPORTEDLIB}")
    IF (${PREFIXED_LIB}_INCLUDE_DIRS)
      MESSAGE(FATAL_ERROR "ERROR: Lib '${IMPORTEDLIB}' being passed through"
      " IMPORTEDLIBS is not allowed to be a TESTONLY lib!"
      "  Use TESTONLYLIBS instead!" )
    ENDIF()
    LIST(FIND ${PACKAGE_NAME}_LIBRARIES ${PREFIXED_LIB} FOUND_IDX)
    IF (NOT FOUND_IDX EQUAL -1)
      MESSAGE(FATAL_ERROR "ERROR: Lib '${IMPORTEDLIB}' in IMPORTEDLIBS is in"
      " this SE package and is *not* an external lib!"
      "  TriBITS takes care of linking against libs the current"
      " SE package automatically.  Please remove '${IMPORTEDLIB}' from IMPORTEDLIBS!")
    ELSEIF (TARGET ${PREFIXED_LIB})
      MESSAGE(FATAL_ERROR "ERROR: Lib '${IMPORTEDLIB}' being passed through"
      " IMPORTEDLIBS is *not* an external library but instead is a library"
      " defined in this CMake project!"
      "  TriBITS takes care of linking against libraries in dependent upstream"
      " SE packages.  If you want to link to a library in an upstream SE"
      " package then add the SE package name for that library to the appropriate"
      " list in this SE package's dependencies file"
      " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    ENDIF()
  ENDFOREACH()

  # Convert from old DEPLIBS to TESTONLYLIBS and IMPORTEDLIBS
  FOREACH(DEPLIB ${PARSE_DEPLIBS})
    SET(PREFIXED_LIB "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${DEPLIB}")
    IF (${PREFIXED_LIB}_INCLUDE_DIRS)
      MESSAGE(WARNING "WARNING: Passing TESTONLY lib '${DEPLIB}' through DEPLIBS"
        " is deprecated!  Instead, please pass through TESTONLYLIBS instead!"
        "  DEPLIBS is deprecated!")
      LIST(APPEND PARSE_TESTONLYLIBS ${DEPLIB})
    ELSEIF (TARGET ${PREFIXED_LIB})
      MESSAGE(WARNING "WARNING: Passing non-TESTONLY lib '${DEPLIB}' through DEPLIBS"
      " is deprecated!  The library '${DEPLIB}' appears to be a"
      " library defined in this CMake project."
      "  TriBITS takes care of linking against libraries in dependent upstream"
      " SE packages.  Therefore, please remove '${DEPLIB}' from this list."
      "   If you want to link to a library from an upstream SE"
      " package, then add the SE package name to the appropriate category"
      " in this SE package's dependencies file: "
      " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
    ELSE()
      MESSAGE(WARNING "WARNING: Passing external lib '${DEPLIB}' through"
        " DEPLIBS is deprecated!  Instead, pass through IMPORTEDLIBS!"
        "  DEPLIBS is deprecated!"
        "  Please note that only external libs are allowed to be passed through"
        " IMPORTEDLIBS.")
      LIST(APPEND PARSE_IMPORTEDLIBS ${DEPLIB})
    ENDIF()
  ENDFOREACH()

  FOREACH(TESTONLYLIB_IN ${PARSE_TESTONLYLIBS})
    SET(TESTONLYLIB "${LIBRARY_NAME_PREFIX}${TESTONLYLIB_IN}")
    IF (${TESTONLYLIB}_INCLUDE_DIRS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Adding include directories ${TESTONLYLIB}_INCLUDE_DIRS ...")
      ENDIF()
      INCLUDE_DIRECTORIES(${${TESTONLYLIB}_INCLUDE_DIRS})
    ENDIF()
  ENDFOREACH()

  IF (PARSE_DEFINES)
    MESSAGE(WARNING "WARNING: Passing extra defines through 'DEFINES' ${PARSE_DEFINES}"
      " is deprecated.  Instead, pass them through 'TARGET_DEFINES'.  The 'DEFINES'"
      " argument was incorrectly implemented by calling ADD_DEFINITIONS() which has"
      " directory scope and not function scope as was documented.  This resulted in"
      " confusing behavior.  If one wishes to set defines at the directly level,"
      " just call ADD_DEFINITIONS() directly.")
    ADD_DEFINITIONS(${PARSE_DEFINES})
  ENDIF()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TRIBITS_ADD_EXECUTABLE: ADD_EXECUTABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})")
  ENDIF()
  ADD_EXECUTABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})
  APPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${EXE_BINARY_NAME})

  IF(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    SET(${PARSE_ADDED_EXE_TARGET_NAME_OUT} ${EXE_BINARY_NAME} PARENT_SCOPE)
  ENDIF()

  IF (PARSE_TARGET_DEFINES)
    TARGET_COMPILE_DEFINITIONS(${EXE_BINARY_NAME} PUBLIC ${PARSE_TARGET_DEFINES})
  ENDIF()

  IF(PARSE_NOEXESUFFIX AND NOT WIN32)
    SET_TARGET_PROPERTIES(${EXE_BINARY_NAME} PROPERTIES SUFFIX "")
  ELSE()
    SET_TARGET_PROPERTIES(${EXE_BINARY_NAME} PROPERTIES SUFFIX
      ${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX})
  ENDIF()

  TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG( ${EXE_BINARY_NAME}
    "${PARSE_LINKER_LANGUAGE}" )

  SET(LINK_LIBS)

  # First, add in the passed in TESTONLY dependent libraries
  IF (PARSE_TESTONLYLIBS)
    FOREACH(LIB ${PARSE_TESTONLYLIBS})
      LIST(APPEND LINK_LIBS "${LIBRARY_NAME_PREFIX}${LIB}")
    ENDFOREACH()
  ENDIF()

  # Second, add the package's own regular libraries
  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    LIST(APPEND LINK_LIBS ${${PACKAGE_NAME}_LIBRARIES})
  ELSE()
    LIST(APPEND LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})
  ENDIF()

  # Third, add the IMPORTEDLIBS
  IF (PARSE_IMPORTEDLIBS)
    LIST(APPEND LINK_LIBS ${PARSE_IMPORTEDLIBS})
  ENDIF()

  # Call INCLUDE_DIRECTORIES() and LINK_DIRECTORIES(...) for upstream
  # dependent Packages and TPLs and accumulate the list of libraries that will
  # need to be linked to.

  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    AND NOT ${PACKAGE_NAME}_INCLUDE_DIRS
    )
    # No libraries have been added for this package so
    # add the upstream package and TPL includes and libraries
    TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)
    TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)
  ENDIF()

  TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
    ${PACKAGE_NAME}  TEST  LINK_LIBS)

  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  TEST  LINK_LIBS)
  ELSE()
    LIST(APPEND LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_TPL_LIBRARIES})
  ENDIF()

  # Last, add last_lib to get extra link options on the link line
  IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    LIST(APPEND LINK_LIBS last_lib)
  ENDIF()

  IF (${PROJECT_NAME}_DUMP_LINK_LIBS)
      MESSAGE("-- ${EXE_NAME}:LINK_LIBS='${LINK_LIBS}'")
  ENDIF()

  TARGET_LINK_LIBRARIES(${EXE_BINARY_NAME} ${LINK_LIBS})

  ASSERT_DEFINED(${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
  IF (${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
    #MESSAGE("${EXE_BINARY_NAME}: Adding property LINK_SEARCH_START_STATIC")
    SET_PROPERTY(TARGET ${EXE_BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
  ENDIF()

  IF(PARSE_DIRECTORY)
    SET_TARGET_PROPERTIES( ${EXE_BINARY_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${PARSE_DIRECTORY} )
  ENDIF()

  SET_PROPERTY(TARGET ${EXE_BINARY_NAME} APPEND PROPERTY
    LABELS ${PACKAGE_NAME}Exes ${PARENT_PACKAGE_NAME}Exes)

  IF(${PROJECT_NAME}_INSTALL_EXECUTABLES AND PARSE_INSTALLABLE)
    INSTALL(
      TARGETS ${EXE_BINARY_NAME}
      EXPORT ${PROJECT_NAME}
        DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}
      COMPONENT ${PACKAGE_NAME}
    )
  ENDIF()
ENDFUNCTION()


#
# Setup include directories and library dependencies
#

#IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
#  MESSAGE("TribitsAddExecutable.cmake")
#  PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
#  PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
#ENDIF()
#
#IF (NOT TRIBITS_ADD_EXECUTABLE_UNIT_TESTING)
#  INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
#    ${${PACKAGE_NAME}_INCLUDE_DIRS})
#  SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
#    ${${PACKAGE_NAME}_LIBRARY_DIRS})
#ENDIF()
