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
INCLUDE(ParseVariableArguments)

###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of the file
### TribitsPackageMacros.cmake!
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
#     [DIRECTORY <dir>]
#     [DEPLIBS <lib0> <lib1> ...]
#     [COMM [serial] [mpi]]
#     [LINKER_LANGUAGE (C|CXX|Fortran)]
#     [DEFINES -D<define0> -D<define1> ...]
#     [INSTALLABLE]
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
#     The root name of the exectuable (and CMake target) (see `Executable and
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
#     not added to the end of the executable name (see `Executable and
#     Target Name (TRIBITS_ADD_EXECUTABLE())`_).
#
#   ``ADD_DIR_TO_NAME``
#
#     If passed in, the directory path relative to the package's base
#     directory (with "/" replaced by "_") is added to the executable name
#     (see `Executable and Target Name (TRIBITS_ADD_EXECUTABLE())`_).  This
#     provides a simple way to create unique test exectuable names inside of a
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
#     ``ADD_EXECUTABLE(<fullExeName> ... )``.  After calling this function,
#     the properties of the source files can be altered using the built-in
#     CMake command ``SET_SOURCE_FILE_PROPERTIES()``.
#
#   ``DIRECTORY <dir>``
#
#     If specified, then the sources for the executable listed in ``SOURCES
#     <src0> <src1> ...`` are assumed to be in the relative or absolute
#     directory ``<dir>`` instead of the current source directory.  This
#     directory path is prepended to each source file name ``<srci>`` unless
#     ``<srci>`` is an absolute path.
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
#   ``DEPLIBS <lib0> <lib1> ...``
#
#     Specifies extra libraries that will be linked to the executable using
#     ``TARGET_LINK_LIBRARY()``.  Note that regular libraries (i.e. not
#     ``TESTONLY``) defined in the current SE package or any upstream SE
#     packages do **NOT** need to be listed!  TriBITS automatically links non
#     ``TESTONLY`` libraries in this package and upstream packages to the
#     executable.  The only libraries that should be listed in this argument
#     are either ``TESTONLY`` libraries, or other libraries that are built
#     external from this CMake project and are not provided through a proper
#     `TriBITS TPL`_.  The latter usage of passing in external libraries is
#     not recommended.  External libraries should be handled as declared
#     `TriBITS TPLs`_.  For a ``TESTONLY`` library, the include directories
#     will automatically be added using::
#
#       INCLUDE_DIRECTORIES(${<libi>_INCLUDE_DIRS})
#
#     where ``<libi>_INCLUDE_DIRS`` was set by::
#
#       TRIBITS_ADD_LIBRARY(<libi> ... TESTONLY ...)
#
#     Therefore, to link to a defined ``TESTONLY`` library in any upstream
#     enabled package, one just needs to pass in the library name through
#     ``DEPLIBS ... <libi> ...`` and that is it!
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
#     CMake target property ``LINKER_LANGUAGE``.  By default, CMake chooses the
#     compiler to be used as the linker based on file extensions.  The most
#     typical use case for this option is when Fortran-only or C-only sources
#     are passed in through ``SOURCES`` but a C++ linker is needed because
#     there are upstream C++ libraries.
#
#   ``DEFINES -D<define0> -D<define1> ...``
#
#     Add the listed defines using ``ADD_DEFINITIONS()``.  These should only
#     affect the listed sources for the built executable and not other
#     compiles in this directory due to the FUNCTION scoping.
#
#   ``INSTALLABLE``
#
#     If passed in, then an install target will be added to install the built
#     executable into the ``${CMAKE_INSTALL_PREFIX}/bin/`` directory (see
#     `Install Target (TRIBITS_ADD_EXECUTABLE())`_).
#
# .. _Executable and Target Name (TRIBITS_ADD_EXECUTABLE()):
#
# **Executable and Target Name (TRIBITS_ADD_EXECUTABLE())**
#
# By default, the full name of the executable and target name
# is::
#
#   <fullExecName> = ${PACKAGE_NAME}_<exeRootName>
#
# If ``ADD_DIR_TO_NAME`` is set, then the directory path relative to the
# package base directory (with "/" replaced with "_"), or ``<relDirName>``, is
# added to the executable name to form::
#
#   <fullExecName> = ${PACKAGE_NAME}_<relDirName>_<exeRootName>
#
# If the option ``NOEXEPREFIX`` is passed in, then the prefix
# ``${PACKAGE_NAME}_`` is removed.
#
# The executable suffix ``${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}`` will be
# added to the actual executable file name if the option ``NOEXESUFFIX`` is
# *not* passed in but this suffix is never added to the target name.
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
# Once ``ADD_EXECUTABLE(<fullExeName> ... )`` is called and this function
# exists, one can set and change properties on the ``<fullExeName>``
# executable target using the built-in ``SET_TARGET_PROPERTIES()`` command as
# well as properties on any of the source files listed in ``SOURCES`` using
# the built-in ``SET_SOURCE_FILE_PROPERTIES()`` command just like in any CMake
# project.
#
# .. _Install Target (TRIBITS_ADD_EXECUTABLE()):
#
# **Install Target (TRIBITS_ADD_EXECUTABLE())**
#
# If ``INSTALLABLE`` is passed in, then an install target using the built-in
# CMake command ``INSTALL(TARGETS <fullExeName> ...)`` is added to install the
# built executable into the ``${CMAKE_INSTALL_PREFIX}/bin/`` directory (actual
# install directory path is determined by
# ``${PROJECT_NAME}_INSTALL_RUNTIME_DIR``, see `Setting the install prefix at
# configure time`_) .
# 
FUNCTION(TRIBITS_ADD_EXECUTABLE EXE_NAME)

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_ADD_EXECUTABLE: ${EXE_NAME} ${ARGN}")
  ENDIF()
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    "SOURCES;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;DIRECTORY;DEPLIBS;COMM;LINKER_LANGUAGE;DEFINES"
    #options
    "NOEXEPREFIX;NOEXESUFFIX;ADD_DIR_TO_NAME;INSTALLABLE"
    ${ARGN}
    )

  #
  # B) Exclude building the test executable based on some several criteria
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

  FOREACH(DEPLIB ${PARSE_DEPLIBS})
    IF (${DEPLIB}_INCLUDE_DIRS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Adding include directories ${DEPLIB}_INCLUDE_DIRS ...")
        #PRINT_VAR(${DEPLIB}_INCLUDE_DIRS)
      ENDIF()
      INCLUDE_DIRECTORIES(${${DEPLIB}_INCLUDE_DIRS})
    ENDIF()
  ENDFOREACH()

  IF (PARSE_DEFINES)
    ADD_DEFINITIONS(${PARSE_DEFINES})
  ENDIF()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TRIBITS_ADD_EXECUTABLE: ADD_EXECTUABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})")
  ENDIF()
  ADD_EXECUTABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})
  APPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${EXE_BINARY_NAME})

  IF(PARSE_NOEXESUFFIX AND NOT WIN32)
    SET_TARGET_PROPERTIES(${EXE_BINARY_NAME} PROPERTIES SUFFIX "")
  ELSE()
    SET_TARGET_PROPERTIES(${EXE_BINARY_NAME} PROPERTIES SUFFIX
      ${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX})
  ENDIF()

  TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG( ${EXE_BINARY_NAME}
    "${PARSE_LINKER_LANGUAGE}" )

  SET(LINK_LIBS)

  # First, add in the passed in dependent libraries
  IF (PARSE_DEPLIBS)
    SET(LIBRARY_NAME_PREFIX "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}")
    SET(PREFIXED_DEPLIBS)
    FOREACH(LIB ${PARSE_DEPLIBS})
      LIST(APPEND PREFIXED_DEPLIBS "${LIBRARY_NAME_PREFIX}${LIB}")
    ENDFOREACH()
    APPEND_SET(LINK_LIBS ${PREFIXED_DEPLIBS})
  ENDIF()
  # 2009/01/09: rabartl: Above, I moved the list of dependent
  # libraries first to get around a problem with test-only libraries
  # creating multiple duplicate libraries on the link line with
  # CMake.

  # Second, add the package's own regular libraries
  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    APPEND_SET(LINK_LIBS ${${PACKAGE_NAME}_LIBRARIES})
  ELSE()
    APPEND_SET(LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})
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
    APPEND_SET(LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_TPL_LIBRARIES})
  ENDIF()

  # Last, add last_lib to get extra link options on the link line
  IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    APPEND_SET(LINK_LIBS last_lib)
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(LINK_LIBS)
  ENDIF()

  TARGET_LINK_LIBRARIES(${EXE_BINARY_NAME} ${LINK_LIBS})

  IF ("${CMAKE_VERSION}" VERSION_GREATER "2.8.4")
    ASSERT_DEFINED(${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
    IF (${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
      #MESSAGE("${EXE_BINARY_NAME}: Adding property LINK_SEARCH_START_STATIC")
      SET_PROPERTY(TARGET ${EXE_BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
    ENDIF()
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
