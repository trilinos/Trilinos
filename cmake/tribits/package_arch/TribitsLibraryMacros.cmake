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

INCLUDE(TribitsCreateClientTemplateHeaders)
INCLUDE(ParseVariableArguments)
INCLUDE(GlobalSet)
INCLUDE(AppendSet)
INCLUDE(AppendGlob)
INCLUDE(AppendGlobalSet)
INCLUDE(AppendStringVar)
INCLUDE(PrependGlobalSet)
INCLUDE(RemoveGlobalDuplicates)
INCLUDE(TribitsGeneralMacros)
INCLUDE(SetAndIncDirs)

###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of the file
### TribitsPackageMacros.cmake!
###


#
# Macro that configures the package's main config.h file
#
FUNCTION(TRIBITS_ADD_CONFIG_DEFINE DEFINE)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("-- " "Package ${PARENT_PACKAGE_NAME}: adding compiler"
      " define to config file: ${DEFINE}")
  ENDIF()
  GLOBAL_SET(${PARENT_PACKAGE_NAME}_CONFIG_DEFINES
    "${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}\n#define ${DEFINE}")
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("-- ${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}")
  ENDIF()
ENDFUNCTION()


#
# @FUNCTION: TRIBITS_CONFIGURE_FILE()
#
# Macro that configures the package's main configured header file (typically
# called ``${PACKAGE_NAME}_config.h`` but any name can be used).
#
# Usage::
#
#   TRIBITS_CONFIGURE_FILE(<packageConfigFile>)
#
# This function requires the file::
#
#    ${PACKAGE_SOURCE_DIR}/cmake/<packageConfigFile>.in
#
# exists and it creates the file::
# 
#   ${CMAKE_CURRENT_BINARY_DIR}/<packageConfigFile>
#
# by calling the built-in ``CONFIGURE_FILE()`` command::
#
#   CONFIGURE_FILE(
#     ${PACKAGE_SOURCE_DIR}/cmake/<packageConfigFile>.in
#     ${CMAKE_CURRENT_BINARY_DIR}/<packageConfigFile>
#     )
#
# which does basic substitution of CMake variables (see documentation for
# built-in CMake ``CONFIGURE_FILE()`` command for rules on how it performs
# substitutions).
#
# In addition to just calling ``CONFIGURE_FILE()``, this function also aids in
# creating configured header files adding macros for deprecating code as
# described below.
#
# **Deprecated Code Macros**
#
# If ``${PARENT_PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS`` is ``TRUE`` (see
# `TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()`_), then the local CMake
# variable ``${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS`` is set which
# adds a define ``<PARENT_PACKAGE_NAME_UC>_DEPRECATED`` (where
# ``<PARENT_PACKAGE_NAME_UC>`` is the package name in all upper-case letters)
# which adds a compiler-specific deprecated warning for an entity.  To take
# advantage of this, just add the line::
#
#   @<PARENT_PACKAGE_NAME_UC>_DEPRECATED_DECLARATIONS@
#
# to the ``<packageConfigFile>.in`` file and it will be expanded at configure
# time.
#
# Then C/C++ code can use this macro to deprecate functions, variables,
# classes, etc., for example, using::
#
#   <PARENT_PACKAGE_NAME_UC>_DEPRECATED class SomeDepreatedClass { ... }.
#
# If the particular compiler does not support deprecated warnings, then this
# macro is defined to be empty.  See `Regulated Backward Compatibility and
# Deprecated Code`_ for more details.
#
FUNCTION(TRIBITS_CONFIGURE_FILE  PACKAGE_NAME_CONFIG_FILE)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_CONFIGURE_FILE: ${PACKAGE_NAME_CONFIG_FILE}")
  ENDIF()

  # Set up the deprecated attribute if showing deprecated warnings
  IF (${PARENT_PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS)
    MULTILINE_SET(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#ifndef ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED  __attribute__((__deprecated__))\n"
      "#  else\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  endif\n"
      "#endif\n"
      )
  ELSE()
    SET(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED")
  ENDIF()

  IF (${PARENT_PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
    APPEND_STRING_VAR(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "\n#define ${PARENT_PACKAGE_NAME_UC}_HIDE_DEPRECATED_CODE")
  ENDIF()

  # Set up the macro to create the define for time monitor
  SET(TIME_MONITOR_DEFINE_NAME ${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR)
  SET(FUNC_TIME_MONITOR_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR)
  SET(FUNC_TIME_MONITOR_DIFF_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR_DIFF)
  IF (${PARENT_PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR)
    MULTILINE_SET(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#ifndef ${FUNC_TIME_MONITOR_MACRO_NAME}\n"
      "#  define ${TIME_MONITOR_DEFINE_NAME}\n"
      "#  define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME) \\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, ${PARENT_PACKAGE_NAME_UC})\n"
      "#  define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF) \\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF)\n"
      "#endif\n"
      )
  ELSE()
    MULTILINE_SET(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME)\n"
      "#define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF)\n"
      )
  ENDIF()

  # Configure the file
  CONFIGURE_FILE(
    ${PACKAGE_SOURCE_DIR}/cmake/${PACKAGE_NAME_CONFIG_FILE}.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    )

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_ADD_LIBRARY()
# 
# Function used to add a CMake library and target using ``ADD_LIBRARY()``.
#
# Usage::
#
#   TRIBITS_ADD_LIBRARY(
#     <libName>
#     [HEADERS <h0> <h1> ...]
#     [NOINSTALLHEADERS <nih0> <hih1> ...]
#     [SOURCES <src0> <src1> ...]
#     [DEPLIBS <deplib0> <deplib1> ...]
#     [IMPORTEDLIBS <ideplib0> <ideplib1> ...]
#     [TESTONLY]
#     [NO_INSTALL_LIB_OR_HEADERS]
#     [CUDALIBRARY]
#     )
#
# *Sections:*
#
# * `Formal Arguments (TRIBITS_ADD_LIBRARY())`_
# * `Include Directories (TRIBITS_ADD_LIBRARY())`_
# * `Install Targets (TRIBITS_ADD_LIBRARY())`_
# * `Additional Library and Source File Properties (TRIBITS_ADD_LIBRARY())`_
# * `Miscellaneous Notes (TRIBITS_ADD_LIBRARY())`_
#
# .. _Formal Arguments (TRIBITS_ADD_LIBRARY()):
#
# **Formal Arguments (TRIBITS_ADD_LIBRARY())**
#
#   ``<libName>``
#
#     Required name of the library.  This is the name passed to
#     ``ADD_LIBRARY(<libName> ...)``.  The name is *not* prefixed by the
#     package name.  CMake will of course add any standard prefix or post-fix
#     to the library file name appropriate for the platform and if this is a
#     static or shared library build (see documentation for the built-in CMake
#     command ``ADD_LIBRARY()``.
#
#   ``HEADERS <h0> <h1> ...``
#
#     List of public header files for using this library.  By default, these
#     header files are assumed to be in the current source directory.  They
#     can also contain the relative path or absolute path to the files if they
#     are not in the current source directory.  This list of headers is passed
#     into ``ADD_LIBRARY(...)`` as well (which is not strictly needed but is
#     helpful for some build tools, like MS Visual Studio).  By default, these
#     headers will be installed (see `Install Targets
#     (TRIBITS_ADD_LIBRARY())`_).
#
#   ``NOINSTALLHEADERS <nih0> <hih1> ...``
#
#     List of private header files which are used by this library. These
#     headers are not installed and do not needed to be passed in for any
#     purpose other than to pass them into ``ADD_LIBRARY()`` as some build
#     tools like to have these listed (e.g. MS Visual Studio).
#
#   ``SOURCES <src0> <src1> ...``
#
#     List of source files passed into ``ADD_LIBRARY()`` that are compiled
#     into header files and included in the library.  The compiler used to
#     compile the files is determined automatically based on the file
#     extension (see CMake documentation for ``ADD_LIBRARY()``).
#
#   ``DEPLIBS <deplib0> <deplib1> ...``
#
#     List of dependent libraries that are built in the current SE package
#     that this library is dependent on.  These libraries are passed into
#     ``TARGET_LINK_LIBRARIES(<libName> ...)`` so that CMake knows about the
#     dependency structure of the libraries within the package.  **NOTE:** One
#     must **not** list libraries in other upstream SE packages or libraries
#     built externally from this TriBITS CMake project.  The TriBITS system
#     automatically handles linking to libraries in upstream TriBITS SE
#     packages.  External libraries need to be listed in the ``IMPORTEDLIBS``
#     argument instead.
#
#   ``IMPORTEDLIBS <ideplib0> <ideplib1> ...``
#
#     List of dependent libraries built externally from this TriBITS CMake
#     project.  These libraries are passed into
#     ``TARGET_LINK_LIBRARIES(<libName> ...)`` so that CMake knows about the
#     dependency.  These libraries are added to the
#     ``${PACKAGE_NAME}_LIBRARIES`` variable so that downstream SE packages
#     will also pick up these libraries and these libraries will show up in
#     the generated ``Makefile.export.${PACKAGE_NAME}`` and
#     ``${PACKAGE_NAME}Config.cmake`` files (if they are generated).  However,
#     not that external libraries are often better handled as `TriBITS TPLs`_.
#     A well constructed TriBITS package and library should never have to use
#     this option.
#
#   ``TESTONLY``
#
#     If passed in, then ``<libName>`` will **not** be added to
#     ``${PACKAGE_NAME}_LIBRARIES`` and an install target for the library will
#     not be added.  In this case, the current include directories will be set
#     in the global variable ``<libName>_INCLUDE_DIR`` which will be used in
#     `TRIBITS_ADD_EXECUTABLE()`_ when a test-only library is linked in
#     through its ``DEPLIBS`` argument.
#
#   ``NO_INSTALL_LIB_OR_HEADERS``
#
#     If specified, then no install targets will be added for the library
#     ``<libName>`` or the header files listed in ``HEADERS``.
#
#   ``CUDALIBRARY``
#
#     If specified then ``CUDA_ADD_LIBRARY()`` is used instead of
#     ``ADD_LIBRARY()`` where ``CUDA_ADD_LIBRARY()`` is assumed to be defined
#     by the standard ``FindCUDA.cmake`` module as processed using the
#     standard TriBITS ``FindTPLCUDA.cmake`` file (see `Standard TriBITS
#     TPLs`_).  For this option to work, this SE package must have an enabled
#     direct or indirect dependency on the TriBITS CUDA TPL or a
#     configure-time error may occur about not knowing about
#     ``CUDA_ALL_LIBRARY()``.
#
# .. _Include Directories (TRIBITS_ADD_LIBRARY()):
#
# **Include Directories (TRIBITS_ADD_LIBRARY())**
#
# Any base directories for the header files listed in the arguments
# ``HEADERS`` or ``NOINSTALLHEADERS`` should be passed into the standard CMake
# command ``INCLUDE_DIRECTORIES()`` *before* calling this function.  These
# include directories will then be added to current packages list of include
# directories ``${PACKAGE_NAME}_INCLUDE_DIRS`` which is then exported to
# downstream SE packages..
#
# .. _Install Targets (TRIBITS_ADD_LIBRARY()):
#
# **Install Targets (TRIBITS_ADD_LIBRARY())**
#
# By default, an install target for the library is created using
# ``INSTALL(TARGETS <libName> ...)`` to install into the directory
# ``${CMAKE_INSTALL_PREFIX}/lib/`` (actual install directory is given by
# ``${PROJECT}_INSTALL_LIB_DIR``, see `Setting the install prefix at configure
# time`_).  However, this install target will not get created if
# ``${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE`` and
# ``BUILD_SHARD_LIBS=OFF``.  But when ``BUILD_SHARD_LIBS=ON``, the install
# target will get created.  Also, this install target will *not* get created
# if ``TESTONLY`` or ``NO_INSTALL_LIB_OR_HEADERS`` are passed in.
#
# By default, an install target for the headers listed in ``HEADERS`` will get
# created using ``INSTALL(FILES <h0> <h1> ...)``, but only if ``TESTONLY`` and
# ``NO_INSTALL_LIB_OR_HEADERS`` are not passed in as well.  These headers get
# installed into the flat directory ``${CMAKE_INSTALL_PREFIX}/include/`` (the
# actual install directory is given by
# ``${PROJECT_NAME}_INSTALL_INCLUDE_DIR``, see `Setting the install prefix at
# configure time`_).  Note that an install target will *not* get created for
# the headers listed in ``NOINSTALLHEADERS``.
#
# .. _Additional Library and Source File Properties (TRIBITS_ADD_LIBRARY()):
#
# **Additional Library and Source File Properties (TRIBITS_ADD_LIBRARY())**
#
# Once ``ADD_LIBRARY(<libName> ... <src0> <src1> ...)`` is called, one can set
# and change properties on the ``<libName>`` library target using the built-in
# CMake command ``SET_TARGET_PROPERTIES()`` as well as set and change
# properties on any of the source files listed in ``SOURCES`` using the
# built-in CMake command ``SET_SOURCE_FILE_PROPERTIES()`` just like in any
# CMake project.
#
# .. _Miscellaneous Notes (TRIBITS_ADD_LIBRARY()):
#
# **Miscellaneous Notes (TRIBITS_ADD_LIBRARY())**
#
# **WARNING:** Do **NOT** use the built-in CMake command ``ADD_DEFINITIONS()``
# to add defines ``-D<someDefine>`` to the compile command line that will
# affect any of the header files in the package!  These CMake-added defines
# are only set locally in this directory and child directories.  These defines
# will **NOT** be set when code in peer directories (e.g. a downstream TriBITS
# packages) compiles that may include these header files.  To add defines that
# affect header files, please use a configured header file (see
# `TRIBITS_CONFIGURE_FILE()`_).
#
FUNCTION(TRIBITS_ADD_LIBRARY LIBRARY_NAME_IN)

  SET(LIBRARY_NAME_PREFIX "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}")

  SET(LIBRARY_NAME ${LIBRARY_NAME_PREFIX}${LIBRARY_NAME_IN})

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nTRIBITS_ADD_LIBRARY: ${LIBRARY_NAME}")
    IF(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
      MESSAGE("\n${PACKAGE_NAME}_LIBRARIES In installation testing mode,"
        " libraries will be found instead of created.")
    ENDIF()
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

  PARSE_ARGUMENTS(
    PARSE #prefix
    "HEADERS;NOINSTALLHEADERS;SOURCES;DEPLIBS;IMPORTEDLIBS;DEFINES" # Lists
    "TESTONLY;NO_INSTALL_LIB_OR_HEADERS;CUDALIBRARY" #Options
    ${ARGN} # Remaining arguments passed in
    )

  IF(PARSE_HEADERS)
    LIST(REMOVE_DUPLICATES PARSE_HEADERS)
  ENDIF()
  IF(PARSE_SOURCES)
    LIST(REMOVE_DUPLICATES PARSE_SOURCES)
  ENDIF()

  # ToDo: Deprecate and remove the usage of DEFINES!  People should be putting
  # defines into configured header files, not adding -D<macroName> to the
  # compile lines!

  IF (NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_TESTONLY)

    # Add the link directory for this library.

    SET_PROPERTY(DIRECTORY  APPEND  PROPERTY  PACKAGE_LIBRARY_DIRS
      ${CMAKE_CURRENT_BINARY_DIR})

    # NOTE: Above, this link path not really used here for anything.
    # Instead it is just added to the other set link library directories
    # that are already set.  These link directories are then extracted
    # and stored into stored in ${PACKAGE_NAME}_LIBRARY_DIRS.

    # Add whatever include directories have been defined so far

    INCLUDE_DIRECTORIES(AFTER ${${PACKAGE_NAME}_INCLUDE_DIRS})

    # Add whatever link directories have been added so far

    SET_PROPERTY(DIRECTORY  APPEND  PROPERTY  PACKAGE_LIBRARY_DIRS
      ${${PACKAGE_NAME}_LIBRARY_DIRS})

    # Local varaible to hold all of the libraries that will be directly linked
    # to this library.
    SET(LINK_LIBS)

    # Add dependent libraries passed directly in

    IF (PARSE_DEPLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "DEPLIBS = ${PARSE_DEPLIBS}")
    ENDIF()
    IF (PARSE_IMPORTEDLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "IMPORTEDLIBS = ${PARSE_IMPORTEDLIBS}")
    ENDIF()

    # Prepend DEPLIBS with LIBRARY_NAME_PREFIX.
    IF (PARSE_DEPLIBS)
      SET(PREFIXED_DEPLIBS)
      FOREACH(LIB ${PARSE_DEPLIBS})
        LIST(APPEND PREFIXED_DEPLIBS "${LIBRARY_NAME_PREFIX}${LIB}")
      ENDFOREACH()
      APPEND_SET(LINK_LIBS ${PREFIXED_DEPLIBS})
    ENDIF()
    IF (PARSE_IMPORTEDLIBS)
      APPEND_SET(LINK_LIBS ${PARSE_IMPORTEDLIBS})
    ENDIF()

    #
    # We only want to link to the dependent package and TPL libraries when we need
    # to.  We only need to link to these dependent libraries when this is the first
    # library being created for this package or if this library does not depend
    # on other libraries created for this package.  Otherwise, we don't need to
    # add the include directories or link libraries because a dependent lib
    # specified in PARSE_DEPLIBS already has everything that we need.
    #
    # We also need to make special considerations for test libraries since
    # things need to be handled a little bit differently (but not much).  In the
    # case of test libaries, we need to also pull the test-only dependencies.
    # In this case, we will always assume that we will add in the test
    # libraries.
    #

    SET(ADD_DEP_PACKAGE_AND_TPL_LIBS TRUE)

    IF (PARSE_DEPLIBS AND NOT PARSE_TESTONLY)
      FOREACH(DEPLIB ${PARSE_DEPLIBS})
        LIST(FIND ${PACKAGE_NAME}_LIBRARIES ${DEPLIB} DEPLIB_IDX)
        IF (NOT DEPLIB_IDX EQUAL -1)
          # The library being created here is dependent on another of this
          # package's libraries so there is no need to add in this package's
          # dependent package and TPL libraries.
          SET(ADD_DEP_PACKAGE_AND_TPL_LIBS FALSE)
        ENDIF()
      ENDFOREACH()
    ELSE()
      # If there are no dependent libs passed in, then this library can not
      # possiblly depend on the package's other libraries so we must link to
      # the dependent libraries in dependent libraries and TPLs.
    ENDIF()

    IF (ADD_DEP_PACKAGE_AND_TPL_LIBS)

      IF (NOT PARSE_TESTONLY)
        SET(LIB_OR_TEST_ARG LIB)
      ELSE()
        SET(LIB_OR_TEST_ARG TEST)
      ENDIF()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Pulling in header and libraries dependencies"
          " for ${LIB_OR_TEST_ARG} ...")
      ENDIF()

      #
      # Call INCLUDE_DIRECTORIES() and LINK_DIRECTORIES(...) for all upstream
      # dependent Packages and TPLs and accumulate the libraries to link against.
      #
      # NOTE: Adding these directories serves two purposes.  First, so that the includes
      # get added the the sources that get built for this library.  Second, so
      # that list full list of include directories can be extracted as a
      # propery and set on ${PACKAGE_NAME}_INCLUDE_DIRS
      #

      TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
        ${PACKAGE_NAME}  ${LIB_OR_TEST_ARG}  LINK_LIBS)

      TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
        ${PACKAGE_NAME}  ${LIB_OR_TEST_ARG}  LINK_LIBS)

    ENDIF()

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(LINK_LIBS)
    ENDIF()

    # Add the library and all the dependencies

    IF (PARSE_DEFINES)
      ADD_DEFINITIONS(${PARSE_DEFINES})
    ENDIF()

    IF (NOT PARSE_CUDALIBRARY)
      ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES})
    ELSE()
      CUDA_ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES})
    ENDIF()

    SET_PROPERTY(TARGET ${LIBRARY_NAME} APPEND PROPERTY
      LABELS ${PACKAGE_NAME}Libs ${PARENT_PACKAGE_NAME}Libs)

    SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES
      VERSION ${${PROJECT_NAME}_VERSION}
      SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION})

    PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
    PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})

    TARGET_LINK_LIBRARIES(${LIBRARY_NAME}  ${LINK_LIBS})

    # Add to the install target

    SET(INSTALL_LIB ON)
    SET(INSTALL_HEADERS ON)
    SET(APPEND_LIB_AND_HEADERS_TO_PACKAGE ON)

    IF (PARSE_TESTONLY)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Skipping installation hooks for this library"
          " because 'TESTONLY' was passed in ...")
      ENDIF()
      SET(INSTALL_LIB OFF)
      SET(INSTALL_HEADERS OFF)
      SET(APPEND_LIB_AND_HEADERS_TO_PACKAGE OFF)
    ELSEIF (PARSE_NO_INSTALL_LIB_OR_HEADERS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Skipping installation hooks for this library"
          " because 'NO_INSTALL_LIB_OR_HEADERS' was passed in ...")
      ENDIF()
      SET(INSTALL_LIB OFF)
      SET(INSTALL_HEADERS OFF)
    ELSEIF (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND NOT BUILD_SHARED_LIBS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Skipping installation of headers and libraries"
          " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and BUILD_SHARED_LIBS=FALSE  ...")
      ENDIF()
      SET(INSTALL_LIB OFF)
      SET(INSTALL_HEADERS OFF)
    ELSEIF (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND BUILD_SHARED_LIBS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Skipping installation of headers but installing libraries"
          " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and BUILD_SHARED_LIBS=TRUE  ...")
      ENDIF()
      SET(INSTALL_HEADERS OFF)
    ENDIF()

    IF (INSTALL_LIB OR INSTALL_HEADERS)
      SET_PROPERTY(GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS ON)
      SET_PROPERTY(GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS ON)
    ENDIF()

    IF (INSTALL_LIB)
      SET_PROPERTY(TARGET ${LIBRARY_NAME} PROPERTY INSTALL_RPATH
        "${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}")
      INSTALL(
        TARGETS ${LIBRARY_NAME}
        EXPORT ${PROJECT_NAME}
          RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}"
          LIBRARY DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
          ARCHIVE DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
          COMPONENT ${PACKAGE_NAME}
        )
    ENDIF()

    IF (INSTALL_HEADERS)
      INSTALL(
        FILES ${PARSE_HEADERS}
        DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
        COMPONENT ${PACKAGE_NAME}
        )
    ENDIF()

    # Append the new include dirs, library dirs, and libraries to this package's lists

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT  PACKAGE_LIBRARY_DIRS)

    IF (APPEND_LIB_AND_HEADERS_TO_PACKAGE)

      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCLUDE_DIRS_CURRENT})
      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS  ${LIBRARY_DIRS_CURRENT})
      IF (PARSE_IMPORTEDLIBS)
        PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES  ${PARSE_IMPORTEDLIBS})
      ENDIF()
      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES  ${LIBRARY_NAME})

      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_INCLUDE_DIRS)
      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARY_DIRS)
      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARIES)

      GLOBAL_SET(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES TRUE)

    ELSE()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Skipping augmentation of package's lists of include"
          " directories and libraries! ...")
      ENDIF()

      GLOBAL_SET(${LIBRARY_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${LIBRARY_NAME}_INCLUDE_DIRS)
      ENDIF()

    ENDIF()
  ENDIF() #if not in installation testing mode

  IF (${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)

    LIST(FIND ${PROJECT_NAME}_INSTALLATION_PACKAGE_LIST ${PACKAGE_NAME}
      ${PACKAGE_NAME}_WAS_INSTALLED)
    IF(${${PACKAGE_NAME}_WAS_INSTALLED} EQUAL -1)
      MESSAGE(FATAL_ERROR
        "The package ${PACKAGE_NAME} was not installed with ${PROJECT_NAME}!"
        "  Please disable package ${PACKAGE_NAME} or install it.")
    ENDIF()

    INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING  BEFORE
       ${${PACKAGE_NAME}_INSTALLATION_INCLUDE_DIRS}
       ${${TRIBITS_PACKAGE}_INSTALLATION_TPL_INCLUDE_DIRS})
    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
      ${${PACKAGE_NAME}_INSTALLATION_LIBRARY_DIRS})

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT PACKAGE_LIBRARY_DIRS)

    GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS ${LIBRARY_DIRS_CURRENT})
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES    ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})

  ENDIF() #instalation testing mode

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

ENDFUNCTION()
