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
INCLUDE(CMakeParseArguments)
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
# built-in CMake `CONFIGURE_FILE()`_ command for rules on how it performs
# substitutions).  This command is typically used to configure the package's
# main `<packageDir>/cmake/<packageName>_config.h.in`_ file.
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
      "\n"
      "#ifndef ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG\n"
      "#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))\n"
      "#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))\n"
      "#  else\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)\n"
      "#  endif\n"
      "#endif\n"
      )
  ELSE()
    MULTILINE_SET(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)\n"
      )
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
#     <libBaseName>
#     [HEADERS <h0> <h1> ...]
#     [HEADERS_INSTALL_SUBDIR <headerssubdir>]
#     [NOINSTALLHEADERS <nih0> <hih1> ...]
#     [SOURCES <src0> <src1> ...]
#     [DEPLIBS <deplib0> <deplib1> ...]
#     [IMPORTEDLIBS <ideplib0> <ideplib1> ...]
#     [STATIC|SHARED]
#     [TESTONLY]
#     [NO_INSTALL_LIB_OR_HEADERS]
#     [CUDALIBRARY]
#     [ADDED_LIB_TARGET_NAME_OUT <libTargetName>]
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
#   ``<libBaseName>``
#
#     Required base name of the library.  The name of the actual library name
#     will be prefixed by ``${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}`` to
#     produce::
#     
#       <libTargetName> = ${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}<libBaseName>
#
#     This is the name passed to ``ADD_LIBRARY(<libTargetName> ...)``.  The
#     name is *not* prefixed by the package name.  CMake will of course add
#     any standard prefix or post-fix to the library file name appropriate for
#     the platform and if this is a static or shared library build (e.g. on
#     Linux prefix = ``'lib'``, postfix = ``'.so'`` for shared lib and postfix
#     = ``'.a'`` static lib) (see documentation for the built-in CMake command
#     ``ADD_LIBRARY()``.
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
#   ``HEADERS_INSTALL_SUBDIR <headerssubdir>``
#
#     Optional subdirectory that the headers will be installed under the
#     standard installation directory.  If ``<headerssubdir>!=""``, then the
#     headers will be installed under
#     ``${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/<headerssubdir>``.  Otherwise,
#     they will be installed under ``${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/``.
#     `Install Targets (TRIBITS_ADD_LIBRARY())`_.
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
#     ``TARGET_LINK_LIBRARIES(<libTargetName> ...)`` so that CMake knows about
#     the dependency structure of the libraries within this SE package.
#     **NOTE:** One must **not** list libraries in other upstream `TriBITS SE
#     Packages`_ or libraries built externally from this TriBITS CMake project
#     in ``DEPLIBS``.  The TriBITS system automatically handles linking to
#     libraries in upstream TriBITS SE packages.  External libraries need to
#     be listed in the ``IMPORTEDLIBS`` argument instead if they are not
#     already specified automatically using a `TriBITS TPL`_.
#
#   ``IMPORTEDLIBS <ideplib0> <ideplib1> ...``
#
#     List of dependent libraries built externally from this TriBITS CMake
#     project.  These libraries are passed into
#     ``TARGET_LINK_LIBRARIES(<libTargetName> ...)`` so that CMake knows about
#     the dependency.  However, note that external libraries are often better
#     handled as `TriBITS TPLs`_.  A well constructed TriBITS package and
#     library should never have to use this option!  So far, the only case
#     where ``IMPORTEDLIBS`` has been shown to be necessary is to pass in the
#     standard C math library ``m``.  In every other case, a TriBITS TPL
#     should be used instead.
#
#   ``STATIC`` or ``SHARED``
#
#     If ``STATIC`` is passed in, then a static library will be created
#     independent of the value of ``BUILD_SHARED_LIBS``.  If ``SHARED`` is
#     passed in, then a shared library will be created independent of the
#     value of ``BUILD_SHARED_LIBS``.  If neither ``STATIC`` or ``SHARED`` are
#     passed in, then a shared library will be created if
#     ``BUILD_SHARED_LIBS`` evaluates to true, otherwise and a static library
#     will be created.  If both ``STATIC`` and ``SHARED`` are passed in (which
#     is obviously a mistake), then a shared library will be created.
#     WARNING: Once you mark a library with ``STATIC``, then all of the
#     downstream libraries in the current SE package and all downstream SE
#     packages must also be also be marked with ``STATIC``.  That is because,
#     generally, one can not link a link a static lib against a downstream
#     shared lib since that is not portable (but can be done on some platforms
#     if, for example, ``-fPIC`` is specified).  So be careful to use
#     ``STATIC`` in all downstream libraries!
#
#   ``TESTONLY``
#
#     If passed in, then ``<libTargetName>`` will **not** be added to
#     ``${PACKAGE_NAME}_LIBRARIES`` and an install target for the library will
#     not be added.  In this case, the current include directories will be set
#     in the global variable ``<libTargetName>_INCLUDE_DIR`` which will be
#     used in `TRIBITS_ADD_EXECUTABLE()`_ when a test-only library is linked
#     in through its ``DEPLIBS`` argument.
#
#   ``NO_INSTALL_LIB_OR_HEADERS``
#
#     If specified, then no install targets will be added for the library
#     ``<libTargetName>`` or the header files listed in ``HEADERS``.
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
#   ``ADDED_LIB_TARGET_NAME_OUT <libTargetName>``
#
#     If specified, then on output the variable ``<libTargetName>`` will be
#     set with the name of the library passed to ``ADD_LIBRARY()``.  Having
#     this name allows the calling ``CMakeLists.txt`` file access and set
#     additional target properties (see `Additional Library and Source File
#     Properties (TRIBITS_ADD_LIBRARY())`_).
#
# .. _Include Directories (TRIBITS_ADD_LIBRARY()):
#
# **Include Directories (TRIBITS_ADD_LIBRARY())**
#
# Any base directories for the header files listed in the arguments
# ``HEADERS`` or ``NOINSTALLHEADERS`` should be passed into the standard CMake
# command ``INCLUDE_DIRECTORIES()`` **before** calling this function.  For
# example, a CMakeLists.txt file will look like::
#
#   ...
#
#   TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
#   CONFIGURE_FILE(...)
#
#   ...
#
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
#
#   ...
#
#   TRIBITS_ADD_LIBRARY( <libName>
#     SOURCES
#       <src0>.c
#       <subdir0>/<src1>.cpp
#       <subdir1>/<src2>.F90
#       ...
#     HEADERS
#        <header0>.h
#        <subdir0>/<header1>.hpp
#         ...
#     NONINSTALLHEADERS <header2>.hpp <header3>.hpp ...
#     ...
#     )
#
# The include of ``${CMAKE_CURRENT_BINARY_DIR}`` is needed for any generated
# header files (e.g. using raw ``CONFIGURE_FILE()`` or
# `TRIBITS_CONFIGURE_FILE()`_) or any generated Fortran ``*.mod`` module files
# generated as a byproduct of compiling F90+ source files (that contain one or
# more Fortran module declarations).
#
# The function ``TRIBITS_ADD_LIBRARY()`` will grab the list of all of the
# include directories in scope from prior calls to ``INCLUDE_DIRECTORIES()``
# and will append these to the variable ``${PACKAGE_NAME}_INCLUDE_DIRS``.
# This list of include directories is exported to downstream SE packages so
# they appear on the compile lines of all downstream object file compiles.
# This is a critical part of the "glue" that allows TriBITS packages to link
# up automatically (just by clearing dependencies in
# `<packageDir>/cmake/Dependencies.cmake`_ files).
#
# .. _Install Targets (TRIBITS_ADD_LIBRARY()):
#
# **Install Targets (TRIBITS_ADD_LIBRARY())**
#
# By default, an install target for the library is created using
# ``INSTALL(TARGETS <libTargetName> ...)`` to install into the directory
# ``${CMAKE_INSTALL_PREFIX}/lib/`` (actual install directory is given by
# ``${PROJECT}_INSTALL_LIB_DIR``, see `Setting the install prefix at configure
# time`_).  However, this install target will not get created if
# `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_ is ``FASLE`` and
# ``BUILD_SHARD_LIBS=OFF``.  But when ``BUILD_SHARD_LIBS=ON``, the install
# target will get added.  Also, this install target will *not* get added if
# ``TESTONLY`` or ``NO_INSTALL_LIB_OR_HEADERS`` are passed in.
#
# By default, an install target for the headers listed in ``HEADERS`` will get
# added using ``INSTALL(FILES <h0> <h1> ...)``, but only if ``TESTONLY`` and
# ``NO_INSTALL_LIB_OR_HEADERS`` are not passed in as well.  Also, the install
# target for the headers will not get added if
# `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_ is ``FASLE``.  If this
# install target is added, then the headers get installed into the flat
# directory ``${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/`` (default is
# ``${CMAKE_INSTALL_PREFIX}/include/``, see `Setting the install prefix at
# configure time`_).  If ``HEADERS_INSTALL_SUBDIR`` is set, then the headers
# will be installed under
# ``${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/<headerssubdir>/``.
#
# Note that an install target will *not* get created for the headers listed in
# ``NOINSTALLHEADERS``.
#
# .. _Additional Library and Source File Properties (TRIBITS_ADD_LIBRARY()):
#
# **Additional Library and Source File Properties (TRIBITS_ADD_LIBRARY())**
#
# Once ``ADD_LIBRARY(<libTargetName> ... <src0> <src1> ...)`` is called, one
# can set and change properties on the ``<libTargetName>`` library target
# using the built-in CMake command ``SET_TARGET_PROPERTIES()`` as well as set
# and change properties on any of the source files listed in ``SOURCES`` using
# the built-in CMake command ``SET_SOURCE_FILE_PROPERTIES()`` just like in any
# CMake project.  For example::
#
#   TRIBITS_ADD_LIBRARY( somelib ...
#     ADDED_LIB_TARGET_NAME_OUT  somelib_TARGET_NAME )
#
#   SET_TARGET_PROPERTIES( ${somelib_TARGET_NAME}
#     PROPERTIES  LINKER_LANGUAGE  CXX )
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

  IF(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_CALLED)
    MESSAGE(FATAL_ERROR "Must call TRIBITS_PACKAGE() before TRIBITS_ADD_LIBRARY() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
  ENDIF()

  IF(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
    MESSAGE(FATAL_ERROR "Must call TRIBITS_ADD_LIBRARY() before TRIBITS_PACKAGE_POSTPROCESS() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
  ENDIF()


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

  CMAKE_PARSE_ARGUMENTS(
    #prefix
    PARSE
    #Options    
    "STATIC;SHARED;TESTONLY;NO_INSTALL_LIB_OR_HEADERS;CUDALIBRARY"
    #one_value_keywords
    ""
    #mulit_value_keywords
    "HEADERS;HEADERS_INSTALL_SUBDIR;NOINSTALLHEADERS;SOURCES;DEPLIBS;IMPORTEDLIBS;DEFINES;ADDED_LIB_TARGET_NAME_OUT"
    ${ARGN}
    )

  TRIBITS_CHECK_FOR_UNPARSED_ARGUMENTS()

  # ToDo: Assert that HEADERS_INSTALL_SUBDIR has 0 or 1 entries!
  # ToDo: Assert that ADDED_LIB_TARGET_NAME_OUT as 0 or 1 entries!

  IF(PARSE_HEADERS)
    LIST(REMOVE_DUPLICATES PARSE_HEADERS)
  ENDIF()
  IF(PARSE_SOURCES)
    LIST(REMOVE_DUPLICATES PARSE_SOURCES)
  ENDIF()

  IF(PARSE_ADDED_LIB_TARGET_NAME_OUT)
    SET(${PARSE_ADDED_LIB_TARGET_NAME_OUT} PARENT_SCOPE)
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

    # Local variable to hold all of the libraries that will be directly linked
    # to this library.
    SET(LINK_LIBS)

    # Add dependent libraries passed directly in

    IF (PARSE_DEPLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "DEPLIBS = ${PARSE_DEPLIBS}")
    ENDIF()
    IF (PARSE_IMPORTEDLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "IMPORTEDLIBS = ${PARSE_IMPORTEDLIBS}")
    ENDIF()

    #
    # Add the DEPLIBS to the LINK_LIBS, assert correct usage of DEPLIBS, and
    # see if we need to link in the upstream SE package and TPL libs.
    #
    # We only want to link to the upstream dependent SE package and TPL
    # libraries if needed.  We only need to link to these upstream dependent
    # libraries when this is the first library being created for this SE
    # package or if this library does not depend on other libraries created
    # for this package.  Otherwise, we don't need to add the include
    # directories or link libraries because a dependent lib specified in
    # PARSE_DEPLIBS already has everything that we need.
    #
    # We also need to make special considerations for test libraries since
    # things need to be handled a little bit differently (but not much).  In the
    # case of test libraries, we need to also pull the test-only dependencies.
    # In this case, we will always assume that we will add in the test
    # libraries.
    #
    # ToDo: Turn the below deprecated WARNING messages to FATAL_ERROR once we
    # give enough time for people to clean up their codes.
    #

    SET(ADD_DEP_PACKAGE_AND_TPL_LIBS TRUE)

    SET(PREFIXED_DEPLIBS)

    FOREACH(LIB ${PARSE_DEPLIBS})

      SET(PREFIXED_LIB "${LIBRARY_NAME_PREFIX}${LIB}")

      # LIB_IN_SE_PKG?
      LIST(FIND ${PACKAGE_NAME}_LIBRARIES ${PREFIXED_LIB} FOUND_IDX)
      IF (FOUND_IDX GREATER -1)
        SET(LIB_IN_SE_PKG TRUE)
      ELSE()
        SET(LIB_IN_SE_PKG FALSE)
      ENDIF()

      # PREFIXED_LIB_IS_TESTONLY?
      IF (${PREFIXED_LIB}_INCLUDE_DIRS)
        SET(LIB_TESTONLY TRUE)
      ELSE()
        SET(LIB_TESTONLY FALSE)
      ENDIF()

      # Check for valid usage (sorted by most common to least common)
      IF (LIB_IN_SE_PKG AND NOT LIB_TESTONLY) #PARSE_TESTONLY=TRUE/FASLE
        # The library being created here is a library dependent on a regular
        # (non-TESTONLY) lib in this SE package.  This is valid usage of
        # DEPLIBS.  There is no need to link this new lib to the SE package's
        # upstream dependent SE package and TPL libraries because thse are
        # already linked into the lib ${LIB}.
        SET(ADD_DEP_PACKAGE_AND_TPL_LIBS FALSE)
      ELSEIF (PARSE_TESTONLY AND LIB_IN_SE_PKG AND NOT LIB_TESTONLY)
        # The library being created here is TESTONLY library and is
        # dependent on a regular (non-TESTONLY) lib.  This is valid usage of
        # DEPLIBS.  In the case of test-only libraries, we always link in
        # the upstream libs.
      ELSEIF (PARSE_TESTONLY AND LIB_TESTONLY) # LIB_IN_SE_PKG=TRUE/FASLE
        # The library being created here is TESTONLY library and is dependent
        # on another TESTONLY library.  This is valid usage of DEPLIBS.  In
        # this case we just hope that this SE package correctly specified a
        # TEST dependency on the upstream SE package that owns this upstream
        # TESTONLY library.
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE("-- "
            "Adding include directories for TESTONLY ${PREFIXED_LIB}_INCLUDE_DIRS ...")
        ENDIF()
        INCLUDE_DIRECTORIES(${${PREFIXED_LIB}_INCLUDE_DIRS})
      ELSEIF (NOT PARSE_TESTONLY AND LIB_TESTONLY) # LIB_IN_SE_PKG=TRUE/FASLE
        MESSAGE(WARNING "WARNING: '${LIB}' in DEPLIBS is a TESTONLY lib"
          " and it is illegal to link to this non-TESTONLY library '${LIBRARY_NAME}'."
          "  Such usage is deprecated (and this warning will soon become an error)!"
          "  If this is a regular library in this SE package or in an dependent upstream SE"
          " package then TriBITS will link automatically to it.  If you remove this and it"
          " does not link, then you need to add a new SE package dependency to"
          " this SE package's dependencies file"
          " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
        # ToDo: Turn the above to FATAL_ERROR after dropping deprecated code
      ELSEIF (NOT LIB_IN_SE_PKG AND TARGET ${PREFIXED_LIB} ) # PARSE_TESTONLY=TRUE/FALSE
        MESSAGE(WARNING "WARNING: '${LIB}' in DEPSLIBS is not"
          " a lib in this SE package but is a library defined in the current"
          " cmake project!  Such usage is  deprecated (and"
          " will result in a configure error soon).  If this is a library in"
          " a dependent upstream SE package, then simply remove it from this list."
          "  TriBITS automatically links in libraries in upstream SE packages."
          "  If you remove '${LIB}' from DEPLIBS and your code does"
          " not link, then you need to add a new SE package dependency to"
          " this SE package's dependencies file"
          " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
      ELSEIF (NOT LIB_IN_SE_PKG AND NOT TARGET ${PREFIXED_LIB} )
        MESSAGE(WARNING "WARNING: '${LIB}' in DEPSLIBS is not"
          " a lib defined in the current cmake project!  Such usage is deprecated (and"
          " will result in a configure error soon).  If this is an external"
          " lib you are trying to link in, it should likely be handled as a TriBITS"
          " TPL.  Otherwise, it should be passed in through IMPORTEDLIBS.  However,"
          " the only case we have found where IMPORTEDLIBS had to be used instead of"
          " through a proper TriBITS TPL is the C math library 'm'.")
      ELSE()
        MESSAGE(WARNING "WARNING: The case PARSE_TESTONLY=${PARSE_TESTONLY},"
          " LIB_IN_SE_PKG=${LIB_IN_SE_PKG}, LIB_TESTONLY=${LIB_TESTONLY}, has"
          " not yet been handled!")
      ENDIF()

      LIST(APPEND PREFIXED_DEPLIBS "${LIBRARY_NAME_PREFIX}${LIB}")

    ENDFOREACH()

    APPEND_SET(LINK_LIBS ${PREFIXED_DEPLIBS})

    #
    # Check IMPORTEDLIBS
    #

    FOREACH(IMPORTEDLIB ${PARSE_IMPORTEDLIBS})
      SET(PREFIXED_LIB "${LIBRARY_NAME_PREFIX}${IMPORTEDLIB}")
      LIST(FIND ${PACKAGE_NAME}_LIBRARIES ${PREFIXED_LIB} FOUND_IDX)
      IF (${PREFIXED_LIB}_INCLUDE_DIRS)
        MESSAGE(WARNING "WARNING: '${IMPORTEDLIB}' in IMPORTEDLIBS is a TESTONLY lib"
          " and it is illegal to pass in through IMPORTEDLIBS!"
          "  Such usage is deprecated (and this warning will soon become an error)!"
          "  Should '${IMPORTEDLIB}' instead be passed through DEPLIBS?")
        # ToDo: Turn the above to FATAL_ERROR after dropping deprecated code
      ELSEIF (FOUND_IDX GREATER -1)
        MESSAGE(WARNING "WARNING: Lib '${IMPORTEDLIB}' in IMPORTEDLIBS is in"
        " this SE package and is *not* an external lib!"
        "  TriBITS takes care of linking against libs the current"
        " SE package automatically.  Please remove it from IMPORTEDLIBS!")
      ELSEIF (TARGET ${PREFIXED_LIB})
        MESSAGE(WARNING "WARNING: Lib '${IMPORTEDLIB}' being passed through"
        " IMPORTEDLIBS is *not* an external library but instead is a library"
        " defined in this CMake project!"
        "  TriBITS takes care of linking against libraries in dependent upstream"
        " SE packages.  If you want to link to a library in an upstream SE"
        " package then add the SE package name to the appropriate category"
        " in this SE package's depencencies file: "
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
      ENDIF()
      # ToDo: Assert that this is not a test-only lib
      LIST(APPEND LINK_LIBS ${IMPORTEDLIB})
    ENDFOREACH()

    #
    # Link in the upstream TEST SE package and TPL libs
    #
    # We link these before those in the LIB SE package and TPL libs because
    # the TEST dependencies tend to be higher in the dependency tree.  It
    # should not really matter but it looks better on the link line.
    #

    IF (PARSE_TESTONLY)

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Pulling in header and libraries dependencies"
          " for TEST dependencies ...")
      ENDIF()

      TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
        ${PACKAGE_NAME}  TEST  LINK_LIBS)

      TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
        ${PACKAGE_NAME}  TEST  LINK_LIBS)

    ENDIF()

    #
    # Add the dependent LIB SE package and TPL libs
    #

    IF (ADD_DEP_PACKAGE_AND_TPL_LIBS)

      # If there are no dependent libs passed in, then this library can not
      # possibly depend on the package's other libraries so we must link to
      # the dependent libraries in dependent libraries and TPLs.

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "Pulling in header and libraries dependencies"
          " for LIB dependencies ...")
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
        ${PACKAGE_NAME}  LIB  LINK_LIBS)

      TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
        ${PACKAGE_NAME}  LIB  LINK_LIBS)

    ENDIF()

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(LINK_LIBS)
    ENDIF()

    # Add the library and all the dependencies

    IF (PARSE_DEFINES)
      ADD_DEFINITIONS(${PARSE_DEFINES})
    ENDIF()

    If (PARSE_STATIC)
      SET(STATIC_KEYWORD "STATIC")
    ELSE()
      SET(STATIC_KEYWORD)
    ENDIF()

    If (PARSE_SHARED)
      SET(SHARED_KEYWORD "SHARED")
    ELSE()
      SET(SHARED_KEYWORD)
    ENDIF()

    IF (NOT PARSE_CUDALIBRARY)
      ADD_LIBRARY(
        ${LIBRARY_NAME}
        ${STATIC_KEYWORD}
        ${SHARED_KEYWORD}
        ${PARSE_HEADERS}
        ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES}
        )
    ELSE()
      CUDA_ADD_LIBRARY(
        ${LIBRARY_NAME}
        ${PARSE_HEADERS}
        ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES}
        )
    ENDIF()

    IF(PARSE_ADDED_LIB_TARGET_NAME_OUT)
      SET(${PARSE_ADDED_LIB_TARGET_NAME_OUT} ${LIBRARY_NAME} PARENT_SCOPE)
    ENDIF()

    SET_PROPERTY(
      TARGET ${LIBRARY_NAME}
      APPEND PROPERTY
      LABELS ${PACKAGE_NAME}Libs ${PARENT_PACKAGE_NAME}Libs
      )

    SET_TARGET_PROPERTIES(
      ${LIBRARY_NAME}
      PROPERTIES
      VERSION ${${PROJECT_NAME}_VERSION}
      SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION}
      )

    PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
    PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})

    IF (${PROJECT_NAME}_DUMP_LINK_LIBS)
      MESSAGE("-- ${LIBRARY_NAME_IN}:LINK_LIBS='${LINK_LIBS}'")
    ENDIF()

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
      INSTALL(
        TARGETS ${LIBRARY_NAME}
        EXPORT ${PACKAGE_NAME}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}"
        LIBRARY DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
        ARCHIVE DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
        COMPONENT ${PACKAGE_NAME}
        )
    ENDIF()

    IF (INSTALL_HEADERS)
      TRIBITS_INSTALL_HEADERS(
        HEADERS  ${PARSE_HEADERS}
        INSTALL_SUBDIR  ${PARSE_HEADERS_INSTALL_SUBDIR}
        COMPONENT  ${PACKAGE_NAME}
        )
    ENDIF()

    # Append the new include dirs, library dirs, and libraries to this package's lists

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT  PACKAGE_LIBRARY_DIRS)

    IF (APPEND_LIB_AND_HEADERS_TO_PACKAGE)

      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCLUDE_DIRS_CURRENT})
      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS  ${LIBRARY_DIRS_CURRENT})
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

      LIST(REMOVE_DUPLICATES INCLUDE_DIRS_CURRENT)
      GLOBAL_SET(${LIBRARY_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${LIBRARY_NAME}_INCLUDE_DIRS)
      ENDIF()

    ENDIF()
  ENDIF() #if not in installation testing mode

  #
  # Adjust for installation testing
  #

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

  ENDIF() #installation testing mode

  #
  # Print the updates to the linkage variables
  #

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

ENDFUNCTION()
