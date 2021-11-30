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

include(TribitsCreateClientTemplateHeaders)
include(CMakeParseArguments)
include(GlobalSet)
include(AppendSet)
include(AppendGlob)
include(AppendGlobalSet)
include(AppendStringVar)
include(PrependGlobalSet)
include(RemoveGlobalDuplicates)
include(TribitsGeneralMacros)
include(TribitsReportInvalidTribitsUsage)
include(SetAndIncDirs)

###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of the file
### TribitsPackageMacros.cmake!
###


#
# Macro that configures the package's main config.h file
#
function(tribits_add_config_define DEFINE)
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("-- " "Package ${PARENT_PACKAGE_NAME}: adding compiler"
      " define to config file: ${DEFINE}")
  endif()
  global_set(${PARENT_PACKAGE_NAME}_CONFIG_DEFINES
    "${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}\n#define ${DEFINE}")
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("-- ${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}")
  endif()
endfunction()


#
# @FUNCTION: tribits_configure_file()
#
# Macro that configures the package's main configured header file (typically
# called ``${PACKAGE_NAME}_config.h`` but any name can be used).
#
# Usage::
#
#   tribits_configure_file(<packageConfigFile>)
#
# This function requires the file::
#
#    ${PACKAGE_SOURCE_DIR}/cmake/<packageConfigFile>.in
#
# exists and it creates the file::
#
#   ${CMAKE_CURRENT_BINARY_DIR}/<packageConfigFile>
#
# by calling the built-in ``configure_file()`` command::
#
#   configure_file(
#     ${PACKAGE_SOURCE_DIR}/cmake/<packageConfigFile>.in
#     ${CMAKE_CURRENT_BINARY_DIR}/<packageConfigFile>
#     )
#
# which does basic substitution of CMake variables (see documentation for
# built-in CMake `configure_file()`_ command for rules on how it performs
# substitutions).  This command is typically used to configure the package's
# main `<packageDir>/cmake/<packageName>_config.h.in`_ file.
#
# In addition to just calling ``configure_file()``, this function also aids in
# creating configured header files adding macros for deprecating code as
# described below.
#
# **Deprecated Code Macros**
#
# If ``${PARENT_PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS`` is ``TRUE`` (see
# `tribits_add_show_deprecated_warnings_option()`_), then the local CMake
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
function(tribits_configure_file  PACKAGE_NAME_CONFIG_FILE)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nPACKAGE_CONFIGURE_FILE: ${PACKAGE_NAME_CONFIG_FILE}")
  endif()

  # Set up the deprecated attribute if showing deprecated warnings
  if (${PARENT_PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS)
    multiline_set(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
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
  else()
    multiline_set(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)\n"
      )
  endif()

  if (${PARENT_PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
    append_string_var(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "\n#define ${PARENT_PACKAGE_NAME_UC}_HIDE_DEPRECATED_CODE")
  endif()

  # Set up the macro to create the define for time monitor
  set(TIME_MONITOR_DEFINE_NAME ${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR)
  set(FUNC_TIME_MONITOR_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR)
  set(FUNC_TIME_MONITOR_DIFF_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR_DIFF)
  if (${PARENT_PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR)
    multiline_set(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#ifndef ${FUNC_TIME_MONITOR_MACRO_NAME}\n"
      "#  define ${TIME_MONITOR_DEFINE_NAME}\n"
      "#  define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME) \\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, ${PARENT_PACKAGE_NAME_UC})\n"
      "#  define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF) \\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF)\n"
      "#endif\n"
      )
  else()
    multiline_set(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME)\n"
      "#define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF)\n"
      )
  endif()

  # Configure the file
  configure_file(
    ${PACKAGE_SOURCE_DIR}/cmake/${PACKAGE_NAME_CONFIG_FILE}.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    )

endfunction()


#
# @FUNCTION: tribits_add_library()
#
# Function used to add a CMake library and target using ``add_library()``.
#
# Usage::
#
#   tribits_add_library(
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
# * `Formal Arguments (tribits_add_library())`_
# * `Include Directories (tribits_add_library())`_
# * `Install Targets (tribits_add_library())`_
# * `Additional Library and Source File Properties (tribits_add_library())`_
# * `Miscellaneous Notes (tribits_add_library())`_
#
# .. _Formal Arguments (tribits_add_library()):
#
# **Formal Arguments (tribits_add_library())**
#
#   ``<libBaseName>``
#
#     Required base name of the library.  The name of the actual library name
#     will be prefixed by ``${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}`` to
#     produce::
#     
#       <libTargetName> = ${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}<libBaseName>
#
#     This is the name passed to ``add_library(<libTargetName> ...)``.  The
#     name is *not* prefixed by the package name.  CMake will of course add
#     any standard prefix or post-fix to the library file name appropriate for
#     the platform and if this is a static or shared library build (e.g. on
#     Linux prefix = ``'lib'``, postfix = ``'.so'`` for shared lib and postfix
#     = ``'.a'`` static lib) (see documentation for the built-in CMake command
#     ``add_library()``.
#
#   ``HEADERS <h0> <h1> ...``
#
#     List of public header files for using this library.  By default, these
#     header files are assumed to be in the current source directory.  They
#     can also contain the relative path or absolute path to the files if they
#     are not in the current source directory.  This list of headers is passed
#     into ``add_library(...)`` as well (which is not strictly needed but is
#     helpful for some build tools, like MS Visual Studio).  By default, these
#     headers will be installed (see `Install Targets
#     (tribits_add_library())`_).
#
#   ``HEADERS_INSTALL_SUBDIR <headerssubdir>``
#
#     Optional subdirectory that the headers will be installed under the
#     standard installation directory.  If ``<headerssubdir>!=""``, then the
#     headers will be installed under
#     ``${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/<headerssubdir>``.  Otherwise,
#     they will be installed under ``${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/``.
#     `Install Targets (tribits_add_library())`_.
#
#   ``NOINSTALLHEADERS <nih0> <hih1> ...``
#
#     List of private header files which are used by this library. These
#     headers are not installed and do not needed to be passed in for any
#     purpose other than to pass them into ``add_library()`` as some build
#     tools like to have these listed (e.g. MS Visual Studio).
#
#   ``SOURCES <src0> <src1> ...``
#
#     List of source files passed into ``add_library()`` that are compiled
#     into header files and included in the library.  The compiler used to
#     compile the files is determined automatically based on the file
#     extension (see CMake documentation for ``add_library()``).
#
#   ``DEPLIBS <deplib0> <deplib1> ...``
#
#     List of dependent libraries that are built in the current SE package
#     that this library is dependent on.  These libraries are passed into
#     ``target_link_libraries(<libTargetName> ...)`` so that CMake knows about
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
#     ``target_link_libraries(<libTargetName> ...)`` so that CMake knows about
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
#     used in `tribits_add_executable()`_ when a test-only library is linked
#     in through its ``DEPLIBS`` argument.
#
#   ``NO_INSTALL_LIB_OR_HEADERS``
#
#     If specified, then no install targets will be added for the library
#     ``<libTargetName>`` or the header files listed in ``HEADERS``.
#
#   ``CUDALIBRARY``
#
#     If specified then ``cuda_add_library()`` is used instead of
#     ``add_library()`` where ``cuda_add_library()`` is assumed to be defined
#     by the standard ``FindCUDA.cmake`` module as processed using the
#     standard TriBITS ``FindTPLCUDA.cmake`` file (see `Standard TriBITS
#     TPLs`_).  For this option to work, this SE package must have an enabled
#     direct or indirect dependency on the TriBITS CUDA TPL or a
#     configure-time error may occur about not knowing about
#     ``cuda_all_library()``.
#
#   ``ADDED_LIB_TARGET_NAME_OUT <libTargetName>``
#
#     If specified, then on output the variable ``<libTargetName>`` will be
#     set with the name of the library passed to ``add_library()``.  Having
#     this name allows the calling ``CMakeLists.txt`` file access and set
#     additional target properties (see `Additional Library and Source File
#     Properties (tribits_add_library())`_).
#
# .. _Include Directories (tribits_add_library()):
#
# **Include Directories (tribits_add_library())**
#
# Any base directories for the header files listed in the arguments
# ``HEADERS`` or ``NOINSTALLHEADERS`` should be passed into the standard CMake
# command ``include_directories()`` **before** calling this function.  For
# example, a CMakeLists.txt file will look like::
#
#   ...
#
#   tribits_configure_file(${PACKAGE_NAME}_config.h)
#   configure_file(...)
#
#   ...
#
#   include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#   include_directories(${CMAKE_CURRENT_BINARY_DIR})
#
#   ...
#
#   tribits_add_library( <libName>
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
# header files (e.g. using raw ``configure_file()`` or
# `tribits_configure_file()`_) or any generated Fortran ``*.mod`` module files
# generated as a byproduct of compiling F90+ source files (that contain one or
# more Fortran module declarations).
#
# The function ``tribits_add_library()`` will grab the list of all of the
# include directories in scope from prior calls to ``include_directories()``
# and will append these to the variable ``${PACKAGE_NAME}_INCLUDE_DIRS``.
# This list of include directories is exported to downstream SE packages so
# they appear on the compile lines of all downstream object file compiles.
# This is a critical part of the "glue" that allows TriBITS packages to link
# up automatically (just by clearing dependencies in
# `<packageDir>/cmake/Dependencies.cmake`_ files).
#
# .. _Install Targets (tribits_add_library()):
#
# **Install Targets (tribits_add_library())**
#
# By default, an install target for the library is created using
# ``install(TARGETS <libTargetName> ...)`` to install into the directory
# ``${CMAKE_INSTALL_PREFIX}/lib/`` (actual install directory is given by
# ``${PROJECT}_INSTALL_LIB_DIR``, see `Setting the install prefix`_).
# However, this install target will not get created if
# `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_ is ``FASLE`` and
# ``BUILD_SHARD_LIBS=OFF``.  But when ``BUILD_SHARD_LIBS=ON``, the install
# target will get added.  Also, this install target will *not* get added if
# ``TESTONLY`` or ``NO_INSTALL_LIB_OR_HEADERS`` are passed in.
#
# By default, an install target for the headers listed in ``HEADERS`` will get
# added using ``install(FILES <h0> <h1> ...)``, but only if ``TESTONLY`` and
# ``NO_INSTALL_LIB_OR_HEADERS`` are not passed in as well.  Also, the install
# target for the headers will not get added if
# `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_ is ``FASLE``.  If this
# install target is added, then the headers get installed into the flat
# directory ``${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/`` (default is
# ``${CMAKE_INSTALL_PREFIX}/include/``, see `Setting the install prefix`_).
# If ``HEADERS_INSTALL_SUBDIR`` is set, then the headers will be installed
# under ``${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/<headerssubdir>/``.
#
# Note that an install target will *not* get created for the headers listed in
# ``NOINSTALLHEADERS``.
#
# .. _Additional Library and Source File Properties (tribits_add_library()):
#
# **Additional Library and Source File Properties (tribits_add_library())**
#
# Once ``add_library(<libTargetName> ... <src0> <src1> ...)`` is called, one
# can set and change properties on the ``<libTargetName>`` library target
# using the built-in CMake command ``set_target_properties()`` as well as set
# and change properties on any of the source files listed in ``SOURCES`` using
# the built-in CMake command ``set_source_file_properties()`` just like in any
# CMake project.  For example::
#
#   tribits_add_library( somelib ...
#     ADDED_LIB_TARGET_NAME_OUT  somelib_TARGET_NAME )
#
#   set_target_properties( ${somelib_TARGET_NAME}
#     PROPERTIES  LINKER_LANGUAGE  CXX )
#
# .. _Miscellaneous Notes (tribits_add_library()):
#
# **Miscellaneous Notes (tribits_add_library())**
#
# When the file ``Version.cmake`` exists and the CMake variables
# ``${PROJECT_NAME}_VERSION`` and ``${PROJECT_NAME}_MAJOR_VERSION`` are
# defined, then produced shared libraries will be given the standard SOVERSION
# symlinks (see `<projectDir>/Version.cmake`_).
#
# **WARNING:** Do **NOT** use the built-in CMake command ``add_definitions()``
# to add defines ``-D<someDefine>`` to the compile command line that will
# affect any of the header files in the package!  These CMake-added defines
# are only set locally in this directory and child directories.  These defines
# will **NOT** be set when code in peer directories (e.g. a downstream TriBITS
# packages) compiles that may include these header files.  To add defines that
# affect header files, please use a configured header file (see
# `tribits_configure_file()`_).
#
function(tribits_add_library LIBRARY_NAME_IN)

  #
  # Confirm that package and subpackage macros/functions have been called in the correct order
  #

  if (CURRENTLY_PROCESSING_SUBPACKAGE)

    # This is a subpackage being processed

    if(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_subpackage() before tribits_add_library()"
        " in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

    if(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_library() before "
        " tribits_subpackage_postprocess() in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

  else()

    # This is a package being processed

    if(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_package() or tribits_package_decl() before"
        " tribits_add_library() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

    if(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_library() before "
        " tribits_package_postprocess() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

  endif()

  set(LIBRARY_NAME_PREFIX "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}")
  set(LIBRARY_NAME ${LIBRARY_NAME_PREFIX}${LIBRARY_NAME_IN})

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_ADD_LIBRARY: ${LIBRARY_NAME}")
    if(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
      message("\n${PACKAGE_NAME}_LIBRARIES In installation testing mode,"
        " libraries will be found instead of created.")
    endif()
  endif()

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(${PACKAGE_NAME}_INCLUDE_DIRS)
    print_var(${PACKAGE_NAME}_LIBRARY_DIRS)
    print_var(${PACKAGE_NAME}_LIBRARIES)
  endif()

  cmake_parse_arguments(
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

  tribits_check_for_unparsed_arguments()

  # ToDo: Assert that HEADERS_INSTALL_SUBDIR has 0 or 1 entries!
  # ToDo: Assert that ADDED_LIB_TARGET_NAME_OUT as 0 or 1 entries!

  if(PARSE_HEADERS)
    list(REMOVE_DUPLICATES PARSE_HEADERS)
  endif()
  if(PARSE_SOURCES)
    list(REMOVE_DUPLICATES PARSE_SOURCES)
  endif()

  if(PARSE_ADDED_LIB_TARGET_NAME_OUT)
    set(${PARSE_ADDED_LIB_TARGET_NAME_OUT} PARENT_SCOPE)
  endif()

  # ToDo: Deprecate and remove the usage of DEFINES!  People should be putting
  # defines into configured header files, not adding -D<macroName> to the
  # compile lines!

  if (NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_TESTONLY)

    # Add the link directory for this library.

    set_property(DIRECTORY  APPEND  PROPERTY  PACKAGE_LIBRARY_DIRS
      ${CMAKE_CURRENT_BINARY_DIR})

    # NOTE: Above, this link path not really used here for anything.
    # Instead it is just added to the other set link library directories
    # that are already set.  These link directories are then extracted
    # and stored into stored in ${PACKAGE_NAME}_LIBRARY_DIRS.

    # Add whatever include directories have been defined so far

    include_directories(AFTER ${${PACKAGE_NAME}_INCLUDE_DIRS})

    # Add whatever link directories have been added so far

    set_property(DIRECTORY  APPEND  PROPERTY  PACKAGE_LIBRARY_DIRS
      ${${PACKAGE_NAME}_LIBRARY_DIRS})

    # Local variable to hold all of the libraries that will be directly linked
    # to this library.
    set(LINK_LIBS)

    # Add dependent libraries passed directly in

    if (PARSE_DEPLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "DEPLIBS = ${PARSE_DEPLIBS}")
    endif()
    if (PARSE_IMPORTEDLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "IMPORTEDLIBS = ${PARSE_IMPORTEDLIBS}")
    endif()

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

    set(ADD_DEP_PACKAGE_AND_TPL_LIBS TRUE)

    set(PREFIXED_DEPLIBS)

    foreach(LIB ${PARSE_DEPLIBS})

      set(PREFIXED_LIB "${LIBRARY_NAME_PREFIX}${LIB}")

      # LIB_IN_SE_PKG?
      list(FIND ${PACKAGE_NAME}_LIBRARIES ${PREFIXED_LIB} FOUND_IDX)
      if (FOUND_IDX GREATER -1)
        set(LIB_IN_SE_PKG TRUE)
      else()
        set(LIB_IN_SE_PKG FALSE)
      endif()

      # PREFIXED_LIB_IS_TESTONLY?
      if (${PREFIXED_LIB}_INCLUDE_DIRS)
        set(LIB_TESTONLY TRUE)
      else()
        set(LIB_TESTONLY FALSE)
      endif()

      # Check for valid usage (sorted by most common to least common)
      if (LIB_IN_SE_PKG AND NOT LIB_TESTONLY) #PARSE_TESTONLY=TRUE/FASLE
        # The library being created here is a library dependent on a regular
        # (non-TESTONLY) lib in this SE package.  This is valid usage of
        # DEPLIBS.  There is no need to link this new lib to the SE package's
        # upstream dependent SE package and TPL libraries because thse are
        # already linked into the lib ${LIB}.
        set(ADD_DEP_PACKAGE_AND_TPL_LIBS FALSE)
      elseif (PARSE_TESTONLY AND LIB_IN_SE_PKG AND NOT LIB_TESTONLY)
        # The library being created here is TESTONLY library and is
        # dependent on a regular (non-TESTONLY) lib.  This is valid usage of
        # DEPLIBS.  In the case of test-only libraries, we always link in
        # the upstream libs.
      elseif (PARSE_TESTONLY AND LIB_TESTONLY) # LIB_IN_SE_PKG=TRUE/FASLE
        # The library being created here is TESTONLY library and is dependent
        # on another TESTONLY library.  This is valid usage of DEPLIBS.  In
        # this case we just hope that this SE package correctly specified a
        # TEST dependency on the upstream SE package that owns this upstream
        # TESTONLY library.
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message("-- "
            "Adding include directories for TESTONLY ${PREFIXED_LIB}_INCLUDE_DIRS ...")
        endif()
        include_directories(${${PREFIXED_LIB}_INCLUDE_DIRS})
      elseif (NOT PARSE_TESTONLY AND LIB_TESTONLY) # LIB_IN_SE_PKG=TRUE/FASLE
        message(WARNING "WARNING: '${LIB}' in DEPLIBS is a TESTONLY lib"
          " and it is illegal to link to this non-TESTONLY library '${LIBRARY_NAME}'."
          "  Such usage is deprecated (and this warning will soon become an error)!"
          "  If this is a regular library in this SE package or in an dependent upstream SE"
          " package then TriBITS will link automatically to it.  If you remove this and it"
          " does not link, then you need to add a new SE package dependency to"
          " this SE package's dependencies file"
          " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
        # ToDo: Turn the above to FATAL_ERROR after dropping deprecated code
      elseif (NOT LIB_IN_SE_PKG AND TARGET ${PREFIXED_LIB} ) # PARSE_TESTONLY=TRUE/FALSE
        message(WARNING "WARNING: '${LIB}' in DEPLIBS is not"
          " a lib in this SE package but is a library defined in the current"
          " cmake project!  Such usage is  deprecated (and"
          " will result in a configure error soon).  If this is a library in"
          " a dependent upstream SE package, then simply remove it from this list."
          "  TriBITS automatically links in libraries in upstream SE packages."
          "  If you remove '${LIB}' from DEPLIBS and your code does"
          " not link, then you need to add a new SE package dependency to"
          " this SE package's dependencies file"
          " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
      elseif (NOT LIB_IN_SE_PKG AND NOT TARGET ${PREFIXED_LIB} )
        message(WARNING "WARNING: '${LIB}' in DEPLIBS is not"
          " a lib defined in the current cmake project!  Such usage is deprecated (and"
          " will result in a configure error soon).  If this is an external"
          " lib you are trying to link in, it should likely be handled as a TriBITS"
          " TPL.  Otherwise, it should be passed in through IMPORTEDLIBS.  However,"
          " the only case we have found where IMPORTEDLIBS had to be used instead of"
          " through a proper TriBITS TPL is the C math library 'm'.")
      else()
        message(WARNING "WARNING: The case PARSE_TESTONLY=${PARSE_TESTONLY},"
          " LIB_IN_SE_PKG=${LIB_IN_SE_PKG}, LIB_TESTONLY=${LIB_TESTONLY}, has"
          " not yet been handled!")
      endif()

      list(APPEND PREFIXED_DEPLIBS "${LIBRARY_NAME_PREFIX}${LIB}")

    endforeach()

    append_set(LINK_LIBS ${PREFIXED_DEPLIBS})

    #
    # Check IMPORTEDLIBS
    #

    foreach(IMPORTEDLIB ${PARSE_IMPORTEDLIBS})
      set(PREFIXED_LIB "${LIBRARY_NAME_PREFIX}${IMPORTEDLIB}")
      list(FIND ${PACKAGE_NAME}_LIBRARIES ${PREFIXED_LIB} FOUND_IDX)
      if (${PREFIXED_LIB}_INCLUDE_DIRS)
        message(WARNING "WARNING: '${IMPORTEDLIB}' in IMPORTEDLIBS is a TESTONLY lib"
          " and it is illegal to pass in through IMPORTEDLIBS!"
          "  Such usage is deprecated (and this warning will soon become an error)!"
          "  Should '${IMPORTEDLIB}' instead be passed through DEPLIBS?")
        # ToDo: Turn the above to FATAL_ERROR after dropping deprecated code
      elseif (FOUND_IDX GREATER -1)
        message(WARNING "WARNING: Lib '${IMPORTEDLIB}' in IMPORTEDLIBS is in"
        " this SE package and is *not* an external lib!"
        "  TriBITS takes care of linking against libs the current"
        " SE package automatically.  Please remove it from IMPORTEDLIBS!")
      elseif (TARGET ${PREFIXED_LIB})
        message(WARNING "WARNING: Lib '${IMPORTEDLIB}' being passed through"
        " IMPORTEDLIBS is *not* an external library but instead is a library"
        " defined in this CMake project!"
        "  TriBITS takes care of linking against libraries in dependent upstream"
        " SE packages.  If you want to link to a library in an upstream SE"
        " package then add the SE package name to the appropriate category"
        " in this SE package's dependencies file: "
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
      endif()
      # ToDo: Assert that this is not a test-only lib
      list(APPEND LINK_LIBS ${IMPORTEDLIB})
    endforeach()

    #
    # Link in the upstream TEST SE package and TPL libs
    #
    # We link these before those in the LIB SE package and TPL libs because
    # the TEST dependencies tend to be higher in the dependency tree.  It
    # should not really matter but it looks better on the link line.
    #

    if (PARSE_TESTONLY)

      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Pulling in header and libraries dependencies"
          " for TEST dependencies ...")
      endif()

      tribits_sort_and_append_package_include_and_link_dirs_and_libs(
        ${PACKAGE_NAME}  TEST  LINK_LIBS)

      tribits_sort_and_append_tpl_include_and_link_dirs_and_libs(
        ${PACKAGE_NAME}  TEST  LINK_LIBS)

    endif()

    #
    # Add the dependent LIB SE package and TPL libs
    #

    if (ADD_DEP_PACKAGE_AND_TPL_LIBS)

      # If there are no dependent libs passed in, then this library can not
      # possibly depend on the package's other libraries so we must link to
      # the dependent libraries in dependent libraries and TPLs.

      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Pulling in header and libraries dependencies"
          " for LIB dependencies ...")
      endif()

      #
      # Call include_directories() and link_directories(...) for all upstream
      # dependent Packages and TPLs and accumulate the libraries to link against.
      #
      # NOTE: Adding these directories serves two purposes.  First, so that the includes
      # get added the the sources that get built for this library.  Second, so
      # that list full list of include directories can be extracted as a
      # property and set on ${PACKAGE_NAME}_INCLUDE_DIRS
      #

      tribits_sort_and_append_package_include_and_link_dirs_and_libs(
        ${PACKAGE_NAME}  LIB  LINK_LIBS)

      tribits_sort_and_append_tpl_include_and_link_dirs_and_libs(
        ${PACKAGE_NAME}  LIB  LINK_LIBS)

    endif()

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(LINK_LIBS)
    endif()

    # Add the library and all the dependencies

    if (PARSE_DEFINES)
      add_definitions(${PARSE_DEFINES})
    endif()

    if (PARSE_STATIC)
      set(STATIC_KEYWORD "STATIC")
    else()
      set(STATIC_KEYWORD)
    endif()

    if (PARSE_SHARED)
      set(SHARED_KEYWORD "SHARED")
    else()
      set(SHARED_KEYWORD)
    endif()

    if (NOT PARSE_CUDALIBRARY)
      add_library(
        ${LIBRARY_NAME}
        ${STATIC_KEYWORD}
        ${SHARED_KEYWORD}
        ${PARSE_HEADERS}
        ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES}
        )
    else()
      cuda_add_library(
        ${LIBRARY_NAME}
        ${PARSE_HEADERS}
        ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES}
        )
    endif()

    if(PARSE_ADDED_LIB_TARGET_NAME_OUT)
      set(${PARSE_ADDED_LIB_TARGET_NAME_OUT} ${LIBRARY_NAME} PARENT_SCOPE)
    endif()

    set_property(
      TARGET ${LIBRARY_NAME}
      APPEND PROPERTY
      LABELS ${PACKAGE_NAME}Libs ${PARENT_PACKAGE_NAME}Libs
      )

    if (NOT "${${PROJECT_NAME}_VERSION}" STREQUAL "" AND
      NOT "${${PROJECT_NAME}_MAJOR_VERSION}" STREQUAL ""
      )
      set_target_properties(
        ${LIBRARY_NAME}
        PROPERTIES
        VERSION ${${PROJECT_NAME}_VERSION}
        SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION}
        )
    endif()

    prepend_global_set(${PARENT_PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
    prepend_global_set(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})

    if (${PROJECT_NAME}_DUMP_LINK_LIBS)
      message("-- ${LIBRARY_NAME_IN}:LINK_LIBS='${LINK_LIBS}'")
    endif()

    target_link_libraries(${LIBRARY_NAME} PUBLIC ${LINK_LIBS})

    if (${PROJECT_NAME}_CXX_STANDARD_FEATURE)
      target_compile_features(${LIBRARY_NAME} PUBLIC "${${PROJECT_NAME}_CXX_STANDARD_FEATURE}")
    ENDIF ()

    # Add to the install target

    set(INSTALL_LIB ON)
    set(INSTALL_HEADERS ON)
    set(APPEND_LIB_AND_HEADERS_TO_PACKAGE ON)

    if (PARSE_TESTONLY)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Skipping installation hooks for this library"
          " because 'TESTONLY' was passed in ...")
      endif()
      set(INSTALL_LIB OFF)
      set(INSTALL_HEADERS OFF)
      set(APPEND_LIB_AND_HEADERS_TO_PACKAGE OFF)
    elseif (PARSE_NO_INSTALL_LIB_OR_HEADERS)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Skipping installation hooks for this library"
          " because 'NO_INSTALL_LIB_OR_HEADERS' was passed in ...")
      endif()
      set(INSTALL_LIB OFF)
      set(INSTALL_HEADERS OFF)
    elseif (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND NOT BUILD_SHARED_LIBS)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Skipping installation of headers and libraries"
          " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and BUILD_SHARED_LIBS=FALSE  ...")
      endif()
      set(INSTALL_LIB OFF)
      set(INSTALL_HEADERS OFF)
    elseif (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND BUILD_SHARED_LIBS)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Skipping installation of headers but installing libraries"
          " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and BUILD_SHARED_LIBS=TRUE  ...")
      endif()
      set(INSTALL_HEADERS OFF)
    endif()

    if (INSTALL_LIB OR INSTALL_HEADERS)
      set_property(GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS ON)
      set_property(GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS ON)
    endif()

    if (INSTALL_LIB)
      install(
        TARGETS ${LIBRARY_NAME}
        EXPORT ${PACKAGE_NAME}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}"
        LIBRARY DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
        ARCHIVE DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
        COMPONENT ${PACKAGE_NAME}
        )
    endif()

    if (INSTALL_HEADERS)
      tribits_install_headers(
        HEADERS  ${PARSE_HEADERS}
        INSTALL_SUBDIR  ${PARSE_HEADERS_INSTALL_SUBDIR}
        COMPONENT  ${PACKAGE_NAME}
        )
    endif()

    # Append the new include dirs, library dirs, and libraries to this package's lists

    get_directory_property(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    get_directory_property(LIBRARY_DIRS_CURRENT  PACKAGE_LIBRARY_DIRS)

    if (APPEND_LIB_AND_HEADERS_TO_PACKAGE)

      prepend_global_set(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCLUDE_DIRS_CURRENT})
      prepend_global_set(${PACKAGE_NAME}_LIBRARY_DIRS  ${LIBRARY_DIRS_CURRENT})
      prepend_global_set(${PACKAGE_NAME}_LIBRARIES  ${LIBRARY_NAME})

      remove_global_duplicates(${PACKAGE_NAME}_INCLUDE_DIRS)
      remove_global_duplicates(${PACKAGE_NAME}_LIBRARY_DIRS)
      remove_global_duplicates(${PACKAGE_NAME}_LIBRARIES)

      if (INSTALL_LIB)
        global_set(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES_TO_INSTALL TRUE)
      endif()

    else()

      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Skipping augmentation of package's lists of include"
          " directories and libraries! ...")
      endif()

      list(REMOVE_DUPLICATES INCLUDE_DIRS_CURRENT)
      global_set(${LIBRARY_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})

      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        print_var(${LIBRARY_NAME}_INCLUDE_DIRS)
      endif()

    endif()
  endif() #if not in installation testing mode

  #
  # Adjust for installation testing
  #

  if (${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)

    list(FIND ${PROJECT_NAME}_INSTALLATION_PACKAGE_LIST ${PACKAGE_NAME}
      ${PACKAGE_NAME}_WAS_INSTALLED)
    if(${${PACKAGE_NAME}_WAS_INSTALLED} EQUAL -1)
      message(FATAL_ERROR
        "The package ${PACKAGE_NAME} was not installed with ${PROJECT_NAME}!"
        "  Please disable package ${PACKAGE_NAME} or install it.")
    endif()

    include_directories(REQUIRED_DURING_INSTALLATION_TESTING  BEFORE
       ${${PACKAGE_NAME}_INSTALLATION_INCLUDE_DIRS}
       ${${TRIBITS_PACKAGE}_INSTALLATION_TPL_INCLUDE_DIRS})
    set_property(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
      ${${PACKAGE_NAME}_INSTALLATION_LIBRARY_DIRS})

    get_directory_property(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
    get_directory_property(LIBRARY_DIRS_CURRENT PACKAGE_LIBRARY_DIRS)

    global_set(${PACKAGE_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})
    global_set(${PACKAGE_NAME}_LIBRARY_DIRS ${LIBRARY_DIRS_CURRENT})
    global_set(${PACKAGE_NAME}_LIBRARIES    ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})

  endif() #installation testing mode

  #
  # Print the updates to the linkage variables
  #

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(${PACKAGE_NAME}_INCLUDE_DIRS)
    print_var(${PACKAGE_NAME}_LIBRARY_DIRS)
    print_var(${PACKAGE_NAME}_LIBRARIES)
  endif()

endfunction()
