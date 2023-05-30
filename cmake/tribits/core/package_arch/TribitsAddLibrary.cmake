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

include(TribitsLibIsTestOnly)
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
include(TribitsDeprecatedHelpers)
include(TribitsSetAndIncDirs)


# @FUNCTION: tribits_add_library()
#
# Function used to add a CMake library and target using ``add_library()`` and
# also the ALIAS target ``${PACKAGE_NAME}::<libname>`` (where ``<libname>`` is
# the full CMake target name as returned from ``${<libTargetName>}``).
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
#     List of dependent libraries that are built in the current package that
#     this library is dependent on.  These libraries are passed into
#     ``target_link_libraries(<libTargetName> ...)`` so that CMake knows about
#     the dependency structure of the libraries within this package.
#     **NOTE:** One must **not** list libraries in other upstream `TriBITS
#     Packages`_ or libraries built externally from this TriBITS CMake project
#     in ``DEPLIBS``.  The TriBITS system automatically handles linking to
#     libraries in upstream TriBITS packages.  External libraries need to be
#     listed in the ``IMPORTEDLIBS`` argument instead if they are not already
#     specified automatically using a `TriBITS TPL`_.
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
#     downstream libraries in the current package and all downstream packages
#     must also be also be marked with ``STATIC``.  That is because,
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
#     in through its ``DEPLIBS`` argument.  Also, the custom property
#     ``TRIBITS_TESTONLY_LIB`` will be set to ``TRUE`` which will ensure that
#     this library will not be added to the ``${PACKAGE_NAME}::all_libs``
#     target.
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
#     TPLs`_).  For this option to work, this package must have an enabled
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
# and will add this to the generated library target using
# ``target_link_libraries()`` so that they get propagated downstream as well.
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
# `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_ is ``FALSE`` and
# ``BUILD_SHARD_LIBS=OFF``.  But when ``BUILD_SHARD_LIBS=ON``, the install
# target will get added.  Also, this install target will *not* get added if
# ``TESTONLY`` or ``NO_INSTALL_LIB_OR_HEADERS`` are passed in.
#
# By default, an install target for the headers listed in ``HEADERS`` will get
# added using ``install(FILES <h0> <h1> ...)``, but only if ``TESTONLY`` and
# ``NO_INSTALL_LIB_OR_HEADERS`` are not passed in as well.  Also, the install
# target for the headers will not get added if
# `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_ is ``FALSE``.  If this
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
function(tribits_add_library  LIBRARY_NAME_IN)

  tribits_add_library_assert_correct_call_context()

  # Set library prefix and name

  set(LIBRARY_NAME_PREFIX "${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}")
  set(LIBRARY_NAME ${LIBRARY_NAME_PREFIX}${LIBRARY_NAME_IN})

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_ADD_LIBRARY: ${LIBRARY_NAME}")
    if(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
      message("\n${PACKAGE_NAME}_LIBRARIES In installation testing mode,"
        " libraries will be found instead of created.")
    endif()
    print_var(${PACKAGE_NAME}_LIBRARIES)
  endif()

  # Parse input args

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

  # Library not added by default
  if(PARSE_ADDED_LIB_TARGET_NAME_OUT)
    set(${PARSE_ADDED_LIB_TARGET_NAME_OUT} "" PARENT_SCOPE)
  endif()

  #
  # Create library target if not doing installation testing or if this is a
  # TESTONLY library.
  #

  if (NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_TESTONLY)

    if (PARSE_DEPLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "DEPLIBS = ${PARSE_DEPLIBS}")
    endif()
    if (PARSE_IMPORTEDLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "IMPORTEDLIBS = ${PARSE_IMPORTEDLIBS}")
    endif()

    # Assert DEPLIBS and IMPORTEDLIBS

    tribits_add_library_assert_deplibs()
    tribits_add_library_assert_importedlibs()

    # Add the library and all the dependencies

    if (PARSE_DEFINES)
      add_definitions(${PARSE_DEFINES})
    endif()

    if (PARSE_STATIC)
      set(STATIC_KEYWORD "STATIC")
    else()
      set(STATIC_KEYWORD "")
    endif()

    if (PARSE_SHARED)
      set(SHARED_KEYWORD "SHARED")
    else()
      set(SHARED_KEYWORD "")
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

    if (PARSE_TESTONLY)
      tribits_set_lib_is_testonly(${LIBRARY_NAME})
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

    #
    # Link ${LIBRARY_NAME} to direct upstream libraries
    #

    # DEPLIBS
    foreach(depLib ${PARSE_DEPLIBS})
      target_link_libraries(${LIBRARY_NAME} PUBLIC "${LIBRARY_NAME_PREFIX}${depLib}")
    endforeach()
    # ${PACKAGE_NAME}_LIBRARIES
    target_link_libraries(${LIBRARY_NAME} PUBLIC ${${PACKAGE_NAME}_LIBRARIES})
    # ${PACKAGE_NAME}_LIB_ENABLED_DEPENDENCIES
    foreach(depPkg IN LISTS ${PACKAGE_NAME}_LIB_ENABLED_DEPENDENCIES)
      target_link_libraries(${LIBRARY_NAME} PUBLIC ${depPkg}::all_libs)
    endforeach()
    # ${PACKAGE_NAME}_TEST_ENABLED_DEPENDENCIES (TESTONLY lib)
    if (PARSE_TESTONLY)
      foreach(depPkg IN LISTS ${PACKAGE_NAME}_TEST_ENABLED_DEPENDENCIES)
        target_link_libraries(${LIBRARY_NAME} PUBLIC ${depPkg}::all_libs)
      endforeach()
    endif()
    # IMPORTEDLIBS
    foreach(importedLib ${PARSE_IMPORTEDLIBS})
      target_link_libraries(${LIBRARY_NAME} PUBLIC "${importedLib}")
    endforeach()

    # ToDo: #63: Above, allow for other link visibilities other than 'PUBLIC'!

    if (${PROJECT_NAME}_CXX_STANDARD_FEATURE)
      target_compile_features(${LIBRARY_NAME} PUBLIC
        "${${PROJECT_NAME}_CXX_STANDARD_FEATURE}")
    endif()

    # Add to the install target

    tribits_add_library_determine_install_lib_and_or_headers(
      installLib  installHeaders  appendLibAndHeadersToPackageVars)

    if (installLib OR installHeaders)
      set_property(GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS ON)
      set_property(GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS ON)
    endif()

    if (installLib)
      install(
        TARGETS ${LIBRARY_NAME}
        EXPORT ${PACKAGE_NAME}
        INCLUDES DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}"
        LIBRARY DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
        ARCHIVE DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
        COMPONENT ${PACKAGE_NAME}
        )
    endif()

    if (installHeaders)
      tribits_install_headers(
        HEADERS  ${PARSE_HEADERS}
        INSTALL_SUBDIR  ${PARSE_HEADERS_INSTALL_SUBDIR}
        COMPONENT  ${PACKAGE_NAME}
        )
    endif()

    # Append the new libraries to this package's lists

    if (appendLibAndHeadersToPackageVars)
      prepend_global_set(${PACKAGE_NAME}_LIBRARIES  ${PACKAGE_NAME}::${LIBRARY_NAME})
      remove_global_duplicates(${PACKAGE_NAME}_LIBRARIES)
      if (installLib)
        global_set(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES_TO_INSTALL TRUE)
      endif()
    else()
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "Skipping augmentation of package's lists of libraries! ...")
      endif()
    endif()

    # Set INTERFACE_INCLUDE_DIRECTORIES property for added library and must
    # only do for the build interface (not the install interface).
    get_directory_property(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    set(buildInterfaceIncludeDirs)
    foreach (includeDir IN LISTS INCLUDE_DIRS_CURRENT)
      list(APPEND buildInterfaceIncludeDirs "$<BUILD_INTERFACE:${includeDir}>")
    endforeach()
    target_include_directories( ${LIBRARY_NAME} PUBLIC ${buildInterfaceIncludeDirs} )

    # Add ALIAS library <PackageName>::<libname>
    add_library(${PACKAGE_NAME}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

    # Optionally Set IMPORTED_NO_SYSTEM
    if (${PROJECT_NAME}_IMPORTED_NO_SYSTEM)
      set_target_properties(${LIBRARY_NAME} PROPERTIES IMPORTED_NO_SYSTEM TRUE)
    endif()

  endif() #if not in installation testing mode

  # Adjust for installation testing

  if (${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    list(FIND ${PROJECT_NAME}_INSTALLATION_PACKAGE_LIST ${PACKAGE_NAME}
      ${PACKAGE_NAME}_WAS_INSTALLED)
    if(${${PACKAGE_NAME}_WAS_INSTALLED} EQUAL -1)
      message(FATAL_ERROR
        "The package ${PACKAGE_NAME} was not installed with ${PROJECT_NAME}!"
        "  Please disable package ${PACKAGE_NAME} or install it.")
    endif()
    global_set(${PACKAGE_NAME}_LIBRARIES  ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})
  endif()

  # Print the updates to the linkage variables

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(${PACKAGE_NAME}_LIBRARIES)
  endif()

endfunction()
#
# ToDo:, above Deprecate and remove the usage of DEFINES!  People should be
# putting defines into configured header files, not adding -D<macroName> to
# the compile lines!


# Function that asserts that tribits_add_library() is called in the correct
# context
#
function(tribits_add_library_assert_correct_call_context)

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

endfunction()


# Assert correct DEPSLIB passed to tribits_add_library()
#
# NOTE: This accesses vars from the enclosed calling function
# tribits_add_library() but does not set any variables in that scope!
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
function(tribits_add_library_assert_deplibs)

  foreach(depLib ${PARSE_DEPLIBS})

    set(prefixedDepLib "${LIBRARY_NAME_PREFIX}${depLib}")

    # Is this lib already listed in ${PACKAGE_NAME}_LIBS?
    list(FIND ${PACKAGE_NAME}_LIBRARIES "${PACKAGE_NAME}::${prefixedDepLib}" FOUND_IDX)
    if (FOUND_IDX GREATER -1)
      set(depLibAlreadyInPkgLibs TRUE)
    else()
      set(depLibAlreadyInPkgLibs FALSE)
    endif()

    # ${PREFIXED_LIB" is TESTONLY?
    tribits_lib_is_testonly(${prefixedDepLib} depLibIsTestOnlyLib)

    # Check for valid usage (sorted by most common to least common)
    if (depLibAlreadyInPkgLibs AND NOT depLibIsTestOnlyLib) # PARSE_TESTONLY=any
      # The library being created here is a (regular or testonly) library
      # dependent on a regular (non-TESTONLY) lib in this package.  This is
      # valid usage of DEPLIBS.  There is no need to link this new lib to
      # the package's upstream dependent package and TPL libraries because
      # these are already linked into one of the of the package's own
      # upstream libs.
    elseif (PARSE_TESTONLY AND depLibAlreadyInPkgLibs AND NOT depLibIsTestOnlyLib)
      # The library being created here is TESTONLY library and is
      # dependent on a regular (non-TESTONLY) lib.  This is valid usage of
      # DEPLIBS.  In the case of test-only libraries, we always link in
      # the upstream libs.
    elseif (PARSE_TESTONLY AND depLibIsTestOnlyLib) # any depLibAlreadyInPkgLibs
      # The library being created here is TESTONLY library and is dependent
      # on another TESTONLY library.  This is valid usage of DEPLIBS.  In
      # this case we just hope that this package correctly specified a TEST
      # dependency on the upstream package that owns this upstream TESTONLY
      # library if it comes from an upstream package.
    elseif (NOT PARSE_TESTONLY AND depLibIsTestOnlyLib) # any depLibAlreadyInPkgLibs
      tribits_deprecated("'${depLib}' in DEPLIBS is a TESTONLY lib"
        " and it is illegal to link to this non-TESTONLY library '${LIBRARY_NAME}'."
        "  Such usage is deprecated (and this warning will soon become an error)!"
        "  If this is a regular library in this package or in an dependent upstream"
        " package then TriBITS will link automatically to it.  If you remove this and it"
        " does not link, then you need to add a new package dependency to"
        " this package's dependencies file"
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake")
      # ToDo: Turn the above to FATAL_ERROR after dropping deprecated code
    elseif (NOT depLibAlreadyInPkgLibs AND TARGET ${prefixedDepLib}) # any PARSE_TESTONLY
      tribits_deprecated("'${depLib}' in DEPLIBS is not"
        " a lib in this package but is a library defined in the current"
        " cmake project!  Such usage is deprecated (and"
        " will result in a configure error soon).  If this is a library in"
        " a dependent upstream package, then simply remove '${depLib}' from this list."
        "  TriBITS automatically links in libraries in upstream packages."
        "  If you remove '${depLib}' from DEPLIBS and your code does"
        " not link, then you need to add a new package dependency to"
        " this package's dependencies file"
        " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake.")
    elseif (NOT depLibAlreadyInPkgLibs AND NOT TARGET ${prefixedDepLib} )
      tribits_deprecated("'${depLib}' in DEPLIBS is not"
        " a lib defined in the current cmake project!  Such usage is deprecated (and"
        " will result in a configure error soon).  If this is an external"
        " lib you are trying to link in, it should likely be handled as a TriBITS"
        " TPL.  Otherwise, it should be passed in through IMPORTEDLIBS.  However,"
        " the only case we have found where IMPORTEDLIBS had to be used instead of"
        " through a proper TriBITS TPL is the C math library 'm'.")
    else()
      message(WARNING "WARNING: The case PARSE_TESTONLY=${PARSE_TESTONLY},"
        " depLibAlreadyInPkgLibs=${depLibAlreadyInPkgLibs},"
          " depLibIsTestOnlyLib=${depLibIsTestOnlyLib}, has"
        " not yet been handled!")
    endif()

  endforeach()

endfunction()


# Assert correct IMPORTEDLIBS passed to tribits_add_library()
#
# NOTE: This accesses vars from the enclosed calling function
# tribits_add_library() but does not set any variables in that scope!
#
# ToDo: Turn the below deprecated WARNING messages to FATAL_ERROR once we
# give enough time for people to clean up their codes.
#
function(tribits_add_library_assert_importedlibs)
  foreach(importedLib ${PARSE_IMPORTEDLIBS})
    set(prefixedImportedLib "${LIBRARY_NAME_PREFIX}${importedLib}")
    list(FIND ${PACKAGE_NAME}_LIBRARIES "${PACKAGE_NAME}::${prefixedImportedLib}"
      FOUND_IMPORTEDLIB_IN_LIBRARIES_IDX)
    tribits_lib_is_testonly(${prefixedImportedLib}  importedLibIsTestOnlyLib)
    if (importedLibIsTestOnlyLib)
      tribits_deprecated("'${importedLib}' in IMPORTEDLIBS is a TESTONLY lib"
        " and it is illegal to pass in through IMPORTEDLIBS!"
        "  Such usage is deprecated (and this warning will soon become an error)!"
        "  Should '${importedLib}' instead be passed through DEPLIBS?")
      # ToDo: Turn the above to FATAL_ERROR after dropping deprecated code
    elseif (FOUND_IMPORTEDLIB_IN_LIBRARIES_IDX GREATER -1)
      message(WARNING "WARNING: Lib '${importedLib}' in IMPORTEDLIBS is in"
      " this package and is *not* an external lib!"
      "  Please move '${importedLib}' from the list IMPORTEDLIBS to DEPLIBS.")
    elseif (TARGET ${prefixedImportedLib})
      message(WARNING "WARNING: Lib '${importedLib}' being passed through"
      " IMPORTEDLIBS is *not* an external library but instead is a library"
      " defined in this CMake project!"
      "  TriBITS takes care of linking against libraries in dependent upstream"
      " packages.  If you want to link to a library in an upstream"
      " package then add the package name to the appropriate category"
      " in this package's dependencies file: "
      " ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/Dependencies.cmake.")
    endif()
  endforeach()
endfunction()


# Determine lib and/or headers should be installed and if the package vars
# should be updated
#
# NOTE: This reads the parsed arguments (prefixed with ``PARSE_``) from the
# calling tribits_add_library() function from the enclosing scope.
#
function(tribits_add_library_determine_install_lib_and_or_headers
    installLibOut  installHeadersOut  appendLibAndHeadersToPackageVarsOut
  )

  set(installLib ON)
  set(installHeaders ON)
  set(appendLibAndHeadersToPackageVars ON)

  if (PARSE_TESTONLY)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "Skipping installation hooks for this library"
        " because 'TESTONLY' was passed in ...")
    endif()
    set(installLib OFF)
    set(installHeaders OFF)
    set(appendLibAndHeadersToPackageVars OFF)
  elseif (PARSE_NO_INSTALL_LIB_OR_HEADERS)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "Skipping installation hooks for this library"
        " because 'NO_INSTALL_LIB_OR_HEADERS' was passed in ...")
    endif()
    set(installLib OFF)
    set(installHeaders OFF)
  elseif (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND NOT BUILD_SHARED_LIBS)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "Skipping installation of headers and libraries"
        " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and"
          " BUILD_SHARED_LIBS=FALSE ...")
    endif()
    set(installLib OFF)
    set(installHeaders OFF)
  elseif (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND BUILD_SHARED_LIBS)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "Skipping installation of headers but installing libraries"
        " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and"
        " BUILD_SHARED_LIBS=TRUE ...")
    endif()
    set(installHeaders OFF)
  endif()

  set(${installLibOut} ${installLib} PARENT_SCOPE)
  set(${installHeadersOut} ${installHeaders} PARENT_SCOPE)
  set(${appendLibAndHeadersToPackageVarsOut} ${appendLibAndHeadersToPackageVars}
    PARENT_SCOPE)

endfunction()
