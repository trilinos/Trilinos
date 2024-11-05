# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

include(TribitsExternalPackageWriteConfigFile)

include(AdvancedSet)
include(AppendSet)
include(AssertDefined)
include(DualScopeSet)
include(GlobalNullSet)
include(GlobalSet)
include(MultilineSet)
include(CMakeParseArguments)
include(SetNotFound)
include(Split)


# @FUNCTION: tribits_tpl_allow_pre_find_package()
#
# Function that determines if a TriBITS find module file
# ``FindTPL<tplName>.cmake`` is allowed to call ``find_package(<tplName>
# ...)`` before calling `tribits_tpl_find_include_dirs_and_libraries()`_.
#
# Usage::
#
#   tribits_tpl_allow_pre_find_package( <tplName>
#     <allowPackagePrefindOut> )
#
# The required arguments are:
#
#   ``<tplName>`` : The input name of the TriBITS TPL (e.g. ``HDF5``).
#
#   ``<allowPackagePrefindOut>`` : Name of a variable which will be set to
#   ``TRUE`` on output if ``find_package(<tplName> ...)`` should be called to
#   find the TPL ``<tplName>`` or ``FALSE`` if it should not be called.
#
# This function will set ``<allowPackagePrefindOut>`` to ``FALSE`` if any of
# the variables ``TPL_<tplName>_INCLUDE_DIRS``, ``${TPL_<tplName>_LIBRARIES``,
# or ``TPL_<tplName>_LIBRARY_DIRS`` are set.  This allows the user to override
# the search for the library components and just specify the absolute
# locations.  The function will also set ``<allowPackagePrefindOut>`` to
# ``FALSE`` if ``<tplName>_INCLUDE_DIRS``, ``<tplName>_LIBRARY_NAMES``, or
# ``<tplName>_LIBRARY_DIRS`` is set and ``<tplName>_FORCE_PRE_FIND_PACKAGE``
# is set to ``FALSE``.  Otherwise, if ``<tplName>_FORCE_PRE_FIND_PACKAGE`` is
# set to ``TRUE``, the function will not return ``FALSE`` for
# ``<allowPackagePrefindOut>`` no matter what the values of
# ``<tplName>_INCLUDE_DIRS``, ``<tplName>_LIBRARY_NAMES``, or
# ``<tplName>_LIBRARY_DIRS``.  Finally, ``<allowPackagePrefindOut>`` is set to
# ``FALSE`` if ``<tplName>_ALLOW_PACKAGE_PREFIND=OFF`` is set in the cache.
#
# The variable ``<tplName>_FORCE_PRE_FIND_PACKAGE`` is needed to allow users
# (or the ``FindTPL<tplName>.cmake`` module itself) to avoid name clashes with
# the variables ``<tplName>_INCLUDE_DIRS`` or ``<tplName>_LIBRARY_DIRS`` in
# the usage of ``find_package(<tplName> ...)`` because a lot of default
# ``Find<tplName>.cmake`` modules also use these variables.  This function
# sets ``<tplName>_FORCE_PRE_FIND_PACKAGE`` as a cache variable with default
# value ``FALSE`` to maintain backward compatibility with existing
# ``FindTPL<tplName>.cmake`` modules.
#
# The cache variable ``<tplName>_ALLOW_PACKAGE_PREFIND`` is to allow the user
# to disable the prefind call to ``find_package()`` even if it would be
# allowed otherwise.
#
# See `Creating FindTPL<tplName>.cmake using find_package() without IMPORTED
# targets`_ for details in how to use this function to create a
# ``FindTPL<tplName>.cmake`` module file.
#
function(tribits_tpl_allow_pre_find_package  TPL_NAME  ALLOW_PACKAGE_PREFIND_OUT)

  if (TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE_DEBUG)
    message("TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE: '${TPL_NAME}'  '${ALLOW_PACKAGE_PREFIND_OUT}'")
    print_var(${TPL_NAME}_INCLUDE_DIRS)
    print_var(${TPL_NAME}_LIBRARY_NAMES)
    print_var(${TPL_NAME}_LIBRARY_DIRS)
    print_var(${TPL_NAME}_FORCE_PRE_FIND_PACKAGE)
  endif()

  advanced_set(${TPL_NAME}_FORCE_PRE_FIND_PACKAGE  FALSE
    CACHE BOOL
    "Determines if the variables ${TPL_NAME}_[INCLUDE_DIRS,LIBRARY_NAMES,LIBRARY_DIRS] should be ignored and the pre-find find_package(${TPL_NAME} should be performed anyway.  But this will *not* do the pre-find if any of the TPL_${TPL_NAME}_[INCLUDE_DIRS,LIBRARY_NAMES,LIBRARY_DIRS] vars are set." )

  # Start out with TRUE and set to FALSE in logic below
  set(${TPL_NAME}_ALLOW_PACKAGE_PREFIND TRUE CACHE BOOL
     "Set to FALSE to skip find_package() prefind for the TriBITS TPL '${TPL_NAME}'")

  # Start out with TRUE and set to FALSE in logic below
  set(ALLOW_PACKAGE_PREFIND TRUE)

  if (NOT ${TPL_NAME}_ALLOW_PACKAGE_PREFIND)
    set(ALLOW_PACKAGE_PREFIND FALSE)
  elseif (
    (NOT "${TPL_${TPL_NAME}_INCLUDE_DIRS}" STREQUAL "")
    OR (NOT "${TPL_${TPL_NAME}_LIBRARIES}" STREQUAL "")
    OR (NOT "${TPL_${TPL_NAME}_LIBRARY_DIRS}" STREQUAL "")
    )
    # The user has selected one or more of the final vars so skip calling
    # find_package(${TPL_NAME} ...) ...
    set(ALLOW_PACKAGE_PREFIND FALSE)
  elseif (
    (NOT "${${TPL_NAME}_INCLUDE_DIRS}" STREQUAL "")
    OR (NOT "${${TPL_NAME}_LIBRARY_NAMES}" STREQUAL "")
    OR (NOT "${${TPL_NAME}_LIBRARY_DIRS}" STREQUAL "")
    )
    # One or more of the ${TPL_NAME}_XXX variables are set
    if (${TPL_NAME}_FORCE_PRE_FIND_PACKAGE)
      # Even with one or more of the ${TPL_NAME}_XXX vars set, we still want
      # to do the find_package(${TPL_NAME} ...) search and ignore this
      # override.
    else()
      # We will not ignore the override of these variables and will instead go
      # ahead and skip the pre-find.
      set(ALLOW_PACKAGE_PREFIND FALSE)
    endif()
  endif()

  set(${ALLOW_PACKAGE_PREFIND_OUT} ${ALLOW_PACKAGE_PREFIND} PARENT_SCOPE)

endfunction()


# @FUNCTION: tribits_tpl_find_include_dirs_and_libraries()
#
# This function reads (cache) variables that specify where to find a `TriBITS
# TPL`_'s headers and libraries and then creates IMPORTED targets, the
# ``<tplName>::all_libs`` target, and writes the file
# ``<tplName>Config.cmake`` into the standard location in the build directory.
# This function is typically called inside of a ``FindTPL<tplName>.cmake``
# module file (see `${TPL_NAME}_FINDMOD`_).
#
# Usage::
#
#   tribits_tpl_find_include_dirs_and_libraries(
#     <tplName>
#     [REQUIRED_HEADERS <header1> <header2> ...]
#     [MUST_FIND_ALL_HEADERS]
#     [REQUIRED_LIBS_NAMES <libname1> <libname2> ...]
#     [MUST_FIND_ALL_LIBS]
#     [NO_PRINT_ENABLE_SUCCESS_FAIL]
#     )
#
# This function can be called to specify/require header files and include
# directories and/or a list of libraries.
#
# The input arguments to this function are:
#
#   ``<tplName>``
#
#     Name of the TPL that is listed in a `<repoDir>/TPLsList.cmake`_ file.
#
#   ``REQUIRED_HEADERS``
#
#     List of header files that are searched in order to find the TPL's
#     include directories files using ``find_path()``. 
#
#   ``MUST_FIND_ALL_HEADERS``
#
#     If set, then all of the header files listed in ``REQUIRED_HEADERS`` must
#     be found (unless ``TPL_<tplName>_INCLUDE_DIRS`` is already set).
#
#   ``REQUIRED_LIBS_NAMES``
#
#     List of libraries that are searched for when looking for the TPL's
#     libraries using ``find_library()``.  This list can be overridden by the
#     user by setting ``<tplName>_LIBRARY_NAMES`` (see below).
#
#   ``MUST_FIND_ALL_LIBS``
#
#     If set, then all of the library files listed in ``REQUIRED_LIBS_NAMES``
#     must be found or the TPL is considered not found (unless
#     ``TPL_<tplName>_LIBRARIES`` is already set).  If the global cache var
#     ``<Project>_MUST_FIND_ALL_TPL_LIBS`` is set to ``TRUE``, then this is
#     turned on as well.  WARNING: The default is not to require finding all
#     of the listed libs.  (This is to maintain backward compatibility with
#     some older ``FindTPL<tplName>.cmake`` modules.)
#
#   ``NO_PRINT_ENABLE_SUCCESS_FAIL``
#
#      If set, then the final success/fail will not be printed
#
# This function implements the TPL find behavior described in `Enabling
# support for an optional Third-Party Library (TPL)`_.
#
# The following (cache) variables, if set, will be used by this function:
#
#   ``<tplName>_INCLUDE_DIRS`` (type ``PATH``)
#
#     List of paths to search first for header files defined in
#     ``REQUIRED_HEADERS <header1> <header2> ...``.
#
#   ``<tplName>_LIBRARY_DIRS`` (type ``PATH``)
#
#     The list of directories to search first for libraries defined in
#     ``REQUIRED_LIBS_NAMES <libname1> <libname2> ...``.  If, for some reason,
#     no libraries should be linked in for this particular configuration, then
#     setting ``<tplName>_LIBRARY_DIRS=OFF`` or is empty will no special paths
#     will be searched.
#
#   ``<tplName>_LIBRARY_NAMES`` (type ``STRING``)
#
#     List of library names to be looked for instead of what is specified in
#     ``REQUIRED_LIBS_NAMES <libname1> <libname2> ...``.
#
#   ``<tplName>_LIB_ENABLED_DEPENDENCIES``
#
#     List of direct upstream external package/TPL dependencies that also
#     define ``<upstreamTplName>::all_libs`` targets.
#
# An addition, the function will avoid calling the find operations if the
# following (cache) variables are set on input:
#
#   ``TPL_<tplName>_INCLUDE_DIRS`` (type ``PATH``)
#
#     A list of common-separated full directory paths that contain the TPL's
#     header files.
#
#   ``TPL_<tplName>_LIBRARIES`` (type ``FILEPATH``)
#
#     A list of commons-separated full library names (i.e. output from
#     ``find_library()``) for all of the libraries for the TPL.
#
# This function produces the following:
#
#   ``TPL_<tplName>_NOT_FOUND`` (type ``BOOL``)
#
#     Will be set to ``ON`` if all of the parts of the TPL could not be found.
#
#   ``<tplName>::<libname>``
#
#     Namespaced IMPORTED target for every library found or specified in
#     ``TPL_<tplName>_LIBRARIES``.  These IMPORTED targets will have the
#     ``<upstreamTplName>::all_libs`` for the upstream external packages/TPLs
#     listed in ``<tplName>_LIB_ENABLED_DEPENDENCIES``.
#
#   ``<tplName>::all_libs``
#
#     INTERFACE target that depends on all of the created IMPORTED targets.
#
#   ``<buildDir>/external_packages/<tplName>/<tplName>Config.cmake``
#
#     A package configure file that contains all of the generated IMPORTED
#     targets ``<tplName>::<libname>`` and the ``<tplName>::all_libs`` target.
#     This fill will also call ``find_dependency()`` to pull in
#     ``<upstreamTplName>Config.cmake`` files for upstream TPLs that are
#     listed in ``<tplName>_LIB_ENABLED_DEPENDENCIES``.  (For more
#     information, see `tribits_extpkg_write_config_file()`_.)
#
# Note, if ``TPL_TENTATIVE_ENABLE_<tplName>=ON``, then if all of the parts of
# the TPL can't be found, then ``TPL_ENABLE_<tplName>`` will be (forced) set
# to ``OFF`` in the cache.  See `tribits_tpl_tentatively_enable()`_.
#
function(tribits_tpl_find_include_dirs_and_libraries TPL_NAME)

  # Make sure the right name is used
  assert_defined(TPL_ENABLE_${TPL_NAME})

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     "MUST_FIND_ALL_LIBS;MUST_FIND_ALL_HEADERS;NO_PRINT_ENABLE_SUCCESS_FAIL"
     #one_value_keywords
     ""
     #multi_value_keywords
     "REQUIRED_HEADERS;REQUIRED_LIBS_NAMES"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    set(TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE TRUE)
  endif()

  if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
    message("TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES: ${TPL_NAME}")
    print_var(PARSE_REQUIRED_HEADERS)
    print_var(PARSE_REQUIRED_LIBS_NAMES)
    print_var(TPL_${TPL_NAME}_INCLUDE_DIRS)
    print_var(TPL_${TPL_NAME}_LIBRARIES)
  endif()

  if (TPL_TENTATIVE_ENABLE_${TPL_NAME})
    message("-- Attempting to tentatively enable TPL '${TPL_NAME}' ...")
    set(ERROR_MSG_MODE)
  else()
    set(ERROR_MSG_MODE SEND_ERROR)
  endif()

  # Allow user override for finding libraries even if REQUIRED_LIBS_NAMES is
  # empty on input

  multiline_set(DOCSTR
    "List of semi-colon separated names of libraries needed to link to for"
    " the TPL ${TPL_NAME}.  This list of libraries will be search for in"
    " find_library(...) calls along with the directories specified with"
    " ${TPL_NAME}_LIBRARY_DIRS.  NOTE: This is not the final list of libraries"
    " used for linking.  That is specified by TPL_${TPL_NAME}_LIBRARIES!"
    )
  advanced_set(${TPL_NAME}_LIBRARY_NAMES ${PARSE_REQUIRED_LIBS_NAMES}
    CACHE STRING ${DOCSTR})
  split("${${TPL_NAME}_LIBRARY_NAMES}" "," ${TPL_NAME}_LIBRARY_NAMES)
  print_var(${TPL_NAME}_LIBRARY_NAMES)

  # Let the user override what the names of the libraries which might
  # actually mean that no libraries are searched for.
  set(REQUIRED_LIBS_NAMES ${${TPL_NAME}_LIBRARY_NAMES})

  #
  # User options
  #

  if (REQUIRED_LIBS_NAMES)

    # Library directories

    multiline_set(DOCSTR
      "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
      " libraries.  This list of paths will be passed into a find_library(...)"
      " command to find the libraries listed in ${TPL_NAME}_LIBRARY_NAMES."
      "  Note that this set of paths is also the default value used for"
      " ${TPL_NAME}_LIBRARY_DIRS.  Therefore, if the headers exist in the"
      " same directories as the library, you do not need to set"
      " ${TPL_NAME}_LIBRARY_DIRS."
      )

    advanced_set(${TPL_NAME}_LIBRARY_DIRS "" CACHE PATH ${DOCSTR})
    if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
      print_var(${TPL_NAME}_LIBRARY_DIRS)
    endif()

    # Libraries

    if (${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS)
      set(MUST_FIND_ALL_LIBS TRUE) 
    else()
      set(MUST_FIND_ALL_LIBS ${PARSE_MUST_FIND_ALL_LIBS}) 
    endif()

    if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
      print_var(${TPL_NAME}_LIBRARY_NAMES)
      print_var(REQUIRED_LIBS_NAMES)
    endif()

  else()

    set(${TPL_NAME}_LIBRARY_DIRS "") # Just to ignore below!

  endif()

  # Include directories

  if (PARSE_REQUIRED_HEADERS)

    multiline_set(DOCSTR
      "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
      " headers.  This list of paths will be passed into a find_path(...)"
      " command to find the headers for ${TPL_NAME} (which are known in advance)."
      )
    advanced_set(${TPL_NAME}_INCLUDE_DIRS ${${TPL_NAME}_LIBRARY_DIRS}
      CACHE PATH ${DOCSTR})

    if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
      print_var(${TPL_NAME}_LIBRARY_DIRS)
      print_var(${TPL_NAME}_INCLUDE_DIRS)
      print_var(PARSE_REQUIRED_HEADERS)
    endif()

  endif()

  #
  # Set the lib extensions to find
  #

  # Save the default the first time through
  if (NOT CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT)
   set(TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT ${CMAKE_FIND_LIBRARY_SUFFIXES})
   #print_var(TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT)
  endif()

  #print_var(TPL_FIND_SHARED_LIBS)
  #print_var(CMAKE_FIND_LIBRARY_SUFFIXES)
  # Set libraries to find
  if (TPL_FIND_SHARED_LIBS)
    # The default should be to find shared libs first
    set(TPL_CMAKE_FIND_LIBRARY_SUFFIXES ${TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT})
  else()
    if (WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a)
    else()
      set(CMAKE_FIND_LIBRARY_SUFFIXES .a )
    endif()
  endif()
  #print_var(CMAKE_FIND_LIBRARY_SUFFIXES)

  #
  # Direct build options
  #

  set(_${TPL_NAME}_ENABLE_SUCCESS TRUE)

  if (REQUIRED_LIBS_NAMES)

    # Libraries

    if (MUST_FIND_ALL_LIBS)
      set(LIB_NOT_FOUND_MSG_PREFIX "ERROR:")
    else()
      set(LIB_NOT_FOUND_MSG_PREFIX "NOTE:")
    endif()

    if (NOT TPL_${TPL_NAME}_LIBRARIES)

      if (MUST_FIND_ALL_LIBS)
        message("-- Must find at least one lib in each of the"
          " lib sets \"${REQUIRED_LIBS_NAMES}\"")
      endif()

      message( "-- Searching for libs in ${TPL_NAME}_LIBRARY_DIRS='${${TPL_NAME}_LIBRARY_DIRS}'")

      set(LIBRARIES_FOUND)

      foreach(LIBNAME_SET ${REQUIRED_LIBS_NAMES})

        message("-- Searching for a lib in the set \"${LIBNAME_SET}\":")

        if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
          print_var(LIBNAME_SET)
        endif()

        set(LIBNAME_LIST ${LIBNAME_SET})
        separate_arguments(LIBNAME_LIST)

        set(LIBNAME_SET_LIB)

        foreach(LIBNAME ${LIBNAME_LIST})

          message("--   Searching for lib '${LIBNAME}' ...")

          if (${TPL_NAME}_LIBRARY_DIRS)
            set(PATHS_ARG PATHS ${${TPL_NAME}_LIBRARY_DIRS})
          else()
            set(PATHS_ARG PATHS)
          endif()

          set_notfound(_${TPL_NAME}_${LIBNAME}_LIBRARY)
          find_library( _${TPL_NAME}_${LIBNAME}_LIBRARY
            NAMES ${LIBNAME}
            ${PATHS_ARG} NO_DEFAULT_PATH )
          find_library( _${TPL_NAME}_${LIBNAME}_LIBRARY
            NAMES ${LIBNAME} )
          mark_as_advanced(_${TPL_NAME}_${LIBNAME}_LIBRARY)

          if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
            print_var(_${TPL_NAME}_${LIBNAME}_LIBRARY)
          endif()

          if (_${TPL_NAME}_${LIBNAME}_LIBRARY)
            message("--     Found lib '${_${TPL_NAME}_${LIBNAME}_LIBRARY}'")
            set(LIBNAME_SET_LIB ${_${TPL_NAME}_${LIBNAME}_LIBRARY})
            break()
          endif()

        endforeach()

        if (NOT LIBNAME_SET_LIB)
          message(
            "-- ${LIB_NOT_FOUND_MSG_PREFIX} Did not find a lib in the lib set \"${LIBNAME_SET}\""
             " for the TPL '${TPL_NAME}'!")
          if (MUST_FIND_ALL_LIBS)
            set(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
          else()
            break()
          endif()
        endif()

        append_set(LIBRARIES_FOUND ${LIBNAME_SET_LIB})

      endforeach()

      multiline_set(DOCSTR
        "List of semi-colon separated full paths to the libraries for the TPL"
        " ${TPL_NAME}.  This is the final variable that is used in the link"
        " commands.  The user variable ${TPL_NAME}_LIBRARY_DIRS is used to look"
        " for the know library names but but is just a suggestion."
        " This variable, however, is the final value and will not be touched."
        )
      advanced_set( TPL_${TPL_NAME}_LIBRARIES ${LIBRARIES_FOUND}
        CACHE FILEPATH ${DOCSTR} FORCE)
      # Above, we have to force the set in case the find failed the last
      # configure in which case this cache var will be empty.  NOTE: If the
      # user specified a non-empty TPL_${TPL_NAME}_LIBRARIES, then we would
      # never get here in the first place!

      if (NOT TPL_${TPL_NAME}_LIBRARIES OR NOT _${TPL_NAME}_ENABLE_SUCCESS)
        message(
          "-- ERROR: Could not find the libraries for the TPL '${TPL_NAME}'!")
        message(
          "-- TIP: If the TPL '${TPL_NAME}' is on your system then you can set:\n"
          "     -D${TPL_NAME}_LIBRARY_DIRS='<dir0>;<dir1>;...'\n"
          "   to point to the directories where these libraries may be found.\n"
          "   Or, just set:\n"
          "     -DTPL_${TPL_NAME}_LIBRARIES='<path-to-libs0>;<path-to-libs1>;...'\n"
          "   to point to the full paths for the libraries which will\n"
          "   bypass any search for libraries and these libraries will be used without\n"
          "   question in the build.  (But this will result in a build-time error\n"
          "   if not all of the necessary symbols are found.)")
        tribits_tpl_find_include_dirs_and_libraries_handle_fail()
      endif()

    endif()

    # Print the final value to be used *always*
    message("-- TPL_${TPL_NAME}_LIBRARIES='${TPL_${TPL_NAME}_LIBRARIES}'")

  else()

    # There are no libraries so set the libraries to null but don't
    # change the cache which should not even have this variable in it.
    # This set command is only to follow the standards for the package
    # support CMake code.
    global_null_set(TPL_${TPL_NAME}_LIBRARIES)

  endif()

  # Include directories

  if (PARSE_REQUIRED_HEADERS)

    if (NOT TPL_${TPL_NAME}_INCLUDE_DIRS)

      if (PARSE_MUST_FIND_ALL_HEADERS)
        message("-- Must find at least one header in each of the"
          " header sets \"${PARSE_REQUIRED_HEADERS}\"")
      endif()

      message( "-- Searching for headers in ${TPL_NAME}_INCLUDE_DIRS='${${TPL_NAME}_INCLUDE_DIRS}'")

      foreach(INCLUDE_FILE_SET ${PARSE_REQUIRED_HEADERS})

        message("-- Searching for a header file in the set \"${INCLUDE_FILE_SET}\":")

        set(INCLUDE_FILE_LIST ${INCLUDE_FILE_SET})
        separate_arguments(INCLUDE_FILE_LIST)
        set(INCLUDE_FILE_SET_PATH) # Start out as empty list

        foreach(INCLUDE_FILE ${INCLUDE_FILE_LIST})

          message("--   Searching for header '${INCLUDE_FILE}' ...")

          if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
            print_var(INCLUDE_FILE)
          endif()

          set_notfound(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
          find_path( _${TPL_NAME}_${INCLUDE_FILE}_PATH
            NAMES ${INCLUDE_FILE}
            PATHS ${${TPL_NAME}_INCLUDE_DIRS}
            NO_DEFAULT_PATH)
          find_path( _${TPL_NAME}_${INCLUDE_FILE}_PATH
            NAMES ${INCLUDE_FILE} )
          mark_as_advanced(_${TPL_NAME}_${INCLUDE_FILE}_PATH)

          if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
            print_var(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
          endif()

          if(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
            message( "--     Found header '${_${TPL_NAME}_${INCLUDE_FILE}_PATH}/${INCLUDE_FILE}'")
            append_set(INCLUDE_FILE_SET_PATH ${_${TPL_NAME}_${INCLUDE_FILE}_PATH})
            break()
          endif()

        endforeach()

        if(NOT INCLUDE_FILE_SET_PATH)
          message("-- ERROR: Could not find a header file in"
            " the set \"${INCLUDE_FILE_SET}\"")
          if(PARSE_MUST_FIND_ALL_HEADERS)
            set(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
          endif()
        endif()

        append_set(INCLUDES_FOUND ${INCLUDE_FILE_SET_PATH})

      endforeach(INCLUDE_FILE_SET ${PARSE_REQUIRED_HEADERS})

      if (INCLUDES_FOUND)
        list(REMOVE_DUPLICATES INCLUDES_FOUND)
      endif()

      multiline_set(DOCSTR
        "List of semi-colon separated paths to append to the compile invocations"
        " to find the headers for the TPL ${TPL_NAME}.  This is the final variable"
        " that is used in the build commands.  The user variable ${TPL_NAME}_INCLUDE_DIRS"
        " is used to look for the given headers first but is just a suggestion."
        " This variable, however, is the final value and will not be touched."
        )

      advanced_set(TPL_${TPL_NAME}_INCLUDE_DIRS ${INCLUDES_FOUND}
        CACHE PATH ${DOCSTR} FORCE)
      # Above, we have to force the set in case the find failed the last
      # configure in which case this cache var will be empty.  NOTE: If the
      # user specified a non-empty TPL_${TPL_NAME}_INCLUDE_DIRS, then we would
      # never get here in the first place!

      if (NOT TPL_${TPL_NAME}_INCLUDE_DIRS OR NOT _${TPL_NAME}_ENABLE_SUCCESS)
        message(
          "-- ERROR: Could not find the include directories for TPL '${TPL_NAME}'!")
        message(
          "-- TIP: If the TPL '${TPL_NAME}' is on your system then you can set:\n"
          "     -D${TPL_NAME}_INCLUDE_DIRS='<dir0>;<dir1>;...'\n"
          "   to point to directories where these header files may be found.\n"
          "   Or, just set:\n"
          "     -DTPL_${TPL_NAME}_INCLUDE_DIRS='<dir0>;<dir1>;...'\n"
          "   to point to the include directories which will bypass any search for\n"
          "   header files and these include directories will be used without\n"
          "   question in the build.  (But this will result in a build-time error\n"
          "   obviously if the necessary header files are not found in these\n"
          "   include directories.)")
        tribits_tpl_find_include_dirs_and_libraries_handle_fail()
      endif()

      if (TPL_${TPL_NAME}_INCLUDE_DIRS)
        message("-- Found TPL '${TPL_NAME}' include dirs '${TPL_${TPL_NAME}_INCLUDE_DIRS}'")
      endif()
    else()

      # TPL_${TPL_NAME}_INCLUDE_DIRS is already in the cache so leave it alone!

    endif()

    # Print the final value to be used *always*
    message("-- TPL_${TPL_NAME}_INCLUDE_DIRS='${TPL_${TPL_NAME}_INCLUDE_DIRS}'")

  else()

    if (${TPL_NAME}_INCLUDE_DIRS)
      advanced_set(TPL_${TPL_NAME}_INCLUDE_DIRS ${${TPL_NAME}_INCLUDE_DIRS}
        CACHE PATH "User provided include dirs in the absence of include files.")
    else()
      if ("${TPL_${TPL_NAME}_INCLUDE_DIRS}" STREQUAL "")
        # Library has no header files, no user override, so just set them to
        # null (unless the user has already set this).
        global_null_set(TPL_${TPL_NAME}_INCLUDE_DIRS)
      endif()
    endif()

  endif()

  # Set library directories to null always.  We do this because
  # the package support code expects this variable and it is used
  # for package dependencies.  Therefore, we need it to allow
  # TPLs and internal packages to be treated in the same way.
  global_null_set(TPL_${TPL_NAME}_LIBRARY_DIRS)

  if (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
    print_var(TPL_${TPL_NAME}_LIBRARY_DIRS)
  endif()
  # 2011/05/09: rabartl: ToDo: Remove this above variable from everywhere!

  #print_var(TPL_TENTATIVE_ENABLE_${TPL_NAME})
  #print_var(_${TPL_NAME}_ENABLE_SUCCESS)
  if (TPL_TENTATIVE_ENABLE_${TPL_NAME})
    if (_${TPL_NAME}_ENABLE_SUCCESS)
      if (NOT PARSE_NO_PRINT_ENABLE_SUCCESS_FAIL)
        message("-- Attempt to tentatively enable TPL '${TPL_NAME}' passed!")
      endif()
    else()
      if (NOT PARSE_NO_PRINT_ENABLE_SUCCESS_FAIL)
        message("-- Attempt to tentatively enable TPL '${TPL_NAME}' failed!"
          "  Setting TPL_ENABLE_${TPL_NAME}=OFF")
      endif()
      set(TPL_ENABLE_${TPL_NAME} OFF CACHE STRING
        "Forced off since tentative enable failed!"  FORCE)
    endif()
  endif()

  if (_${TPL_NAME}_ENABLE_SUCCESS)
    global_set(TPL_${TPL_NAME}_NOT_FOUND FALSE)
  endif()

  set(buildDirExternalPkgsDir
    "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}")
  set(tplConfigFileBaseDir "${buildDirExternalPkgsDir}/${TPL_NAME}")
  set(tplConfigFile "${tplConfigFileBaseDir}/${TPL_NAME}Config.cmake")
  tribits_extpkg_write_config_file(${TPL_NAME} "${tplConfigFile}")
  if (NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    include("${tplConfigFile}")
  endif()
  # NOTE: The file <tplName>ConfigVersion.cmake will get created elsewhere as
  # will the install targets for the files <tplName>Config and
  # <tplName>ConfigVersion.cmake.

endfunction()


# @FUNCTION: tribits_tpl_tentatively_enable()
#
# Function that sets up for an optionally enabled TPL that is attempted to be
# enabled but will be disabled if all of the parts are not found.
#
# Usage::
#
#   tribits_tpl_tentatively_enable(<tplName>)
# 
# This function can be called from any CMakeLists.txt file to put a TPL in
# tentative enable mode.  But typically, it is called from an Package's
# `<packageDir>/cmake/Dependencies.cmake`_ file (see `How to tentatively
# enable an external package/TPL`_).
#
# This should only be used for optional TPLs.  It will not work correctly for
# required TPLs because any enabled packages that require this TPL will not be
# disabled and instead will fail to configure or fail to build.
#
# All this function does is to force set ``TPL_ENABLE_<tplName>=ON`` if it has
# not already been set, and sets ``TPL_TENTATIVE_ENABLE_<tplName>=ON`` in the
# cache.
#
# NOTE: This function will only tentatively enable a TPL if its enable has not
# be explicitly set on input, i.e. if ``-D TPL_ENABLE_<tplName>=""``.  If the
# TPL has been explicitly enabled (i.e. ``-D TPL_ENABLE_<tplName>=ON``) or
# disabled (i.e. ``-D TPL_ENABLE_<tplName>=OFF``), then this function has no
# effect and the TPL will be unconditionally enabled or disabled.
#
function(tribits_tpl_tentatively_enable  TPL_NAME)

  if ("${TPL_ENABLE_${TPL_NAME}}" STREQUAL "")
    # The TPL's enable status has not been set so tentatively enable it.
    set(TPL_ENABLE_${TPL_NAME} ON CACHE STRING
      "Set by tribits_tpl_tentatively_enable()" FORCE)
    advanced_set(TPL_TENTATIVE_ENABLE_${TPL_NAME} ON CACHE STRING
      "Set by tribits_tpl_tentatively_enable()" FORCE)
  else()
    # The TPL's enable status has already be hard set to be ON or OFF so we
    # will leave it alone.
  endif()

endfunction()


# Set find error and print error message
#
macro(tribits_tpl_find_include_dirs_and_libraries_handle_fail) 
  set(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
  global_set(TPL_${TPL_NAME}_NOT_FOUND TRUE)
  message(
    "-- ERROR: Failed finding all of the parts of TPL '${TPL_NAME}' (see above), Aborting!\n" )
  if ("${ERROR_MSG_MODE}" STREQUAL "SEND_ERROR")
    #message("ERROR_MSG_MODE=SEND_ERROR, Aborting")
    return()
  endif()
endmacro()
