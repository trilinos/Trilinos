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

IF (TribitsTplFindIncludeDirsAndLibraries_INCLUDED)
  RETURN()
ELSE()
  SET(TribitsTplFindIncludeDirsAndLibraries_INCLUDED TRUE)
ENDIF()

INCLUDE(AdvancedSet)
INCLUDE(AppendSet)
INCLUDE(AssertDefined)
INCLUDE(DualScopeSet)
INCLUDE(GlobalNullSet)
INCLUDE(GlobalSet)
INCLUDE(MultilineSet)
INCLUDE(ParseVariableArguments)
INCLUDE(SetNotFound)

#
# @FUNCTION: TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE()
#
# Function that determines if a TriBITS find module file
# ``FindTPL<tplName>.cmake`` is allowed to call ``FIND_PACKAGE(<tplName>
# ...)`` before calling `TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`_.
#
# Usage::
#
#   TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE( <tplName>
#     <allowPackagePrefindOut> )
#
# The required arguments are:
#
#   ``<tplName>`` : The input name of the TriBITS TPL (e.g. ``HDF5``).
#
#   ``<allowPackagePrefindOut>`` : Name of a variable which will be set to
#   ``TRUE`` on output if ``FIND_PACKAGE(<tplName> ...)`` should be called to
#   find the TPL ``<tplName>`` or ``FALSE`` if it should not be called.
#
# This function will set ``<allowPackagePrefindOut>`` to ``FALSE`` if any of
# the variables ``TPL_<tplName>_INCLUDE_DIRS``, ``${TPL_<tplName>_LIBRARIES``,
# or ``TPL_<tplName>_LIBRARY_DIRS`` are set.  This allows the user to override
# the search for the library components and just specify the absolute
# locations.  The function will also set ``<allowPackagePrefindOut>`` to
# ``FALSE`` if ``<tplName>_INCLUDE_DIRS``, ``<tplName>_LIBRARY_NAMES``, or
# ``<tplName>_LIBRARY_DIRS`` is set and ``<tplName>_FORCE_PRE_FIND_PACKAGE`` is
# set to ``FALSE``.  Otherwise, if ``<tplName>_FORCE_PRE_FIND_PACKAGE`` is set
# to ``TRUE``, the function will not return ``FALSE`` for
# ``<allowPackagePrefindOut>`` no matter what the values of
# ``<tplName>_INCLUDE_DIRS``, ``<tplName>_LIBRARY_NAMES``, or
# ``<tplName>_LIBRARY_DIRS``.
#
# The variable ``<tplName>_FORCE_PRE_FIND_PACKAGE`` is needed to allow users
# (or the ``FindTPL<tplName>.cmake`` module itself) to avoid name clashes with
# the variables ``<tplName>_INCLUDE_DIRS`` or ``<tplName>_LIBRARY_DIRS`` in
# the usage of ``FIND_PACKAGE(<tplName> ...)`` because a lot of default
# ``Find<tplName>.cmake`` modules also use these variables.  This function
# sets ``<tplName>_FORCE_PRE_FIND_PACKAGE`` as a cache variable with default
# value ``FALSE`` to maintain backward compatibility with existing
# ``FindTPL<tplName>.cmake`` modules.
#
# See `How to use FIND_PACKAGE() for a TriBITS TPL`_ for details in how to use
# this function to create a ``FindTPL<tplName>.cmake`` module file.
#
FUNCTION(TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE  TPL_NAME  ALLOW_PACAKGE_PREFIND_OUT)

  IF (TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE_DEBUG)
    MESSAGE("TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE: '${TPL_NAME}'  '${ALLOW_PACAKGE_PREFIND_OUT}'")
    PRINT_VAR(${TPL_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${TPL_NAME}_LIBRARY_NAMES)
    PRINT_VAR(${TPL_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${TPL_NAME}_FORCE_PRE_FIND_PACKAGE)
  ENDIF()

  ADVANCED_SET(${TPL_NAME}_FORCE_PRE_FIND_PACKAGE  FALSE
    CACHE BOOL
    "Determines if the variables ${TPL_NAME}_[INCLUDE_DIRS,LIBRARY_NAMES,LIBRARY_DIRS] should be ignored and the pre-find FIND_PACKAGE(${TPL_NAME} should be performed anyway.  But this will *not* do the pre-find if any of the TPL_${TPL_NAME}_[INCLUDE_DIRS,LIBRARY_NAMES,LIBRARY_DIRS] vars are set." )

  # Start out with TRUE and set to FALSE in logic below
  SET(ALLOW_PACAKGE_PREFIND TRUE)

  IF (
    (NOT "${TPL_${TPL_NAME}_INCLUDE_DIRS}" STREQUAL "")
    OR (NOT "${TPL_${TPL_NAME}_LIBRARIES}" STREQUAL "")
    OR (NOT "${TPL_${TPL_NAME}_LIBRARY_DIRS}" STREQUAL "")
    )
    # The user has selected one or more of the final vars so skip calling
    # FIND_PACKAGE(${TPL_NAME} ...) ...
    SET(ALLOW_PACAKGE_PREFIND FALSE)
  ELSEIF (
    (NOT "${${TPL_NAME}_INCLUDE_DIRS}" STREQUAL "")
    OR (NOT "${${TPL_NAME}_LIBRARY_NAMES}" STREQUAL "")
    OR (NOT "${${TPL_NAME}_LIBRARY_DIRS}" STREQUAL "")
    )
    # One ore more of the ${TPL_NAME}_XXX variables are set
    IF (${TPL_NAME}_FORCE_PRE_FIND_PACKAGE)
      # Even with one or more of the ${TPL_NAME}_XXX vars set, we still want
      # to do the FIND_PACKAGE(${TPL_NAME} ...) search and ignore this
      # override.
    ELSE()
      # We will not ignore the override of these variables and will instead go
      # ahead and skip the pre-find.
      SET(ALLOW_PACAKGE_PREFIND FALSE)
    ENDIF()
  ENDIF()

  SET(${ALLOW_PACAKGE_PREFIND_OUT} ${ALLOW_PACAKGE_PREFIND} PARENT_SCOPE)

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()
#
# Function that sets up cache variables for users to specify where to find a
# `TriBITS TPL`_'s headers and libraries.  This function is typically called
# inside of a ``FindTPL<tplName>.cmake`` moulde file (see
# `${TPL_NAME}_FINDMOD`_).
#
# Usage::
#
#   TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES(
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
#     include directories files using ``FIND_PATH()``. 
#
#   ``MUST_FIND_ALL_HEADERS``
#
#     If set, then all of the header files listed in ``REQUIRED_HEADERS`` must
#     be found in order for ``TPL_<tplName>_INCLUDE_DIRS`` to be defined.
#
#   ``REQUIRED_LIBS_NAMES``
#
#     List of libraries that are searched for when looking for the TPL's
#     libraries using ``FIND_LIBRARY()``.  This list can be overridden by the
#     user by setting ``<tplName>_LIBRARY_DIRS`` (see below).
#
#   ``MUST_FIND_ALL_LIBS``
#
#     If set, then all of the library files listed in ``REQUIRED_LIBS_NAMES``
#     must be found or the TPL is considered not found!
#
#   ``NO_PRINT_ENABLE_SUCCESS_FAIL``
#
#      If set, then the final success/fail will not be printed
#
# This function implements the TPL find behavior described in `Enabling
# support for an optional Third-Party Library (TPL)`_.
#
# The following (cache) variables, if set, will be used by that this function:
#
#   ``<tplName>_INCLUDE_DIRS`` (type ``PATH``)
#
#     List of paths to search first for header files defined in
#     ``REQUIRED_HEADERS``.
#
#   ``<tplName>_LIBRARY_DIRS`` (type ``PATH``)
#
#     The list of directories to search first for libraries defined in
#     ``REQUIRED_LIBS_NAMES``.  If, for some reason, no libraries should be
#     linked in for this particular configuration, then setting
#     ``<tplName>_LIBRARY_DIRS=OFF`` will 
#
#   ``<tplName>_LIBRARY_NAMES`` (type ``STRING``)
#
#     List of library names to be looked for instead of what is specified in
#     ``REQUIRED_LIBS_NAMES``.
#
# This function sets global variables to return state so it can be called from
# anywhere in the call stack.  The following cache variables are defined that
# are intended for the user to set and/or use:
#
#   ``TPL_<tplName>_INCLUDE_DIRS`` (type ``PATH``)
#
#     A list of common-separated full directory paths that contain the TPL's
#     header files.  If this variable is set before calling this function,
#     then no headers are searched for and this variable will be assumed to
#     have the correct list of header paths.
#
#   ``TPL_<tplName>_LIBRARIES`` (type ``FILEPATH``)
#
#     A list of commons-separated full library names (i.e. output from
#     ``FIND_LIBRARY()``) for all of the libraries found for the TPL.  If this
#     variable is set before calling this function, then no libraries are
#     searched for and this variable will be assumed to have the correct list
#     of libraries to link to.
#
#   ``TPL_<tplName>_NOT_FOUND`` (type ``BOOL``)
#
#     Will be set to ``ON`` if all of the parts of the TPL could not be found.
#
# ToDo: Document the behavior of this function for finding headers and
# libraries and when a find is successful and when it is not.
#
# Note, if ``TPL_TENTATIVE_ENABLE_<tplName>=ON``, then if all of the parts of
# the TPL can't be found, then ``TPL_ENABLE_<tplName>`` will be (forced) set
# to ``OFF`` in the cache.  See `TRIBITS_TPL_TENTATIVELY_ENABLE()`_.
#
FUNCTION(TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES TPL_NAME)

  # Make sure the right name is used
  ASSERT_DEFINED(TPL_ENABLE_${TPL_NAME})

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "REQUIRED_HEADERS;REQUIRED_LIBS_NAMES"
     #options
     "MUST_FIND_ALL_LIBS;MUST_FIND_ALL_HEADERS;NO_PRINT_ENABLE_SUCCESS_FAIL"
     ${ARGN}
     )

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    SET(TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE TRUE)
  ENDIF()

  IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
    MESSAGE("TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES: ${TPL_NAME}")
    PRINT_VAR(PARSE_REQUIRED_HEADERS)
    PRINT_VAR(PARSE_REQUIRED_LIBS_NAMES)
    PRINT_VAR(TPL_${TPL_NAME}_INCLUDE_DIRS)
    PRINT_VAR(TPL_${TPL_NAME}_LIBRARIES)
  ENDIF()

  IF (TPL_TENTATIVE_ENABLE_${TPL_NAME})
    MESSAGE("-- Attempting to tentatively enable TPL '${TPL_NAME}' ...")
    SET(ERROR_MSG_MODE)
  ELSE()
    SET(ERROR_MSG_MODE SEND_ERROR)
  ENDIF()

  #
  # User options
  #

  IF (PARSE_REQUIRED_LIBS_NAMES)

    # Library directories

    MULTILINE_SET(DOCSTR
      "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
      " libraries.  This list of paths will be passed into a FIND_LIBRARY(...)"
      " command to find the libraries listed in ${TPL_NAME}_LIBRARY_NAMES."
      "  Note that this set of paths is also the default value used for"
      " ${TPL_NAME}_LIBRARY_DIRS.  Therefore, if the headers exist in the"
      " same directories as the library, you do not need to set"
      " ${TPL_NAME}_LIBRARY_DIRS."
      )

    ADVANCED_SET(${TPL_NAME}_LIBRARY_DIRS "" CACHE PATH ${DOCSTR})
    IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
      PRINT_VAR(${TPL_NAME}_LIBRARY_DIRS)
    ENDIF()

    # Libraries

    MULTILINE_SET(DOCSTR
      "List of semi-colon separated names of libraries needed to link to for"
      " the TPL ${TPL_NAME}.  This list of libraries will be search for in"
      " FIND_LIBRARY(...) calls along with the directories specified with"
      " ${TPL_NAME}_LIBRARY_DIRS.  NOTE: This is not the final list of libraries"
      " used for linking.  That is specified by TPL_${TPL_NAME}_LIBRARIES!"
      )
    ADVANCED_SET(${TPL_NAME}_LIBRARY_NAMES ${PARSE_REQUIRED_LIBS_NAMES}
      CACHE STRING ${DOCSTR})

    # Let the user override what the names of the libraries which might
    # actually mean that no libraries are searched for.
    SET(PARSE_REQUIRED_LIBS_NAMES ${${TPL_NAME}_LIBRARY_NAMES})

    IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
      PRINT_VAR(${TPL_NAME}_LIBRARY_NAMES)
      PRINT_VAR(PARSE_REQUIRED_LIBS_NAMES)
    ENDIF()

  ELSE()

    SET(${TPL_NAME}_LIBRARY_DIRS) # Just to ignore below!

  ENDIF()

  # Include directories

  IF (PARSE_REQUIRED_HEADERS)

    MULTILINE_SET(DOCSTR
      "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
      " headers.  This list of paths will be passed into a FIND_PATH(...)"
      " command to find the headers for ${TPL_NAME} (which are known in advance)."
      )
    ADVANCED_SET(${TPL_NAME}_INCLUDE_DIRS ${${TPL_NAME}_LIBRARY_DIRS}
      CACHE PATH ${DOCSTR})

    IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
      PRINT_VAR(${TPL_NAME}_LIBRARY_DIRS)
      PRINT_VAR(${TPL_NAME}_INCLUDE_DIRS)
      PRINT_VAR(PARSE_REQUIRED_HEADERS)
    ENDIF()

  ENDIF()

  #
  # Set the lib extensions to find
  #

  # Save the default the first time through
  IF (NOT CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT)
   SET(TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT ${CMAKE_FIND_LIBRARY_SUFFIXES})
   #PRINT_VAR(TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT)
  ENDIF()

  #PRINT_VAR(TPL_FIND_SHARED_LIBS)
  #PRINT_VAR(CMAKE_FIND_LIBRARY_SUFFIXES)
  # Set libraries to find
  IF (TPL_FIND_SHARED_LIBS)
    # The default should be to find shared libs first
    SET(TPL_CMAKE_FIND_LIBRARY_SUFFIXES ${TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT})
  ELSE()
    IF (WIN32)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a)
    ELSE()
      SET(CMAKE_FIND_LIBRARY_SUFFIXES .a )
    ENDIF()
  ENDIF()
  #PRINT_VAR(CMAKE_FIND_LIBRARY_SUFFIXES)

  #
  # Direct build options
  #

  SET(_${TPL_NAME}_ENABLE_SUCCESS TRUE)

  IF (PARSE_REQUIRED_LIBS_NAMES)

    # Libraries

    IF (NOT TPL_${TPL_NAME}_LIBRARIES)

      IF (PARSE_MUST_FIND_ALL_LIBS)
        MESSAGE("-- Must find at least one lib in each of the"
          " lib sets \"${PARSE_REQUIRED_LIBS_NAMES}\"")
      ENDIF()

      MESSAGE( "-- Searching for libs in ${TPL_NAME}_LIBRARY_DIRS='${${TPL_NAME}_LIBRARY_DIRS}'")

      SET(LIBRARIES_FOUND)

      FOREACH(LIBNAME_SET ${${TPL_NAME}_LIBRARY_NAMES})

        MESSAGE("-- Searching for a lib in the set \"${LIBNAME_SET}\":")

        IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
          PRINT_VAR(LIBNAME_SET)
        ENDIF()

        SET(LIBNAME_LIST ${LIBNAME_SET})
        SEPARATE_ARGUMENTS(LIBNAME_LIST)

        SET(LIBNAME_SET_LIB)

        FOREACH(LIBNAME ${LIBNAME_LIST})

          MESSAGE("--   Searching for lib '${LIBNAME}' ...")

          IF (${TPL_NAME}_LIBRARY_DIRS)
            SET(PATHS_ARG PATHS ${${TPL_NAME}_LIBRARY_DIRS})
          ELSE()
            SET(PATHS_ARG PATHS)
          ENDIF()

          SET_NOTFOUND(_${TPL_NAME}_${LIBNAME}_LIBRARY)
          FIND_LIBRARY( _${TPL_NAME}_${LIBNAME}_LIBRARY
            NAMES ${LIBNAME}
            ${PATHS_ARG} NO_DEFAULT_PATH )
          FIND_LIBRARY( _${TPL_NAME}_${LIBNAME}_LIBRARY
            NAMES ${LIBNAME} )
          MARK_AS_ADVANCED(_${TPL_NAME}_${LIBNAME}_LIBRARY)

          IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
            PRINT_VAR(_${TPL_NAME}_${LIBNAME}_LIBRARY)
          ENDIF()

          IF (_${TPL_NAME}_${LIBNAME}_LIBRARY)
            MESSAGE("--     Found lib '${_${TPL_NAME}_${LIBNAME}_LIBRARY}'")
            SET(LIBNAME_SET_LIB ${_${TPL_NAME}_${LIBNAME}_LIBRARY})
            BREAK()
          ENDIF()

        ENDFOREACH()

        IF (NOT LIBNAME_SET_LIB)
          MESSAGE(
            "-- ERROR: Did not find a lib in the lib set \"${LIBNAME_SET}\""
             " for the TPL '${TPL_NAME}'!")
          IF (PARSE_MUST_FIND_ALL_LIBS)
	    SET(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
          ELSE()
            BREAK()
          ENDIF()
        ENDIF()

        APPEND_SET(LIBRARIES_FOUND ${LIBNAME_SET_LIB})

      ENDFOREACH()

      MULTILINE_SET(DOCSTR
        "List of semi-colon separated full paths to the libraries for the TPL"
        " ${TPL_NAME}.  This is the final variable that is used in the link"
        " commands.  The user variable ${TPL_NAME}_LIBRARY_DIRS is used to look"
        " for the know library names but but is just a suggestion."
        " This variable, however, is the final value and will not be touched."
        )
      ADVANCED_SET( TPL_${TPL_NAME}_LIBRARIES ${LIBRARIES_FOUND}
        CACHE FILEPATH ${DOCSTR} FORCE)
      # Above, we have to force the set in case the find failed the last
      # configure in which case this cache var will be empty.  NOTE: If the
      # user specified a non-empty TPL_${TPL_NAME}_LIBRARIES, then we would
      # never get here in the first place!

      IF (NOT TPL_${TPL_NAME}_LIBRARIES OR NOT _${TPL_NAME}_ENABLE_SUCCESS)
        MESSAGE(
          "-- ERROR: Could not find the libraries for the TPL '${TPL_NAME}'!")
        MESSAGE(
          "-- TIP: If the TPL '${TPL_NAME}' is on your system then you can set:\n"
          "     -D${TPL_NAME}_LIBRARY_DIRS='<dir0>;<dir1>;...'\n"
          "   to point to the directories where these libraries may be found.\n"
          "   Or, just set:\n"
          "     -DTPL_${TPL_NAME}_LIBRARIES='<path-to-libs0>;<path-to-libs1>;...'\n"
          "   to point to the full paths for the libraries which will\n"
          "   bypass any search for libraries and these libraries will be used without\n"
          "   question in the build.  (But this will result in a build-time error\n"
          "   if not all of the necessary symbols are found.)")
        TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_HANDLE_FAIL()
      ENDIF()

    ENDIF()

    # Print the final value to be used *always*
    MESSAGE("-- TPL_${TPL_NAME}_LIBRARIES='${TPL_${TPL_NAME}_LIBRARIES}'")

  ELSE()

    # There are no libraries so set the libraries to null but don't
    # change the cache which should not even have this variable in it.
    # This set command is only to follow the standards for the package
    # support CMake code.
    GLOBAL_NULL_SET(TPL_${TPL_NAME}_LIBRARIES)

  ENDIF()

  # Include directories

  IF (PARSE_REQUIRED_HEADERS)

    IF (NOT TPL_${TPL_NAME}_INCLUDE_DIRS)

      IF (PARSE_MUST_FIND_ALL_HEADERS)
        MESSAGE("-- Must find at least one header in each of the"
          " header sets \"${PARSE_REQUIRED_HEADERS}\"")
      ENDIF()

      MESSAGE( "-- Searching for headers in ${TPL_NAME}_INCLUDE_DIRS='${${TPL_NAME}_INCLUDE_DIRS}'")

      FOREACH(INCLUDE_FILE_SET ${PARSE_REQUIRED_HEADERS})

        MESSAGE("-- Searching for a header file in the set \"${INCLUDE_FILE_SET}\":")

        SET(INCLUDE_FILE_LIST ${INCLUDE_FILE_SET})
        SEPARATE_ARGUMENTS(INCLUDE_FILE_LIST)
        SET(INCLUDE_FILE_SET_PATH) # Start out as empty list

        FOREACH(INCLUDE_FILE ${INCLUDE_FILE_LIST})

          MESSAGE("--   Searching for header '${INCLUDE_FILE}' ...")

          IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
            PRINT_VAR(INCLUDE_FILE)
          ENDIF()

          SET_NOTFOUND(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
          FIND_PATH( _${TPL_NAME}_${INCLUDE_FILE}_PATH
            NAMES ${INCLUDE_FILE}
            PATHS ${${TPL_NAME}_INCLUDE_DIRS}
            NO_DEFAULT_PATH)
          FIND_PATH( _${TPL_NAME}_${INCLUDE_FILE}_PATH
            NAMES ${INCLUDE_FILE} )
          MARK_AS_ADVANCED(_${TPL_NAME}_${INCLUDE_FILE}_PATH)

          IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
            PRINT_VAR(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
          ENDIF()

          IF(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
            MESSAGE( "--     Found header '${_${TPL_NAME}_${INCLUDE_FILE}_PATH}/${INCLUDE_FILE}'")
            APPEND_SET(INCLUDE_FILE_SET_PATH ${_${TPL_NAME}_${INCLUDE_FILE}_PATH})
            BREAK()
          ENDIF()

        ENDFOREACH()

        IF(NOT INCLUDE_FILE_SET_PATH)
          MESSAGE("-- ERROR: Could not find a header file in"
            " the set \"${INCLUDE_FILE_SET}\"")
          IF(PARSE_MUST_FIND_ALL_HEADERS)
            SET(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
          ENDIF()
        ENDIF()

        APPEND_SET(INCLUDES_FOUND ${INCLUDE_FILE_SET_PATH})

      ENDFOREACH(INCLUDE_FILE_SET ${PARSE_REQUIRED_HEADERS})

      IF (INCLUDES_FOUND)
        LIST(REMOVE_DUPLICATES INCLUDES_FOUND)
      ENDIF()

      MULTILINE_SET(DOCSTR
        "List of semi-colon separated paths to append to the compile invocations"
        " to find the headers for the TPL ${TPL_NAME}.  This is the final variable"
        " that is used in the build commands.  The user variable ${TPL_NAME}_INCLUDE_DIRS"
        " is used to look for the given headers first but is just a suggestion."
        " This variable, however, is the final value and will not be touched."
        )

      ADVANCED_SET(TPL_${TPL_NAME}_INCLUDE_DIRS ${INCLUDES_FOUND}
        CACHE PATH ${DOCSTR} FORCE)
      # Above, we have to force the set in case the find failed the last
      # configure in which case this cache var will be empty.  NOTE: If the
      # user specified a non-empty TPL_${TPL_NAME}_INCLUDE_DIRS, then we would
      # never get here in the first place!

      IF (NOT TPL_${TPL_NAME}_INCLUDE_DIRS OR NOT _${TPL_NAME}_ENABLE_SUCCESS)
        MESSAGE(
          "-- ERROR: Could not find the include directories for TPL '${TPL_NAME}'!")
        MESSAGE(
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
        TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_HANDLE_FAIL()
      ENDIF()

      IF (TPL_${TPL_NAME}_INCLUDE_DIRS)
        MESSAGE("-- Found TPL '${TPL_NAME}' include dirs '${TPL_${TPL_NAME}_INCLUDE_DIRS}'")
      ENDIF()
    ELSE()

      # TPL_${TPL_NAME}_INCLUDE_DIRS is already in the cache so leave it alone!

    ENDIF()

    # Print the final value to be used *always*
    MESSAGE("-- TPL_${TPL_NAME}_INCLUDE_DIRS='${TPL_${TPL_NAME}_INCLUDE_DIRS}'")

  ELSE()

    IF (${TPL_NAME}_INCLUDE_DIRS)
      ADVANCED_SET(TPL_${TPL_NAME}_INCLUDE_DIRS ${${TPL_NAME}_INCLUDE_DIRS}
        CACHE PATH "User provided include dirs in the absence of include files.")
    ELSE()
      # Library has no header files, no user override, so just set them to null
      GLOBAL_NULL_SET(TPL_${TPL_NAME}_INCLUDE_DIRS)
    ENDIF()

  ENDIF()

  # Set library directories to null always.  We do this because
  # the package support code expects this variable and it is used
  # for package dependencies.  Therefore, we need it to allow
  # TPLs and internal packages to be treated in the same way.
  GLOBAL_NULL_SET(TPL_${TPL_NAME}_LIBRARY_DIRS)

  IF (TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE)
    PRINT_VAR(TPL_${TPL_NAME}_LIBRARY_DIRS)
  ENDIF()
  # 2011/05/09: rabartl: ToDo: Remove this above variable from everywhere!

  #PRINT_VAR(TPL_TENTATIVE_ENABLE_${TPL_NAME})
  #PRINT_VAR(_${TPL_NAME}_ENABLE_SUCCESS)
  IF (TPL_TENTATIVE_ENABLE_${TPL_NAME})
    IF (_${TPL_NAME}_ENABLE_SUCCESS)
      IF (NOT PARSE_NO_PRINT_ENABLE_SUCCESS_FAIL)
        MESSAGE("-- Attempt to tentatively enable TPL '${TPL_NAME}' passed!")
      ENDIF()
    ELSE()
      IF (NOT PARSE_NO_PRINT_ENABLE_SUCCESS_FAIL)
        MESSAGE("-- Attempt to tentatively enable TPL '${TPL_NAME}' failed!"
          "  Setting TPL_ENABLE_${TPL_NAME}=OFF")
      ENDIF()
      SET(TPL_ENABLE_${TPL_NAME} OFF CACHE STRING
        "Forced off since tentative enable failed!"  FORCE)
    ENDIF()
  ENDIF()

  IF (_${TPL_NAME}_ENABLE_SUCCESS)
    GLOBAL_SET(TPL_${TPL_NAME}_NOT_FOUND FALSE)
  ENDIF()

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_TPL_TENTATIVELY_ENABLE()
#
# Function that sets up for an optionally enabled TPL that is attempted to be
# enabled but will be disabled if all of the parts are not found.
#
# Usage::
#
#   TRIBITS_TPL_TENTATIVELY_ENABLE(<tplName>)
# 
# This function can be called from any CMakeLists.txt file to put a TPL in
# tentative enable mode.
#
# This should only be used for optional TPLs.  It will not work correctly for
# required TPLs because any enabled packages that require this TPL will not be
# disabled and instead will fail to configure or fail to build.
#
# All this function does is to force set ``TPL_ENABLE_<tplName>=ON`` if it has
# not already been set. and ``TPL_TENTATIVE_ENABLE_<tplName>`` in the cache.
#
# NOTE: This function will only tentatively enable a TPL it its enable has not
# be explicitly set on input, i.e. if ``TPL_ENABLE_<tplName>=""``.
#
FUNCTION(TRIBITS_TPL_TENTATIVELY_ENABLE  TPL_NAME)

  IF ("${TPL_ENABLE_${TPL_NAME}}" STREQUAL "")
    # The TPL's enable status has not been set so tentatively enable it.
    SET(TPL_ENABLE_${TPL_NAME} ON CACHE STRING "autoset" FORCE)
    ADVANCED_SET(TPL_TENTATIVE_ENABLE_${TPL_NAME} ON CACHE STRING "autoset" FORCE)
  ELSE()
    # The TPL's enable status has already be hard set to be ON or OFF so we
    # will leave it alone.
  ENDIF()

ENDFUNCTION()


#
# Utility macro
#

MACRO(TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_HANDLE_FAIL) 
  SET(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
  GLOBAL_SET(TPL_${TPL_NAME}_NOT_FOUND TRUE)
  MESSAGE(
    "-- ERROR: Failed finding all of the parts of TPL '${TPL_NAME}' (see above), Aborting!\n" )
  IF ("${ERROR_MSG_MODE}" STREQUAL "SEND_ERROR")
    #MESSAGE("ERROR_MSG_MODE=SEND_ERROR, Aborting")
    RETURN()
  ENDIF()
ENDMACRO()
