#############################################################################
#
# Function to generate an export makefile, simplifying the building of
# clients of Trilinos. This function will also build an example of a 
# client's makefile that includes the export makefile.. 
# 
# This function can be called any time after the processing of 
# dependency-based enabling of packages and TPLs. 
# If called before that point, it won't work correctly as it will 
# only process packages and TPLs that have been explicitly enabled. 
#
# The function works as follows: 
# (1) initialize empty lists of libraries, TPL libraries, and TPL 
#     include paths
# (2) loop over all packages in order from most to least dependent.
#     (a) Skip packages that aren't enabled. 
#     (b) Append this package's libraries to the library list. Each package
#         lists its libraries from most to least dependent, i.e., the order
#         in which they'd appear on a linker command line. 
#     (c) Append the TPLs used by this library to the TPL library list. Skip
#         any optional TPLs that have not been enabled.
#     (d) Append the include paths for TPLs used by this library to 
#         the list of include paths.
# (3) Prepend the command-line prefix (e.g., "-l") to each library name.  
# (4) Form a command-line string from the library names, e.g., 
#         "-lfoo -lbar -lsnafu"
#     If necessary, write an RPATH specification.
# (5) Reverse the TPL lib list order, remove duplicates, and reverse again. 
#     This step is needed because several packages may list the same TPLs. 
#     The removals must proceed in reverse order on order to preserve the
#     correct dependency order.
# (6) Form command-line strings for the TPL libs and includes.
# (7) Use PACKAGE_CONFIGURE_FILE to splice everything into the export makefile
#     and client makefile.
#     
# Possible issues for the future:
# (1) I've had to assume "-I" for the include path prefix. This might not
#     be portable to all platforms. 
# (2) Package libraries are listed without suffices as "foo" and "bar", 
#     and so need to be prefixed, e.g., "-lfoo -lbar". TPL libraries are
#     given with suffices and paths, and so can *not* be prefixed with 
#     "-l". This as not a problem, provided this pattern is consistent for
#     all TPLs. 
# (3) Multiple TPLs for a package are assumed to be listed in reverse
#     dependency order (e.g, blas before lapack). This seems to be the
#     case, but I don't know whether it's been formally specified in the
#     TPL setup process. 
#
# Kevin Long, Texas Tech University
# kevin.long@ttu.edu
#
#############################################################################


FUNCTION(PACKAGE_WRITE_EXPORT_MAKEFILE Proj PACKAGE EMFName CMFName)
# Input arguments:
#   (*) Proj    -- name of project, e.g., "Trilinos"
#   (*) PACKAGE -- name of the package calling this function, e.g., "Epetra"
#   (*) EMFName -- name of export makefile, e.g., "Makefile.export"
#   (*) CMFName -- name of client makefile, e.g., "Makefile.client" 


# Get the list of all packages, assumed to be in order of dependency
# (most-dependent packages last). The Trilinos cmake system specifies 
# that the master list of packages be in this order, so this assumption
# should be valid unless there's an error in the Trilinos package file.
  SET(PACK_LIST ${${Proj}_PACKAGES})

# Remove from the package list all packages after ${PACKAGE} in the dependency
# chain. This way each package can create its own minimalist export makefile
# with no upstream libraries. 
  SET(TMP_PACK_LIST)
  SET(SKIP FALSE)
  FOREACH(PACK_NAME ${PACK_LIST})
    IF(NOT SKIP)
      LIST(APPEND TMP_PACK_LIST ${PACK_NAME})
      IF(PACKAGE STREQUAL ${PACK_NAME})
        SET(SKIP TRUE)
      ENDIF()
    ENDIF()
  ENDFOREACH()

  SET(PACK_LIST ${TMP_PACK_LIST})


# Reverse the order of the package list, letting us loop 
# from most-dependent to least-dependent. 
  LIST(REVERSE PACK_LIST)

# Create empty variables for the lists of libraries, TPL libs, and include
# paths. 
  SET(LIB_LIST)
  SET(TPL_LIB_LIST)
  SET(TPL_INCLUDE_LIST)

# ------------------------------------
# --------- Main loop over packages. 
# ------------------------------------
  FOREACH(PACK_NAME ${PACK_LIST})

    # Skip packages that haven't been enabled. 
    IF(Trilinos_ENABLE_${PACK_NAME})
  
      # ------- library handling -------------

      # Append this package's libraries to the list.
      LIST(APPEND LIB_LIST ${${PACK_NAME}_LIBRARIES})

      # ------ required TPL handling ---------

      # Get the required TPLs for this package.
      SET(REQ_TPLS ${${PACK_NAME}_LIB_REQUIRED_DEP_TPLS})
      # The TPL list is in reverse order for some reason. 
      IF(REQ_TPLS)
        LIST(REVERSE REQ_TPLS)
      ENDIF()
      # Append each required TPL lib to the main TPL lib list. Likewise the
      # include paths. 
      FOREACH(TPL_NAME ${REQ_TPLS})
        LIST(APPEND TPL_LIB_LIST ${TPL_${TPL_NAME}_LIBRARIES})
        LIST(APPEND TPL_INCLUDE_LIST ${TPL_${TPL_NAME}_INCLUDE_DIRS})
      ENDFOREACH()

      # ------ optional TPL handling ---------

      # Get the optional TPLs for this package.
      SET(OPT_TPLS ${${PACK_NAME}_LIB_OPTIONAL_DEP_TPLS})
      # The TPL list is in reverse order for some reason. 
      IF(OPT_TPLS)
        LIST(REVERSE OPT_TPLS)
      ENDIF()

      # Append each enabled optional TPL lib to the main TPL lib list. 
      # Likewise the include paths. 
      FOREACH(TPL_NAME ${OPT_TPLS})
        # Skip TPLs that haven't been anabled.
        IF(TPL_ENABLE_${TPL_NAME}) 
          LIST(APPEND TPL_LIB_LIST ${TPL_${TPL_NAME}_LIBRARIES})
          LIST(APPEND TPL_INCLUDE_LIST ${TPL_${TPL_NAME}_INCLUDE_DIRS})
        ENDIF()
      ENDFOREACH()
    ENDIF()
  ENDFOREACH()

# ------------------------------------
# --------- End main loop over packages. At this point, we need
# --------- to set up strings to be written to the command lines
# --------- in the export makefile
# ------------------------------------

  # Prepend the library command-line flag (e.g., "-l" on unix systems) to
  # each library.
  SET(LIB_LIST_COPY)
  FOREACH(LIB ${LIB_LIST})
    LIST(APPEND LIB_LIST_COPY ${CMAKE_LINK_LIBRARY_FLAG}${LIB})
  ENDFOREACH()
  SET(LIB_LIST ${LIB_LIST_COPY})

  # Gather libs into a single string as would appear on the 
  # linker command line. 
  SET(LIB_STR "")
  FOREACH(LIB ${LIB_LIST})
    SET(LIB_STR "${LIB_STR} ${LIB}")
  ENDFOREACH()

  # Write the specification of the rpath if necessary. This is only needed
  # if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/lib)
  ENDIF()

  # Remove duplicates from the TPL list. To preserve dependency order the
  # last entry in the list must be the one that remains. For example, if the
  # list contains [lapack,fred,blas,joe,bob,joe,lapack,blas] the stripped
  # list will be [fred,bob,joe,lapack,blas].
  IF (TPL_LIB_LIST)
    LIST(REVERSE TPL_LIB_LIST)
    LIST(REMOVE_DUPLICATES TPL_LIB_LIST)
    LIST(REVERSE TPL_LIB_LIST)
  ENDIF()


  # Gather TPL libs into a single string as would appear on the 
  # linker command line. 
  SET(TPL_LIB_STR "")
  FOREACH(LIB ${TPL_LIB_LIST})
    SET(TPL_LIB_STR "${TPL_LIB_STR} ${LIB}")
  ENDFOREACH()

 

  # Gather TPL includes into a single string as would appear on the 
  # linker command line. 
  SET(TPL_INCLUDE_STR "")
  FOREACH(INC_DIR ${TPL_INCLUDE_LIST})
    ################################################################
    # WARNING: I've hardwired the "-I" flag here. I don't know 
    # whether that's portable, but I couldn't find a CMake variable
    # for the include prefix so I had to hardwire "-I." 
    # Perhaps most compilers use "-I"? 
    # - KL 16 Mar 2009. 
    ################################################################
    SET(TPL_INCLUDE_STR "${TPL_INCLUDE_STR} -I${INC_DIR}")
  ENDFOREACH()

  # Print diagnostic info if verbose
  IF (${Proj}_VERBOSE_CONFIGURE)
    PRINT_VAR(LIB_STR)
    PRINT_VAR(TPL_LIB_STR)
    PRINT_VAR(TPL_INCLUDE_STR)
  ENDIF()

  # Get compiler flags
  PACKAGE_ARCH_DEFINE_STANDARD_COMPILE_FLAGS_VARS()

  # Find the directory name associated with a package. Try the package name
  # and an all-lower-case version of the package name. 
  STRING(TOLOWER ${PACKAGE} LOWER_CASE_PACKAGE)
  IF(EXISTS ${PROJECT_SOURCE_DIR}/packages/${PACKAGE})
    SET(PACKAGE_DIR_NAME ${PACKAGE})
  ELSEIF(EXISTS ${PROJECT_SOURCE_DIR}/packages/${LOWER_CASE_PACKAGE} )
    SET(PACKAGE_DIR_NAME ${LOWER_CASE_PACKAGE})
  ELSE()
    MESSAGE(SEND_ERROR "PackageWriteExportMakefile could not find ${PROJECT_SOURCE_DIR}/packages/${PACKAGE} or ${PROJECT_SOURCE_DIR}/packages/${LOWER_CASE_PACKAGE}")
  ENDIF()


  # Generate an all-caps version of the package name. This will be used
  # in package-specific variables inside the export makefile, in keeping
  # with usual makefile naming style.
  STRING(TOUPPER ${PACKAGE} CAPS_PACKAGE_NAME)

  # Generate a note discouraging editing of the export makefile
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")

  # Splice the library, include, and compiler information into the export
  # makefile. 
  ### WARNING: hardwired "/" separator here.
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/cmake/${EMFName}.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${EMFName}
    )

  # Splice the path and package name information into the client makefile
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/cmake/${CMFName}.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${CMFName}
    )



ENDFUNCTION()