set(REQUIRED_HEADERS umfpack.h amd.h )
set(REQUIRED_LIBS_NAMES umfpack amd )
set(IMPORTED_TARGETS_FOR_ALL_LIBS SuiteSparse::UMFPACK )

tribits_tpl_allow_pre_find_package(UMFPACK  UMFPACK_ALLOW_PREFIND)

if (UMFPACK_ALLOW_PREFIND)
  message("-- Using find_package(UMFPACK ...) ...")
  find_package(UMFPACK)
  if (UMFPACK_FOUND)
    message("-- Found UMFPACK_DIR='${UMFPACK_DIR}'")
    message("-- Generating UMFPACK::all_libs and UMFPACKConfig.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(UMFPACK
      INNER_FIND_PACKAGE_NAME  UMFPACK
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET UMFPACK::all_libs)

  set(UMFPACK_INCLUDE_DIRS_DEFAULT "/usr/include/suitesparse")
  set(UMFPACK_INCLUDE_DIRS "${UMFPACK_INCLUDE_DIRS_DEFAULT}" CACHE PATH
    "Default path to find UMFPACK include files")

  tribits_tpl_find_include_dirs_and_libraries( UMFPACK
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} )
endif()


# Amesos2 has Umfpack wrappers which depend on being able to
# send complex as Real-Im as a single array.
# Old versions of Umfpack don't support this so for now we are
# just disabling the complex tests for Umfpack version < 5.

include(CheckCSourceCompiles)
FUNCTION(CHECK_UMFPACK_HAS_VERSION_5  VARNAME)
  SET(SOURCE
  "
  #include <stdio.h>
  #include <umfpack.h>
  int main()
  {
    // this got added in Umfpack version 4.5
    #if UMFPACK_MAIN_VERSION >= 5
      return 0;
    #else
      umfpack_version_failure
    #endif
  }

  "
  )
  SET(CMAKE_REQUIRED_INCLUDES ${TPL_UMFPACK_INCLUDE_DIRS})
  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_UMFPACK_LIBRARIES})
  SET(CMAKE_REQUIRED_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

IF(TPL_ENABLE_UMFPACK)
  CHECK_UMFPACK_HAS_VERSION_5(HAVE_UMFPACK_VERSION_5)
ENDIF()
