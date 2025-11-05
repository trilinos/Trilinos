# Tertiary stable code since amd is also built as part of UMFPACK.

set(REQUIRED_HEADERS amd.h )
set(REQUIRED_LIBS_NAMES amd )
set(IMPORTED_TARGETS_FOR_ALL_LIBS SuiteSparse::AMD )

tribits_tpl_allow_pre_find_package(AMD  AMD_ALLOW_PREFIND)

if (AMD_ALLOW_PREFIND)
  message("-- Using find_package(AMD ...) ...")
  find_package(AMD)
  if (AMD_FOUND)
    message("-- Found AMD_DIR='${AMD_DIR}'")
    message("-- Generating AMD::all_libs and AMDConfig.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(AMD
      INNER_FIND_PACKAGE_NAME  AMD
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET AMD::all_libs)

  set(AMD_INCLUDE_DIRS_DEFAULT "/usr/include/suitesparse")
  set(AMD_INCLUDE_DIRS "${AMD_INCLUDE_DIRS_DEFAULT}" CACHE PATH
    "Default path to find AMD include files")

  tribits_tpl_find_include_dirs_and_libraries( AMD
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} )
endif()
