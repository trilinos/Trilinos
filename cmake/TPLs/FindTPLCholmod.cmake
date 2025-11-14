set(REQUIRED_HEADERS cholmod.h cholmod_core.h )
set(REQUIRED_LIBS_NAMES cholmod amd colamd ccolamd camd suitesparseconfig )
set(IMPORTED_TARGETS_FOR_ALL_LIBS SuiteSparse::CHOLMOD )

tribits_tpl_allow_pre_find_package(Cholmod  Cholmod_ALLOW_PREFIND)

if (Cholmod_ALLOW_PREFIND)
  message("-- Using find_package(CHOLMOD ...) ...")
  find_package(CHOLMOD)
  if (CHOLMOD_FOUND)
    message("-- Found Cholmod_DIR='${CHOLMOD_DIR}'")
    message("-- Generating Cholmod::all_libs and CholmodConfig.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(Cholmod
      INNER_FIND_PACKAGE_NAME  CHOLMOD
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET Cholmod::all_libs)
  set(Cholmod_INCLUDE_DIRS_DEFAULT "/usr/include/suitesparse")
  set(Cholmod_INCLUDE_DIRS "${Cholmod_INCLUDE_DIRS_DEFAULT}" CACHE PATH
    "Default path to find Cholmod include files")

  tribits_tpl_find_include_dirs_and_libraries( Cholmod
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} )
endif()
