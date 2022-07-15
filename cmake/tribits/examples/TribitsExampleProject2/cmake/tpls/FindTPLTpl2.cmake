set(REQUIRED_HEADERS  Tpl2a.hpp) # Only look for one header file to find include dir
set(REQUIRED_LIBS_NAMES  tpl2b tpl2a)
set(IMPORTED_TARGETS_FOR_ALL_LIBS  tpl2::tpl2a tpl2::tpl2b)

tribits_tpl_allow_pre_find_package(Tpl2  Tpl2_ALLOW_PREFIND)

if (Tpl2_ALLOW_PREFIND)
  message("-- Using find_package(Tpl2 ...) ...")
  find_package(Tpl2)
  if (Tpl2_FOUND)
    message("-- Found Tpl2_DIR='${Tpl2_DIR}'")
    message("-- Generating Tpl2::all_libs and Tpl2Config.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(Tpl2
      INNER_FIND_PACKAGE_NAME  Tpl2
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET Tpl2::all_libs)
  tribits_tpl_find_include_dirs_and_libraries( Tpl2
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} )
endif()
