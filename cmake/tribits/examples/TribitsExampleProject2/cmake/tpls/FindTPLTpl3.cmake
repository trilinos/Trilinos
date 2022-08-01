set(REQUIRED_HEADERS  Tpl3.hpp)
set(REQUIRED_LIBS_NAMES  tpl3)
set(IMPORTED_TARGETS_FOR_ALL_LIBS  tpl3::tpl3)

tribits_tpl_allow_pre_find_package(Tpl3  Tpl3_ALLOW_PREFIND)

if (Tpl3_ALLOW_PREFIND)
  message("-- Using find_package(Tpl3 ...) ...")
  find_package(Tpl3)
  if (Tpl3_FOUND)
    message("-- Found Tpl3_DIR='${Tpl3_DIR}'")
    message("-- Generating Tpl3::all_libs and Tpl3Config.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(Tpl3
      INNER_FIND_PACKAGE_NAME  Tpl3
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET Tpl3::all_libs)
  tribits_tpl_find_include_dirs_and_libraries( Tpl3
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} )
endif()
