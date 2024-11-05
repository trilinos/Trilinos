set(REQUIRED_HEADERS  Tpl4.hpp)
set(IMPORTED_TARGETS_FOR_ALL_LIBS  tpl4::tpl4)

tribits_tpl_allow_pre_find_package(Tpl4  Tpl4_ALLOW_PREFIND)

if (Tpl4_ALLOW_PREFIND)
  message("-- Using find_package(Tpl4 ...) ...")
  find_package(Tpl4)
  if (Tpl4_FOUND)
    message("-- Found Tpl4_DIR='${Tpl4_DIR}'")
    message("-- Generating Tpl4::all_libs and Tpl4Config.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(Tpl4
      INNER_FIND_PACKAGE_NAME  Tpl4
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET Tpl4::all_libs)
  tribits_tpl_find_include_dirs_and_libraries( Tpl4
    REQUIRED_HEADERS ${REQUIRED_HEADERS} )
endif()
