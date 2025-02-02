set(REQUIRED_HEADERS boost/version.hpp boost/mpl/at.hpp )
set(REQUIRED_LIBS_NAMES boost_program_options boost_system )
set(IMPORTED_TARGETS_FOR_ALL_LIBS Boost )

tribits_tpl_allow_pre_find_package(Boost  Boost_ALLOW_PREFIND)

if (Boost_ALLOW_PREFIND)
  message("-- Using find_package(BOOST ...) ...")
  find_package(BOOST COMPONENTS program_options system)
  if (BOOST_FOUND)
    message("-- Found Boost_DIR='${BOOST_DIR}'")
    message("-- Generating Boost::all_libs and BoostConfig.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(BoostLib
      INNER_FIND_PACKAGE_NAME  Boost
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET BoostLib::all_libs)
  tribits_tpl_find_include_dirs_and_libraries( BoostLib
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} )
endif()
