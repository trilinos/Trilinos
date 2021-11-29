include(TribitsGetImportedLocationProperty)


#
# Functions
#

function(tpl1_write_config_file tpl1Dir)
  set(configFileStr "")
  string(APPEND configFileStr
    "include(CMakeFindDependencyMacro)\n"
    "set(Tpl1_DIR \"${tpl1Dir}\")\n"
    "find_dependency(Tpl1)\n"
    "add_library(Tpl1::all_libs INTERFACE IMPORTED GLOBAL)\n"
    "target_link_libraries(Tpl1::all_libs INTERFACE tpl1::tpl1)\n"
    )
  set(buildDirExternalPkgsDir
    "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}")
  set(tplConfigFile
    "${buildDirExternalPkgsDir}/${TPL_NAME}/${TPL_NAME}Config.cmake")
  file(WRITE "${tplConfigFile}" "${configFileStr}")

endfunction()


#
# Executable Code
#

set(REQUIRED_HEADERS Tpl1.hpp)
set(REQUIRED_LIBS_NAMES tpl1)

tribits_tpl_allow_pre_find_package(Tpl1  Tpl1_ALLOW_PREFIND)

if (Tpl1_ALLOW_PREFIND)
  message("-- Using find_package(Tpl1 ...) ...")
  find_package(Tpl1)
  if (Tpl1_FOUND)
    message("-- Found Tpl1_DIR='${Tpl1_DIR}'")
    if (Tpl1_EXTRACT_INFO_AFTER_FIND_PACKAGE)
      message("-- Extracting include dirs and libraries from target tpl1::tpl1")
      get_target_property(inclDirs tpl1::tpl1 INTERFACE_INCLUDE_DIRECTORIES)
      tribits_get_imported_location_property(tpl1::tpl1 libfile)
      set(TPL_Tpl1_INCLUDE_DIRS "${inclDirs}" CACHE PATH "Include dirs for Tpl1")
      set(TPL_Tpl1_LIBRARIES "${libfile}" CACHE PATH "Libraries for Tpl1")
    else()
      # Create imported target Tpl1::all_libs
      add_library(Tpl1::all_libs INTERFACE IMPORTED GLOBAL)
      target_link_libraries(Tpl1::all_libs INTERFACE tpl1::tpl1)
      set(TPL_Tpl1_LIBRARIES Tpl1::all_libs CACHE STRING
        "Set in ${CMAKE_CURRENT_LIST_FILE}")
      set(TPL_Tpl1_INCLUDE_DIRS "" CACHE STRING
        "Set in ${CMAKE_CURRENT_LIST_FILE}")
      set(TPL_Tpl1_LIBRARY_DIRS "" CACHE STRING
        "Set in ${CMAKE_CURRENT_LIST_FILE}")
      print_var(TPL_Tpl1_LIBRARIES)
      print_var(TPL_Tpl1_INCLUDE_DIRS)
      # Write a specialized Tpl1Config.cmake file
      tpl1_write_config_file("${Tpl1_DIR}")
    endif()
  endif()
endif()

if (NOT TARGET Tpl1::all_libs)
  tribits_tpl_find_include_dirs_and_libraries( Tpl1
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
    )
endif()
