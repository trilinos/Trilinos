find_package(ROCSPARSE)
if(ROCSPARSE_FOUND)

  get_target_property(kk_rocsparse_include_dir_list roc::rocsparse INTERFACE_INCLUDE_DIRECTORIES)
  list(GET kk_rocsparse_include_dir_list 0 ROCSPARSE_INCLUDE_DIRS)

  set(kk_rocsparse_config "")
  get_target_property(kk_rocsparse_configs roc::rocsparse IMPORTED_CONFIGURATIONS)
  list(FIND kk_rocsparse_configs RELEASE rocsparse_release_config)
  list(FIND kk_rocsparse_configs RELWITHDEBINFO rocsparse_relwithdebinfo_config)
  list(FIND kk_rocsparse_configs DEBUG rocsparse_debug_config)

  if(${rocsparse_release_config} GREATER_EQUAL "0")
    list(GET kk_rocsparse_configs ${rocsparse_release_config} kk_rocsparse_config)
  elseif(${rocsparse_relwithdebinfo_config} GREATER_EQUAL "0")
    list(GET kk_rocsparse_configs ${rocsparse_relwithdebinfo_config} kk_rocsparse_config)
  elseif(${rocsparse_debug_config} GREATER_EQUAL "0")
    list(GET kk_rocsparse_configs ${rocsparse_debug_config} kk_rocsparse_config)
  endif()

  get_target_property(kk_rocsparse_library roc::rocsparse IMPORTED_LOCATION_${kk_rocsparse_config})
  get_filename_component(ROCSPARSE_LIBRARIES ${kk_rocsparse_library} DIRECTORY)


  kokkoskernels_find_imported(ROCSPARSE INTERFACE
    LIBRARIES rocsparse
    LIBRARY_PATHS ${ROCSPARSE_LIBRARIES}
    HEADER_PATHS "${ROCSPARSE_INCLUDE_DIRS};${ROCSPARSE_INCLUDE_DIRS}/rocsparse"
    HEADER rocsparse.h
  )

elseif(ROCSPARSE_LIBRARIES AND ROCSPARSE_LIBRARY_DIRS AND ROCSPARSE_INCLUDE_DIRS)
  message(WARNING "ROCSPARSE_LIBRARIES AND ROCSPARSE_LIBRARY_DIRS AND ROCSPARSE_INCLUDE_DIRS are deprecated please use ROCSPARSE_DIR or ROCSPARSE_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSPARSE INTERFACE
    LIBRARIES ${ROCSPARSE_LIBRARIES}
    LIBRARY_PATHS ${ROCSPARSE_LIBRARY_DIRS}
    HEADER_PATHS ${ROCSPARSE_INCLUDE_DIRS}
  )
elseif(ROCSPARSE_LIBRARIES AND ROCSPARSE_LIBRARY_DIRS)
  message(WARNING "ROCSPARSE_LIBRARIES AND ROCSPARSE_LIBRARY_DIRS are deprecated please use ROCSPARSE_DIR or ROCSPARSE_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSPARSE INTERFACE
    LIBRARIES ${ROCSPARSE_LIBRARIES}
    LIBRARY_PATHS ${ROCSPARSE_LIBRARY_DIRS}
    HEADER rocsparse.h
  )
elseif(ROCSPARSE_LIBRARIES)
  message(WARNING "ROCSPARSE_LIBRARIES are deprecated please use ROCSPARSE_DIR or ROCSPARSE_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSPARSE INTERFACE
    LIBRARIES ${ROCSPARSE_LIBRARIES}
    HEADER rocsparse.h
  )
elseif(ROCSPARSE_LIBRARY_DIRS)
  message(WARNING "ROCSPARSE_LIBRARY_DIRS are deprecated please use ROCSPARSE_DIR or ROCSPARSE_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSPARSE INTERFACE
    LIBRARIES rocsparse
    LIBRARY_PATHS ${ROCSPARSE_LIBRARY_DIRS}
    HEADER rocsparse.h
  )
elseif(ROCSPARSE_ROOT OR KokkosKernels_ROCSPARSE_ROOT) # nothing specific provided, just ROOT
  set(ROCSPARSE_ROOT_DIR ${ROCSPARSE_ROOT})
  if(not ROCSPARSE_ROOT_DIR)
    set(ROCSPARSE_ROOT_DIR ${KokkosKernels_ROCSPARSE_ROOT})
  endif()
  kokkoskernels_find_imported(ROCSPARSE INTERFACE
    LIBRARIES rocsparse
    LIBRARY_PATHS ${ROCSPARSE_ROOT_DIR}/lib
    HEADER rocsparse.h
    HEADERS_PATHS "${ROCSPARSE_ROOT_DIR}/include;${ROCSPARSE_ROOT_DIR}/include/rocsparse"
  )
endif()
