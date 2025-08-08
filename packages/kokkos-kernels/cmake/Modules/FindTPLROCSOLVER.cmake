find_package(ROCSOLVER)
if(ROCSOLVER_FOUND)

  get_target_property(kk_rocsolver_include_dir_list roc::rocsolver INTERFACE_INCLUDE_DIRECTORIES)
  list(GET kk_rocsolver_include_dir_list 0 ROCSOLVER_INCLUDE_DIRS)

  set(kk_rocsolver_config "")
  get_target_property(kk_rocsolver_configs roc::rocsolver IMPORTED_CONFIGURATIONS)
  list(FIND kk_rocsolver_configs RELEASE rocsolver_release_config)
  list(FIND kk_rocsolver_configs RELWITHDEBINFO rocsolver_relwithdebinfo_config)
  list(FIND kk_rocsolver_configs DEBUG rocsolver_debug_config)

  if(${rocsolver_release_config} GREATER_EQUAL "0")
    list(GET kk_rocsolver_configs ${rocsolver_release_config} kk_rocsolver_config)
  elseif(${rocsolver_relwithdebinfo_config} GREATER_EQUAL "0")
    list(GET kk_rocsolver_configs ${rocsolver_relwithdebinfo_config} kk_rocsolver_config)
  elseif(${rocsolver_debug_config} GREATER_EQUAL "0")
    list(GET kk_rocsolver_configs ${rocsolver_debug_config} kk_rocsolver_config)
  endif()

  get_target_property(kk_rocsolver_library roc::rocsolver IMPORTED_LOCATION_${kk_rocsolver_config})
  get_filename_component(ROCSOLVER_LIBRARIES ${kk_rocsolver_library} DIRECTORY)


  kokkoskernels_find_imported(ROCSOLVER INTERFACE
    LIBRARIES rocsolver
    LIBRARY_PATHS ${ROCSOLVER_LIBRARIES}
    HEADER_PATHS "${ROCSOLVER_INCLUDE_DIRS};${ROCSOLVER_INCLUDE_DIRS}/rocsolver"
    HEADER rocsolver.h
  )

elseif(ROCSOLVER_LIBRARIES AND ROCSOLVER_LIBRARY_DIRS AND ROCSOLVER_INCLUDE_DIRS)
  message(WARNING "ROCSOLVER_LIBRARIES AND ROCSOLVER_LIBRARY_DIRS AND ROCSOLVER_INCLUDE_DIRS are deprecated please use ROCSOLVER_DIR or ROCSOLVER_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSOLVER INTERFACE
    LIBRARIES ${ROCSOLVER_LIBRARIES}
    LIBRARY_PATHS ${ROCSOLVER_LIBRARY_DIRS}
    HEADER_PATHS ${ROCSOLVER_INCLUDE_DIRS}
  )
elseif(ROCSOLVER_LIBRARIES AND ROCSOLVER_LIBRARY_DIRS)
  message(WARNING "ROCSOLVER_LIBRARIES AND ROCSOLVER_LIBRARY_DIRS are deprecated please use ROCSOLVER_DIR or ROCSOLVER_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSOLVER INTERFACE
    LIBRARIES ${ROCSOLVER_LIBRARIES}
    LIBRARY_PATHS ${ROCSOLVER_LIBRARY_DIRS}
    HEADER rocsolver.h
  )
elseif(ROCSOLVER_LIBRARIES)
  message(WARNING "ROCSOLVER_LIBRARIES are deprecated please use ROCSOLVER_DIR or ROCSOLVER_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSOLVER INTERFACE
    LIBRARIES ${ROCSOLVER_LIBRARIES}
    HEADER rocsolver.h
  )
elseif(ROCSOLVER_LIBRARY_DIRS)
  message(WARNING "ROCSOLVER_LIBRARY_DIRS are deprecated please use ROCSOLVER_DIR or ROCSOLVER_ROOT so find_packge logic can be used instead.")
  kokkoskernels_find_imported(ROCSOLVER INTERFACE
    LIBRARIES rocsolver
    LIBRARY_PATHS ${ROCSOLVER_LIBRARY_DIRS}
    HEADER rocsolver.h
  )
elseif(ROCSOLVER_ROOT OR KokkosKernels_ROCSOLVER_ROOT) # nothing specific provided, just ROOT
  set(ROCSOLVER_ROOT_DIR ${ROCSOLVER_ROOT})
  if(not ROCSOLVER_ROOT_DIR)
    set(ROCSOLVER_ROOT_DIR ${KokkosKernels_ROCSOLVER_ROOT})
  endif()
  kokkoskernels_find_imported(ROCSOLVER INTERFACE
    LIBRARIES rocsolver
    LIBRARY_PATHS ${ROCSOLVER_ROOT_DIR}/lib
    HEADER rocsolver.h
    HEADER_PATHS "${ROCSOLVER_ROOT_DIR}/include;${ROCSOLVER_ROOT_DIR}/include/rocsolver"
  )
endif()

