find_package(ROCBLAS)
if(ROCBLAS_FOUND)

  get_target_property(kk_rocblas_include_dir_list roc::rocblas INTERFACE_INCLUDE_DIRECTORIES)
  list(GET kk_rocblas_include_dir_list 0 ROCBLAS_INCLUDE_DIRS)

  set(kk_rocblas_config "")
  get_target_property(kk_rocblas_configs roc::rocblas IMPORTED_CONFIGURATIONS)
  list(FIND kk_rocblas_configs RELEASE rocblas_release_config)
  list(FIND kk_rocblas_configs RELWITHDEBINFO rocblas_relwithdebinfo_config)
  list(FIND kk_rocblas_configs DEBUG rocblas_debug_config)

  if(${rocblas_release_config} GREATER_EQUAL "0")
    list(GET kk_rocblas_configs ${rocblas_release_config} kk_rocblas_config)
  elseif(${rocblas_relwithdebinfo_config} GREATER_EQUAL "0")
    list(GET kk_rocblas_configs ${rocblas_relwithdebinfo_config} kk_rocblas_config)
  elseif(${rocblas_debug_config} GREATER_EQUAL "0")
    list(GET kk_rocblas_configs ${rocblas_debug_config} kk_rocblas_config)
  endif()

  get_target_property(kk_rocblas_library roc::rocblas IMPORTED_LOCATION_${kk_rocblas_config})
  get_filename_component(ROCBLAS_LIBRARIES ${kk_rocblas_library} DIRECTORY)

  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES rocblas
    LIBRARY_PATHS ${ROCBLAS_LIBRARIES}
    HEADER rocblas.h
    HEADER_PATHS "${ROCBLAS_INCLUDE_DIRS};${ROCBLAS_INCLUDE_DIRS}/rocblas"
  )

elseif(ROCBLAS_LIBRARIES AND ROCBLAS_LIBRARY_DIRS AND ROCBLAS_INCLUDE_DIRS)
  message(WARNING "ROCBLAS_LIBRARIES AND ROCBLAS_LIBRARY_DIRS AND ROCBLAS_INCLUDE_DIRS are deprecated please use ROCBLAS_DIR or ROCBLAS_ROOT so find_package logic can be used instead.")
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES ${ROCBLAS_LIBRARIES}
    LIBRARY_PATHS ${ROCBLAS_LIBRARY_DIRS}
    HEADER_PATHS ${ROCBLAS_INCLUDE_DIRS}
  )
elseif(ROCBLAS_LIBRARIES AND ROCBLAS_LIBRARY_DIRS)
  message(WARNING "ROCBLAS_LIBRARIES AND ROCBLAS_LIBRARY_DIRS are deprecated please use ROCBLAS_DIR or ROCBLAS_ROOT so find_package logic can be used instead.")
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES ${ROCBLAS_LIBRARIES}
    LIBRARY_PATHS ${ROCBLAS_LIBRARY_DIRS}
    HEADER rocblas.h
  )
elseif(ROCBLAS_LIBRARIES)
  message(WARNING "ROCBLAS_LIBRARIES are deprecated please use ROCBLAS_DIR or ROCBLAS_ROOT so find_package logic can be used instead.")
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES ${ROCBLAS_LIBRARIES}
    HEADER rocblas.h
  )
elseif(ROCBLAS_LIBRARY_DIRS)
  message(WARNING "ROCBLAS_LIBRARY_DIRS are deprecated please use ROCBLAS_DIR or ROCBLAS_ROOT so find_package logic can be used instead.")
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES rocblas
    LIBRARY_PATHS ${ROCBLAS_LIBRARY_DIRS}
    HEADER rocblas.h
  )
else(ROCBLAS_ROOT OR KokkosKernels_ROCBLAS_ROOT) # nothing specific provided, just ROOT
  set(ROCBLAS_ROOT_DIR ${ROCBLAS_ROOT})
  IF(NOT ROCBLAS_ROOT_DIR)
    SET(ROCBLAS_ROOT_DIR ${KokkosKernels_ROCBLAS_ROOT})
  ENDIF()
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES rocblas
    LIBRARY_PATHS ${ROCBLAS_ROOT_DIR}/lib
    HEADER rocblas.h
    HEADER_PATHS "${ROCBLAS_ROOT_DIR}/include;${ROCBLAS_ROOT_DIR}/include/rocblas"
  )
endif()
