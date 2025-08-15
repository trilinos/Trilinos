find_package(CUDAToolkit)

if(CUDAToolkit_FOUND)
  get_target_property(kk_cusolver_include_dir_list CUDA::cusolver INTERFACE_INCLUDE_DIRECTORIES)
  list(GET kk_cusolver_include_dir_list 0 kk_cusolver_include_dir)
  get_target_property(kk_cusolver_library CUDA::cusolver IMPORTED_LOCATION)
  get_filename_component(kk_cusolver_library_dir ${kk_cusolver_library} DIRECTORY)
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES cusolver
    LIBRARY_PATHS ${kk_cusolver_library_dir}
    HEADER cusolverDn.h
    HEADER_PATHS ${kk_cusolver_include_dir}
  )
elseif(CUSOLVER_LIBRARIES AND CUSOLVER_LIBRARY_DIRS AND CUSOLVER_INCLUDE_DIRS)
  message(WARNING "CUSOLVER_LIBRARIES and CUSOLVER_LIBRARY_DIRS and CUSOLVER_INCLUDE_DIRS are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES ${CUSOLVER_LIBRARIES}
    LIBRARY_PATHS ${CUSOLVER_LIBRARY_DIRS}
    HEADER_PATHS ${CUSOLVER_INCLUDE_DIRS}
  )
elseif(CUSOLVER_LIBRARIES AND CUSOLVER_LIBRARY_DIRS)
  message(WARNING "CUSOLVER_LIBRARIES and CUSOLVER_LIBRARY_DIRS are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES ${CUSOLVER_LIBRARIES}
    LIBRARY_PATHS ${CUSOLVER_LIBRARY_DIRS}
    HEADER cusolverDn.h
  )
elseif(CUSOLVER_LIBRARIES)
  message(WARNING "CUSOLVER_LIBRARIES is deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES ${CUSOLVER_LIBRARIES}
    HEADER cusolverDn.h
  )
elseif(CUSOLVER_LIBRARY_DIRS)
  message(WARNING "CUSOLVER_LIBRARY_DIRS is deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES cusolver
    LIBRARY_PATHS ${CUSOLVER_LIBRARY_DIRS}
    HEADER cusolverDn.h
  )
elseif(CUSOLVER_ROOT OR KokkosKernels_CUSOLVER_ROOT) # nothing specific provided, just ROOT
  message(WARNING "CUSOLVER_ROOT and KokkosKernels_CUSOLVER_ROOT are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES cusolver
    HEADER cusolverDn.h
  )
elseif(CMAKE_VERSION VERSION_LESS "3.27")
  # backwards compatible way using FIND_PACKAGE(CUDA) (removed in 3.27)
  FIND_PACKAGE(CUDA)
  INCLUDE(FindPackageHandleStandardArgs)
  IF (NOT CUDA_FOUND)
    #Important note here: this find Module is named TPLCUSOLVER
    #The eventual target is named CUSOLVER. To avoid naming conflicts
    #the find module is called TPLCUSOLVER. This call will cause
    #the find_package call to fail in a "standard" CMake way
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLCUSOLVER REQUIRED_VARS CUDA_FOUND)
  ELSE()
    #The libraries might be empty - OR they might explicitly be not found
    IF("${CUDA_cusolver_LIBRARY}" MATCHES "NOTFOUND")
      FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLCUSOLVER REQUIRED_VARS CUDA_cusolver_LIBRARY)
    ELSE()
      KOKKOSKERNELS_CREATE_IMPORTED_TPL(CUSOLVER INTERFACE LINK_LIBRARIES "${CUDA_cusolver_LIBRARY}")
    ENDIF()
  ENDIF()
endif()
