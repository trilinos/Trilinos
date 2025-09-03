find_package(CUDAToolkit)

if(CUDAToolkit_FOUND)
  get_target_property(kk_cublas_include_dir_list CUDA::cublas INTERFACE_INCLUDE_DIRECTORIES)
  list(GET kk_cublas_include_dir_list 0 kk_cublas_include_dir)
  get_target_property(kk_cublas_library CUDA::cublas IMPORTED_LOCATION)
  get_filename_component(kk_cublas_library_dir ${kk_cublas_library} DIRECTORY)
  kokkoskernels_find_imported(CUBLAS INTERFACE
    LIBRARIES cublas
    LIBRARY_PATHS ${kk_cublas_library_dir}
    HEADER cublas.h
    HEADER_PATHS ${kk_cublas_include_dir}
  )
elseif(CUBLAS_LIBRARIES AND CUBLAS_LIBRARY_DIRS AND CUBLAS_INCLUDE_DIRS)
  message(WARNING "CUBLAS_LIBRARIES and CUBLAS_LIBRARY_DIRS and CUBLAS_INCLUDE_DIRS are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUBLAS INTERFACE
    LIBRARIES ${CUBLAS_LIBRARIES}
    LIBRARY_PATHS ${CUBLAS_LIBRARY_DIRS}
    HEADER_PATHS ${CUBLAS_INCLUDE_DIRS}
  )
elseif(CUBLAS_LIBRARIES AND CUBLAS_LIBRARY_DIRS)
  message(WARNING "CUBLAS_LIBRARIES and CUBLAS_LIBRARY_DIRS are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUBLAS INTERFACE
    LIBRARIES ${CUBLAS_LIBRARIES}
    LIBRARY_PATHS ${CUBLAS_LIBRARY_DIRS}
    HEADER cublas.h
  )
elseif(CUBLAS_LIBRARIES)
  message(WARNING "CUBLAS_LIBRARIES is deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUBLAS INTERFACE
    LIBRARIES ${CUBLAS_LIBRARIES}
    HEADER cublas.h
  )
elseif(CUBLAS_LIBRARY_DIRS)
  message(WARNING "CUBLAS_LIBRARY_DIRS is deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUBLAS INTERFACE
    LIBRARIES cublas
    LIBRARY_PATHS ${CUBLAS_LIBRARY_DIRS}
    HEADER cublas.h
  )
elseif(CUBLAS_ROOT OR KokkosKernels_CUBLAS_ROOT) # nothing specific provided, just ROOT
  message(WARNING "CUBLAS_ROOT and KokkosKernels_CUBLAS_ROOT are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUBLAS INTERFACE
    LIBRARIES cublas
    HEADER cublas.h
  )
elseif(CMAKE_VERSION VERSION_LESS "3.27")
  # backwards compatible way using FIND_PACKAGE(CUDA) (removed in 3.27)
  FIND_PACKAGE(CUDA)
  INCLUDE(FindPackageHandleStandardArgs)
  IF (NOT CUDA_FOUND)
    #Important note here: this find Module is named TPLCUBLAS
    #The eventual target is named CUBLAS. To avoid naming conflicts
    #the find module is called TPLCUBLAS. This call will cause
    #the find_package call to fail in a "standard" CMake way
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLCUBLAS REQUIRED_VARS CUDA_FOUND)
  ELSE()
    #The libraries might be empty - OR they might explicitly be not found
    IF("${CUDA_CUBLAS_LIBRARIES}" MATCHES "NOTFOUND")
      FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLCUBLAS REQUIRED_VARS CUDA_CUBLAS_LIBRARIES)
    ELSE()
      KOKKOSKERNELS_CREATE_IMPORTED_TPL(CUBLAS INTERFACE
        LINK_LIBRARIES "${CUDA_CUBLAS_LIBRARIES}")
    ENDIF()
  ENDIF()
endif()
