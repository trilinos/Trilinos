find_package(CUDAToolkit)

if(CUDAToolkit_FOUND)
  get_target_property(kk_cusparse_include_dir_list CUDA::cusparse INTERFACE_INCLUDE_DIRECTORIES)
  list(GET kk_cusparse_include_dir_list 0 kk_cusparse_include_dir)
  get_target_property(kk_cusparse_library CUDA::cusparse IMPORTED_LOCATION)
  get_filename_component(kk_cusparse_library_dir ${kk_cusparse_library} DIRECTORY)
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES cusparse
    LIBRARY_PATHS ${kk_cusparse_library_dir}
    HEADER cusparse.h
    HEADER_PATHS ${kk_cusparse_include_dir}
  )
elseif(CUSPARSE_LIBRARIES AND CUSPARSE_LIBRARY_DIRS AND CUSPARSE_INCLUDE_DIRS)
  message(WARNING "CUSPARSE_LIBRARIES and CUSPARSE_LIBRARY_DIRS and CUSPARSE_INCLUDE_DIRS are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES ${CUSPARSE_LIBRARIES}
    LIBRARY_PATHS ${CUSPARSE_LIBRARY_DIRS}
    HEADER_PATHS ${CUSPARSE_INCLUDE_DIRS}
  )
elseif(CUSPARSE_LIBRARIES AND CUSPARSE_LIBRARY_DIRS)
  message(WARNING "CUSPARSE_LIBRARIES and CUSPARSE_LIBRARY_DIRS are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES ${CUSPARSE_LIBRARIES}
    LIBRARY_PATHS ${CUSPARSE_LIBRARY_DIRS}
    HEADER cusparse.h
  )
elseif(CUSPARSE_LIBRARIES)
  message(WARNING "CUSPARSE_LIBRARIES is deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES ${CUSPARSE_LIBRARIES}
    HEADER cusparse.h
  )
elseif(CUSPARSE_LIBRARY_DIRS)
  message(WARNING "CUSPARSE_LIBRARY_DIRS is deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES cusparse
    LIBRARY_PATHS ${CUSPARSE_LIBRARY_DIRS}
    HEADER cusparse.h
  )
elseif(CUSPARSE_ROOT OR KokkosKernels_CUSPARSE_ROOT) # nothing specific provided, just ROOT
  message(WARNING "CUSPARSE_ROOT and KokkosKernels_CUSPARSE_ROOT are deprecated. Use CUDAToolkit_ROOT to guide CMake's built-in search.")
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES cusparse
    HEADER cusparse.h
  )
elseif(CMAKE_VERSION VERSION_LESS "3.27")
  # backwards compatible way using FIND_PACKAGE(CUDA) (removed in 3.27)
  FIND_PACKAGE(CUDA)
  INCLUDE(FindPackageHandleStandardArgs)
  IF (NOT CUDA_FOUND)
    #Important note here: this find Module is named TPLCUSPARSE
    #The eventual target is named CUSPARSE. To avoid naming conflicts
    #the find module is called TPLCUSPARSE. This call will cause
    #the find_package call to fail in a "standard" CMake way
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLCUSPARSE REQUIRED_VARS CUDA_FOUND)
  ELSE()
    #The libraries might be empty - OR they might explicitly be not found
    IF("${CUDA_cusparse_LIBRARY}" MATCHES "NOTFOUND")
      FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLCUSPARSE REQUIRED_VARS CUDA_cusparse_LIBRARY)
    ELSE()
      KOKKOSKERNELS_CREATE_IMPORTED_TPL(CUSPARSE INTERFACE LINK_LIBRARIES "${CUDA_cusparse_LIBRARY}")
    ENDIF()
  ENDIF()
endif()
