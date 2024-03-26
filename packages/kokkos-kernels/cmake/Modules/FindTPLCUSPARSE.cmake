if(CUSPARSE_LIBRARIES AND CUSPARSE_LIBRARY_DIRS AND CUSPARSE_INCLUDE_DIRS)
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES ${CUSPARSE_LIBRARIES}
    LIBRARY_PATHS ${CUSPARSE_LIBRARY_DIRS}
    HEADER_PATHS ${CUSPARSE_INCLUDE_DIRS}
  )
elseif(CUSPARSE_LIBRARIES AND CUSPARSE_LIBRARY_DIRS)
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES ${CUSPARSE_LIBRARIES}
    LIBRARY_PATHS ${CUSPARSE_LIBRARY_DIRS}
    HEADER cusparse.h
  )
elseif(CUSPARSE_LIBRARIES)
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES ${CUSPARSE_LIBRARIES}
    HEADER cusparse.h
  )
elseif(CUSPARSE_LIBRARY_DIRS)
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES cusparse
    LIBRARY_PATHS ${CUSPARSE_LIBRARY_DIRS}
    HEADER cusparse.h
  )
elseif(CUSPARSE_ROOT OR KokkosKernels_CUSPARSE_ROOT) # nothing specific provided, just ROOT
  kokkoskernels_find_imported(CUSPARSE INTERFACE
    LIBRARIES cusparse
    HEADER cusparse.h
  )
else() # backwards-compatible way
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
