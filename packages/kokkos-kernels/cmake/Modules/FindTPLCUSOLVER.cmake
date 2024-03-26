if(CUSOLVER_LIBRARIES AND CUSOLVER_LIBRARY_DIRS AND CUSOLVER_INCLUDE_DIRS)
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES ${CUSOLVER_LIBRARIES}
    LIBRARY_PATHS ${CUSOLVER_LIBRARY_DIRS}
    HEADER_PATHS ${CUSOLVER_INCLUDE_DIRS}
  )
elseif(CUSOLVER_LIBRARIES AND CUSOLVER_LIBRARY_DIRS)
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES ${CUSOLVER_LIBRARIES}
    LIBRARY_PATHS ${CUSOLVER_LIBRARY_DIRS}
    HEADER cusolverDn.h
  )
elseif(CUSOLVER_LIBRARIES)
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES ${CUSOLVER_LIBRARIES}
    HEADER cusolverDn.h
  )
elseif(CUSOLVER_LIBRARY_DIRS)
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES cusolver
    LIBRARY_PATHS ${CUSOLVER_LIBRARY_DIRS}
    HEADER cusolverDn.h
  )
elseif(CUSOLVER_ROOT OR KokkosKernels_CUSOLVER_ROOT) # nothing specific provided, just ROOT
  kokkoskernels_find_imported(CUSOLVER INTERFACE
    LIBRARIES cusolver
    HEADER cusolverDn.h
  )
else() # backwards-compatible way
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
