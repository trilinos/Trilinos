IF(ROCBLAS_LIBRARIES AND ROCBLAS_LIBRARY_DIRS AND ROCBLAS_INCLUDE_DIRS)
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES ${ROCBLAS_LIBRARIES}
    LIBRARY_PATHS ${ROCBLAS_LIBRARY_DIRS}
    HEADER_PATHS ${ROCBLAS_INCLUDE_DIRS}
  )
ELSEIF(ROCBLAS_LIBRARIES AND ROCBLAS_LIBRARY_DIRS)
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES ${ROCBLAS_LIBRARIES}
    LIBRARY_PATHS ${ROCBLAS_LIBRARY_DIRS}
    HEADER rocblas.h
  )
ELSEIF(ROCBLAS_LIBRARIES)
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES ${ROCBLAS_LIBRARIES}
    HEADER rocblas.h
  )
ELSEIF(ROCBLAS_LIBRARY_DIRS)
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES rocblas
    LIBRARY_PATHS ${ROCBLAS_LIBRARY_DIRS}
    HEADER rocblas.h
  )
ELSEIF(ROCBLAS_ROOT OR KokkosKernels_ROCBLAS_ROOT) # nothing specific provided, just ROOT
  kokkoskernels_find_imported(ROCBLAS INTERFACE
    LIBRARIES rocblas
    HEADER rocblas.h
  )
ELSE() # backwards-compatible way
  FIND_PACKAGE(ROCBLAS)
  INCLUDE(FindPackageHandleStandardArgs)
  IF (NOT ROCBLAS_FOUND)
    #Important note here: this find Module is named TPLROCBLAS
    #The eventual target is named ROCBLAS. To avoid naming conflicts
    #the find module is called TPLROCBLAS. This call will cause
    #the find_package call to fail in a "standard" CMake way
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLROCBLAS REQUIRED_VARS ROCBLAS_FOUND)
  ELSE()
    #The libraries might be empty - OR they might explicitly be not found
    IF("${ROCBLAS_LIBRARIES}" MATCHES "NOTFOUND")
      FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLROCBLAS REQUIRED_VARS ROCBLAS_LIBRARIES)
    ELSE()
      KOKKOSKERNELS_CREATE_IMPORTED_TPL(ROCBLAS INTERFACE
        LINK_LIBRARIES "${ROCBLAS_LIBRARIES}")
    ENDIF()
  ENDIF()
ENDIF()
