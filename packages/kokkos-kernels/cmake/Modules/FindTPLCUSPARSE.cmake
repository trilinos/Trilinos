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
     KOKKOSKERNELS_CREATE_IMPORTED_TPL(CUSPARSE LIBRARY ${CUDA_cusparse_LIBRARY})
  ENDIF()
ENDIF()
