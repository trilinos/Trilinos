
IF (NOT TPL_ENABLE_CUDA)
  MESSAGE(FATAL_ERROR "\nCUSOLVER: This TPL requires CUDA")
ELSE()
  find_library(CUDA_cusolver_LIBRARY
    cusolver
    HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib
  )
  IF(CUDA_cusolver_LIBRARY STREQUAL "CUDA_cusolver_LIBRARY-NOTFOUND") 
    MESSAGE(FATAL_ERROR "\nCUSOLVER: could not find cusolver library.")
  ENDIF()
  SET(TPL_CUSOLVER_LIBRARIES ${CUDA_cusolver_LIBRARY})
ENDIF()

tribits_tpl_find_include_dirs_and_libraries(CUSOLVER  REQUIRED_LIBS_NAMES  cusparse)

unset(TPL_CUSOLVER_LIBRARIES)
