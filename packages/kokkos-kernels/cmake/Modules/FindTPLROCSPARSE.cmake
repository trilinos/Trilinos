# MPL: 05/01/2023: This file follows the partern of FindTPLROCBLAS.cmake
FIND_PACKAGE(ROCSPARSE)
if(TARGET roc::rocsparse)
  SET(TPL_ROCSPARSE_IMPORTED_NAME roc::rocsparse)
  SET(TPL_IMPORTED_NAME roc::rocsparse)
  ADD_LIBRARY(KokkosKernels::ROCSPARSE ALIAS roc::rocsparse)
ELSE()
  MESSAGE(FATAL_ERROR "Package ROCSPARSE requested but not found")
ENDIF()
