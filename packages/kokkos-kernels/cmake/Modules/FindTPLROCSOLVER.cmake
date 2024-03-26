# LBV: 11/08/2023: This file follows the partern of FindTPLROCBLAS.cmake/FindTPLROCSPARSE.cmake
FIND_PACKAGE(ROCSOLVER)
if(TARGET roc::rocsolver)
  SET(TPL_ROCSOLVER_IMPORTED_NAME roc::rocsolver)
  SET(TPL_IMPORTED_NAME roc::rocsolver)
  ADD_LIBRARY(KokkosKernels::ROCSOLVER ALIAS roc::rocsolver)
ELSE()
  MESSAGE(FATAL_ERROR "Package ROCSOLVER requested but not found")
ENDIF()
