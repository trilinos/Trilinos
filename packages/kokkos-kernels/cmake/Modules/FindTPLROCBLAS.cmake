# MPL: 12/29/2022: CMake regular way to find a package
FIND_PACKAGE(ROCBLAS)
if(TARGET roc::rocblas)
## MPL: 12/29/2022: Variable TPL_ROCBLAS_IMPORTED_NAME follows the requested convention
## of KokkosKernel (method kokkoskernels_import_tpl of kokkoskernels_tpls.cmake)
  SET(TPL_ROCBLAS_IMPORTED_NAME roc::rocblas)
  SET(TPL_IMPORTED_NAME roc::rocblas)
## MPL: 12/29/2022: A target comming from a TPL must follows the requested convention
## of KokkosKernel (method kokkoskernels_link_tpl of kokkoskernels_tpls.cmake)
  ADD_LIBRARY(KokkosKernels::ROCBLAS ALIAS roc::rocblas)
ELSE()
  MESSAGE(FATAL_ERROR "Package ROCBLAS requested but not found")
ENDIF()
