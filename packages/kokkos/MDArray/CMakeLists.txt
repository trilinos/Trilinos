INCLUDE(PackageLibraryMacros)
INCLUDE(AddOptionAndDefine)

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_MDArray_ThreadPool_Host
  HAVE_KOKKOS_MDARRAY_THREADPOOL_HOST
  "Enable mdarrays using threadpool."
  ON
  )

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_MDArray_Cuda
  HAVE_KOKKOS_MDARRAY_CUDA
  "Enable mdarrays using cuda."
  ${Kokkos_ENABLE_CUDA}
  )

# No longer needed; Trilinos does this at the top level.
# IF (Kokkos_ENABLE_MDArray_Cuda)
#   FIND_PACKAGE(CUDA)
# ENDIF()

#PACKAGE_ADD_TEST_DIRECTORIES(unit_tests)
