if ("${TPL_ENABLE_CUDA}" STREQUAL "ON")
  set(COMPADRE_REQUIRED_TPLS CUSOLVER CUBLAS)
else()
  set(COMPADRE_REQUIRED_TPLS LAPACK BLAS)
endif()

tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES
    KokkosCore KokkosContainers KokkosAlgorithms
  LIB_REQUIRED_TPLS
    ${COMPADRE_REQUIRED_TPLS}
  LIB_OPTIONAL_TPLS
    MPI CUDA
  )
