TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
        LIB_REQUIRED_PACKAGES KokkosCore KokkosContainers KokkosAlgorithms
        LIB_OPTIONAL_TPLS quadmath MKL BLAS LAPACK CUSPARSE MAGMA METIS SuperLU Cholmod LAPACKE CBLAS ARMPL ROCBLAS ROCSPARSE
        TEST_OPTIONAL_TPLS yaml-cpp
)
# NOTE: If you update names in LIB_OPTIONAL_TPLS above, make sure to map those names in
# the macro 'KOKKOSKERNELS_ADD_TPL_OPTION' that resides in cmake/kokkoskernels_tpls.cmake.