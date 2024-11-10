TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
        LIB_REQUIRED_PACKAGES Kokkos
        LIB_OPTIONAL_TPLS quadmath MKL BLAS LAPACK METIS SuperLU Cholmod CUBLAS CUSPARSE CUSOLVER ROCBLAS ROCSPARSE ROCSOLVER
        TEST_OPTIONAL_TPLS yamlcpp
)
# NOTE: If you update names in LIB_OPTIONAL_TPLS above, make sure to map those names in
# the macro 'KOKKOSKERNELS_ADD_TPL_OPTION' that resides in cmake/kokkoskernels_tpls.cmake.

if (TPL_ENABLE_CUDA)
  tribits_tpl_tentatively_enable(CUBLAS)
endif()

