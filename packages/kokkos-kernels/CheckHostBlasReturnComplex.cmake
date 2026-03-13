include(CheckCXXSourceRuns)

function(check_host_blas_return_complex VARNAME)

  if(KOKKOSKERNELS_HAS_TRILINOS)
    set(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})
  else()
    # For TPLs, just pull out the required libraries from the target properies.
    if(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
      get_target_property(CMAKE_REQUIRED_LIBRARIES KokkosKernels::ARMPL INTERFACE_LINK_LIBRARIES)
    else()
      set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
    endif()
  endif()

  set(SOURCE
    "
#include <complex>

#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

extern \"C\" {
  void F77_BLAS_MANGLE(zdotc,ZDOTC)(
    std::complex<double>* result, const int* n,
    const std::complex<double> x[], const int* incx,
    const std::complex<double> y[], const int* incy);
}

int main() {
  const int NUM=2;
  const int INC=1;
  std::complex<double> f[NUM];
  const std::complex<double>
    ONE = std::complex<double>(0.0,1.0),
    TWO = std::complex<double>(0.0,2.0);
  f[0] =  ONE;
  f[1] =  TWO;
  std::complex<double> ret;
  F77_BLAS_MANGLE(zdotc,ZDOTC)(&ret, &NUM, f, &INC, f, &INC);
  return (ret.real() == double(5.0) ? 0 : 1);
}
  ")

  # Test whether the above program, which assumes BLAS can give back complex results
  # via pointer arguments, compiles and runs correctly.
  # If it does, assume that we don't need to get complex results as direct return values,
  # which causes -Wreturn-type-c-linkage warnings.
  check_cxx_source_runs("${SOURCE}" KK_BLAS_RESULT_AS_POINTER_ARG)

  if(${KK_BLAS_RESULT_AS_POINTER_ARG})
    set(${VARNAME} OFF PARENT_SCOPE)
  else()
    set(${VARNAME} ON PARENT_SCOPE)
  endif()

endfunction()
