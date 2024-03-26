INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_HOST_BLAS_RETURN_COMPLEX VARNAME)

  IF (KOKKOSKERNELS_HAS_TRILINOS)
    SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})
  ELSE()
    # For TPLs, just pull out the required libraries from the target properies.
    IF (KOKKOSKERNELS_ENABLE_TPL_ARMPL)
      GET_TARGET_PROPERTY(CMAKE_REQUIRED_LIBRARIES KokkosKernels::ARMPL INTERFACE_LINK_LIBRARIES)
    ELSE()
      SET(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
    ENDIF()
  ENDIF()

  SET(SOURCE
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
  "
  )

# Test whether the above program, which assumes BLAS can give back complex results
# via pointer arguments, compiles and runs correctly.
# If it does, assume that we don't need to get complex results as direct return values,
# which causes -Wreturn-type-c-linkage warnings.
CHECK_CXX_SOURCE_RUNS("${SOURCE}" KK_BLAS_RESULT_AS_POINTER_ARG)

IF(${KK_BLAS_RESULT_AS_POINTER_ARG})
  SET(${VARNAME} OFF PARENT_SCOPE)
ELSE()
  SET(${VARNAME} ON PARENT_SCOPE)
ENDIF()

ENDFUNCTION()
