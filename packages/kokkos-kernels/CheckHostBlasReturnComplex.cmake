INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_HOST_BLAS_RETURN_COMPLEX VARNAME)

  IF (KOKKOSKERNELS_HAS_TRILINOS)
    SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})
  ELSE()
    SET(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
  ENDIF()

  SET(SOURCE
  "
#include <complex>

#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

extern \"C\" {
  std::complex<double> F77_BLAS_MANGLE(zdotc,ZDOTC)(
    const int* n, 
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
  std::complex<double> ret
   = F77_BLAS_MANGLE(zdotc,ZDOTC)(&NUM, f, &INC, f, &INC);
  return (ret.real() == double(5.0) ? 0 : 1);
}
  "
  )

  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})

ENDFUNCTION()
