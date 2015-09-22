INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_CXX_COMPLEX_BLAS_PROBLEM VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})

  SET(SOURCE
  "
#include <complex>

#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

extern \"C\" {
  std::complex<float> F77_BLAS_MANGLE(cdotc,CDOTC)(
    const int* n, const  std::complex<float> x[], const int* incx,
    const std::complex<float> y[], const int* incy); 
}

int main() {
  const int NUM=2;
  const int INC=1;
  std::complex<float> f[NUM];
  const std::complex<float>
    ONE = std::complex<float>(1.0,0.0),
    TWO = std::complex<float>(2.0,0.0);
  f[0] =  ONE;
  f[1] =  ONE;
  const std::complex<float> ret =
    F77_BLAS_MANGLE(cdotc,CDOTC)(&NUM, f, &INC, f, &INC);
  return (ret == TWO ? 0 : 1);
}
  "
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" CXX_COMPLEX_BLAS_WORKS)

  IF (CXX_COMPLEX_BLAS_WORKS)
    GLOBAL_SET(${VARNAME} FALSE)
  ELSE()
    GLOBAL_SET(${VARNAME} TRUE)
  ENDIF()
  
ENDFUNCTION()
