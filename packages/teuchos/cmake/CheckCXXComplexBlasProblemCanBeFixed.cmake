INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_CXX_COMPLEX_BLAS_PROBLEM_CAN_BE_FIXED VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})

  SET(SOURCE
  "
#include <complex>

#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

extern \"C\" {
  std::complex<float> F77_BLAS_MANGLE(cdotc,CDOTC)(
    std::complex<float> *ret, const int* n, const  std::complex<float> x[],
    const int* incx, const std::complex<float> y[], const int* incy); 
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
  std::complex<float> ret(0.0,0.0);
  F77_BLAS_MANGLE(cdotc,CDOTC)(&ret, &NUM, f, &INC, f, &INC);
  return (ret == TWO ? 0 : 1);
}
  "
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()

FUNCTION(CHECK_COMPLEX_BLAS_VECLIB_OKAY VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})

  SET(SOURCE
  "
#include <complex>
#include <vecLib/cblas.h>

int main() {
  const int NUM=2;
  const int INC=1;
  std::complex<float> f[NUM];
  const std::complex<float>
    ONE = std::complex<float>(1.0,0.0),
    TWO = std::complex<float>(2.0,0.0);
  f[0] =  ONE;
  f[1] =  ONE;
  std::complex<float> ret(0.0,0.0);
  cblas_cdotc_sub(NUM,f,INC,f,INC,&ret);
  return (ret == TWO ? 0 : 1);
}
  "
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()
