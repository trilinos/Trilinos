INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_BLAS_FLOAT_APPLE_VECLIB_BUGFIX VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})

  SET(SOURCE
  "
#include <complex>
#include <vecLib/cblas.h>

int main()
{
  const int n = 3;
  const float x[n] = { -1.0, -2.0, -2.0};
  const std::complex<float> xc[n] = { std::complex<float>(3.0,-4.0),
                                      std::complex<float>(0.0,0.0),
                                      std::complex<float>(0.0,0.0) };
  const int incx = 1;
  float val;

  val = cblas_sdot(n, x, incx, x, incx);
  if (val != 9.0) return 1;

  val = cblas_sasum(n, x, incx);
  if (val != 5.0) return 1;

  val = cblas_snrm2(n, x, incx);
  if (val != 3.0) return 1;

  val = cblas_scasum(n, xc, incx);
  if (val != 7.0) return 1;

  val = cblas_scnrm2(n, xc, incx);
  if (val != 5.0) return 1;

  return 0;
}
  "
  )

  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})

ENDFUNCTION()
