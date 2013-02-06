INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_BLAS_FLOAT_DOUBLE_RETURN VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})

  SET(SOURCE
  "
#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

#include <complex>

#define SDOT_F77    F77_BLAS_MANGLE(sdot,   SDOT)
#define SASUM_F77   F77_BLAS_MANGLE(sasum,  SASUM)
#define SCASUM_F77  F77_BLAS_MANGLE(scasum, SCASUM)
#define SNRM2_F77   F77_BLAS_MANGLE(snrm2,  SNRM2)
#define SCNRM2_F77  F77_BLAS_MANGLE(scnrm2, SCNRM2)

extern \"C\" {
  double SDOT_F77(const int* n, const float x[], const int* incx, const float y[], const int* incy);
  double SASUM_F77(const int* n, const float x[], const int* incx);
  double SCASUM_F77(const int* n, const std::complex<float> x[], const int* incx);
  double SNRM2_F77(const int* n, const float x[], const int* incx);
  double SCNRM2_F77(const int* n, const std::complex<float> x[], const int* incx);
}

int main()
{
  const int n = 3;
  const float x[n] = { -1.0, -2.0, -2.0};
  const std::complex<float> xc[n] = { std::complex<float>(3.0,-4.0),
                                      std::complex<float>(0.0,0.0),
                                      std::complex<float>(0.0,0.0) };
  const int incx = 1;
  double val;

  val = SDOT_F77(&n, x, &incx, x, &incx);
  if (val != 9.0) return 1;

  val = SASUM_F77(&n, x, &incx);
  if (val != 5.0) return 1;

  val = SNRM2_F77(&n, x, &incx);
  if (val != 3.0) return 1;

  val = SCASUM_F77(&n, xc, &incx);
  if (val != 7.0) return 1;

  val = SCNRM2_F77(&n, xc, &incx);
  if (val != 5.0) return 1;

  return 0;
}
  "
  )

  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})

ENDFUNCTION()
