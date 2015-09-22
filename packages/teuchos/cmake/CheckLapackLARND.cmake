INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_LAPACK_LARND VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_LAPACK_LIBRARIES})

  SET(SOURCE
  "
#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

#define DLARND_F77   F77_BLAS_MANGLE(dlarnd,DLARND)

extern \"C\" { double DLARND_F77(const int* idist, int* seed); }

int main()
{

  const int idist = 1;
  int seed[4] = { 0.0, 0.0, 0.0, 1.0 };

  double val = DLARND_F77(&idist, seed);

  return (val < 0.0 ? 1 : 0);

}
  "
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()
