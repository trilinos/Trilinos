INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_BLAS_FLOAT VARNAME)

  CONFIGURE_FILE(
    ${PACKAGE_SOURCE_DIR}/cmake/FortranMangle.h.in
    ${PACKAGE_BINARY_DIR}/cmake/FortranMangle.h
    )

  SET(CMAKE_REQUIRED_INCLUDES ${PACKAGE_BINARY_DIR}/cmake)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BLAS_LIBRARIES})

  SET(SOURCE
  "
#include \"FortranMangle.h\"

#define SASUM_F77   F77_BLAS_MANGLE(sasum,SASUM)

extern \"C\" { float SASUM_F77(const int* n, const float x[], const int* incx); }

int main()
{

  const int n = 1;
  const float x[n] = { -1.0 };
  const int incx = 1;

  float val = SASUM_F77(&n, x, &incx);

  return (val == 1.0 ? 0 : 1);

}
  "
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()
