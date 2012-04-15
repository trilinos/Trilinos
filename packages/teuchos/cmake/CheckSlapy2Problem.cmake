INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_SLAPY2_PROBLEM VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES  ${TPL_LAPACK_LIBRARIES} ${TPL_BLAS_LIBRARIES})

  SET(SOURCE
  "
#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

extern \"C\" {
  float F77_BLAS_MANGLE(slapy2,SLAPY2)(
    const float* x, const float* y);
}

int main() {
  float fx = 3.0f, fy = 4.0f;
  float fret = F77_BLAS_MANGLE(slapy2,SLAPY2)(&fx, &fy);
  return (fret == 5.0f ? 0 : 1);
}
  "
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" LAPACK_SLAPY2_WORKS)

  IF (LAPACK_SLAPY2_WORKS)
    SET(${VARNAME} FALSE PARENT_SCOPE)
  ELSE()
    SET(${VARNAME} TRUE PARENT_SCOPE)
  ENDIF()
  
ENDFUNCTION()
