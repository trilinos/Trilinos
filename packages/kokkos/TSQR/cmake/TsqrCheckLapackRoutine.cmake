INCLUDE(CheckCXXSourceRuns)
INCLUDE(TsqrBlasAndLapack)

FUNCTION(TSQR_CHECK_LAPACK_ROUTINE ROUTINE VARNAME)

  IF (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
    message (STATUS "Checking whether your LAPACK library includes ${ROUTINE}")
  ENDIF()

  # Remember the current CMAKE_REQUIRED_LIBRARIES so we can restore it
  # after we are done testing.
  set (CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
  set (CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${TPL_LAPACK_LIBRARIES} ${TPL_BLAS_LIBRARIES})

  #message (STATUS "CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")

  # Define lowercase and uppercase versions of the routine name.
  # C(++) macros don't give you the full power of the language; their
  # compile-time computational abilities are limited.  We have to
  # compute lowercase and uppercase string manipulations for them,
  # which is what we are doing here.
  string (TOLOWER ${ROUTINE} LOWERCASE_ROUTINE)
  string (TOUPPER ${ROUTINE} UPPERCASE_ROUTINE)

  # We enclose the call to the routine in "if (argc > 1)" so that it
  # won't be invoked.  This is really just a link test analogous to
  # CHECK_FUNCTION_EXISTS, but we have to implement this way because:
  #
  # 1. We can't assume that Trilinos is being built with Fortran
  #    support, so we can't call CHECK_FORTRAN_FUNCTION_EXISTS.
  #
  # 2. That means we need to access the Fortran function through C++,
  #    which means we need to use the F77_BLAS_MANGLE macro.  This
  #    means we can't just call CHECK_FUNCTION_EXISTS.
  SET(SOURCE
  "
#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

#define ROUTINE_F77   F77_BLAS_MANGLE(${LOWERCASE_ROUTINE},${UPPERCASE_ROUTINE})

extern \"C\" { void ROUTINE_F77 (); }

int main(int argc, char* argv[]) {
  if (argc > 1)
    ROUTINE_F77();

  return 0;
}
"
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})
  
  # Restore the original value of CMAKE_REQUIRED_LIBRARIES.
  set (CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVE})

ENDFUNCTION()
