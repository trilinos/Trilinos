KOKKOSKERNELS_ADD_OPTION(
  INST_DOUBLE
  ${KOKKOSKERNELS_ADD_DEFAULT_ETI}
  BOOL
  "Whether to pre instantiate kernels for the scalar type double.  This option is KokkosKernels_INST_DOUBLE=ON by default.  Disabling this may increase build times. Default: ON"
  )

IF (Trilinos_ENABLE_FLOAT)
  SET(KOKKOSKERNELS_INST_FLOAT_DEFAULT  ON)
ELSE()
  SET(KOKKOSKERNELS_INST_FLOAT_DEFAULT  OFF)
ENDIF()

KOKKOSKERNELS_ADD_OPTION(
  INST_FLOAT
  ${KOKKOSKERNELS_INST_FLOAT_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the scalar type float.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_FLOAT."
  )

SET(FLOATS
  FLOAT
  DOUBLE
  COMPLEX_FLOAT
  COMPLEX_DOUBLE)
SET(DOUBLE_CPP_TYPE "double")
SET(FLOAT_CPP_TYPE "float")
SET(COMPLEX_FLOAT_CPP_TYPE "Kokkos::complex<float>")
SET(COMPLEX_DOUBLE_CPP_TYPE "Kokkos::complex<double>")

#Just leave the Trilinos logic here alone
#Outside of Trilinos, let the defaults be always OFF
IF (KOKKOSKERNELS_INST_DOUBLE AND Trilinos_ENABLE_COMPLEX_DOUBLE)
  SET(KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT ON)
ELSE()
  SET(KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT OFF)
ENDIF()
IF (KOKKOSKERNELS_INST_FLOAT AND Trilinos_ENABLE_COMPLEX_FLOAT)
  SET(KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT ON)
ELSE()
  SET(KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT OFF)
ENDIF()

KOKKOSKERNELS_ADD_OPTION(
  INST_COMPLEX_DOUBLE
  ${KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the scalar type complex<double>.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_DOUBLE."
  )

KOKKOSKERNELS_ADD_OPTION(
  INST_COMPLEX_FLOAT
  ${KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the scalar type complex<float>.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_FLOAT."
  )

IF (KOKKOSKERNELS_INST_DOUBLE)
  LIST(APPEND SCALAR_LIST "double")
ENDIF()

IF (KOKKOSKERNELS_INST_FLOAT)
  LIST(APPEND SCALAR_LIST "float")
ENDIF()

IF (KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
  LIST(APPEND SCALAR_LIST "complex<double>")
ENDIF()

IF (KOKKOSKERNELS_INST_COMPLEX_FLOAT)
  LIST(APPEND SCALAR_LIST "complex<float>")
ENDIF()
