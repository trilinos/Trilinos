kokkoskernels_add_option("INST_DOUBLE" ${KOKKOSKERNELS_ADD_DEFAULT_ETI} BOOL
  "Whether to pre instantiate kernels for the scalar type double.  This option is KokkosKernels_INST_DOUBLE=ON by default.  Disabling this may increase build times. Default: ON")

if(Trilinos_ENABLE_FLOAT)
  set(KOKKOSKERNELS_INST_FLOAT_DEFAULT ON)
else()
  set(KOKKOSKERNELS_INST_FLOAT_DEFAULT OFF)
endif()

kokkoskernels_add_option("INST_FLOAT" ${KOKKOSKERNELS_INST_FLOAT_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the scalar type float.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_FLOAT.")

kokkoskernels_add_option("INST_HALF" OFF BOOL
  "Whether to pre instantiate kernels for the scalar type Kokkos::Experimental::half_t.  Disabling this may increase build times. Default: OFF")

kokkoskernels_add_option("INST_BHALF" OFF BOOL
  "Whether to pre instantiate kernels for the scalar type Kokkos::Experimental::bhalf_t.  Disabling this may increase build times. Default: OFF")

set(REAL_FLOATS              FLOAT DOUBLE)
set(FLOATS                   FLOAT DOUBLE COMPLEX_FLOAT COMPLEX_DOUBLE)
set(DOUBLE_CPP_TYPE          "double")
set(FLOAT_CPP_TYPE           "float")
set(HALF_CPP_TYPE            "Kokkos::Experimental::half_t")
set(BHALF_CPP_TYPE           "Kokkos::Experimental::bhalf_t")
set(COMPLEX_FLOAT_CPP_TYPE   "Kokkos::complex<float>")
set(COMPLEX_DOUBLE_CPP_TYPE  "Kokkos::complex<double>")

#Just leave the Trilinos logic here alone
#Outside of Trilinos, let the defaults be always OFF
if(KOKKOSKERNELS_INST_DOUBLE AND Trilinos_ENABLE_COMPLEX_DOUBLE)
  set(KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT ON)
else()
  set(KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT OFF)
endif()
if(KOKKOSKERNELS_INST_FLOAT AND Trilinos_ENABLE_COMPLEX_FLOAT)
  set(KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT ON)
else()
  set(KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT OFF)
endif()

kokkoskernels_add_option("INST_COMPLEX_DOUBLE" ${KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the scalar type complex<double>.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_DOUBLE.")

kokkoskernels_add_option("INST_COMPLEX_FLOAT" ${KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the scalar type complex<float>.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_FLOAT.")

if(KOKKOSKERNELS_INST_DOUBLE)
  list(APPEND SCALAR_LIST "double")
endif()

if(KOKKOSKERNELS_INST_FLOAT)
  list(APPEND SCALAR_LIST "float")
endif()

# TODO: Fix build errors in kokkos when half_t is used in ETI
#IF (KOKKOSKERNELS_INST_HALF)
#  LIST(APPEND SCALAR_LIST "Kokkos::Experimental::half_t")
#ENDIF()

if(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
  list(APPEND SCALAR_LIST "complex<double>")
endif()

if(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
  list(APPEND SCALAR_LIST "complex<float>")
endif()
