
#
# A) Define the package
#

TRIBITS_PACKAGE(Amesos2)

SET(${PACKAGE_NAME_UC}_VERSION "0.3d")

#
# B) Set up package-specific options
#

IF (${PACKAGE_NAME}_ENABLE_MUMPS AND NOT Tpetra_INST_INT_INT)
   MESSAGE(WARNING "*****Amesos2 will not provide MUMPS Support with Tpetra_INST_INT_INT OFF.*****")
  SET(HAVE_${PACKAGE_NAME_UC}_MUMPS OFF)
ENDIF ()

IF (Amesos2_ENABLE_UMFPACK)
  IF (Trilinos_ENABLE_COMPLEX_FLOAT AND NOT Trilinos_ENABLE_COMPLEX_DOUBLE)
    MESSAGE(WARNING "*****Amesos2 is building Umfpack without a complex<double> type. Umfpack natively supports only complex<double> types, not complex<float> types, and this will lead to type conversions. *****")
  ENDIF()
  # This may not happen anyways… right now it seems we can’t compile this package with Trilinos_ENABLE_FLOAT and TPetra_INST_DOUBLE turned off - need to resolve.
  IF (NOT Tpetra_INST_DOUBLE)
    MESSAGE(WARNING "*****Amesos2 is building Umfpack without a double type. Umfpack natively supports only double types, not float types, and this will lead to type conversions. *****")
  ENDIF()
ENDIF()

IF (Amesos2_ENABLE_CSS_MKL)
  IF (NOT ${PACKAGE_NAME}_ENABLE_PARDISO_MKL OR NOT TPL_ENABLE_MPI)
    MESSAGE(WARNING "*****Amesos2_ENABLE_CSS_MKL requires PARDISO_MKL and MPI. Turning off CSS_MKL. *****")
    SET(Amesos2_ENABLE_CSS_MKL OFF)
    SET(HAVE_AMESOS2_CSS_MKL OFF)
  ENDIF()
ENDIF()

# if using SuperLU_5.0, must use the proper API
IF(${PACKAGE_NAME}_ENABLE_SuperLU)
  IF(NOT ${PROJECT_NAME}_ENABLE_SuperLU5_API)

   SET(HAVE_AMESOS2_SUPERLU5_API OFF)

  ELSE()

   SET(HAVE_AMESOS2_SUPERLU5_API ON)

  ENDIF()
ENDIF()


TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_VERBOSE_DEBUG
  HAVE_${PACKAGE_NAME_UC}_VERBOSE_DEBUG
  "Enable verbose output during debug-mode"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TIMERS
  HAVE_${PACKAGE_NAME_UC}_TIMERS
  "Enable performance timers"
  OFF )

#
# TPL Solver interfaces are enabled automatically based on Trilinos
# TPL configuration.
#
# Config options such as Amesos2_ENABLE_XXXXXX are automatically
# generated in case a user would wish to explicitely disable Amesos2
# support for an otherwise enabled third-part solver library
#


TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KLU2
  HAVE_AMESOS2_KLU2
  "Enable KLU2 in Amesos2"
  ON
)

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Basker
   HAVE_AMESOS2_BASKER
   "Enable Basker in Amesos2"
   OFF
)

IF (NOT Kokkos_ENABLE_OPENMP OR Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP)
  tribits_disable_optional_dependency(ShyLU_NodeBasker "NOTE: Disabling ${PACKAGE_NAME}_ENABLE_ShyLU_NodeBasker because CUDA or HIP is enabled in Kokkos or OpenMP is not enabled in Kokkos")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_AMESOS2_EXPERIMENTAL
  "Enable building and installation of experimental Amesos2 features."
  NO 
)

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

#leaving in place for future excludes
#TRIBITS_EXCLUDE_FILES(
#  )

TRIBITS_PACKAGE_POSTPROCESS()
