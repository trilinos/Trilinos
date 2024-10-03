

#
# A) Define the package
#

TRIBITS_PACKAGE(Amesos)

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

ADVANCED_SET(Amesos_Allow_SuperLUDist_Without_ParMETIS OFF
  CACHE BOOL "Allow Amesos to use the SuperLUDist TPL without the enable of the ParMETIS TPL")
if (NOT Amesos_Allow_SuperLUDist_Without_ParMETIS)
  IF(${PACKAGE_NAME}_ENABLE_SuperLUDist AND NOT ${PACKAGE_NAME}_ENABLE_ParMETIS)
    MESSAGE(FATAL_ERROR "The Amesos support for the SuperLUDist TPL requires the ParMETIS TPL.  Either disable Amesos SuperLUDist support or enable the ParMETIS TPL.")
  ENDIF()
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_PARAKLETE)
  MESSAGE(FATAL_ERROR "Support for the PARAKLETE package has been removed from Amesos.")
ENDIF()

# if using MUMPS, must be parallel
IF(${PACKAGE_NAME}_ENABLE_MUMPS)
IF(NOT TPL_ENABLE_MPI)
   
   MESSAGE(WARNING "\n\n*** Amesos does not support MUMPS without MPI ** \n\n")
   
   SET(${PACKAGE_NAME}_ENABLE_MUMPS OFF)

ENDIF()
ENDIF()

# if using CSC, must be parallel and needs Pardiso_MKL
IF (${PACKAGE_NAME}_ENABLE_CSS_MKL)
  IF (NOT TPL_ENABLE_MPI)
    MESSAGE(WARNING "*****Amesos_ENABLE_CSS_MKL requires MPI. Turning off CSS_MKL. *****")
    SET(${PACKAGE_NAME}_ENABLE_CSS_MKL OFF)
    SET(HAVE_AMESOS_CSS_MKL OFF)
  ENDIF()
ENDIF()

# if using parallel MUMPS, must also link in scalapack and BLACS

IF(${PACKAGE_NAME}_ENABLE_MUMPS AND ${PACKAGE_NAME}_ENABLE_MPI AND NOT ${PACKAGE_NAME}_ENABLE_SCALAPACK)
  MESSAGE(FATAL_ERROR "The Amesos support for parallel MUMPS requires SCALAPACK.")
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_MUMPS AND ${PACKAGE_NAME}_ENABLE_MPI AND NOT ${PACKAGE_NAME}_ENABLE_BLACS)
  MESSAGE(FATAL_ERROR "The Amesos support for parallel MUMPS requires BLACS.")
ENDIF()


# if using SuperLU_5.0, must use the proper API
IF(${PACKAGE_NAME}_ENABLE_SuperLU)
  IF(NOT ${PROJECT_NAME}_ENABLE_SuperLU5_API)

   SET(HAVE_AMESOS_SUPERLU5_API OFF)

  ELSE()

   SET(HAVE_AMESOS_SUPERLU5_API ON)

  ENDIF()
ENDIF()

## Most of the options below are bogus.  The handling of TPLs is done
## through the dependencies file.  This needs to be cleaned up. I'm
## leaving them for now so that we have a list of all options.  They
## will be removed as correct TPL support is added.

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DSCPACK
  HAVE_AMESOS_DSCPACK
  "Enable DSCPACK functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MC64
  HAVE_AMESOS_MC64
  "Enable MC64 functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KLU
  HAVE_AMESOS_KLU
  "Enable KLU (Tim Davis's Files)."
  ON )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PARAKLETE
  HAVE_AMESOS_PARAKLETE
  "Enable PARAKLETE functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_LAPACK
  HAVE_AMESOS_LAPACK
  "Enable LAPACK functionality."
  ON )
  
# ...

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_FILES(
  doc/AmesosOverview
  doc/PARA06
  example/RunParaklete.cpp
  example/Thyra_AmesosLinearOpWithSolveFactory.cpp
  example/pk.h
  example/run_pk.c
  example/simpleStratimikosSolve.cpp
  example/simpleStratimikosSolve.hpp
  example/stratimikos_example.cpp
  src/Amesos_BTF.h
  src/Amesos_Component.h
  src/Amesos_Merikos.h
  src/Amesos_BTF.h
  src/src-repository
  src/stamp-h.in
  test/TestOptions/Dummy
  test/Test_Basic/NotQuiteDense.triU
  test/Test_Performance/In_Test_UmfpackPerformance.csh
  test/scripts/daily/mpi/TestBasic.sh
  test/scripts/daily/serial/TestAmesos.sh
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
