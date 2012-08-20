
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(Amesos)

#
# B) Set up package-specific options
#

# if using SuperLUDist, must also link in ParMETIS for some reason
IF(${PACKAGE_NAME}_ENABLE_SuperLUDist AND NOT ${PACKAGE_NAME}_ENABLE_ParMETIS)
  MESSAGE(FATAL_ERROR "The Amesos support for the SuperLUIDist TPL requires the ParMETIS TPL.  Either disable Amesos SuperLUDist support or enable the ParMETIS TPL.")
ENDIF()

# if using parallel MUMPS, must also link in scalapack and BLACS

IF(${PACKAGE_NAME}_ENABLE_MUMPS AND ${PACKAGE_NAME}_ENABLE_MPI AND NOT ${PACKAGE_NAME}_ENABLE_SCALAPACK)
  MESSAGE(FATAL_ERROR "The Amesos support for parallel MUMPS requires SCALAPACK.")
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_MUMPS AND ${PACKAGE_NAME}_ENABLE_MPI AND NOT ${PACKAGE_NAME}_ENABLE_BLACS)
  MESSAGE(FATAL_ERROR "The Amesos support for parallel MUMPS requires BLACS.")
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

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PARDISO
  HAVE_AMESOS_PARDISO
  "Enable PARDISO functionality."
  OFF )

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

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

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
  src/SuiteSparse/CAMD/Source/amesos_camd_1.c
  src/SuiteSparse/CAMD/Source/amesos_camd_aat.c
  src/SuiteSparse/CAMD/Source/amesos_camd_control.c
  src/SuiteSparse/CAMD/Source/amesos_camd_defaults.c
  src/SuiteSparse/CAMD/Source/amesos_camd_dump.c
  src/SuiteSparse/CAMD/Source/amesos_camd_info.c
  src/SuiteSparse/CAMD/Source/amesos_camd_order.c
  src/SuiteSparse/CAMD/Source/amesos_camd_preprocess.c
  src/SuiteSparse/CAMD/Source/amesos_camd_valid.c
  src/SuiteSparse/CCOLAMD/Source/amesos_ccolamd.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_amd.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_analyze.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_colamd.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_etree.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_factorize.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_postorder.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_rcond.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_resymbol.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_rowcolcounts.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_rowfac.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_solve.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_spsolve.c
  src/SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_rowfac.c
  src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_aat.c
  src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_add.c
  src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_band.c
  src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_change_factor.c
  src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_copy.c
  src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_dense.c
  src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_factor.c
  src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_camd.c
  src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_ccolamd.c
  src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_csymamd.c
  src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_metis.c
  src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_nesdis.c
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
