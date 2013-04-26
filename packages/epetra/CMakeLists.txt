
INCLUDE(TribitsPackageMacros)

TRIBITS_PACKAGE(Epetra)

SET(${PACKAGE_NAME}_VERSION "3.7d")

# 
#   Exception added to deal with fixed form Fortran code on XL
#
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "XL|VisualAge")
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed=72 -qxlines")
   endif()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ABC
  HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  "Enable array bounds checking on Epetra objects."
  ${${PROJECT_NAME}_ENABLE_DEBUG} )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_FORMAT_IO
  HAVE_FORMAT_IO
  "Enable format IO capability."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_WARNING_MESSAGES
  HAVE_WARNING_MESSAGES
  "Enable building with warning messages."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_FATAL_MESSAGES
  HAVE_FATAL_MESSAGES
  "Enable printing of fatal error messages."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_THREADS
  HAVE_THREADS
  "Enable threading support."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Fortran
  HAVE_FORTRAN_SUPPORT
  "Enable Fortran support in Epetra."
  ${${PROJECT_NAME}_ENABLE_Fortran} )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_DISABLE_READY_SEND_IN_DO_POSTS
  EPETRA_NO_READY_SEND_IN_DO_POSTS
  "Replace (barrier+MPI_Rsend) with MPI_Send in Epetra_MpiDistributor::DoPosts (EXPERIMENTAL)."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(Trilinos_NO_32BIT_GLOBAL_INDICES
  EPETRA_NO_32BIT_GLOBAL_INDICES
  "If option is enabled, 32 bit global indices are disabled.  One must use long long global indices only for maps."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(Trilinos_NO_64BIT_GLOBAL_INDICES
  EPETRA_NO_64BIT_GLOBAL_INDICES
  "If option is enabled, 64 bit global indices are disabled.  One must use int global indices only for maps."
  OFF )

IF(TPL_ENABLE_MKL)
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MKL_SPARSE
    Epetra_ENABLE_MKL_SPARSE
    "If option is enabled, MKL functions will be used for sparse operations in Epetra."
    OFF )

  IF(Epetra_ENABLE_MKL_SPARSE)
    INCLUDE_DIRECTORIES(${TPL_MKL_INCLUDE_DIRS})

    TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_DISABLE_MKL_SPARSE_MM
      Epetra_DISABLE_MKL_SPARSE_MM
      "If option is turned on, MKL matrix-matrix functions will not be used for sparse operations in Epetra.  This does not affect the usage of MKL matrix-vector functions."
      OFF )

  ENDIF()

ELSE()
  IF(Epetra_ENABLE_MKL_SPARSE)
    MESSAGE(FATAL_ERROR "Epetra_ENABLE_MKL_SPARSE requires that MKL TPL be enabled.  Either enable the TPL or do not enable Epetra_ENABLE_MKL_SPARSE.")
  ENDIF()
ENDIF()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  doc
  example/C_wrappers
  example/Fortran
  example/ImportExport
  example/InverseIteration
  example/MapColoring
  example/ReducedLinearProblem
  example/petra_howle
  example/petra_nonlinear
  example/petra_transpose
  src/Epetra_FastCrsMatrix.cpp
  src/Epetra_FastCrsMatrix.h
  src/Epetra_InvOperator.cpp
  src/Epetra_LinearProblemRedistor.cpp
  src/Epetra_LinearProblemRedistor.h
  src/Epetra_MpiSmpComm.*
  src/stamp-h.in
  src/xxemacs
  test/BasicPerfTest/runSummary
  test/Comm/simple_mpi.cpp
  test/Comm/threaded_Makefile
  test/Comm/threaded_main.cpp
  test/EpetraBenchmarkTest
  test/LinearProblemRedistor
  test/Makefile.template
  test/Map/c_main.c
  test/MultiVector/Makefile.purify
  test/OSKI
  test/VbrMatrix/Suppressions.in
  test/Vector/Makefile.purify
  test/testAll.*
  )

TRIBITS_PACKAGE_POSTPROCESS()


