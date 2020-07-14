

TRIBITS_ADD_EXECUTABLE(
  Solver_Test
  SOURCES Solver_Test
  )


##### cuSOLVER Tests #####

IF (${PACKAGE_NAME}_ENABLE_CUSPARSE)
IF (${PACKAGE_NAME}_ENABLE_CUSOLVER)
IF (${PACKAGE_NAME}_ENABLE_METIS)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestcuSolverFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES cuSOLVER_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME cuSOLVER_Solver_Test
    ARGS "--xml-params=cuSOLVER_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

ENDIF()
ENDIF()
ENDIF()


##### KLU2 Tests #####

##IF (${PACKAGE_NAME}_ENABLE_KLU2 AND NOT ${PROJECT_NAME}_ENABLE_DEBUG AND NOT Kokkos_ENABLE_DEBUG)
IF (${PACKAGE_NAME}_ENABLE_KLU2)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
    KLU2_UnitTests
#    EXCLUDE_IF_NOT_TRUE ${PROJECT_NAME}_ENABLE_DEBUG
    SOURCES
      KLU2_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestKLU2Files
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES klu2_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME KLU2_Solver_Test
    ARGS "--xml-params=klu2_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    COMM serial mpi
    )

ENDIF()

##### SuperLu Tests #####


IF (${PACKAGE_NAME}_ENABLE_SuperLU)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopySuperLUFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES superlu_test.xml
    EXEDEPS Solver_Test
    )

    
  TRIBITS_ADD_TEST(
    Solver_Test
    NAME SuperLU_Solver_Test
    ARGS "--xml-params=superlu_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )

  IF (KokkosKernels_ENABLE_SUPERNODAL_SPTRSV AND KokkosKernels_ENABLE_TPL_SUPERLU)

    TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestSuperLUTRSVFiles
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      SOURCE_FILES superlu_trsv_test.xml
      EXEDEPS Solver_Test
      )

    TRIBITS_ADD_TEST(
      Solver_Test
      NAME SuperLU_TriangularSolve_Solver_Test
      ARGS "--xml-params=superlu_trsv_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

  ENDIF()

IF (NOT ${PROJECT_NAME}_ENABLE_DEBUG)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Superlu_UnitTests
#    EXCLUDE_IF_NOT_TRUE ${PROJECT_NAME}_ENABLE_DEBUG
    SOURCES
      Superlu_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )
ENDIF()

ENDIF()


##### Umfpack Tests #####

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_UMFPACK)
IF (${PACKAGE_NAME}_ENABLE_UMFPACK)

  IF (HAVE_UMFPACK_VERSION_5)
    SET(umpack_test_name "Umfpack_Solver_Test")
    SET(umfpack_test_xml_file "umfpack_test.xml")
  ELSE()
    MESSAGE(STATUS "Disabling Amesos2 Umfpack Complex Test because Umfpack version is not >= 5. This Umfpack does not support Real-Im as a single array which is currently used in the wrapper.")
    SET(umpack_test_name "Umfpack_Solver_Test_No_Complex")
    SET(umfpack_test_xml_file "umfpack_test_no_complex.xml")
  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyUmfpackFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${umfpack_test_xml_file}
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME ${umpack_test_name}
    ARGS "--xml-params=${umfpack_test_xml_file} --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )

ENDIF()

##### Tacho Tests #####

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ShyLU_NodeTacho)
IF (${PACKAGE_NAME}_ENABLE_ShyLU_NodeTacho)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyTachoFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES tacho_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME Tacho_Solver_Test
    ARGS "--xml-params=tacho_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1 # Need to determine if we want more MPI when not running OpenMP with threads
    COMM serial mpi
    )

ENDIF()

##### CHOLMOD Tests #####

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Cholmod)
IF (${PACKAGE_NAME}_ENABLE_Cholmod)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCholmodFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES cholmod_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME Cholmod_Solver_Test
    ARGS "--xml-params=cholmod_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  IF (KokkosKernels_ENABLE_SUPERNODAL_SPTRSV AND KokkosKernels_ENABLE_TPL_CHOLMOD)

    TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCholmodTRSVFiles
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      SOURCE_FILES cholmod_trsv_test.xml
      EXEDEPS Solver_Test
      )

    TRIBITS_ADD_TEST(
      Solver_Test
      NAME Cholmod_TriangularSolve_Solver_Test
      ARGS "--xml-params=cholmod_trsv_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

  ENDIF()

ENDIF()

##### SuperLU_MT Tests #####

IF (${PACKAGE_NAME}_ENABLE_SuperLUMT)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopySuperLUMTFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES superlumt_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME SuperLU_MT_Solver_Test
    ARGS "--xml-params=superlumt_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )
ENDIF()


##### SuperLU_DIST Tests #####


IF (${PACKAGE_NAME}_ENABLE_SuperLUDist)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopySuperLUDistFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES superludist_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME SuperLU_DIST_Solver_Test
    ARGS "--xml-params=superludist_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM mpi
    )
ENDIF()


##### PARDISO_MKL Tests #####

IF (${PACKAGE_NAME}_ENABLE_PARDISO_MKL)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyPardisoMKLFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES pardiso_mkl_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME Pardiso_MKL_Solver_Test
    ARGS "--xml-params=pardiso_mkl_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )
ENDIF()


##### LAPACK Tests #####

IF (${PACKAGE_NAME}_ENABLE_LAPACK)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyLapackFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES lapack_test.xml
    EXEDEPS Solver_Test
    )

  TRIBITS_ADD_TEST(
    Solver_Test
    NAME LAPACK_Solver_Test
    ARGS "--xml-params=lapack_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )
ENDIF()


####   BASKER Tests ####
## IF (${PACKAGE_NAME}_ENABLE_Basker AND NOT ${PROJECT_NAME}_ENABLE_DEBUG)
IF (${PACKAGE_NAME}_ENABLE_Basker)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Basker_UnitTests
    SOURCES
      Basker_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )

TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyBaskerFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES Basker_test.xml
    EXEDEPS Solver_Test
    )

TRIBITS_ADD_TEST(
    Solver_Test
    NAME Basker_Solver_Test
    ARGS "--xml-params=Basker_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    COMM serial mpi
    )
    
ENDIF()

####   SHYLUBASKER Tests ####
IF (${PACKAGE_NAME}_ENABLE_ShyLU_NodeBasker)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ShyLUBasker_UnitTests
    SOURCES
      ShyLUBasker_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )
    
    ##Copy shylubasker_test.xml, but do not test
    ##depends on not included matrices
    TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyShyLUBaskerFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES shylubasker_test.xml
    EXEDEPS Solver_Test
    )

ENDIF()

##### MUMPS Tests ####
#MUMPS requires int/int type
IF (${PACKAGE_NAME}_ENABLE_MUMPS AND Tpetra_INST_INT_INT)
   TRIBITS_ADD_EXECUTABLE_AND_TEST(
   MUMPS_UnitTests
   SOURCES
      MUMPS_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
   COMM serial mpi
   STANDARD_PASS_OUTPUT
   )
ENDIF()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactory_UnitTests
  SOURCES
  SolverFactory
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

