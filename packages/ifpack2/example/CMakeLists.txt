TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_ADD_EXECUTABLE(
  ex1
  SOURCES ex1.cpp
)

IF (Tpetra_INST_DOUBLE)
TRIBITS_ADD_EXECUTABLE(
  RelaxationWithEquilibration
  SOURCES RelaxationWithEquilibration.cpp
)
ENDIF()


TRIBITS_ADD_EXECUTABLE(
  BlockTriDiagonalSolver
  SOURCES BlockTriDi.cpp
  CATEGORIES BASIC PERFORMANCE
)


TRIBITS_ADD_TEST(
  BlockTriDiagonalSolver
  NAME BlockTriDiExample
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED (
  ${PACKAGE_NAME}_ENABLE_Xpetra
  ${PACKAGE_NAME}_ENABLE_Galeri
)

IF(${PACKAGE_NAME}_ENABLE_Xpetra AND ${PACKAGE_NAME}_ENABLE_Galeri)

  set(blockSize 11)
  set(nx 256)
  set(ny 128)
  set(nz 128)

  set(COMMON_ARGS "--withoutOverlapCommAndComp --matrixType=Laplace3D --blockSize=${blockSize} --nx=${nx} --ny=${ny} --nz=${nz} --tol=0 --numIters=5 --withStackedTimer --numRepeats=5")
  set(BTDS_ARGS "${COMMON_ARGS} --sublinesPerLine=1 --sublinesPerLineSchur=1")
  set(BTDS_SCHUR_ARGS "${COMMON_ARGS} --sublinesPerLine=1 --sublinesPerLineSchur=2")
  set(BTDS_DEFAULT_SCHUR_ARGS "${COMMON_ARGS} --sublinesPerLine=1 --sublinesPerLineSchur=-1")
  set(BLOCK_JACOBI_ARGS "${COMMON_ARGS} --sublinesPerLine=-1")

  set(TEST_ARGS ${BTDS_ARGS} ${BTDS_SCHUR_ARGS} ${BTDS_DEFAULT_SCHUR_ARGS} ${BLOCK_JACOBI_ARGS})
  set(TEST_NAMES "BTDS" "BTDS_Schur" "BTDS_Schur_Default" "BTDS_Jacobi")
  set(ARCHS "CUDA" "HIP")

  list(LENGTH TEST_ARGS len1)
  math(EXPR N_TESTS "${len1} - 1")

  ASSERT_DEFINED(
    Tpetra_INST_CUDA
    Tpetra_INST_HIP
  )

  foreach(ARCH ${ARCHS})
    IF (Tpetra_INST_${ARCH})
      MESSAGE(STATUS "Ifpack2: BlockTriDiagonalSolver ${ARCH} test ENABLED")
      set(COUNTER 2)
      set(MAX_COUNT 8)

      WHILE(COUNTER LESS MAX_COUNT)
        math(EXPR NP "1 << ${COUNTER}")
        foreach(I_TEST RANGE ${N_TESTS})
          list(GET TEST_ARGS ${I_TEST} TEST_ARG)
          list(GET TEST_NAMES ${I_TEST} TEST_NAME)      
          TRIBITS_ADD_TEST(
            BlockTriDiagonalSolver
            NAME "${TEST_NAME}_${ARCH}"
            COMM mpi
            ARGS "${TEST_ARG} --problemName=\"${TEST_NAME} ${ARCH} ${blockSize}x${nx}x${ny}x${nz} ${NP} ranks\""
            NUM_MPI_PROCS ${NP}
            RUN_SERIAL
            CATEGORIES PERFORMANCE
          )

        endforeach()

        math(EXPR COUNTER "${COUNTER} + 1")
      ENDWHILE()
    ENDIF()
  endforeach()
ENDIF()