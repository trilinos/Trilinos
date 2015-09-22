

IF(${PACKAGE_NAME}_ENABLE_Triutils)
  SET(${PACKAGE_NAME}_TEST_BASIC_SOURCES
    SparseDirectTimingVars.cpp
    SparseDirectTimingVars.h
    SparseSolverResult.cpp
    SparseSolverResult.h 
    Amesos_TestDriver.cpp
    Amesos_TestSolver.cpp
    Amesos_TestSolver.h
    Amesos_TestMultiSolver.cpp
    Amesos_TestMrhsSolver.cpp    )
ELSE()
  SET(${PACKAGE_NAME}_TEST_BASIC_SOURCES DoNothing.cpp)
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  amesos_test
  SOURCES ${${PACKAGE_NAME}_TEST_BASIC_SOURCES}
  COMM serial mpi
  NOEXEPREFIX
  )

IF(${PACKAGE_NAME}_ENABLE_Epetraext)
  SET(AMESOS_PARAKLETE_TEST AmesosParaklete.csh)
ELSE()
  SET(AMESOS_PARAKLETE_TEST AmesosParaklete_Sans_EpetraExt.csh)
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_TESTBASIC_COPY1
  DEST_FILES 
    TestBasic.csh
    AmesosSuperludist.csh
    AmesosDscpack.csh
    AmesosUmfpack.csh
    AmesosKlu.csh
    AmesosSuperlu.csh
    AmesosScalapack.csh
    AmesosLapack.csh
    AmesosTaucs.csh
    AmesosPardiso.csh
    ${AMESOS_PARAKLETE_TEST}
    AmesosMumps.csh
    AmesosKluSerial.csh
    AmesosSuperluSerial.csh
    AmesosUmfpackSerial.csh
    Suppressions
  SOURCE_PREFIX "In_"
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_TESTBASIC_COPY2
  SOURCE_FILES 662_bus.rsa 
  DEST_FILES 662_bus_out.rsa
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_TESTBASIC_COPY3
  DEST_FILES 
    ImpcolA.rua
    ImpcolB.rua
    ImpcolC.rua
    ImpcolD.rua
    ImpcolE.rua
    fidapm05.rua
    Dense.triU
    SuperLU.triU
    SuperLU.rua
    bcsstk01.mtx
    bcsstk02.mtx
    bcsstk04.mtx
    nos1.mtx
    nos5.mtx
    FourByFour.mtx
    StructurallySingular.mtx
    NumericallySingular.mtx
    Diagonal.mtx
    MissingADiagonal.mtx
    KheadK.mtx
    KheadSorted.mtx
    Khead.triS
    count_coverage
  SOURCE_PREFIX "In_"
  )

#requires dscpack and klu
#not clear if this is being run in autotools either
#TRIBITS_ADD_TEST(
#  TestBasic.csh
#  NAME TestBasic
#  ARGS "KLU"
#  COMM serial mpi
#  NOEXEPREFIX
#  NOEXESUFFIX
#  )

