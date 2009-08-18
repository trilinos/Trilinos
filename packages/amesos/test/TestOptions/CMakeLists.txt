
INCLUDE(PackageAddExecutableAndTest)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
  SET(${PACKAGE_NAME}_TESTOPTION_SOURCES
    TestOptions.cpp
    CreateTridi.cpp
    CreateTridi.h
    TestSuperludist.cpp
    TestScalapack.cpp
    TestScalapack.h
    TestKlu.cpp
    TestKlu.h
    TestSuperludist.h
    NewMatNewMap.cpp
    NewMatNewMap.h
    TestAllClasses.cpp
    TestAllClasses.h
    PartialFactorization.cpp
    PartialFactorization.h
    PerformOneSolveAndTest.cpp
    PerformOneSolveAndTest.h
    TestOtherClasses.cpp
    TestOtherClasses.h
    )
ELSE()
  SET(${PACKAGE_NAME}_TESTOPTION_SOURCES DoNothing.cpp)
ENDIF()


PACKAGE_ADD_EXECUTABLE(
  TestOptions
  SOURCES ${${PACKAGE_NAME}_TESTOPTION_SOURCES}
  COMM serial mpi
  )

# NOTE: This test currently fails in MPI DEBUG mode (bug 4539)
IF ( (NOT TPL_ENABLE_MPI) AND (NOT BUILD_TYPE STREQUAL "DEBUG") )

  PACKAGE_ADD_TEST(
    TestOptions
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()

