INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/utils)

SET(UTILS_SOURCES 
  ${PACKAGE_SOURCE_DIR}/test/utils/Sacado_Random.hpp
  ${PACKAGE_SOURCE_DIR}/test/utils/Sacado_RandomImp.hpp
)

IF (Sacado_ENABLE_Teuchos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FadSerializationTests
    SOURCES Fad_SerializationTests.cpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TaySerializationTests
    SOURCES Tay_SerializationTests.cpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FadCommTests
    SOURCES Fad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ELRFadCommTests
    SOURCES ELRFad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CacheFadCommTests
    SOURCES CacheFad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ELRCacheFadCommTests
    SOURCES ELRCacheFad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TayCommTests
    SOURCES Tay_CommTests.cpp Tay_CommTests.hpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )
ENDIF()
