
INCLUDE(PackageAddExecutableAndTest)

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   BasicPerfTest
#   SOURCES
#     BasicPerfTest
#     ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
#   # ARGS 
#   COMM serial # mpi
#   STANDARD_PASS_OUTPUT
#   )


IF (Tpetra_ENABLE_Triutils)
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  HybridPlatformTest
  SOURCES
    HybridPlatformTest
  # ARGS 
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
  )
ENDIF()
