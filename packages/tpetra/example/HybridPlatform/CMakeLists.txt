INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  HybridPlatformTest
  SOURCES
    HybridPlatformTest
  # ARGS 
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files1
  SOURCE_DIR ${Tpetra_SOURCE_DIR}/example/HybridPlatform
  SOURCE_FILES gpuonly.xml hybrid.xml mpionly.xml threaded.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files2
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)

# IF (Anasazi_ENABLE_Tpetra AND ${PROJECT_NAME}_ENABLE_Ifpack2)
# 
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ifpack2/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../ifpack2/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../anasazi/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../anasazi/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../anasazi/tpetra/src)
# 
#   TRIBITS_ADD_EXECUTABLE_AND_TEST(
#     HybridPlatformAnasazi
#     SOURCES
#       HybridPlatformAnasazi
#     # ARGS 
#     COMM mpi # serial
#     STANDARD_PASS_OUTPUT
#     DEPLIBS tpetra ifpack2 anasazi
#   )
# 
#   INSTALL(TARGETS Tpetra_HybridPlatformAnasazi
#           RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
#            
# ENDIF()
