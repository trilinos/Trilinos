INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  HybridPlatformTest
  SOURCES
    HybridPlatformTest
  # ARGS 
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files1
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/tpetra/example/HybridPlatform
  SOURCE_FILES gpuonly.xml hybrid.xml mpionly.xml threaded.xml
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files2
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)

# IF (Belos_ENABLE_Tpetra AND Trilinos_ENABLE_Ifpack2)
# 
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ifpack2/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../ifpack2/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../belos/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../belos/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../belos/tpetra/src)
# 
#   PACKAGE_ADD_EXECUTABLE_AND_TEST(
#     HybridPlatformBelos
#     SOURCES
#       HybridPlatformBelos
#     # ARGS 
#     COMM mpi # serial
#     STANDARD_PASS_OUTPUT
#     DEPLIBS tpetra ifpack2 belos
#   )
# 
#   INSTALL(TARGETS Tpetra_HybridPlatformBelos
#           RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
#            
# ENDIF()
# 
# IF (Anasazi_ENABLE_Tpetra AND Trilinos_ENABLE_Ifpack2)
# 
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ifpack2/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../ifpack2/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../anasazi/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../anasazi/src)
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../anasazi/tpetra/src)
# 
#   PACKAGE_ADD_EXECUTABLE_AND_TEST(
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
