INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# IF (Anasazi_ENABLE_Tpetra)
# ENDIF()

# IF (Belos_ENABLE_Tpetra)
# ENDIF()

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/../kokkos/LinAlg)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  LocalMatOpExample
  SOURCES 
    LocalMatOpExample
    LocalMatOpExplicitInstant
  # ARGS 
  COMM mpi serial
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
)
