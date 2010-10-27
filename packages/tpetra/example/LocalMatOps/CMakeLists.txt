INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# IF (Anasazi_ENABLE_Tpetra)
# ENDIF()

# IF (Belos_ENABLE_Tpetra)
# ENDIF()


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
