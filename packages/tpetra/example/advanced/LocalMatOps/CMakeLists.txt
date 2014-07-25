INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF (NOT KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode" AND
     NOT KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosOpenMPWrapperNode" AND
     NOT KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosThreadsWrapperNode")
  TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/../kokkos/LinAlg)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LocalMatOpExample
    SOURCES 
      LocalMatOpExample
      LocalMatOpExplicitInstant
    # ARGS 
    COMM mpi serial
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )
     
ENDIF ()

