INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_AmesosFactory_Tridiag
  SOURCES example_AmesosFactory_Tridiag.cpp
  COMM serial mpi
  )

IF(Amesos_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_AmesosFactory
    SOURCES example_AmesosFactory.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_AmesosFactory_HB
    SOURCES example_AmesosFactory_HB.cpp
    COMM serial mpi
    )

  IF(Amesos_ENABLE_Triutils)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      compare_solvers
      SOURCES compare_solvers.cpp
      COMM serial mpi
      )
  ENDIF()
ENDIF()

IF(Amesos_ENABLE_MC64)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_MC64
    SOURCES example_MC64.cpp
    COMM serial mpi
    )
ENDIF()

IF(TPL_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    a_trivial_mpi_test
    SOURCES a_trivial_mpi_test.cpp
    COMM serial mpi
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_EXAMPLE_COPY1
  DEST_FILES 662_bus.rsa
  SOURCE_PREFIX "In_"
  )
