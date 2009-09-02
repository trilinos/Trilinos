
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

IF(Amesos_ENABLE_Galeri)
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Test_Singular
  SOURCES cxx_main.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_TEST_SINGULAR_COPY
  SOURCE_FILES bcsstm05_ns.rua
  )
ENDIF()
