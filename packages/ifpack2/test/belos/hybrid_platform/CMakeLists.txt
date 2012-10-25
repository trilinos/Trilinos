INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Ifpack2Tpetra_HybridPlatformTest_Files1
  SOURCE_DIR ${Tpetra_XML_MACHINE_FILE_DIR}
  SOURCE_FILES gpuonly.xml hybrid.xml mpionly.xml threaded.xml titan.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Ifpack2Tpetra_HybridPlatformTest_Files2
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BelosTpetraHybridPlatformExample
  SOURCES
    BelosTpetraHybridPlatformExample
  # ARGS 
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
)

INSTALL(TARGETS Ifpack2_BelosTpetraHybridPlatformExample
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
