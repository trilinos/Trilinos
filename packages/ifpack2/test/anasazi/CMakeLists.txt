INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../belos/hybrid_platform)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Ifpack2_HybridPlatformAnasaziTest_Files1
  SOURCE_DIR ${Tpetra_MACHINE_XML_FILE_DIR}
  SOURCE_FILES gpuonly.xml hybrid.xml mpionly.xml threaded.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Ifpack2_HybridPlatformAnasaziTest_Files2
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  HybridPlatformAnasazi
  SOURCES
    HybridPlatformAnasazi
  # ARGS 
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
)

INSTALL(TARGETS Ifpack2_HybridPlatformAnasazi
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
