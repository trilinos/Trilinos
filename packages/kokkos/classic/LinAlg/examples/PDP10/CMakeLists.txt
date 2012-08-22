INCLUDE(TribitsAddExecutableAndTest)   
INCLUDE(TribitsLibraryMacros)

IF (KokkosClassic_ENABLE_CUSPARSE AND KokkosClassic_ENABLE_ThreadPool) 
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  set(source_files PDP10_TestOps.cu)
  CUDA_COMPILE(CUDA_FILES PDP10_TestOps.cu)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(  PDP10_driver
    SOURCES 
      PDP10_driver
      ${CUDA_FILES}
      ${source_files}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  ) 
ENDIF()
