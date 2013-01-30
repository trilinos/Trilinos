INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  TestEpetraMatrixFreeApply
  SOURCES TestEpetraMatrixFreeApply.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )
 
ADD_SUBDIRECTORY(KokkosArraySPMVKernels)
