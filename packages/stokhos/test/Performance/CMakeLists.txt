INCLUDE(TribitsAddExecutableAndTest)

IF (Stokhos_ENABLE_EpetraExt)
  TRIBITS_ADD_EXECUTABLE(
    TestEpetraMatrixFreeApply
    SOURCES TestEpetraMatrixFreeApply.cpp
    COMM serial mpi
    DEPLIBS stokhos
    )
ENDIF()

ADD_SUBDIRECTORY(KokkosArraySPMVKernels)
ADD_SUBDIRECTORY(MPVectorKernels)
ADD_SUBDIRECTORY(CuspSpMM)
ADD_SUBDIRECTORY(CuspAMG)
