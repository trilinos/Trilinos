
IF (Stokhos_ENABLE_EpetraExt)
  TRIBITS_ADD_EXECUTABLE(
    TestEpetraMatrixFreeApply
    SOURCES TestEpetraMatrixFreeApply.cpp
    COMM serial mpi
    )
ENDIF()

ADD_SUBDIRECTORY(KokkosArraySPMVKernels)
IF (Stokhos_ENABLE_TpetraKernels)
  IF (Stokhos_ENABLE_Sacado)
    ADD_SUBDIRECTORY(MPVectorKernels)
  ENDIF()
  ADD_SUBDIRECTORY(KokkosSpMM)
ENDIF()
ADD_SUBDIRECTORY(CuspSpMM)
ADD_SUBDIRECTORY(CuspAMG)
IF (Stokhos_ENABLE_Sacado AND Stokhos_ENABLE_TpetraKernels)
  ADD_SUBDIRECTORY(PCEMeanMultiply)
ENDIF()
IF (Stokhos_ENABLE_Sacado AND Stokhos_ENABLE_KokkosContainers)
  ADD_SUBDIRECTORY(MPAssembly)
ENDIF()
