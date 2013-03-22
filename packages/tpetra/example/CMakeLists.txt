INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORY(PowerMethod)
ADD_SUBDIRECTORY(HybridPlatform)

IF (KokkosClassic_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(LocalMatOps)
ENDIF()

IF (Tpetra_ENABLE_RTI)
  ADD_SUBDIRECTORY(RTInterface)
ENDIF()

ADD_SUBDIRECTORY(MultiPrec)
