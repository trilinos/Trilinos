INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

IF (KokkosClassic_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(LocalMatOps)
ENDIF()


ADD_SUBDIRECTORY(HybridPlatform)

IF (Tpetra_ENABLE_RTI)
  ADD_SUBDIRECTORY(MultiPrec)
  ADD_SUBDIRECTORY(RTInterface)
ENDIF()


