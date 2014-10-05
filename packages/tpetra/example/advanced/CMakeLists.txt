INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORY(HybridPlatform)

IF (Tpetra_ENABLE_RTI)
  ADD_SUBDIRECTORY(MultiPrec)
  ADD_SUBDIRECTORY(RTInterface)
ENDIF()


