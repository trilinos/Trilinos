
# FIXME (mfh 09 Jan 2015) Disable Tpetra::HybridPlatform tests for
# now, until I can figure out how to make them build again, given
# recent changes.  HybridPlatform is not a priority for us at the
# moment.
#ADD_SUBDIRECTORY(HybridPlatform)

IF (Tpetra_ENABLE_RTI)
  ADD_SUBDIRECTORY(MultiPrec)
  ADD_SUBDIRECTORY(RTInterface)
ENDIF()


