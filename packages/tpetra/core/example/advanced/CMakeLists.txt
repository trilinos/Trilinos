
IF (${PACKAGE_NAME}_ENABLE_RTI)
  # FIXME (mfh 15 Jan 2015) This stuff is interesting, but in order to
  # work, it wants control of Tpetra's set of explicit template
  # instantiation (ETI) Scalar types.  Better to disable for now to
  # avoid troubles.
  # 
  #ADD_SUBDIRECTORY(MultiPrec)
  
  ADD_SUBDIRECTORY(RTInterface)
ENDIF()


