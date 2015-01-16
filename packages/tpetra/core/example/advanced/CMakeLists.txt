
IF (${PACKAGE_NAME}_ENABLE_RTI)
  # FIXME (mfh 16 Jan 2015) The MultiPrec examples do not link
  # correctly by default, because Tpetra still doesn't do ETI for
  # CrsMatrix<double>::convert<float>().  This should be considered a
  # bug, since Scalar = float and double are both enabled for this
  # build.  Thus, I have not yet enabled the MultiPrec examples.
  # They do build correctly; they just don't link correctly.
  #ADD_SUBDIRECTORY(MultiPrec)
  
  ADD_SUBDIRECTORY(RTInterface)
ENDIF()


