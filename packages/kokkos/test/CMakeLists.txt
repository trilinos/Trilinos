ADD_SUBDIRECTORY(BaseSparseMultiply)
ADD_SUBDIRECTORY(BaseSparseSolve)
ADD_SUBDIRECTORY(DenseMultiVector)
ADD_SUBDIRECTORY(HbMatrix)
ADD_SUBDIRECTORY(Permutation)

IF(${PACKAGE_NAME}_ENABLE_Experimental)
#
#  It looks like these tests require a tpl of some sort
#  called oski. There appears to be no reference to it 
#  in the configure.ac or Makefile.am though.
#
#  ADD_SUBDIRECTORY(OskiMultiVector)
#  ADD_SUBDIRECTORY(OskiMatrix)
#  ADD_SUBDIRECTORY(OskiSparseMultiply)
ENDIF()
