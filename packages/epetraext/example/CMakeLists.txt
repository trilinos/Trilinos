# Note that the HDF5_IO.cpp example was added as a TRIBITS executable since the API changed.
# Other examples in EpetraExt should also be incorporated into the regular example build
# system but it was decided to just do this one for now. I will add an issue for the others.

# Here just checking for ${PACKAGE_NAME}_BUILD_INOUT but below is the code copied from 
# the test folder - as the other examples are developed we may want something more like that.
IF (${PACKAGE_NAME}_BUILD_INOUT)
  ADD_SUBDIRECTORY(inout)
ENDIF()



# For reference - this is how things are setup in the test folder
# Possibly we want the above more like this or this can just be deleted.

#IF (${PACKAGE_NAME}_BUILD_INOUT AND ${PACKAGE_NAME}_ENABLE_Triutils)
#  IF (NOT Trilinos_NO_32BIT_GLOBAL_INDICES)
#    ADD_SUBDIRECTORY(inout)
#  ENDIF()
#
#  IF (NOT Trilinos_NO_64BIT_GLOBAL_INDICES)
#    ADD_SUBDIRECTORY(inout_LL)
#  ENDIF()
#ENDIF()
