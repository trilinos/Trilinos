INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  geometric
  probing
  #utest
  #rebalance_1d_default
  )

IF (TPL_ENABLE_MPI)
  #ADD_SUBDIRECTORY(zoltan)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
  ADD_SUBDIRECTORY(mtx)
ENDIF()
