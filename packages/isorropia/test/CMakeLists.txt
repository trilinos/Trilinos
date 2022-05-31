
ADD_SUBDIRECTORIES(
  geometric
  part_sizes
  #utest
  #rebalance_1d_default
  )

IF (TPL_ENABLE_MPI)
  #ADD_SUBDIRECTORY(zoltan)
ENDIF()

IF ("FALSE") # if (${PACKAGE_NAME}_ENABLE_EpetraExt)  # See #10534
  ADD_SUBDIRECTORY(mtx)
  ADD_SUBDIRECTORY(probing)
ENDIF()

