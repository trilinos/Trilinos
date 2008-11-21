INCLUDE(Add_Subdirectories)

ADD_SUBDIRECTORIES(
  geometric
  utest
  rebalance_1d_default
  )

IF (Trilinos_ENABLE_MPI)
  #ADD_SUBDIRECTORY(zoltan)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
  #ADD_SUBDIRECTORY(mtx)
ENDIF()
