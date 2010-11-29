INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  geometric
  part_sizes
  #probing
  #utest
  #rebalance_1d_default
  )

IF (TPL_ENABLE_MPI)
  #ADD_SUBDIRECTORY(zoltan)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
  ADD_SUBDIRECTORY(mtx)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_TestLevelSolver)
  ADD_SUBDIRECTORY(levelsolver)
ENDIF()
