ADD_SUBDIRECTORY(Factory)

IF(${PACKAGE_NAME}_ENABLE_Experimental)
  ADD_SUBDIRECTORY(ProjectedLeastSquaresSolver)
ENDIF()

ADD_SUBDIRECTORY(RealSolverManager)

