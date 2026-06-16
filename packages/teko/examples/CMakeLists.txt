ADD_SUBDIRECTORY(BuildPreconditioner)
ADD_SUBDIRECTORY(AddMultiplyPrecs)

ADD_SUBDIRECTORY(StridedSolve)

IF (ML_ENABLE_TekoSmoothers)
   ADD_SUBDIRECTORY(ML-Teko-Coupling)
ENDIF()
