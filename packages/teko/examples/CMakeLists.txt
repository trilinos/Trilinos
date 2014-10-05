ADD_SUBDIRECTORY(StridedSolve)
ADD_SUBDIRECTORY(AddMultiplyPrecs)
ADD_SUBDIRECTORY(BuildPreconditioner)
ADD_SUBDIRECTORY(BasicExamples)

IF (ML_ENABLE_TekoSmoothers)
   ADD_SUBDIRECTORY(ML-Teko-Coupling)
ENDIF()

INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_COPY_FILES_TO_BINARY_DIR(example_copyfiles
  SOURCE_FILES nsjac_test.mm
     nslhs_test.mm
     nsrhs_test.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
  EXEDEPS StridedSolve_Example
  )
