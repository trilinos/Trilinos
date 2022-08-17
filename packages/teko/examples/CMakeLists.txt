ADD_SUBDIRECTORY(StridedSolve)
ADD_SUBDIRECTORY(AddMultiplyPrecs)
ADD_SUBDIRECTORY(BuildPreconditioner)
ADD_SUBDIRECTORY(BasicExamples)

IF (ML_ENABLE_TekoSmoothers)
   ADD_SUBDIRECTORY(ML-Teko-Coupling)
ENDIF()

# Consider 'Teko_ENABLE_Epetra' only if '${PROJECT_NAME}_ENABLE_Epetra=ON'
# Use global_set, so TEKO_HAVE_EPETRA is visible in other subdirs
IF(${PROJECT_NAME}_ENABLE_Epetra AND DEFINED Teko_ENABLE_Epetra)
   global_set(TEKO_HAVE_EPETRA ${Teko_ENABLE_Epetra})
ELSE()
   global_set(TEKO_HAVE_EPETRA ${${PROJECT_NAME}_ENABLE_Epetra})
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(example_copyfiles
  SOURCE_FILES nslhs_test.mm
     nsrhs_test.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
  EXEDEPS StridedSolve_Example
  )

IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(example_copyfiles_from_tests
    SOURCE_FILES nsjac.mm
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../tests/data"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
    EXEDEPS StridedSolve_Example
    )
ENDIF()
