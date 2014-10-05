INCLUDE(AddSubdirectories)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt ${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters ${PACKAGE_NAME}_ENABLE_Stratimikos)
IF(${PACKAGE_NAME}_ENABLE_EpetraExt AND 
   ${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters AND
   ${PACKAGE_NAME}_ENABLE_Stratimikos)

  ADD_SUBDIRECTORIES(
    diagonalTransient
    epetra
    timeDiscretizedBackwardEuler
    )

ENDIF()
