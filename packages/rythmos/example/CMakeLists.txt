INCLUDE(AddSubdirectories)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt ${PACKAGE_NAME}_ENABLE_Stratimikos)
IF(${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Stratimikos)

  ADD_SUBDIRECTORIES(
    diagonalTransient
    #timeDiscretizedBackwardEuler
    )

ENDIF()
