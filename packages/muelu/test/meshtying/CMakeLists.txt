INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Amesos2
  )

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND
    ${PACKAGE_NAME}_ENABLE_Ifpack2 AND
    ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MeshTying_cpdata
    SOURCE_FILES Q_mm.txt G_mm.txt GT_mm.txt C_mm.txt f_mm.txt Lagr2Dof.txt
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MeshTyingBlocked_SimpleSmoother
    SOURCES MeshTyingBlocked_SimpleSmoother.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MeshTying_cp
    SOURCE_FILES myXML.xml
    )

ENDIF()
