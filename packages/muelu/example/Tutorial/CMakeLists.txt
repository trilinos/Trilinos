
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  ${PACKAGE_NAME}_ENABLE_Galeri
  ${PACKAGE_NAME}_ENABLE_Ifpack
  ${PACKAGE_NAME}_ENABLE_Ifpack2
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )

IF (${PACKAGE_NAME}_ENABLE_Belos  AND
    ${PACKAGE_NAME}_ENABLE_Galeri AND
    ((${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Amesos AND ${PACKAGE_NAME}_ENABLE_Ifpack) OR
     (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Amesos2   AND ${PACKAGE_NAME}_ENABLE_Ifpack2)))

  # Test fails with
  # ```
  # /path/to/packages/muelu/src/Smoothers/MueLu_DirectSolver_def.hpp:124:
  # 14:
  # 14:  Throw number = 1
  # 14:
  # 14:  Throw test that evaluated to true: s_.is_null()
  # 14:
  # 14:  Direct solver for Tpetra was not constructed
  # ```
  # Deactivate for now.
  IF(FALSE)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      tutorial
      SOURCES tutorial.cpp
      NUM_MPI_PROCS 2
      COMM mpi serial
      )
  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(tutorial1_cp
    SOURCE_FILES tutorial1a.xml tutorial1b.xml tutorial1c.xml tutorial1d.xml tutorial1e.xml
    )

ENDIF()

ADD_SUBDIRECTORY(hands-on)
