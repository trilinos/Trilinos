INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_AztecOO
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  ${PACKAGE_NAME}_ENABLE_Galeri
  ${PACKAGE_NAME}_ENABLE_Ifpack
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos       AND
    ${PACKAGE_NAME}_ENABLE_AztecOO      AND
    ${PACKAGE_NAME}_ENABLE_Epetra       AND
    ${PACKAGE_NAME}_ENABLE_EpetraExt    AND
    ${PACKAGE_NAME}_ENABLE_Galeri       AND
    ${PACKAGE_NAME}_ENABLE_Ifpack
   )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    laplace2d
    SOURCES laplace2d.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(laplace2d_cp
    SOURCE_FILES hands-on.sh hands-on.py hands-on.py-expert MueLu_Agg2VTK.py
    tmpl/MueLu_Agg2VTK.py_TMPL tmpl/muelu.xml_TMPL tmpl/muelu_easy.xml_TMPL
    xml/s2a.xml xml/s2b.xml xml/s4a.xml xml/s4b.xml xml/s4c.xml xml/s5a.xml
    xml/s2_easy.xml xml/s3a_easy.xml xml/s3b_easy.xml xml/s3c_easy.xml
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    recirc2d
    SOURCES recirc2d.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(recirc2d_cp
    SOURCE_FILES xml/s3a.xml xml/s3b.xml xml/s3c.xml
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    recirc2d_api
    SOURCES recirc2d_api.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
    )
ENDIF()
