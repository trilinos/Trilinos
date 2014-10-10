
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
    tutorial_laplace2d
    SOURCES laplace2d.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(tutorial_laplace2d_cp
    SOURCE_FILES hands-on.py MueLu_Agg2VTK.py
    tmpl/MueLu_Agg2VTK.py_TMPL tmpl/muelu.xml_TMPL tmpl/muelu_easy.xml_TMPL
    laplace2d.cpp
    xml/s1_easy.xml xml/s1_easy_3levels_smoothed.xml xml/s1_easy_3levels_unsmoothed.xml
    xml/s1_easy_10levels.xml xml/s1_easy_exercise.xml xml/s1_easy_jacobi.xml xml/s1_easy_jacobi2.xml
    xml/s2_adv_a.xml xml/s2_adv_b.xml xml/s2_adv_c.xml
    xml/s2a.xml xml/s2b.xml
    xml/s4a.xml xml/s4b.xml xml/s4c.xml
    xml/s5a.xml xml/s5b.xml
    xml/s6_dep.xml xml/s6_export.xml
    xml/s2_easy.xml xml/s3a_easy.xml xml/s3b_easy.xml xml/s3c_easy.xml
    xml/n1_easy.xml xml/n1_easy_pg.xml
    xml/n2_easy_agg.xml xml/n2_easy_agg2.xml xml/n2_easy_agg3.xml xml/n2_easy_export.xml
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tutorial_recirc2d
    SOURCES recirc2d.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(tutorial_recirc2d_cp
    SOURCE_FILES xml/s3a.xml xml/s3b.xml xml/s3c.xml xml/s3b1.xml xml/s3b2.xml xml/s3b3.xml
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tutorial_recirc2d_api
    SOURCES recirc2d_api.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(userguide_srcfiles_cp
    SOURCE_FILES ScalingTest.cpp ScalingTestParamList.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Challenge_XML
    SOURCES Challenge.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(userguide_challenge_cp
    SOURCE_FILES challenges/condif2d_A.txt challenges/condif2d_b.txt challenges/condif2d_ns.txt challenges/condif2d_reference.xml
    challenges/stru2d_A.txt challenges/stru2d_b.txt challenges/stru2d_ns.txt challenges/stru2d_reference.xml
    )
ENDIF()
