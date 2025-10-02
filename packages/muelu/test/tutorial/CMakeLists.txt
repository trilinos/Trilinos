INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri ${PACKAGE_NAME}_ENABLE_Isorropia ${PACKAGE_NAME}_ENABLE_Tpetra)

SET(XML_FILES_ALL
    challenges/Elasticity2D_reference.xml
    challenges/Elasticity2D_parameters.xml
    challenges/Recirc2D_reference.xml
    challenges/Recirc2D_parameters.xml
    blocked_rebalancing.xml
    ml_ParameterList.xml
    ml_RebalancedParameterList.xml
    muelu_ParameterList.xml
    n1_easy.xml
    n1_easy_pg.xml
    n2_easy_agg.xml
    n2_easy_agg2.xml
    n2_easy_agg3.xml
    n2_easy_export.xml
    s1_easy.xml
    s1_easy_3levels_smoothed.xml
    s1_easy_3levels_unsmoothed.xml
    s1_easy_exercise.xml
    s1_easy_jacobi.xml
    s1_easy_gaussseidel.xml
    s2a.xml
    s2_adv_a.xml
    s2_adv_b.xml
    s2_adv_c.xml
    s2b.xml
    s2_easy.xml
    s3a.xml
    s3a_easy.xml
    s3b.xml
    s3b1.xml
    s3b2.xml
    s3b3.xml
    s3b_easy.xml
    s3c.xml
    s3c_easy.xml
    s4a.xml
    s4av.xml
    s4b.xml
    s4c.xml
    s5a.xml
    s5b.xml
    s6_dep.xml
    s6_export.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  CopyTutorialXmlFiles
  SOURCE_FILES
    ${XML_FILES_ALL}
    hands-on.py
    MueLu_Agg2VTK.py
    tmpl/muelu.xml_TMPL
    tmpl/MueLu_Agg2VTK.py_TMPL
    tmpl/muelu_easy.xml_TMPL
)

IF(${PACKAGE_NAME}_ENABLE_Galeri)

  # TRIBITS_ADD_EXECUTABLE(
  #   Tutorial_MLParameterList
  #   SOURCES MLParameterList.cpp
  #   COMM serial mpi
  #   )

  # TRIBITS_ADD_TEST(
  #   Tutorial_MLParameterList
  #   NAME "Test_Tutorial_MLParameterList"
  #   ARGS "--xml=ml_ParameterList.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 4
  #   )

  # TRIBITS_ADD_TEST(
  #   Tutorial_MLParameterList
  #   NAME "Test_Tutorial_MLParameterList_Rebalanced"
  #   ARGS "--xml=ml_RebalancedParameterList.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 4
  #   )

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_cppInterface
    SOURCES Tutorial_cppInterface.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_TEST(
    Tutorial_cppInterface
    NAME "Tutorial_cppInterface"
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_EXECUTABLE(
    TutorialDriver
    SOURCES TutorialDriver.cpp
    COMM serial mpi
  )

  SET(XML_FILES_TEST ${XML_FILES_ALL})
  LIST(
    REMOVE_ITEM
    XML_FILES_TEST
    challenges/Elasticity2D_reference.xml
    challenges/Elasticity2D_parameters.xml
    challenges/Recirc2D_reference.xml
    challenges/Recirc2D_parameters.xml
    blocked_rebalancing.xml
    ml_ParameterList.xml
    ml_RebalancedParameterList.xml
    muelu_ParameterList.xml
    n2_easy_agg.xml
    n2_easy_agg2.xml
    n2_easy_agg3.xml
    n2_easy_export.xml
    s4av.xml
    s4c.xml
    s5a.xml
    s5b.xml
    # xmls never mentioned in the tutorial
    n1_easy_pg.xml
    s2_easy.xml
    s2b.xml
    s3a_easy.xml
    s3b_easy.xml
    s3c.xml
    s3c_easy.xml
  )

  FOREACH(xml_file IN LISTS XML_FILES_TEST)
    # without directory and without extension
    GET_FILENAME_COMPONENT(xml_we "${xml_file}" NAME_WE)

    TRIBITS_ADD_TEST(
      TutorialDriver
      NAME "Test_Tutorial_laplace2d_${xml_we}"
      ARGS "--matrixType=Laplace2D --matrixType=Laplace2D --xml=${xml_file}"
      COMM serial mpi
      NUM_MPI_PROCS 2
    )
  ENDFOREACH()

  # foreach() are split so that tests are ordered by matrix type
  FOREACH(xml_file IN LISTS XML_FILES_TEST)
    GET_FILENAME_COMPONENT(xml_we "${xml_file}" NAME_WE)

    TRIBITS_ADD_TEST(
      TutorialDriver
      NAME "Test_Tutorial_recirc2d_${xml_we}"
      ARGS "--matrixType=Recirc2D --matrixType=Recirc2D --xml=${xml_file} --diffusion=1e-5"
      COMM serial mpi
      NUM_MPI_PROCS 4
    )
  ENDFOREACH()

ENDIF()
