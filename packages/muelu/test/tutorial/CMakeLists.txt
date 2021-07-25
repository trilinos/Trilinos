INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Isorropia
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra)

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyXmlFiles
    SOURCE_FILES
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
      s1_easy_jacobi2.xml
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

#C++-Files

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_laplace2d
    SOURCES laplace2d.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d"
    ARGS "--xml=s1_easy.xml --xml=s1_easy_3levels_smoothed.xml --xml=s1_easy_3levels_unsmoothed.xml --xml=s1_easy_exercise.xml --xml=s1_easy_jacobi.xml --xml=s1_easy_jacobi2.xml --xml=s2a.xml --xml=s2_adv_a.xml --xml=s2_adv_b.xml --xml=s2_adv_c.xml --xml=s2b.xml --xml=s2_easy.xml --xml=s3a.xml --xml=s3a_easy.xml --xml=s3b.xml --xml=s3b1.xml --xml=s3b2.xml --xml=s3b3.xml --xml=s3b_easy.xml --xml=s3c.xml--xml=s3c_easy.xml --xml=s4a.xml --xml=s4av.xml --xml=s4b.xml --xml=s4c.xml --xml=s5a.xml --xml=s5b.xml --xml=s6_dep.xml --xml=s6_export.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  IF (${PACKAGE_NAME}_ENABLE_Isorropia)
    TRIBITS_ADD_TEST(
      Tutorial_laplace2d
      NAME "Test_Tutorial_laplace2d_repartitioning"
      ARGS "--xml=s1_easy.xml --xml=s1_easy_3levels_smoothed.xml --xml=s1_easy_3levels_unsmoothed.xml --xml=s1_easy_exercise.xml --xml=s1_easy_jacobi.xml --xml=s1_easy_jacobi2.xml --xml=s2a.xml --xml=s2_adv_a.xml --xml=s2_adv_b.xml --xml=s2_adv_c.xml --xml=s2b.xml --xml=s2_easy.xml --xml=s3a.xml --xml=s3a_easy.xml --xml=s3b.xml --xml=s3b1.xml --xml=s3b2.xml --xml=s3b3.xml --xml=s3b_easy.xml --xml=s3c.xml--xml=s3c_easy.xml --xml=s4a.xml --xml=s4av.xml --xml=s4b.xml --xml=s4c.xml --xml=s5a.xml --xml=s5b.xml --xml=s6_dep.xml --xml=s6_export.xml"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )
  ENDIF()

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_recirc2d
    SOURCES recirc2d.cpp
    COMM serial mpi
    )

    TRIBITS_ADD_TEST(
    Tutorial_recirc2d
    NAME "Test_Tutorial_recirc2d"
    ARGS "--xml=n1_easy.xml --xml=n1_easy_pg.xml --xml=n2_easy_agg.xml --xml=n2_easy_agg2.xml --xml=n2_easy_agg3.xml --xml=n2_easy_export.xml"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_MLParameterList
    SOURCES MLParameterList.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Tutorial_MLParameterList
    NAME "Test_Tutorial_MLParameterList"
    ARGS "--xml=ml_ParameterList.xml --xml=ml_RebalancedParameterList.xml"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  #TRIBITS_ADD_EXECUTABLE(
  #  Tutorial_ScalingTestParamList
  #  SOURCES ScalingTestParamList.cpp
  #  COMM serial mpi
  #  )

  #TRIBITS_ADD_TEST(
  #  Tutorial_ScalingTestParamList
  #  NAME "Test_Tutorial_ScalingTestParamList"
  #  ARGS ""
  #  COMM serial mpi
  #  NUM_MPI_PROCS 4
  #  )

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_ScalingTest
    SOURCES ScalingTest.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Tutorial_ScalingTest
    NAME "Test_Tutorial_ScalingTest"
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

ENDIF()
