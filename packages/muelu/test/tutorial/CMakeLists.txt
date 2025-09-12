INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Galeri
  ${PACKAGE_NAME}_ENABLE_Isorropia
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )
  
set(XML_FILES_ALL
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

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTutorialXmlFiles
  SOURCE_FILES
    ${XML_FILES_ALL}
    hands-on.py
    MueLu_Agg2VTK.py
    tmpl/MueLu_Agg2VTK.py_TMPL
    tmpl/muelu.xml_TMPL
    tmpl/muelu_easy.xml_TMPL
  )

IF (${PACKAGE_NAME}_ENABLE_Galeri)

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
  
  TRIBITS_ADD_EXECUTABLE(
    TutorialDriver
    SOURCES TutorialDriver.cpp
    COMM serial mpi
    )
  #TRIBITS_ADD_EXECUTABLE(
  #  TutorialDriver_api
  #  SOURCES TutorialDriver_api.cpp
  #  COMM serial mpi
  #  )
    
  foreach(xml_file IN LISTS XML_FILES_ALL)
    # xml_name = filename without directory and without extension
    get_filename_component(xml_we "${xml_file}" NAME_WE)

    TRIBITS_ADD_TEST(
      TutorialDriver
      NAME "Test_Tutorial_laplace2d_${xml_we}"
      ARGS "--matrixType=Laplace2D --matrixType=Laplace2D --xml=${xml_file}"
      COMM serial mpi
      NUM_MPI_PROCS 2
    )
  endforeach()
  
  # foreach() are split so that tests are ordered by matrix type
  foreach(xml_file IN LISTS XML_FILES_ALL)
    # xml_name = filename without directory and without extension
    get_filename_component(xml_we "${xml_file}" NAME_WE)

    TRIBITS_ADD_TEST(
      TutorialDriver
      NAME "Test_Tutorial_recirc2d_${xml_we}"
      ARGS "--matrixType=Recirc2D --matrixType=Recirc2D --xml=${xml_file} --diffusion=1e-5"
      COMM serial mpi
      NUM_MPI_PROCS 4
    )
  endforeach()

  IF(FALSE)
  ## These tests are incorrectly named. ScalingTestParamList.cpp always uses Laplace2D. I don't know if the plan is to generalize the different problem types into using the same file? But this is mostly an older file anyways (most edits by Tobias Wiesner and other 10 years ago)
    TRIBITS_ADD_EXECUTABLE(
      Tutorial_ScalingTestParamList
      SOURCES ScalingTestParamList.cpp
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      Tutorial_ScalingTestParamList
      NAME "Test_Tutorial_Sc_recirc2d_n1_easy"
      ARGS "--diffusion=1.0e-5 --convection=1.0 --xml=n1_easy.xml"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Tutorial_ScalingTestParamList
      NAME "Test_Tutorial_Sc_recirc2d_n1_easy_pg"
      ARGS "--diffusion=1.0e-5 --convection=1.0 --xml=n1_easy_pg.xml --notimings"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Tutorial_ScalingTestParamList
      NAME "Test_Tutorial_Sc_recirc2d_n2_easy_agg"
      ARGS "--diffusion=1.0e-5 --convection=1.0 --xml=n2_easy_agg.xml --notimings"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Tutorial_ScalingTestParamList
      NAME "Test_Tutorial_Sc_recirc2d_n2_easy_agg2"
      ARGS "--diffusion=1.0e-5 --convection=1.0 --xml=n2_easy_agg2.xml --notimings"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Tutorial_ScalingTestParamList
      NAME "Test_Tutorial_Sc_recirc2d_n2_easy_agg3"
      ARGS "--diffusion=1.0e-5 --convection=1.0 --xml=n2_easy_agg3.xml --notimings"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Tutorial_ScalingTestParamList
      NAME "Test_Tutorial_Sc_recirc2d_n2_easy_export"
      ARGS "--diffusion=1.0e-5 --convection=1.0 --xml=n2_easy_export.xml --notimings"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )
  ENDIF()

ENDIF()

# IF (${PACKAGE_NAME}_ENABLE_Epetra)

#   TRIBITS_ADD_EXECUTABLE(
#     Tutorial_recirc2d
#     SOURCES recirc2d.cpp
#     COMM serial mpi
#     )

#     TRIBITS_ADD_TEST(
#     Tutorial_recirc2d
#     NAME "Test_Tutorial_recirc2d"
#     ARGS "--xml=n1_easy.xml --xml=n1_easy_pg.xml --xml=n2_easy_agg.xml --xml=n2_easy_agg2.xml --xml=n2_easy_agg3.xml --xml=n2_easy_export.xml"
#     COMM serial mpi
#     NUM_MPI_PROCS 4
#     )

#   #TRIBITS_ADD_EXECUTABLE(
#   #  Tutorial_ScalingTestParamList
#   #  SOURCES ScalingTestParamList.cpp
#   #  COMM serial mpi
#   #  )

#   #TRIBITS_ADD_TEST(
#   #  Tutorial_ScalingTestParamList
#   #  NAME "Test_Tutorial_ScalingTestParamList"
#   #  ARGS ""
#   #  COMM serial mpi
#   #  NUM_MPI_PROCS 4
#   #  )

# ENDIF()
