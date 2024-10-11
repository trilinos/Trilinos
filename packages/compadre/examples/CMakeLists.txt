include_tribits_build()

###############
#
#
# Template Functions
#
#
###############

function(add_exe_w_compadre EXE_NAME CPP_NAME)
  add_executable(${EXE_NAME} ${CPP_NAME})
  target_link_libraries(${EXE_NAME} PRIVATE compadre)
  bob_export_target(${EXE_NAME})
endfunction(add_exe_w_compadre)

###############
#
#
# GMLS Examples / Tests
#
#
###############
if (Compadre_EXAMPLES)

  add_exe_w_compadre(GMLS_Host_Test GMLS_Host.cpp)
  add_exe_w_compadre(GMLS_Device_Test GMLS_Device.cpp)
  add_exe_w_compadre(GMLS_NeumannGradScalar_Test GMLS_NeumannGradScalar.cpp)
  add_exe_w_compadre(GMLS_Vector_Test GMLS_Vector.cpp)
  add_exe_w_compadre(GMLS_Divergence_Test GMLS_DivergenceFree.cpp)
  add_exe_w_compadre(GMLS_SmallBatchReuse_Device_Test GMLS_SmallBatchReuse_Device.cpp)
  add_exe_w_compadre(GMLS_Manifold_Test GMLS_Manifold.cpp)
  add_exe_w_compadre(GMLS_Staggered GMLS_Staggered.cpp)
  add_exe_w_compadre(GMLS_Staggered_Manifold_Test GMLS_Staggered_Manifold.cpp)
  add_exe_w_compadre(GMLS_MultiSite_Test GMLS_Multiple_Evaluation_Sites.cpp)
  add_exe_w_compadre(GMLS_Manifold_MultiSite_Test GMLS_Manifold_Multiple_Evaluation_Sites.cpp)
  add_exe_w_compadre(TestUtility UtilityTest.cpp)
  add_exe_w_compadre(NeighborSearchTest NeighborSearchTest.cpp)

  if (Compadre_TESTS)
    if (Trilinos_PREFIX AND Trilinos_VERSION VERSION_GREATER_EQUAL 13.5)
      IF (GtestID GREATER -1 )
        add_exe_w_compadre(UnitTests Compadre_UnitTests.cpp)
        TARGET_LINK_LIBRARIES(UnitTests PRIVATE Gtest::all_libs)
        ADD_TEST(CompadreGtests UnitTests)
      ELSE()
        ADD_TEST(NAME CompadreGtests COMMAND ${CMAKE_COMMAND} -E echo "MISSING_GTEST")
      ENDIF()
    else()
      add_exe_w_compadre(UnitTests Compadre_UnitTests.cpp)
      TARGET_LINK_LIBRARIES(UnitTests PRIVATE gtest)
      ADD_TEST(CompadreGtests UnitTests)
    endif()
    SET_TESTS_PROPERTIES(CompadreGtests PROPERTIES SKIP_REGULAR_EXPRESSION "MISSING_GTEST" LABELS "Gtest;UnitTest;unittest;Unit;unit")

    # Host views tests for GMLS
    ADD_TEST(NAME GMLS_Host_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Host_Test "--p" "4" "--nt" "200" "--d" "3" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Host_Dim3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 60)
    
    ADD_TEST(NAME GMLS_Host_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Host_Test "--p" "4" "--nt" "200" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Host_Dim2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 20)
    
    ADD_TEST(NAME GMLS_Host_Dim1_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Host_Test "--p" "4" "--nt" "200" "--d" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Host_Dim1_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 20)

    # Device views tests for GMLS
    ADD_TEST(NAME GMLS_Device_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "--p" "4" "--nt" "200" "--d" "3" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Device_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "--p" "4" "--nt" "200" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;batch" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Device_Dim1_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "--p" "4" "--nt" "200" "--d" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim1_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

    # Device views tests for GMLS - LU solver
    ADD_TEST(NAME GMLS_Device_Dim3_LU COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "--p" "4" "--nt" "200" "--d" "3" "--solver" "LU" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim3_LU PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

    # Device views tests with Neumann BC for GMLS - LU solver
    ADD_TEST(NAME GMLS_NeumannGradScalar_Dim3_LU COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_NeumannGradScalar_Test "--p" "3" "--nt" "200" "--d" "3" "--solver" "LU" "--constraint" "NEUMANN_GRAD_SCALAR" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_NeumannGradScalar_Dim3_LU PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

    # Device views tests with Neumann BC for STAGGERED GMLS
    ADD_TEST(NAME GMLS_StaggeredNeumannGradScalar_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered "--p" "3" "--nt" "200" "--d" "3" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_StaggeredNeumannGradScalar_Dim3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;staggered" TIMEOUT 20)

    # Device views tests for GMLS (vector basis)
    ADD_TEST(NAME GMLS_Vector_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "--p" "3" "--nt" "20" "--d" "3" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;vector" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Vector_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "--p" "3" "--nt" "20" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;vector" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Vector_Dim1_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "--p" "3" "--nt" "20" "--d" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim1_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;vector" TIMEOUT 10)

    # Device views tests for GMLS (vector basis) with LU
    ADD_TEST(NAME GMLS_Vector_Dim3_LU COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "--p" "3" "--nt" "200" "--d" "3" "--solver" "LU" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim3_LU PROPERTIES LABELS "IntegrationTest;integration;kokkos;vector" TIMEOUT 10)

    # Device views tests for small batch GMLS, reusing GMLS class object
    ADD_TEST(NAME GMLS_SmallBatchReuse_Device_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_SmallBatchReuse_Device_Test "--p" "4" "--nt" "200" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_SmallBatchReuse_Device_Dim2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 20)
    
    ADD_TEST(NAME GMLS_SmallBatchReuse_Device_Dim1_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_SmallBatchReuse_Device_Test "--p" "4" "--nt" "200" "--d" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_SmallBatchReuse_Device_Dim1_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 20)

    # Multisite test for GMLS
    ADD_TEST(NAME GMLS_MultiSite_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_MultiSite_Test "--p" "4" "--nt" "200" "--d" "3" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_MultiSite_Dim3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

    # Staggered scheme test for GMLS on non-manifold
    # Note: Using even polynomial order may cause this test to fail
    ADD_TEST(NAME GMLS_Staggered_Dim3_QR COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Staggered "--p" "3" "--nt" "100" "--d" "3" "${KOKKOS_THREADS_ARG}=4")
    SET_TESTS_PROPERTIES(GMLS_Staggered_Dim3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;staggered" TIMEOUT 60)

    # Staggered scheme test for GMLS on non-manifold
    # Note: Using even polynomial order may cause this test to fail
    ADD_TEST(NAME GMLS_Staggered_Dim2_QR COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Staggered "--p" "3" "--nt" "200" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Staggered_Dim2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;staggered" TIMEOUT 20)

    if (NOT Compadre_DEBUG)
      # This test is too slow in DEBUG (3x longer than all other tests combined)
      # Multisite manifold test for GMLS
      #ADD_TEST(NAME GMLS_Manifold_MultiSite COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_MultiSite_Test "3" "4" "${KOKKOS_THREADS_ARG}=2")
      #SET_TESTS_PROPERTIES(GMLS_Manifold_MultiSite PROPERTIES LABELS "IntegrationTest;integration;kokkos;multi;manifold" TIMEOUT 10)
      if (PYTHON_EXECUTABLE)
        # Python driven test of a C++ executable (Python changes command line arguments given to executable)
        CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py" @ONLY)
        ADD_TEST(NAME GMLS_Manifold_MultiSite_Refinement_Study COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py" "--porder=3" "--grids=3")
        SET_TESTS_PROPERTIES(GMLS_Manifold_MultiSite_Refinement_Study
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold" 
            TIMEOUT 60
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_MultiSite_Test>
        )
      endif() # PYTHON_EXECUTABLE

      # Divergence-free basis test for GMLS on non-manifold
      ADD_TEST(NAME GMLS_DivergenceFree_Dim3_P3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "--p" "3" "--nt" "200" "--d" "3" "${KOKKOS_THREADS_ARG}=2")
      SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim3_P3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr;batched" TIMEOUT 60)
      #ADD_TEST(NAME GMLS_DivergenceFree_Dim3_P2_QR COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "2" "200" "3" "0" "0" "0" "${KOKKOS_THREADS_ARG}=2")
      #SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim3_P2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr" TIMEOUT 60)

      ADD_TEST(NAME GMLS_DivergenceFree_Dim2_P3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "--p" "3" "--nt" "200" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
      SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim2_P3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr;batched" TIMEOUT 60)
      #ADD_TEST(NAME GMLS_DivergenceFree_Dim2_P2_QR COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "2" "200" "2" "0" "0" "0" "${KOKKOS_THREADS_ARG}=2")
      #SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim2_P2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr" TIMEOUT 60)
    endif()

    #    # Manifold tests for GMLS
    #    ADD_TEST(NAME GMLS_Manifold_Dim3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Test "4" "200" "3" "1" "0" "0" "${KOKKOS_THREADS_ARG}=2")
    #    SET_TESTS_PROPERTIES(GMLS_Manifold_Dim3 PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

    if (PYTHON_EXECUTABLE)
      # Python driven test of a C++ executable (Python changes command line arguments given to executable) - calling QR solver
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" @ONLY)
      ADD_TEST(NAME GMLS_Manifold_Refinement_Study_QR COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=3" "--grids=4" "--solver-type=QR")
      SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_QR
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold" 
            TIMEOUT 20
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_Test>
      )

      # Python driven test of a C++ executable (Python changes command line arguments given to executable) - calling LU solver
      ADD_TEST(NAME GMLS_Manifold_Refinement_Study_LU COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=3" "--grids=4" "--solver-type=LU")
      SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_LU
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold" 
            TIMEOUT 20
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_Test>
      )

      # Python driven test of a C++ executable (Python changes command line arguments given to executable)
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Staggered_Manifold.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py" @ONLY)
      ADD_TEST(NAME GMLS_Staggered_Manifold_Refinement_Study COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py" "--porder=3" "--grids=4")
      SET_TESTS_PROPERTIES(GMLS_Staggered_Manifold_Refinement_Study
        PROPERTIES
          LABELS "ConvergenceTest;convergence;manifold;staggered" 
          TIMEOUT 20
          REQUIRED_FILES
            $<TARGET_FILE:GMLS_Staggered_Manifold_Test>
      )
    endif() # PYTHON_EXECUTABLE

    # Utility test - Filter By ID
    ADD_TEST(NAME Test_Utilities COMMAND ${CMAKE_CURRENT_BINARY_DIR}/TestUtility "200" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(Test_Utilities PROPERTIES LABELS "UtilityTest;utility;kokkos" TIMEOUT 5)

    # Neighbor radius search - 2D
    ADD_TEST(NAME NeighborRadiusSearch2DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "200" "6.5" "0" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch2DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborRadiusSearch2DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "300" "4.5" "0" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch2DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborRadiusSearch2DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "400" "1.8" "${KOKKOS_THREADS_ARG}=2")
    #SET_TESTS_PROPERTIES(NeighborRadiusSearch2DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # Neighbor radius search - 3D
    ADD_TEST(NAME NeighborRadiusSearch3DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "100" "4.5" "0" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch3DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborRadiusSearch3DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "150" "3.5" "0" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch3DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborRadiusSearch3DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "160" "1.8" "${KOKKOS_THREADS_ARG}=2")
    #SET_TESTS_PROPERTIES(NeighborRadiusSearch3DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # Neighbor KNN search - 2D
    ADD_TEST(NAME NeighborKNNSearch2DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "200" "6.5" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch2DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborKNNSearch2DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "300" "4.5" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch2DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborKNNSearch2DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "400" "1.8" "1" "${KOKKOS_THREADS_ARG}=2")
    #SET_TESTS_PROPERTIES(NeighborKNNSearch2DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # Neighbor KNN search - 3D
    ADD_TEST(NAME NeighborKNNSearch3DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "100" "4.5" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborKNNSearch3DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "150" "3.5" "1" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborKNNSearch3DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "160" "1.8" "1" "${KOKKOS_THREADS_ARG}=2")
    #SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
    if (Compadre_USE_PYTHON AND PYTHON_EXECUTABLE)
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/test_pycompadre.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/test_pycompadre_manifold.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/test_tangent_planes_manifold.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/kokkos_test_case.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/Python_3D_Convergence.py.in" "${CMAKE_CURRENT_BINARY_DIR}/Python_3D_Convergence.py" @ONLY)
      ADD_TEST(NAME GMLS_Python_Convergence_Test_3d_Point_Reconstruction COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Python_3D_Convergence.py" "${KOKKOS_THREADS_ARG}=4")
      SET_TESTS_PROPERTIES(GMLS_Python_Convergence_Test_3d_Point_Reconstruction PROPERTIES LABELS "IntegrationTest;integration;python;kokkos" TIMEOUT 10)

      if (Compadre_USE_MATLAB)
        ADD_TEST(NAME GMLS_Matlab_Python_Interface_1d_Point_Reconstruction COMMAND "matlab" "-nodisplay" "-nojvm" "-nosplash" "-nodesktop" "-r \"try, run('${SWIG_PREFIX}/Matlab_1D_Using_Python_Interface.m'), catch, exit(1), end, exit(0);\"")
        SET_TESTS_PROPERTIES(GMLS_Matlab_Python_Interface_1d_Point_Reconstruction PROPERTIES LABELS "IntegrationTest;integration;python;kokkos;matlab" TIMEOUT 10)
      endif()

      ADD_TEST(NAME PyCOMPADRE_Tests COMMAND "${PYTHON_EXECUTABLE}" "-m" "unittest" "discover" "-s" "${PYTHON_LIBRARY_PREFIX}/pycompadre/examples")
      SET_TESTS_PROPERTIES(PyCOMPADRE_Tests PROPERTIES ENVIRONMENT PYTHONPATH=${PYTHON_LIBRARY_PREFIX})
      SET_TESTS_PROPERTIES(PyCOMPADRE_Tests PROPERTIES LABELS "pycompadre" TIMEOUT 30)

    endif()

  endif (Compadre_TESTS)
endif (Compadre_EXAMPLES)

bob_end_subdir()


