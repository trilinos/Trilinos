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
    add_exe_w_compadre(UnitTests Compadre_UnitTests.cpp)
    TARGET_LINK_LIBRARIES(UnitTests PRIVATE compadre_gtest)
    ADD_TEST(CompadreGtests UnitTests LABELS "Gtest;UnitTest;unittest;Unit;unit")

    # KokkosParser initializing or joining Kokkos; Kokkos initializes once per process, so each mode is a test
    add_exe_w_compadre(KokkosParser_Test KokkosParser_Test.cpp)
    foreach(mode args no_args external)
      ADD_TEST(NAME KokkosParser_${mode} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/KokkosParser_Test "${mode}" "${KOKKOS_THREADS_ARG}=2")
      SET_TESTS_PROPERTIES(KokkosParser_${mode} PROPERTIES LABELS "UnitTest;unittest;Unit;unit;kokkos" TIMEOUT 20)
    endforeach()

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

    # Device views test of CardinalCubicBSpline for GMLS
    ADD_TEST(NAME GMLS_Device_Dim3_QR_CardinalCubicBSpline COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "--p" "4" "--nt" "200" "--d" "3" "${KOKKOS_THREADS_ARG}=2" "--kernel-type" "cardinalcubicbspline")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim3_QR_CardinalCubicBSpline PROPERTIES LABELS "IntegrationTest;integration;kokkos;bspline" TIMEOUT 10)

    # Device views tests for GMLS - Cholesky solver
    ADD_TEST(NAME GMLS_Device_Dim3_CHOL COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "--p" "4" "--nt" "200" "--d" "3" "--solver" "CHOL" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim3_CHOL PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

    # Device views tests with Neumann BC for GMLS - Cholesky solver
    ADD_TEST(NAME GMLS_NeumannGradScalar_Dim3_CHOL COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_NeumannGradScalar_Test "--p" "3" "--nt" "200" "--d" "3" "--solver" "CHOL" "--constraint" "NEUMANN_GRAD_SCALAR" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_NeumannGradScalar_Dim3_CHOL PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

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

    # Device views tests for GMLS (vector basis) with Cholesky
    ADD_TEST(NAME GMLS_Vector_Dim3_CHOL COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "--p" "3" "--nt" "200" "--d" "3" "--solver" "CHOL" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim3_CHOL PROPERTIES LABELS "IntegrationTest;integration;kokkos;vector" TIMEOUT 10)

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
    ADD_TEST(NAME GMLS_Staggered_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered "--p" "3" "--nt" "100" "--d" "3" "${KOKKOS_THREADS_ARG}=4")
    SET_TESTS_PROPERTIES(GMLS_Staggered_Dim3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;staggered" TIMEOUT 60)

    # Staggered scheme test for GMLS on non-manifold
    # Note: Using even polynomial order may cause this test to fail
    ADD_TEST(NAME GMLS_Staggered_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered "--p" "3" "--nt" "200" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
    SET_TESTS_PROPERTIES(GMLS_Staggered_Dim2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;staggered" TIMEOUT 20)

    if (NOT Compadre_DEBUG)
      # This test is too slow in DEBUG (3x longer than all other tests combined)
      # Multisite manifold test for GMLS
      #ADD_TEST(NAME GMLS_Manifold_MultiSite COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_MultiSite_Test "3" "4" "${KOKKOS_THREADS_ARG}=2")
      #SET_TESTS_PROPERTIES(GMLS_Manifold_MultiSite PROPERTIES LABELS "IntegrationTest;integration;kokkos;multi;manifold" TIMEOUT 10)
      if (Python_EXECUTABLE)
        # Python driven test of a C++ executable (Python changes command line arguments given to executable)
        CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py" @ONLY)
        ADD_TEST(NAME GMLS_Manifold_MultiSite_Refinement_Study COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py" "--porder=3" "--grids=3")
        SET_TESTS_PROPERTIES(GMLS_Manifold_MultiSite_Refinement_Study
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold" 
            TIMEOUT 60
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_MultiSite_Test>
        )
      endif() # Python_EXECUTABLE

      # Divergence-free basis test for GMLS on non-manifold
      ADD_TEST(NAME GMLS_DivergenceFree_Dim3_P3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "--p" "3" "--nt" "200" "--d" "3" "${KOKKOS_THREADS_ARG}=2")
      SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim3_P3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr;batched" TIMEOUT 60)
      #ADD_TEST(NAME GMLS_DivergenceFree_Dim3_P2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "2" "200" "3" "0" "0" "0" "${KOKKOS_THREADS_ARG}=2")
      #SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim3_P2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr" TIMEOUT 60)

      ADD_TEST(NAME GMLS_DivergenceFree_Dim2_P3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "--p" "3" "--nt" "200" "--d" "2" "${KOKKOS_THREADS_ARG}=2")
      SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim2_P3_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr;batched" TIMEOUT 60)
      #ADD_TEST(NAME GMLS_DivergenceFree_Dim2_P2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "2" "200" "2" "0" "0" "0" "${KOKKOS_THREADS_ARG}=2")
      #SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim2_P2_QR PROPERTIES LABELS "IntegrationTest;integration;kokkos;divergencefree;qr" TIMEOUT 60)
    endif()

    #    # Manifold tests for GMLS
    #    ADD_TEST(NAME GMLS_Manifold_Dim3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Test "4" "200" "3" "1" "0" "0" "${KOKKOS_THREADS_ARG}=2")
    #    SET_TESTS_PROPERTIES(GMLS_Manifold_Dim3 PROPERTIES LABELS "IntegrationTest;integration;kokkos" TIMEOUT 10)

    if (Python_EXECUTABLE)
      # Python driven test of a C++ executable (Python changes command line arguments given to executable) - calling QR solver
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" @ONLY)
      ADD_TEST(NAME GMLS_Manifold_Refinement_Study_QR COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=3" "--grids=4" "--solver-type=QR")
      SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_QR
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold" 
            TIMEOUT 20
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_Test>
      )

      # Python driven test of a C++ executable (Python changes command line arguments given to executable) - calling Cholesky solver
      ADD_TEST(NAME GMLS_Manifold_Refinement_Study_CHOL COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=3" "--grids=4" "--solver-type=CHOL")
      SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_CHOL
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold" 
            TIMEOUT 20
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_Test>
      )

      # Python driven test of a C++ executable - curvature order above the main order, so the
      # curvature basis is larger than the main basis (buffers sized for the larger one) - calling QR solver
      ADD_TEST(NAME GMLS_Manifold_Refinement_Study_Curvature_Order_QR COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=2" "--curvature-porder=3" "--grids=4" "--solver-type=QR")
      SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_Curvature_Order_QR
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold"
            TIMEOUT 20
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_Test>
      )

      # Python driven test of a C++ executable - curvature order above the main order, so the
      # curvature basis is larger than the main basis (buffers sized for the larger one) - calling CHOL solver
      ADD_TEST(NAME GMLS_Manifold_Refinement_Study_Curvature_Order_CHOL COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=2" "--curvature-porder=3" "--grids=4" "--solver-type=CHOL")
      SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_Curvature_Order_CHOL
          PROPERTIES
            LABELS "ConvergenceTest;convergence;manifold"
            TIMEOUT 20
            REQUIRED_FILES
              $<TARGET_FILE:GMLS_Manifold_Test>
      )

      # Python driven test of a C++ executable (Python changes command line arguments given to executable)
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Staggered_Manifold.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py" @ONLY)
      ADD_TEST(NAME GMLS_Staggered_Manifold_Refinement_Study COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py" "--porder=3" "--grids=4")
      SET_TESTS_PROPERTIES(GMLS_Staggered_Manifold_Refinement_Study
        PROPERTIES
          LABELS "ConvergenceTest;convergence;manifold;staggered" 
          TIMEOUT 20
          REQUIRED_FILES
            $<TARGET_FILE:GMLS_Staggered_Manifold_Test>
      )
    endif() # Python_EXECUTABLE

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
    SET_TESTS_PROPERTIES(NeighborRadiusSearch3DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 15)
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
    SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 15)
    #ADD_TEST(NAME NeighborKNNSearch3DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "160" "1.8" "1" "${KOKKOS_THREADS_ARG}=2")
    #SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
    if (Compadre_USE_PYTHON AND Python_EXECUTABLE)
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/test_pycompadre.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/test_pycompadre_manifold.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/test_tangent_planes_manifold.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/kokkos_test_case.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../pycompadre/examples/kokkos_lifetime_tests.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../pycompadre/examples")
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/Python_3D_Convergence.py.in" "${CMAKE_CURRENT_BINARY_DIR}/Python_3D_Convergence.py" @ONLY)
      ADD_TEST(NAME GMLS_Python_Convergence_Test_3d_Point_Reconstruction COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Python_3D_Convergence.py" "${KOKKOS_THREADS_ARG}=4")
      SET_TESTS_PROPERTIES(GMLS_Python_Convergence_Test_3d_Point_Reconstruction PROPERTIES LABELS "IntegrationTest;integration;python;kokkos" TIMEOUT 10)

      if (Compadre_USE_MATLAB)
        ADD_TEST(NAME GMLS_Matlab_Python_Interface_1d_Point_Reconstruction COMMAND "matlab" "-nodisplay" "-nojvm" "-nosplash" "-nodesktop" "-r \"try, run('${SWIG_PREFIX}/Matlab_1D_Using_Python_Interface.m'), catch, exit(1), end, exit(0);\"")
        SET_TESTS_PROPERTIES(GMLS_Matlab_Python_Interface_1d_Point_Reconstruction PROPERTIES LABELS "IntegrationTest;integration;python;kokkos;matlab" TIMEOUT 10)
      endif()

      ADD_TEST(NAME PyCOMPADRE_Tests COMMAND "${Python_EXECUTABLE}" "-m" "unittest" "discover" "-s" "${PYTHON_LIBRARY_PREFIX}/pycompadre/examples")
      SET_TESTS_PROPERTIES(PyCOMPADRE_Tests PROPERTIES ENVIRONMENT "PYTHONPATH=${PYTHON_LIBRARY_PREFIX};OMP_NUM_THREADS=4;OMP_PROC_BIND=false;OPENBLAS_NUM_THREADS=1"
          PROCESSORS 4)
      if (NOT Compadre_DEBUG)
        SET_TESTS_PROPERTIES(PyCOMPADRE_Tests PROPERTIES LABELS "pycompadre" TIMEOUT 60)
      else()
        SET_TESTS_PROPERTIES(PyCOMPADRE_Tests PROPERTIES LABELS "pycompadre" TIMEOUT 60)
      endif()

      # pycompadre.KokkosParser in a Python interpreter embedded in a C++ program
      target_link_libraries(KokkosParser_Test PRIVATE Python::Python)
      target_compile_definitions(KokkosParser_Test PRIVATE COMPADRE_TEST_EMBEDDED_PYTHON)
      foreach(mode embedded embedded_external)
        ADD_TEST(NAME KokkosParser_${mode} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/KokkosParser_Test "${mode}" "${KOKKOS_THREADS_ARG}=2")
        SET_TESTS_PROPERTIES(KokkosParser_${mode} PROPERTIES ENVIRONMENT PYTHONPATH=${PYTHON_LIBRARY_PREFIX} LABELS "pycompadre" TIMEOUT 30)
      endforeach()

      # KokkosParser: Kokkos runtime lifetime and arguments (each scenario runs in its own interpreter)
      ADD_TEST(NAME PyCOMPADRE_Kokkos_Lifetime COMMAND "${Python_EXECUTABLE}" "${PYTHON_LIBRARY_PREFIX}/pycompadre/examples/kokkos_lifetime_tests.py")
      SET_TESTS_PROPERTIES(PyCOMPADRE_Kokkos_Lifetime PROPERTIES ENVIRONMENT PYTHONPATH=${PYTHON_LIBRARY_PREFIX} LABELS "pycompadre" TIMEOUT 120)

    endif()

  endif (Compadre_TESTS)
endif (Compadre_EXAMPLES)

bob_end_subdir()


