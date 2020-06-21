include_tribits_build()

# ERROR MESSAGES IF A TEST IS NOT RUN
SET(OPENBLAS_NUM_THREADS "1")


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

  if(Compadre_USE_LAPACK)
    add_exe_w_compadre(LAPACK_Test UnitTest_ThreadedLapack.cpp) 
  endif()
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

    if(Compadre_USE_LAPACK)
      # Test if LAPACK+BLAS are compatible for use in the toolkit
      ADD_TEST(NAME LAPACK_THREADSAFE COMMAND ${CMAKE_CURRENT_BINARY_DIR}/LAPACK_Test "--kokkos-threads=8" --verbose --output-on-failure)
      SET_TESTS_PROPERTIES(LAPACK_THREADSAFE PROPERTIES LABELS "UnitTest;unit;lapack;Lapack;LAPACK;kokkos" TIMEOUT 5)
      SET_TESTS_PROPERTIES(LAPACK_THREADSAFE PROPERTIES SKIP_RETURN_CODE 77)
    endif()

    # Host views tests for GMLS
    ADD_TEST(NAME GMLS_Host_Dim3_SVD COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Host_Test "4" "200" "3" "0" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Host_Dim3_SVD PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 60)
    
    ADD_TEST(NAME GMLS_Host_Dim2_SVD COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Host_Test "4" "200" "2" "0" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Host_Dim2_SVD PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 20)
    
    ADD_TEST(NAME GMLS_Host_Dim1_SVD COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Host_Test "4" "200" "1" "0" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Host_Dim1_SVD PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 20)

    # Device views tests for GMLS
    ADD_TEST(NAME GMLS_Device_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "4" "200" "3" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim3_QR PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Device_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "4" "200" "2" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim2_QR PROPERTIES LABELS "UnitTest;unit;kokkos;batch" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Device_Dim1_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "4" "200" "1" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim1_QR PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 10)

    # Device views tests for GMLS - LU solver
    ADD_TEST(NAME GMLS_Device_Dim3_LU COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Device_Test "4" "200" "3" "2" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Device_Dim3_LU PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 10)

    # Device views tests with Neumann BC for GMLS - LU solver
    ADD_TEST(NAME GMLS_NeumannGradScalar_Dim3_LU COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_NeumannGradScalar_Test "3" "200" "3" "2" "0" "1" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_NeumannGradScalar_Dim3_LU PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 10)

    # Device views tests with Neumann BC for STAGGERED GMLS - SVD solver
    ADD_TEST(NAME GMLS_StaggeredNeumannGradScalar_Dim3_SVD COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered "3" "200" "3" "0" "0" "1" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_StaggeredNeumannGradScalar_Dim3_SVD PROPERTIES LABELS "UnitTest;unit;kokkos;staggered" TIMEOUT 10)

    # Device views tests for GMLS (vector basis)
    ADD_TEST(NAME GMLS_Vector_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "3" "20" "3" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim3_QR PROPERTIES LABELS "UnitTest;unit;kokkos;vector" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Vector_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "3" "20" "2" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim2_QR PROPERTIES LABELS "UnitTest;unit;kokkos;vector" TIMEOUT 10)
    
    ADD_TEST(NAME GMLS_Vector_Dim1_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "3" "20" "1" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim1_QR PROPERTIES LABELS "UnitTest;unit;kokkos;vector" TIMEOUT 10)

    # Device views tests for GMLS (vector basis) with LU
    ADD_TEST(NAME GMLS_Vector_Dim3_LU COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Vector_Test "3" "200" "3" "2" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Vector_Dim3_LU PROPERTIES LABELS "UnitTest;unit;kokkos;vector" TIMEOUT 10)

    # Device views tests for small batch GMLS, reusing GMLS class object
    ADD_TEST(NAME GMLS_SmallBatchReuse_Device_Dim2_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_SmallBatchReuse_Device_Test "4" "200" "2" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_SmallBatchReuse_Device_Dim2_QR PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 20)
    
    ADD_TEST(NAME GMLS_SmallBatchReuse_Device_Dim1_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_SmallBatchReuse_Device_Test "4" "200" "1" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_SmallBatchReuse_Device_Dim1_QR PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 20)

    # Multisite test for GMLS
    ADD_TEST(NAME GMLS_MultiSite_Dim3_QR COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_MultiSite_Test "4" "200" "3" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_MultiSite_Dim3_QR PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 10)

    # Staggered scheme test for GMLS on non-manifold
    # Note: Using even polynomial order may cause this test to fail
    ADD_TEST(NAME GMLS_Staggered_Dim3_QR COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Staggered "3" "100" "3" "1" "0" "0" "--kokkos-threads=4")
    SET_TESTS_PROPERTIES(GMLS_Staggered_Dim3_QR PROPERTIES LABELS "UnitTest;unit;kokkos;staggered" TIMEOUT 60)

    # Staggered scheme test for GMLS on non-manifold
    # Note: Using even polynomial order may cause this test to fail
    ADD_TEST(NAME GMLS_Staggered_Dim2_QR COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Staggered "3" "200" "2" "1" "0" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(GMLS_Staggered_Dim2_QR PROPERTIES LABELS "UnitTest;unit;kokkos;staggered" TIMEOUT 20)

    if (NOT Compadre_DEBUG)
      # This test is too slow in DEBUG (3x longer than all other tests combined)
      # Multisite manifold test for GMLS
      #ADD_TEST(NAME GMLS_Manifold_MultiSite COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_MultiSite_Test "3" "4" "--kokkos-threads=2")
      #SET_TESTS_PROPERTIES(GMLS_Manifold_MultiSite PROPERTIES LABELS "UnitTest;unit;kokkos;multi;manifold" TIMEOUT 10)
      # Python driven test of a C++ executable (Python changes command line arguments given to executable)
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py" @ONLY)
      ADD_TEST(NAME GMLS_Manifold_MultiSite_Refinement_Study COMMAND "python" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py" "--porder=3" "--grids=3")
      SET_TESTS_PROPERTIES(GMLS_Manifold_MultiSite_Refinement_Study PROPERTIES LABELS "ConvergenceTest;convergence;manifold" TIMEOUT 60)

      # Divergence-free basis test for GMLS on non-manifold
      # Note: SVD is needed to be used here due to the null space introduced
      ADD_TEST(NAME GMLS_DivergenceFree_Dim3_P3_SVD COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "3" "200" "3" "0" "0" "0" "--kokkos-threads=2")
      SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim3_P3_SVD PROPERTIES LABELS "UnitTest;unit;kokkos;divergencefree;svd;batched" TIMEOUT 60)
      #ADD_TEST(NAME GMLS_DivergenceFree_Dim3_P2_SVD COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "2" "200" "3" "0" "0" "0" "--kokkos-threads=2")
      #SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim3_P2_SVD PROPERTIES LABELS "UnitTest;unit;kokkos;divergencefree;svd" TIMEOUT 60)

      ADD_TEST(NAME GMLS_DivergenceFree_Dim2_P3_SVD COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "3" "200" "2" "0" "0" "0" "--kokkos-threads=2")
      SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim2_P3_SVD PROPERTIES LABELS "UnitTest;unit;kokkos;divergencefree;svd;batched" TIMEOUT 60)
      #ADD_TEST(NAME GMLS_DivergenceFree_Dim2_P2_SVD COMMAND ${CMAKD_CURRENT_BINARY_DIR}/GMLS_Divergence_Test "2" "200" "2" "0" "0" "0" "--kokkos-threads=2")
      #SET_TESTS_PROPERTIES(GMLS_DivergenceFree_Dim2_P2_SVD PROPERTIES LABELS "UnitTest;unit;kokkos;divergencefree;svd" TIMEOUT 60)
    endif()

    #    # Manifold tests for GMLS
    #    ADD_TEST(NAME GMLS_Manifold_Dim3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Test "4" "200" "3" "1" "0" "0" "--kokkos-threads=2")
    #    SET_TESTS_PROPERTIES(GMLS_Manifold_Dim3 PROPERTIES LABELS "UnitTest;unit;kokkos" TIMEOUT 10)

    # Python driven test of a C++ executable (Python changes command line arguments given to executable) - calling QR solver
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" @ONLY)
    ADD_TEST(NAME GMLS_Manifold_Refinement_Study_QR COMMAND "python" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=3" "--grids=4" "--solver-type=1")
    SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_QR PROPERTIES LABELS "ConvergenceTest;convergence;manifold" TIMEOUT 20)

    # Python driven test of a C++ executable (Python changes command line arguments given to executable) - calling LU solver
    ADD_TEST(NAME GMLS_Manifold_Refinement_Study_LU COMMAND "python" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py" "--porder=3" "--grids=4" "--solver-type=2")
    SET_TESTS_PROPERTIES(GMLS_Manifold_Refinement_Study_LU PROPERTIES LABELS "ConvergenceTest;convergence;manifold" TIMEOUT 20)

    # Python driven test of a C++ executable (Python changes command line arguments given to executable)
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Staggered_Manifold.py.in" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py" @ONLY)
    ADD_TEST(NAME GMLS_Staggered_Manifold_Refinement_Study COMMAND "python" "${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py" "--porder=3" "--grids=4")
    SET_TESTS_PROPERTIES(GMLS_Staggered_Manifold_Refinement_Study PROPERTIES LABELS "ConvergenceTest;convergence;manifold;staggered" TIMEOUT 20)

    # Utility test - Filter By ID
    ADD_TEST(NAME Test_Utilities COMMAND ${CMAKE_CURRENT_BINARY_DIR}/TestUtility "200" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(Test_Utilities PROPERTIES LABELS "UtilityTest;utility;kokkos" TIMEOUT 5)

    # Neighbor radius search - 2D
    ADD_TEST(NAME NeighborRadiusSearch2DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "200" "6.5" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch2DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborRadiusSearch2DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "300" "4.5" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch2DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborRadiusSearch2DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "400" "1.8" "--kokkos-threads=2")
    #SET_TESTS_PROPERTIES(NeighborRadiusSearch2DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # Neighbor radius search - 3D
    ADD_TEST(NAME NeighborRadiusSearch3DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "100" "4.5" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch3DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborRadiusSearch3DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "150" "3.5" "0" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborRadiusSearch3DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborRadiusSearch3DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "160" "1.8" "--kokkos-threads=2")
    #SET_TESTS_PROPERTIES(NeighborRadiusSearch3DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # Neighbor KNN search - 2D
    ADD_TEST(NAME NeighborKNNSearch2DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "200" "6.5" "1" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch2DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborKNNSearch2DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "300" "4.5" "1" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch2DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborKNNSearch2DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "2" "400" "1.8" "1" "--kokkos-threads=2")
    #SET_TESTS_PROPERTIES(NeighborKNNSearch2DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # Neighbor KNN search - 3D
    ADD_TEST(NAME NeighborKNNSearch3DTest_1 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "100" "4.5" "1" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_1 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    ADD_TEST(NAME NeighborKNNSearch3DTest_2 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "150" "3.5" "1" "--kokkos-threads=2")
    SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_2 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)
    #ADD_TEST(NAME NeighborKNNSearch3DTest_3 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/NeighborSearchTest "3" "160" "1.8" "1" "--kokkos-threads=2")
    #SET_TESTS_PROPERTIES(NeighborKNNSearch3DTest_3 PROPERTIES LABELS "kdtree;nanoflann;" TIMEOUT 5)

    # WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
    if (Compadre_USE_PYTHON)
      CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/Python_3D_Convergence.py.in" "${CMAKE_CURRENT_BINARY_DIR}/Python_3D_Convergence.py" @ONLY)
      ADD_TEST(NAME GMLS_Python_Convergence_Test_3d_Point_Reconstruction COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Python_3D_Convergence.py" "--kokkos-threads=4")
      SET_TESTS_PROPERTIES(GMLS_Python_Convergence_Test_3d_Point_Reconstruction PROPERTIES LABELS "UnitTest;unit;python;kokkos" TIMEOUT 10)

      if (Compadre_USE_MATLAB)
        ADD_TEST(NAME GMLS_Matlab_Python_Interface_1d_Point_Reconstruction COMMAND "matlab" "-nodisplay" "-nojvm" "-nosplash" "-nodesktop" "-r \"try, run('${SWIG_PREFIX}/Matlab_1D_Using_Python_Interface.m'), catch, exit(1), end, exit(0);\"")
        SET_TESTS_PROPERTIES(GMLS_Matlab_Python_Interface_1d_Point_Reconstruction PROPERTIES LABELS "UnitTest;unit;python;kokkos;matlab" TIMEOUT 10)
      endif()
    endif()

  endif (Compadre_TESTS)
endif (Compadre_EXAMPLES)

bob_end_subdir()


