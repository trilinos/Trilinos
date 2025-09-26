tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#tribits_add_executable(
#  UnitTests
#  SOURCES
#    Compadre_UnitTests.cpp
#  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Host_Test
  SOURCES
    GMLS_Host.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Device_Test
  SOURCES
    GMLS_Device.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_NeumannGradScalar_Test
  SOURCES
    GMLS_NeumannGradScalar.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Vector_Test
  SOURCES
    GMLS_Vector.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Divergence_Test
  SOURCES
    GMLS_DivergenceFree.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_SmallBatchReuse_Device_Test
  SOURCES
    GMLS_SmallBatchReuse_Device.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Manifold_Test
  SOURCES
    GMLS_Manifold.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Staggered
  SOURCES
    GMLS_Staggered.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Staggered_Manifold_Test
  SOURCES
    GMLS_Staggered_Manifold.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_MultiSite_Test
  SOURCES
    GMLS_Multiple_Evaluation_Sites.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  GMLS_Manifold_MultiSite_Test
  SOURCES
    GMLS_Manifold_Multiple_Evaluation_Sites.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  TestUtility
  SOURCES
    UtilityTest.cpp
  ) # end tribits_add_executable

tribits_add_executable(
  NeighborSearchTest
  SOURCES
    NeighborSearchTest.cpp
  ) # end tribits_add_executable

#set(testName Compadre_Unit_Tests)
#tribits_add_test(
#  UnitTests
#  NAME
#    ${testName}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
#  ) # end tribits_add_test
#if (${testName}_CREATED)
#  set_tests_properties(
#    ${${testName}_CREATED}
#    PROPERTIES
#      LABELS
#        "UnitTest;unittest;Unit;unit"
#      TIMEOUT
#        60
#    ) # end set_tests_properties
#endif() # test created

# Host views tests for GMLS
set(testName GMLS_Host_Dim3_QR)
tribits_add_test(
  GMLS_Host_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 3 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName GMLS_Host_Dim2_QR)
tribits_add_test(
  GMLS_Host_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 2 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName GMLS_Host_Dim1_QR)
tribits_add_test(
  GMLS_Host_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Device views tests for GMLS
set(testName GMLS_Device_Dim3_QR)
tribits_add_test(
  GMLS_Device_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 3 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName GMLS_Device_Dim2_QR)
tribits_add_test(
  GMLS_Device_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 2 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName GMLS_Device_Dim1_QR)
tribits_add_test(
  GMLS_Device_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Device views tests for GMLS - LU solver
set(testName GMLS_Device_Dim3_LU)
tribits_add_test(
  GMLS_Device_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 3 --solver LU --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Device views tests with Neumann BC for GMLS - LU solver
set(testName GMLS_NeumannGradScalar_Dim3_LU)
tribits_add_test(
  GMLS_NeumannGradScalar_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 3 --nt 200 --d 3 --solver LU --constraint NEUMANN_GRAD_SCALAR --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

## Device views tests with Neumann BC for STAGGERED GMLS - QR solver
#set(testName GMLS_StaggeredNeumannGradScalar_Dim3_QR)
#tribits_add_test(
#  GMLS_Staggered
#  NAME
#    ${testName}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS
#    "3 200 3 0 0 1 --kokkos-num-threads=2"
#  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
#  ) # end tribits_add_test
#if (${testName}_CREATED)
#  set_tests_properties(
#    ${${testName}_CREATED}
#    PROPERTIES
#      LABELS
#        "IntegrationTest;integration;kokkos;staggered"
#      TIMEOUT
#        60
#    ) # end set_tests_properties
#endif() # test created

# Device views tests for GMLS (vector basis)
set(testName GMLS_Vector_Dim3_QR)
tribits_add_test(
  GMLS_Vector_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 3 --nt 20 --d 3 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos;vector"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName GMLS_Vector_Dim2_QR)
tribits_add_test(
  GMLS_Vector_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 3 --nt 20 --d 2 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos;vector"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName GMLS_Vector_Dim1_QR)
tribits_add_test(
  GMLS_Vector_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 3 --nt 20 --d 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos;vector"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Device views tests for GMLS (vector basis) with LU
set(testName GMLS_Vector_Dim3_LU)
tribits_add_test(
  GMLS_Vector_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 3 --nt 200 --d 3 --solver LU --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos;vector"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Device views tests for small batch GMLS, reusing GMLS class object
set(testName GMLS_SmallBatchReuse_Device_Dim2_QR)
tribits_add_test(
  GMLS_SmallBatchReuse_Device_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  COMM serial mpi
  ARGS
    "--p 4 --nt 200 --d 2 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName GMLS_SmallBatchReuse_Device_Dim1_QR)
tribits_add_test(
  GMLS_SmallBatchReuse_Device_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Multisite test for GMLS
set(testName GMLS_MultiSite_Dim3_QR)
tribits_add_test(
  GMLS_MultiSite_Test
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 4 --nt 200 --d 3 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Staggered scheme test for GMLS on non-manifold
# Note: Using even polynomial order may cause this test to fail
set(testName GMLS_Staggered_Dim3_QR)
tribits_add_test(
  GMLS_Staggered
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 3 --nt 100 --d 3 --kokkos-num-threads=4"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos;staggered"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Staggered scheme test for GMLS on non-manifold
# Note: Using even polynomial order may cause this test to fail
set(testName GMLS_Staggered_Dim2_QR)
tribits_add_test(
  GMLS_Staggered
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--p 3 --nt 200 --d 2 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "IntegrationTest;integration;kokkos;staggered"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

if (NOT(Compadre_DEBUG OR Compadre_EXTREME_DEBUG))
  # This test is too slow in DEBUG (3x longer than all other tests
  # combined)

  if (Python3_EXECUTABLE)
    # Python driven test of a C++ executable (Python changes command line
    # arguments given to executable)
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py.in
      ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py
      @ONLY
      ) # end configure_file
    set(testName GMLS_Manifold_Multiple_Evaluation_Sites)
    TRIBITS_ADD_ADVANCED_TEST(
      ${testName}
      TEST_0 CMND ${Python3_EXECUTABLE} ARGS ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold_Multiple_Evaluation_Sites.py --porder=3 --grids=3 --in-trilinos=True
      PASS_REGULAR_EXPRESSION "Passed."
      COMM mpi serial
      ADDED_TEST_NAME_OUT ${testName}_CREATED
    )
    if (${testName}_CREATED)
      set_tests_properties(
        ${${testName}_CREATED}
        PROPERTIES
          LABELS
            "ConvergenceTest;convergence;manifold"
          TIMEOUT
            60
          REQUIRED_FILES
            $<TARGET_FILE:Compadre_GMLS_Manifold_MultiSite_Test>
      )
    endif() # test created
  endif() # Python3_EXECUTABLE

  # Divergence-free basis test for GMLS on non-manifold
  # Note: QR is needed to be used here due to the null space introduced
  set(testName GMLS_DivergenceFree_Dim3_P3_QR)
  tribits_add_test(
    GMLS_Divergence_Test
    NAME
      ${testName}
    NUM_MPI_PROCS 1
    ARGS
      "--p 3 --nt 200 --d 3 --kokkos-num-threads=2"
    ADDED_TESTS_NAMES_OUT ${testName}_CREATED
    ) # end tribits_add_test
  if (${testName}_CREATED)
    set_tests_properties(
      ${${testName}_CREATED}
      PROPERTIES
        LABELS
          "IntegrationTest;integration;kokkos;divergencefree;qr;batched"
        TIMEOUT
          60
      ) # end set_tests_properties
  endif() # test created
  #set(testName GMLS_DivergenceFree_Dim3_P2_QR)
  #tribits_add_test(
  #  GMLS_Divergence_Test
  #  NAME
  #    ${testName}
  #  NUM_MPI_PROCS 1
  #  ARGS
  #    "2 200 3 0 0 0 --kokkos-num-threads=2"
  #  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  #  ) # end tribits_add_test
  #if (${testName}_CREATED)
  #  set_tests_properties(
  #    ${${testName}_CREATED}
  #    PROPERTIES
  #      LABELS
  #        "IntegrationTest;integration;kokkos;divergencefree;qr"
  #      TIMEOUT
  #        60
  #    ) # end set_tests_properties
  #endif() # test created

  # Divergence-free basis test for GMLS on non-manifold
  set(testName GMLS_DivergenceFree_Dim2_P3_QR)
  tribits_add_test(
    GMLS_Divergence_Test
    NAME
      ${testName}
    NUM_MPI_PROCS 1
    ARGS
      "--p 3 --nt 200 --d 2 --kokkos-num-threads=2"
    ADDED_TESTS_NAMES_OUT ${testName}_CREATED
    ) # end tribits_add_test
  if (${testName}_CREATED)
    set_tests_properties(
      ${${testName}_CREATED}
      PROPERTIES
        LABELS
          "IntegrationTest;integration;kokkos;divergencefree;qr;batched"
        TIMEOUT
          60
      ) # end set_tests_properties
  endif() # test created
  #set(testName GMLS_DivergenceFree_Dim2_P2_QR)
  #tribits_add_test(
  #  GMLS_Divergence_Test
  #  NAME
  #    ${testName}
  #  NUM_MPI_PROCS 1
  #  ARGS
  #    "2 200 2 0 0 0 --kokkos-num-threads=2"
  #  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  #  ) # end tribits_add_test
  #if (${testName}_CREATED)
  #  set_tests_properties(
  #    ${${testName}_CREATED}
  #    PROPERTIES
  #      LABELS
  #        "IntegrationTest;integration;kokkos;divergencefree;qr"
  #      TIMEOUT
  #        60
  #    ) # end set_tests_properties
  #endif() # test created
endif() # not debug

if (Python3_EXECUTABLE)
  # Python driven test of a C++ executable (Python changes command line
  # arguments given to executable) - calling QR solver
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Manifold.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py
    @ONLY
    ) # end configure_file
  set(testName GMLS_Manifold_Refinement_Study_QR)
  TRIBITS_ADD_ADVANCED_TEST(
    ${testName}
    TEST_0 CMND ${Python3_EXECUTABLE} ARGS ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py --porder=3 --grids=4 --solver-type=QR --in-trilinos=True
    PASS_REGULAR_EXPRESSION "Passed."
    COMM mpi serial
    ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
  if (${testName}_CREATED)
    set_tests_properties(
      ${${testName}_CREATED}
      PROPERTIES
        LABELS
          "ConvergenceTest;convergence;manifold"
        TIMEOUT
          60
        REQUIRED_FILES
          $<TARGET_FILE:Compadre_GMLS_Manifold_Test>
    )
  endif() # test created
  
  # Python driven test of a C++ executable (Python changes command line
  # arguments given to executable) - calling LU solver
  set(testName GMLS_Manifold_Refinement_Study_LU)
  TRIBITS_ADD_ADVANCED_TEST(
    ${testName}
    TEST_0 CMND ${Python3_EXECUTABLE} ARGS ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Manifold.py --porder=3 --grids=4 --solver-type=LU --in-trilinos=True
    PASS_REGULAR_EXPRESSION "Passed."
    COMM mpi serial
    ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
  if (${testName}_CREATED)
    set_tests_properties(
      ${${testName}_CREATED}
      PROPERTIES
        LABELS
          "ConvergenceTest;convergence;manifold"
        TIMEOUT
          60
        REQUIRED_FILES
          $<TARGET_FILE:Compadre_GMLS_Manifold_Test>
    )
  endif() # test created

  # Python driven test of a C++ executable (Python changes command line
  # arguments given to executable) - calling QR solver
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/GMLS_Staggered_Manifold.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py
    @ONLY
    ) # end configure_file
  # Python driven test of a C++ executable (Python changes command line
  # arguments given to executable)
  set(testName GMLS_Staggered_Manifold_Refinement_Study)
  TRIBITS_ADD_ADVANCED_TEST(
    ${testName}
    TEST_0 CMND ${Python3_EXECUTABLE} ARGS ${CMAKE_CURRENT_BINARY_DIR}/GMLS_Staggered_Manifold.py --porder=3 --grids=4 --in-trilinos=True
    PASS_REGULAR_EXPRESSION "Passed."
    COMM mpi serial
    ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
  if (${testName}_CREATED)
    set_tests_properties(
      ${${testName}_CREATED}
      PROPERTIES
        LABELS
          "ConvergenceTest;convergence;manifold;staggered"
        TIMEOUT
          60
        REQUIRED_FILES
          $<TARGET_FILE:Compadre_GMLS_Staggered_Manifold_Test>
    )
  endif() # test created
endif() # Python3_EXECUTABLE

# Utility test - Filter By ID
set(testName Test_Utilities)
tribits_add_test(
  TestUtility
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "200 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "UtilityTest;utility;kokkos"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

# Neighbor radius search - 2D
set(testName NeighborRadiusSearch2DTest_1)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "2 200 6.5 0 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName NeighborRadiusSearch2DTest_2)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "2 300 4.5 0 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

#set(testName NeighborRadiusSearch2DTest_3)
#tribits_add_test(
#  NeighborSearchTest
#  NAME
#    ${testName}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS
#    "2 400 1.8 --kokkos-num-threads=2"
#  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
#  ) # end tribits_add_test
#if (${testName}_CREATED)
#  set_tests_properties(
#    ${${testName}_CREATED}
#    PROPERTIES
#      LABELS
#        "kdtree;nanoflann"
#      TIMEOUT
#        60
#    ) # end set_tests_properties
#endif() # test created

# Neighbor radius search - 3D
set(testName NeighborRadiusSearch3DTest_1)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "3 100 4.5 0 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName NeighborRadiusSearch3DTest_2)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "3 150 3.5 0 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

#set(testName NeighborRadiusSearch3DTest_3)
#tribits_add_test(
#  NeighborSearchTest
#  NAME
#    ${testName}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS
#    "3 160 1.8 --kokkos-num-threads=2"
#  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
#  ) # end tribits_add_test
#if (${testName}_CREATED)
#  set_tests_properties(
#    ${${testName}_CREATED}
#    PROPERTIES
#      LABELS
#        "kdtree;nanoflann"
#      TIMEOUT
#        60
#    ) # end set_tests_properties
#endif() # test created

# Neighbor KNN search - 2D
set(testName NeighborKNNSearch2DTest_1)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "2 200 6.5 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName NeighborKNNSearch2DTest_2)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "2 300 4.5 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

#set(testName NeighborKNNSearch2DTest_3)
#tribits_add_test(
#  NeighborSearchTest
#  NAME
#    ${testName}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS
#    "2 400 1.8 --kokkos-num-threads=2"
#  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
#  ) # end tribits_add_test
#if (${testName}_CREATED)
#  set_tests_properties(
#    ${${testName}_CREATED}
#    PROPERTIES
#      LABELS
#        "kdtree;nanoflann"
#      TIMEOUT
#        60
#    ) # end set_tests_properties
#endif() # test created

# Neighbor KNN search - 3D
set(testName NeighborKNNSearch3DTest_1)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "3 100 4.5 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

set(testName NeighborKNNSearch3DTest_2)
tribits_add_test(
  NeighborSearchTest
  NAME
    ${testName}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "3 150 3.5 1 --kokkos-num-threads=2"
  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
  ) # end tribits_add_test
if (${testName}_CREATED)
  set_tests_properties(
    ${${testName}_CREATED}
    PROPERTIES
      LABELS
        "kdtree;nanoflann"
      TIMEOUT
        60
    ) # end set_tests_properties
endif() # test created

#set(testName NeighborKNNSearch3DTest_3)
#tribits_add_test(
#  NeighborSearchTest
#  NAME
#    ${testName}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS
#    "3 160 1.8 1 --kokkos-num-threads=2"
#  ADDED_TESTS_NAMES_OUT ${testName}_CREATED
#  ) # end tribits_add_test
#if (${testName}_CREATED)
#  set_tests_properties(
#    ${${testName}_CREATED}
#    PROPERTIES
#      LABELS
#        "kdtree;nanoflann"
#      TIMEOUT
#        60
#    ) # end set_tests_properties
#endif() # test created
