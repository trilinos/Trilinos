
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  MixedPoissonExample
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  )

IF (PANZER_HAVE_EPETRA_STACK)
  TRIBITS_ADD_ADVANCED_TEST(
    MixedPoissonExample
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedPoissonExample
      ARGS --use-epetra --x-elements=5 --y-elements=5 --z-elements=5 --output-filename="base_mixed_poisson_epetra_"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
    TEST_1 EXEC MixedPoissonExample
      ARGS --use-tpetra --x-elements=5 --y-elements=5 --z-elements=5 --output-filename="base_mixed_poisson_tpetra_"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
    COMM mpi
  )
ELSE (PANZER_HAVE_EPETRA_STACK)
  TRIBITS_ADD_ADVANCED_TEST(
    MixedPoissonExample
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedPoissonExample
      ARGS --use-tpetra --x-elements=5 --y-elements=5 --z-elements=5 --output-filename="base_mixed_poisson_tpetra_"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
    COMM mpi
  )
ENDIF (PANZER_HAVE_EPETRA_STACK)

## basis order 1-3
FOREACH( ORDER 1 2 3)
  SET(HDIV_ORDER ${ORDER})
  MATH(EXPR HGRAD_ORDER "${ORDER}")

  TRIBITS_ADD_ADVANCED_TEST(
    MixedPoissonExample-ConvTest-Hex-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedPoissonExample
      ARGS --x-elements=3 --y-elements=3 --z-elements=3 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 3
      OUTPUT_FILE MPE-ConvTest-Hex-p${ORDER}-03
    TEST_1 EXEC MixedPoissonExample
      ARGS --x-elements=6 --y-elements=6 --z-elements=6 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Hex-p${ORDER}-06
    TEST_2 EXEC MixedPoissonExample
      ARGS --x-elements=12 --y-elements=12 --z-elements=12 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Hex-p${ORDER}-12
    TEST_3 CMND python3
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         ${ORDER}
         MPE-ConvTest-Hex-p${ORDER}-
         3
         6
         12
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
  )

ENDFOREACH()

## Variable coefficient and anisotropic tensor variants, basis order 1-2.
## VarCoeff weights the div-div term by three field multipliers, which
## exercises the multi-field-multiplier path of
## Integrator_DivBasisTimesScalar. Tensor uses a full SPD tensor in the
## (K grad phi, grad q) term, which exercises
## Integrator_GradBasisDotTensorTimesVector. Both AD kernels are tested.
FOREACH( VARIANT VarCoeff Tensor)
  IF (VARIANT STREQUAL "VarCoeff")
    SET(VARIANT_ARGS "--variable-coefficient")
  ELSE ()
    SET(VARIANT_ARGS "--anisotropic-tensor")
  ENDIF ()
  STRING(TOLOWER ${VARIANT} VARIANT_LC)

  FOREACH( ORDER 1 2)
    FOREACH( AD_MEM shared no-shared)
      IF (AD_MEM STREQUAL "shared")
        SET(AD_MEM_ARGS "--use-shared-mem-for-ad")
      ELSE ()
        SET(AD_MEM_ARGS "--no-use-shared-mem-for-ad")
      ENDIF ()

      TRIBITS_ADD_ADVANCED_TEST(
        MixedPoissonExample-${VARIANT}-ConvTest-Hex-Order-${ORDER}-${AD_MEM}
        EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
        TEST_0 EXEC MixedPoissonExample
          ARGS --x-elements=3 --y-elements=3 --z-elements=3 --hgrad-basis-order=${ORDER} --hdiv-basis-order=${ORDER} ${VARIANT_ARGS} ${AD_MEM_ARGS} --output-filename="${VARIANT_LC}-${ORDER}-${AD_MEM}-"
          PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
          NUM_MPI_PROCS 3
          OUTPUT_FILE MPE-${VARIANT}-ConvTest-Hex-p${ORDER}-${AD_MEM}-03
        TEST_1 EXEC MixedPoissonExample
          ARGS --x-elements=6 --y-elements=6 --z-elements=6 --hgrad-basis-order=${ORDER} --hdiv-basis-order=${ORDER} ${VARIANT_ARGS} ${AD_MEM_ARGS} --output-filename="${VARIANT_LC}-${ORDER}-${AD_MEM}-"
          PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
          NUM_MPI_PROCS 4
          OUTPUT_FILE MPE-${VARIANT}-ConvTest-Hex-p${ORDER}-${AD_MEM}-06
        TEST_2 EXEC MixedPoissonExample
          ARGS --x-elements=12 --y-elements=12 --z-elements=12 --hgrad-basis-order=${ORDER} --hdiv-basis-order=${ORDER} ${VARIANT_ARGS} ${AD_MEM_ARGS} --output-filename="${VARIANT_LC}-${ORDER}-${AD_MEM}-"
          PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
          NUM_MPI_PROCS 4
          OUTPUT_FILE MPE-${VARIANT}-ConvTest-Hex-p${ORDER}-${AD_MEM}-12
        TEST_3 CMND python3
          ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
             ${ORDER}
             MPE-${VARIANT}-ConvTest-Hex-p${ORDER}-${AD_MEM}-
             3
             6
             12
          PASS_REGULAR_EXPRESSION "Test Passed"
        COMM mpi
      )
    ENDFOREACH()
  ENDFOREACH()
ENDFOREACH()

## Multiblock basis order 1
FOREACH( ORDER 1)
  SET(HDIV_ORDER ${ORDER})
  MATH(EXPR HGRAD_ORDER "${ORDER}")

  TRIBITS_ADD_ADVANCED_TEST(
    MixedPoissonMultiblockExample-ConvTest-Hex-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedPoissonExample
      ARGS --x-blocks=2 --x-elements=6 --y-elements=6 --z-elements=6 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER} --output-filename="multiblock-${ORDER}-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 3
      OUTPUT_FILE MPE-Multiblock-ConvTest-Hex-p${ORDER}-06
    TEST_1 EXEC MixedPoissonExample
      ARGS --x-blocks=2 --x-elements=12 --y-elements=12 --z-elements=12 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER} --output-filename="multiblock-${ORDER}-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-Multiblock-ConvTest-Hex-p${ORDER}-12
    TEST_2 CMND python3
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         ${ORDER}
         MPE-Multiblock-ConvTest-Hex-p${ORDER}-
         6
         12
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
  )
ENDFOREACH()
