
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA AND NOX_ENABLE_LOCA AND
   NOX_ENABLE_Sacado AND NOX_ENABLE_EpetraExt AND Sacado_ENABLE_FEApp)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Thyra_ModelEvaluator_Test
    SOURCES 
    tcubed.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Thyra_ModelEvaluator_TurningPoint_Test
    SOURCES 
    tcubedTP.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()
