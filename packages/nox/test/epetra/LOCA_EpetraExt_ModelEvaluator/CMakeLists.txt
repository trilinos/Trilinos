
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA AND
  NOX_ENABLE_Sacado AND NOX_ENABLE_EpetraExt AND Sacado_ENABLE_FEApp)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_EpetraExt_ModelEvaluator_Test
    SOURCES 
    tcubed.C
    COMM serial mpi
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()
