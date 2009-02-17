
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA AND
   NOX_ENABLE_Sacado AND NOX_ENABLE_EpetraExt AND Sacado_ENABLE_FEApp)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Sacado_FEApp_Brusselator_Example
    SOURCES brusselator.C
    COMM serial mpi
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Sacado_FEApp_Bratu_Example
    SOURCES bratu.C
    COMM serial mpi
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()
