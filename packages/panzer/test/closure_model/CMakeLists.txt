INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

SET(SOURCES
    closure_model.cpp 
    user_app_ClosureModel_Factory.hpp
    user_app_ClosureModel_FactoryT.hpp 
    user_app_ClosureModel_Factory_TemplateBuilder.hpp
    user_app_ConstantModel.hpp
    user_app_ConstantModelT.hpp
    ) 

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  closure_model
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
