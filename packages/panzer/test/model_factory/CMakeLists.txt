
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

SET(SOURCES
    model_factory.cpp 
    user_app_ModelFactory.hpp
    user_app_ModelFactoryT.hpp 
    user_app_ConstantModel.hpp
    user_app_ConstantModelT.hpp
    user_app_ModelFactory_TemplateBuilder.hpp
    ) 

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  model_factory
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
