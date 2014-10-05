INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

SET(SOURCES
    closure_model.cpp 
    user_app_ClosureModel_Factory.hpp
    user_app_ClosureModel_Factory_impl.hpp 
    user_app_ClosureModel_Factory_TemplateBuilder.hpp
    user_app_ConstantModel.hpp
    user_app_ConstantModel_impl.hpp
    ) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  closure_model
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
