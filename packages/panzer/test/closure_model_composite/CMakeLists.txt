INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

SET(SOURCES
    closure_model_composite.cpp 
    user_app_ClosureModel_Factory_Physics1.hpp
    user_app_ClosureModel_Factory_Physics1_impl.hpp 
    user_app_ClosureModel_Factory_Physics1_TemplateBuilder.hpp
    user_app_ClosureModel_Factory_Physics2.hpp
    user_app_ClosureModel_Factory_Physics2_impl.hpp 
    user_app_ClosureModel_Factory_Physics2_TemplateBuilder.hpp
    ) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  closure_model_composite
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
