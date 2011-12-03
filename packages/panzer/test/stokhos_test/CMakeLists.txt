
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/ugi_rig)

SET(RIG_SOURCES
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/UnitTest_UniqueGlobalIndexer.cpp)
SET(RIG_HEADERS
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/UnitTest_UniqueGlobalIndexer.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_BCStrategy_Dirichlet_Constant.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_BCStrategy_Dirichlet_ConstantT.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_BCStrategy_Factory.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_ClosureModel_Factory.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_ClosureModel_FactoryT.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_ClosureModel_Factory_TemplateBuilder.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_ConstantModel.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_ConstantModelT.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_Convection.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_ConvectionT.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_EquationSetFactory.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_EquationSet_Energy.hpp
   ${PACKAGE_SOURCE_DIR}/test/ugi_rig/user_app_EquationSet_EnergyT.hpp)

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StokhosAssembly
  SOURCES 
     assembly_engine.cpp 
     ${RIG_SOURCES}
     ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )
