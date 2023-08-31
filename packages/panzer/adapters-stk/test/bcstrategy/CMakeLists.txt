

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/equation_set)
TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/closure_model)

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

SET(COMMON_HEADERS
    user_app_BCStrategy_Dirichlet_Constant.hpp
    user_app_BCStrategy_Dirichlet_Constant_impl.hpp
    user_app_BCStrategy_Neumann_Constant.hpp
    user_app_BCStrategy_Neumann_Constant_impl.hpp
    user_app_BCStrategy_Factory.hpp
    )



SET(SOURCES_TPETRA
    bcstrategy_tpetra.cpp
    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bcstrategy_tpetra
  SOURCES ${SOURCES_TPETRA} ${COMMON_HEADERS} ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bcstrategy_composite_factory
  SOURCES bcstrategy_composite_factory.cpp ${UNIT_TEST_DRIVER} user_app_BCStrategy_Dirichlet_Constant.hpp user_app_BCStrategy_Dirichlet_Constant_impl.hpp user_app_BCStrategy_Factory_Physics1.hpp user_app_BCStrategy_Factory_Physics2.hpp
  NUM_MPI_PROCS 1
  )

IF (PANZER_HAVE_EPETRA)
    SET(SOURCES_EPETRA
        bcstrategy.cpp
        )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bcstrategy_epetra
      SOURCES ${SOURCES_EPETRA} ${COMMON_HEADERS} ${UNIT_TEST_DRIVER}
      NUM_MPI_PROCS 1
      )
ENDIF(PANZER_HAVE_EPETRA)


