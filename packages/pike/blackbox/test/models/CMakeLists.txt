INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(HEADERS "")
SET(NOINSTALLHEADERS "")
SET(SOURCES "")

APPEND_SET(NOINSTALLHEADERS
  Pike_Mock_ModelEvaluator.hpp
  Pike_Mock_DataTransfer.hpp
  Pike_LinearHeatConduction_ModelEvaluator.hpp
  Pike_LinearHeatConduction_DataTransfer.hpp
  Pike_Mock_UserSolverFactory.hpp
  Pike_Mock_UserStatusTestFactory.hpp
  Pike_VanderPol_ModelEvaluator_Eq1.hpp
  Pike_VanderPol_ModelEvaluator_Eq2.hpp
  Pike_VanderPol_DataTransfer_Eq1ToEq2.hpp
  Pike_VanderPol_DataTransfer_Eq2ToEq1.hpp
  )

APPEND_SET(SOURCES
  Pike_Mock_ModelEvaluator.cpp
  Pike_Mock_DataTransfer.cpp
  Pike_LinearHeatConduction_ModelEvaluator.cpp
  Pike_LinearHeatConduction_DataTransfer.cpp
  Pike_Mock_UserSolverFactory.cpp
  Pike_Mock_UserStatusTestFactory.cpp
  Pike_VanderPol_ModelEvaluator_Eq1.cpp
  Pike_VanderPol_ModelEvaluator_Eq2.cpp
  Pike_VanderPol_DataTransfer_Eq1ToEq2.cpp
  Pike_VanderPol_DataTransfer_Eq2ToEq1.cpp
  )

TRIBITS_ADD_LIBRARY(
  pike-test-apps
  NO_INSTALL_LIB_OR_HEADERS
  NOINSTALLHEADERS ${NOINSTALLHEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  DEPLIBS pike-blackbox
  )

SET(UNIT_TEST_DRIVER ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_models
  SOURCES test_models.cpp ${UNIT_TEST_DRIVER}
  TESTONLYLIBS pike-test-apps
  NUM_MPI_PROCS 1
  )
