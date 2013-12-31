INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(HEADERS "")
SET(NOINSTALLHEADERS "")
SET(SOURCES "")

APPEND_SET(NOINSTALLHEADERS
  Pike_Mock_ModelEvaluator.hpp
  Pike_LinearHeatConduction_ModelEvaluator.hpp
  Pike_LinearHeatConduction_DataTransfer.hpp
  )

APPEND_SET(SOURCES
  Pike_Mock_ModelEvaluator.cpp
  Pike_LinearHeatConduction_ModelEvaluator.cpp
  Pike_LinearHeatConduction_DataTransfer.cpp
  )

TRIBITS_ADD_LIBRARY(
  pike-test-apps
  NO_INSTALL_LIB_OR_HEADERS
  NOINSTALLHEADERS ${NOINSTALLHEADERS}
  SOURCES ${SOURCES}
  DEPLIBS pike-blackbox
  )

SET(UNIT_TEST_DRIVER ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  applications
  SOURCES applications.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(core_tests
#  SOURCE_FILES parameter_list_acceptance_test.xml parameter_list_acceptance_test2.xml
#  EXEDEPS parameter_list_acceptance_test
#  )
