
ROL_ADD_EXECUTABLE(
  PolynomialRegression1D
  SOURCES example_01.cpp
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE(
  RobustGeneralizedPolynomialRegression1D
  SOURCES example_02.cpp
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE(
  LinearRegression1D
  SOURCES example_03.cpp
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE(
  IndependentMomentOperators
  SOURCES example_04.cpp
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE(
  DropoutGeneralizedPolynomialRegression1D
  SOURCES example_05.cpp
  ADD_DIR_TO_NAME
  )

IF(${PROJECT_NAME}_ENABLE_Teuchos)
ROL_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

ROL_ADD_EXECUTABLE(
  MultiObjectivePolynomialRegression1D
  SOURCES example_06.cpp
  ADD_DIR_TO_NAME
  )

ENDIF()

ROL_COPY_FILES_TO_BINARY_DIR(
  OEDExampleDataCopy
  SOURCE_FILES
    input.xml input_02.xml input_03.xml input_06.xml
    plotResults.m pointsGL.txt weightsGL.txt
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
)
