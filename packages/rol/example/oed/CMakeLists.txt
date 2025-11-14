
ROL_ADD_EXECUTABLE(
  PolynomialRegression1D
  SOURCES example_01.cpp
  ADD_DIR_TO_NAME
  )

ROL_COPY_FILES_TO_BINARY_DIR(
  OEDExampleDataCopy
  SOURCE_FILES
    input.xml plotResults.m pointsGL.txt weightsGL.txt
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
)
