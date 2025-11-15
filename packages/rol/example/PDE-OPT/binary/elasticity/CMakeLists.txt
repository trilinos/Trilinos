
IF(${PROJECT_NAME}_ENABLE_Intrepid AND 
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_pebbl    AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/pebbl)

  # Need ROL/PEBBL interface
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/pebbl/src/algorithm)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/pebbl/src/function)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/pebbl/src/vector)

  ROL_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    BinTopoOptDataCopy
    SOURCE_FILES
      input_ex01.xml
      plotresults.m
      bwr.m
      plotdensity.m
      processpics.m
      picme
      process3D.sh
      meshfiles/brick1.txt
      meshfiles/brick1.e
      meshfiles/brick2.txt
      meshfiles/brick2.e
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
