
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  TRIBITS_ADD_EXECUTABLE(
    Burgers
    SOURCES example_Burgers.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    1dPoisson
    SOURCES example_1dPoisson.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    2dPoisson
    SOURCES example_2dPoisson.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    MPAKouriSurowiec2019DataCopy
    SOURCE_FILES
      input_Burgers.xml
      input_1dPoisson.xml
      input_2dPoisson.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
