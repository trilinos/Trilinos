INCLUDE(TribitsAddExecutableAndTest)

IF(${PROJECT_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE(
    EpetraBatchManagerInterface
    SOURCES test_01.cpp
    ADD_DIR_TO_NAME
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
      solEpetraDataCopy
      SOURCE_FILES
        input.xml
      SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
      DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    )
ENDIF()

