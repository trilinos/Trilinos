
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  TRIBITS_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    DynamicNavierStokesDataCopy
    SOURCE_FILES
      input.xml
      plotControl.m
      plotFlow.m
      plotInitialFlow.m
      flowMovie.m
      vorticity.m
      trisurfMovie.m
      plotmesh.m
      parseData.m
      channel.txt
      channel.jou
      process.sh
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
