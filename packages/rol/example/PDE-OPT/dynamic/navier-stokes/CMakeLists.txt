IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  # Need for parallel-in-time.
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/mpi/src/function)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/mpi/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/mpi)

  ROL_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )
  ROL_ADD_EXECUTABLE(
    example_02
    SOURCES example_02.cpp
    ADD_DIR_TO_NAME
  )
  IF(TPL_ENABLE_MPI)
    ROL_ADD_EXECUTABLE(
      AugmentedSystem
      SOURCES AugmentedSystem.cpp ROL_PinTVectorCommunication_StdTpetraComposite.hpp
      ADD_DIR_TO_NAME
    )
  ENDIF()

  ROL_COPY_FILES_TO_BINARY_DIR(
    DynamicNavierStokesDataCopy
    SOURCE_FILES
      input.xml
      input_02.xml
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
