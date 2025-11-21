
IF(${PROJECT_NAME}_ENABLE_Intrepid AND 
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

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

  ROL_COPY_FILES_TO_BINARY_DIR(
    TRSPG_TopoOptDataCopy
    SOURCE_FILES
      input_ex01.xml
      plotresults.m
      bwr.m
      plotdensity.m
      processpics.m
      picme
      process3D.sh
      meshfiles/truss2d.jou
      meshfiles/truss3d.jou
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()

IF(${PROJECT_NAME}_ENABLE_Intrepid2 AND 
   ${PROJECT_NAME}_ENABLE_Amesos2   AND
   ${PROJECT_NAME}_ENABLE_Ifpack2   AND
   ${PROJECT_NAME}_ENABLE_MueLu     AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  ROL_ADD_EXECUTABLE(
    example_01_Kokkos
    SOURCES example_01K.cpp
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE(
    example_02_Kokkos
    SOURCES example_02K.cpp
    ADD_DIR_TO_NAME
  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    TRSPG_TopoOptDataCopyK
    SOURCE_FILES
      input_ex01.xml
      plotresults.m
      bwr.m
      plotdensity.m
      processpics.m
      picme
      process3D.sh
      meshfiles/truss2d.jou
      meshfiles/truss3d.jou
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()
