
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

  ROL_ADD_EXECUTABLE(
    example_03
    SOURCES example_03.cpp
    ADD_DIR_TO_NAME
  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    TopoOptDataCopy
    SOURCE_FILES
      input_ex01.xml
      input_ex02.xml
      input_ex03.xml
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
      meshfiles/Lbeam.jou
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

  ROL_ADD_EXECUTABLE(
    example_03_Kokkos
    SOURCES example_03K.cpp
    ADD_DIR_TO_NAME
  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    TopoOptDataCopyK
    SOURCE_FILES
      input_ex01.xml
      input_ex02.xml
      input_ex03.xml
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
      meshfiles/Lbeam.jou
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()
