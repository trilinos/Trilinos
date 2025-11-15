
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  # Need ROL_TpetraTeuchosBatchManager.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  ROL_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

#  ROL_ADD_EXECUTABLE(
#    example_02
#    SOURCES example_02.cpp
#    ADD_DIR_TO_NAME
#  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    AllenCahnDataCopy
    SOURCE_FILES
      input.xml plotresults.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()

IF(${PROJECT_NAME}_ENABLE_Intrepid2 AND
   ${PROJECT_NAME}_ENABLE_Ifpack2   AND
   ${PROJECT_NAME}_ENABLE_MueLu     AND
   ${PROJECT_NAME}_ENABLE_Amesos2   AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  # Need ROL_TpetraTeuchosBatchManager.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  ROL_ADD_EXECUTABLE(
    example_01_Kokkos
    SOURCES example_01K.cpp
    ADD_DIR_TO_NAME
  )

#  ROL_ADD_EXECUTABLE(
#    example_02
#    SOURCES example_02.cpp
#    ADD_DIR_TO_NAME
#  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    AllenCahnDataCopyK
    SOURCE_FILES
      input.xml plotresults.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()
