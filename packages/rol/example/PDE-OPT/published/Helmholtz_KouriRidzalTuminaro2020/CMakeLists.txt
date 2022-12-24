
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  TRIBITS_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

  FILE(GLOB INP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} inputs/*.xml)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    HelmholtzKouriRidzalTuminaro2020DataCopy
    SOURCE_FILES
      input_ex01.xml
      README.md
      study.m
      run.m
      prdw.m
      imshift.m
      mesh/README.md
      mesh/helmholtz-mesh-0.e
      mesh/helmholtz-mesh-1.e
      mesh/helmholtz-mesh-0.txt
      mesh/helmholtz-mesh-1.txt
      mesh/helmholtz-mesh.jou
      results/README.md
      ${INP_FILES}
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
