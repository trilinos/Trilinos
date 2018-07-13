INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../test/unit_tests)


# IF (${PACKAGE_NAME}_ENABLE_Galeri)

#   IF ((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR
#       (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))

    TRIBITS_ADD_EXECUTABLE(
      RegionsDriver
      SOURCES Driver.cpp
      COMM serial mpi
      )

    # TRIBITS_COPY_FILES_TO_BINARY_DIR(tawiesn_cp
    #   SOURCE_FILES driver.xml
    #   )

#   ENDIF()

# ENDIF() # Galeri
