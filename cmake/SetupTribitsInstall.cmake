# Install TriBITS so that other projects can use it.
ADVANCED_SET(${PROJECT_NAME}_INSTALL_TriBITS ON CACHE BOOL
  "If ture, install TriBITS into <lib-install-dir>/cmake/tribits/")
IF (${PROJECT_NAME}_INSTALL_TriBITS)
  ASSERT_DEFINED(Trilinos_SOURCE_DIR)
  ASSERT_DEFINED(${PROJECT_NAME}_INSTALL_LIB_DIR)
  INSTALL(
    DIRECTORY "${Trilinos_SOURCE_DIR}/cmake/tribits"
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake"
    PATTERN "*.pyc" EXCLUDE
    )
ENDIF()
