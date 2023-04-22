if (${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
  add_subdirectory(epetra)
endif()

if (${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
  add_subdirectory(tpetra)
endif()
