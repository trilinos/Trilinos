TRIBITS_SUBPACKAGE(Core)

# Consider 'Panzer_ENABLE_Epetra' only if '${PROJECT_NAME}_ENABLE_Epetra=ON'
# Use global_set, so PANZER_HAVE_EPETRA is visible in other subdirs
IF(${PROJECT_NAME}_ENABLE_Epetra AND ${PROJECT_NAME}_ENABLE_ThyraEpetraAdapters AND ${PROJECT_NAME}_ENABLE_ThyraEpetraExtAdapters)
  global_set(PANZER_HAVE_EPETRA ON)
ELSE()
   global_set(PANZER_HAVE_EPETRA OFF)
ENDIF()

IF(PANZER_HAVE_EPETRA)
   MESSAGE(STATUS "Global Indexer Epetra Utilities On")
ELSE()
   MESSAGE(STATUS "Global Indexer Epetra Utilities Off")
ENDIF()

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_SUBPACKAGE_POSTPROCESS()
