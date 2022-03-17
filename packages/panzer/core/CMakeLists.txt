TRIBITS_SUBPACKAGE(Core)

# Consider 'Panzer_ENABLE_Epetra' only if 'Trilinos_ENABLE_Epetra=ON'
IF(Trilinos_ENABLE_Epetra AND DEFINED Panzer_ENABLE_Epetra)
   SET(PANZER_HAVE_EPETRA ${Panzer_ENABLE_Epetra})
ELSE()
   SET(PANZER_HAVE_EPETRA ${Trilinos_ENABLE_Epetra})
ENDIF()

IF(PANZER_HAVE_EPETRA)
   MESSAGE(STATUS "Global Indexer Epetra Utilities On")
ELSE()
   MESSAGE(STATUS "Global Indexer Epetra Utilities Off")
ENDIF()

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_SUBPACKAGE_POSTPROCESS()
