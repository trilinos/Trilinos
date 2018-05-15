IF(${PACKAGE_NAME}_ENABLE_Epetra)
ADD_SUBDIRECTORIES(
    TestLaplacian_Epetra
    TestInterfaceSets
)
ENDIF ()

IF(${PACKAGE_NAME}_ENABLE_Tpetra)
ADD_SUBDIRECTORIES(
    TestLaplacian_Tpetra
)
ENDIF ()
