IF(ShyLU_DDFROSch_ENABLE_Epetra)

ADD_SUBDIRECTORIES(
    Import_Tpetra
    InterfacePartitionOfUnity
# Disabled because of Epetra
    # InterfaceSets
# Disabled because of Epetra
    # LocalPartitionOfUnityBasis
    Thyra_Xpetra_Laplace
)

ENDIF()
