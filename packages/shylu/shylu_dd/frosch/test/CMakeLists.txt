IF(ShyLU_DDFROSch_ENABLE_Epetra)
ADD_SUBDIRECTORIES(
    Import_Tpetra
    InterfacePartitionOfUnity
# Disabled because of Epetra
    # InterfaceSets
# Disabled because of Epetra
    # LocalPartitionOfUnityBasis
)
ENDIF()

IF(ShyLU_DDFROSch_ENABLE_Belos AND ShyLU_DDFROSch_ENABLE_Galeri AND ShyLU_DDFROSch_ENABLE_Stratimikos AND ShyLU_DDFROSch_ENABLE_Thyra)
ADD_SUBDIRECTORIES(
    Thyra_Xpetra_Laplace
)
ENDIF()

IF(ShyLU_DDFROSch_ENABLE_Belos AND ShyLU_DDFROSch_ENABLE_Stratimikos AND ShyLU_DDFROSch_ENABLE_Thyra AND ShyLU_DDFROSch_ENABLE_EpetraExt AND ShyLU_DDFROSch_ENABLE_HDF5)
ADD_SUBDIRECTORIES(
    Thyra_Xpetra_Stokes_HDF5
)
ENDIF()
