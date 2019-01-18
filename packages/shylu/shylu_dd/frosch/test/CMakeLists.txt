ADD_SUBDIRECTORIES(
    Import_Tpetra
    InterfacePartitionOfUnity
    InterfaceSets
    LocalPartitionOfUnityBasis
)

IF(ShyLU_DDFROSch_ENABLE_Belos AND ShyLU_DDFROSch_ENABLE_Stratimikos AND ShyLU_DDFROSch_ENABLE_Thyra)
ADD_SUBDIRECTORIES(
    Thyra_Xpetra_Laplace
)
ENDIF()

IF(ShyLU_DDFROSch_ENABLE_Belos AND ShyLU_DDFROSch_ENABLE_Stratimikos AND ShyLU_DDFROSch_ENABLE_Thyra AND ShyLU_DDFROSch_ENABLE_EpetraExt AND ShyLU_DDFROSch_ENABLE_HDF5)
ADD_SUBDIRECTORIES(
    Thyra_Xpetra_Stokes_HDF5
)
ENDIF()
