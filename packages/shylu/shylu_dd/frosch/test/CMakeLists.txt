ADD_SUBDIRECTORIES(
    Import_Tpetra
    InterfacePartitionOfUnity
    InterfaceSets
    LocalPartitionOfUnityBasis
)

IF(ShyLU_DDFROSch_ENABLE_Belos AND Kokkos_ENABLE_SERIAL AND ShyLU_DDFROSch_ENABLE_Stratimikos AND ShyLU_DDFROSch_ENABLE_Thyra)
ADD_SUBDIRECTORIES(
    Gathering
    Overlap
    Thyra_Xpetra_Elasticity
    Thyra_Xpetra_Laplace
)
ENDIF()

IF(ShyLU_DDFROSch_ENABLE_Belos AND Kokkos_ENABLE_SERIAL AND ShyLU_DDFROSch_ENABLE_Stratimikos AND ShyLU_DDFROSch_ENABLE_Thyra AND HAVE_AMESOS2_EPETRA AND HAVE_SHYLU_DDFROSCH_EPETRAEXT AND ShyLU_DDFROSch_ENABLE_HDF5)
ADD_SUBDIRECTORIES(
    Thyra_Xpetra_Stokes_HDF5
)
ENDIF()
