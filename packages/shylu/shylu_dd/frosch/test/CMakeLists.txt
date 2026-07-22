ADD_SUBDIRECTORIES(
    Import_Tpetra
    InterfacePartitionOfUnity
    InterfaceSets
    LocalPartitionOfUnityBasis
    SolverFactory
)

IF(ShyLU_DDFROSch_ENABLE_Belos AND Kokkos_ENABLE_SERIAL AND ShyLU_DDFROSch_ENABLE_Stratimikos AND ShyLU_DDFROSch_ENABLE_Thyra)
ADD_SUBDIRECTORIES(
    Gathering
    Overlap
    Thyra_Xpetra_Elasticity
    Thyra_Xpetra_Laplace
)
ENDIF()

