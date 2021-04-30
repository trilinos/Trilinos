TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverFactoryCopyFiles
DEST_FILES ParameterList_Amesos_Klu.xml
           ParameterList_Amesos2_Klu.xml
           ParameterList_Belos_GMRES.xml
           ParameterList_FROSchPreconditioner_TwoLevelBlockPreconditioner.xml
           ParameterList_FROSchPreconditioner_TwoLevelPreconditioner.xml
           ParameterList_Ifpack2_RILUK.xml
           ParameterList_MueLu.xml
           ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml
           ParameterList_ThyraSolver_Belos_GMRES.xml
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
EXEDEPS solverfactory
)
