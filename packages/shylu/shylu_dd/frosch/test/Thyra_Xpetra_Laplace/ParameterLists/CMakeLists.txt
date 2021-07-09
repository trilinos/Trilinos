TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraXpetraLaplaceCopyFiles
DEST_FILES ParameterList_GDSWPreconditioner.xml
           ParameterList_RGDSWPreconditioner.xml
           ParameterList_TwoLevelBlockPreconditioner_GDSW.xml
           ParameterList_TwoLevelPreconditioner_GDSW.xml
           ParameterList_TwoLevelPreconditioner_GDSW_muelu.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_ML.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW.xml
           ParameterList_TwoLevelPreconditioner_RGDSW.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS thyra_xpetra_laplace
)
