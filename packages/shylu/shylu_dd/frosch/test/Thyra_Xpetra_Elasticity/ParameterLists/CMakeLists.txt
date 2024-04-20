TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraXpetraElasticityCopyFiles
DEST_FILES ParameterList_GDSWPreconditioner_ILU.xml
           ParameterList_GDSWPreconditioner_Klu.xml
           ParameterList_RGDSWPreconditioner_ILU.xml
           ParameterList_RGDSWPreconditioner_Klu.xml
	   ParameterList_RGDSWPreconditioner_Klu_Half.xml
           ParameterList_RGDSWPreconditioner_Klu_Tpetra.xml
           ParameterList_TwoLevelPreconditioner_GDSW_ILU.xml
           ParameterList_TwoLevelPreconditioner_GDSW_Klu.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_ILU.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu_ML.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_ILU.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_Klu.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_ILU.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_2D.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_3D.xml
           ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_ML.xml
           ParameterList_TwoLevelPreconditioner_RGDSW_ILU.xml
           ParameterList_TwoLevelPreconditioner_RGDSW_Klu.xml
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
EXEDEPS thyra_xpetra_elasticity
CATEGORIES BASIC
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraXpetraElasticityCopyFilesPerformance
DEST_FILES ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
EXEDEPS thyra_xpetra_elasticity
CATEGORIES PERFORMANCE
)
