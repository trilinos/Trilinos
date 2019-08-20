TRIBITS_COPY_FILES_TO_BINARY_DIR(OverlapCopyFiles
DEST_FILES ParameterList_OneLevelPreconditioner_Graph.xml
           ParameterList_OneLevelPreconditioner_Matrix.xml
           ParameterList_OneLevelPreconditioner_Old.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS overlap
)
