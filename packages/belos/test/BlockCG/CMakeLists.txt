
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)
INCLUDE(Copy_Files_To_Binary_Dir)

 COPY_FILES_TO_BINARY_DIR(CopyTestBlockCgComplexFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
    SOURCE_FILES mhd1280b.cua mhd1280a.cua
    EXEDEPS bl_gmres_diag
    )
