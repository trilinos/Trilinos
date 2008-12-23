
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

 PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCgComplexFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
    SOURCE_FILES mhd1280b.cua mhd1280a.cua
    EXEDEPS bl_gmres_diag
    )
