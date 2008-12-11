
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)
INCLUDE(Copy_Files_To_Binary_Dir)

#bl_gmres_diag
TRILINOS_ADD_EXECUTABLE(
  bl_gmres_diag
  SOURCES test_bl_gmres_diag.cpp createEpetraProblem.cpp
  )

IF(Trilinos_ENABLE_Triutils)
    
  # bl_gmres_hb
  TRILINOS_ADD_EXECUTABLE(
    bl_gmres_hb
    SOURCES test_bl_gmres_hb.cpp createEpetraProblem.cpp
    )
  # resolve_gmres_hb
  TRILINOS_ADD_EXECUTABLE(
    resolve_gmres_hb
    SOURCES test_resolve_gmres_hb.cpp createEpetraProblem.cpp
    )
  # pseudo_gmres_hb
  TRILINOS_ADD_EXECUTABLE(
    pseudo_gmres_hb
    SOURCES test_pseudo_gmres_hb.cpp createEpetraProblem.cpp
    )

  COPY_FILES_TO_BINARY_DIR(CopyTestBlockGmresFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/example/BlockGmres
    SOURCE_FILES orsirr1.hb fidap036.hb
    )

  IF(Trilinos_ENABLE_Ifpack)

    # bl_pgmres_hb
    TRILINOS_ADD_EXECUTABLE(
      bl_pgmres_hb
      SOURCES test_bl_pgmres_hb.cpp createEpetraProblem.cpp
      )
    # pseudo_pgmres_hb
    TRILINOS_ADD_EXECUTABLE(
      pseudo_pgmres_hb
      SOURCES test_pseudo_pgmres_hb.cpp createEpetraProblem.cpp
      )
    # bl_fgmres_hb
    TRILINOS_ADD_EXECUTABLE(
      bl_fgmres_hb
      SOURCES test_bl_fgmres_hb.cpp createEpetraProblem.cpp
      )

  ENDIF(Trilinos_ENABLE_Ifpack)

  IF(Teuchos_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    # bl_gmres_complex_hb
    TRILINOS_ADD_EXECUTABLE(
      bl_gmres_complex_hb
      SOURCES test_bl_gmres_complex_hb.cpp createEpetraProblem.cpp
      )

    COPY_FILES_TO_BINARY_DIR(CopyTestBlockGmresComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
      SOURCE_FILES mhd1280b.cua mhd1280a.cua
      )

    ## Bad belos/test/definition file - points to BlockCG instead of BlockGmres to get .hb file for first two complex tests.
    COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
      DEST_DIR ${CMAKE_BINARY_DIR}/packages/belos/test/BlockCG
      SOURCE_FILES mhd1280b.cua mhd1280a.cua
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(Trilinos_ENABLE_Triutils)
