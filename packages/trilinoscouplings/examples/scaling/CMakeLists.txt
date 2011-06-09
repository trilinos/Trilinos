INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Amesos AND Trilinos_ENABLE_AztecOO AND
   Trilinos_ENABLE_Intrepid AND Trilinos_ENABLE_ML AND
   Trilinos_ENABLE_Pamgen)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson
    SOURCES example_Poisson.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_DivLSFEM
    SOURCES example_DivLSFEM.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_CurlLSFEM
    SOURCES example_CurlLSFEM.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_StabilizatedADR
    SOURCES example_StabilizedADR.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  PACKAGE_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
  SOURCE_FILES Poisson.xml
               ADR.xml
               CurlLSFEMin.xml
               CurlLSFEMblock_in_block.xml
               DivLSFEMin.xml
               ML_nonsym.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  ) 
ENDIF()

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Amesos AND Trilinos_ENABLE_AztecOO AND
   Trilinos_ENABLE_Intrepid AND Trilinos_ENABLE_ML AND
   Trilinos_ENABLE_STK)
  PACKAGE_ADD_EXECUTABLE(
    Example_Poisson_STK
    SOURCES example_Poisson_stk.cpp
    DEPLIBS epetra epetraext aztecoo amesos ml
    ADD_DIR_TO_NAME
    )
  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyMeshFiles
  SOURCE_FILES unit_cube_10int_hex.exo
               unit_cube_5int_tet.exo
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  ) 
ENDIF()
