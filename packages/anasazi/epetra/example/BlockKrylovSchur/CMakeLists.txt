INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Epetra_BKS_example
  SOURCES BlockKrylovSchurEpetraEx.cpp
  COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE(
  Epetra_BKS_SVD_example
  SOURCES BlockKrylovSchurEpetraExSVD.cpp
  COMM serial mpi
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BKS_File_example
    SOURCES BlockKrylovSchurEpetraExFile.cpp
    COMM serial mpi
  )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)


IF (${PACKAGE_NAME}_ENABLE_Amesos)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BKS_Amesos_example
    SOURCES BlockKrylovSchurEpetraExGenAmesos.cpp
    COMM serial mpi
    )

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BKS_Buckling_example
    SOURCES BlockKrylovSchurEpetraExBuckling.cpp
    COMM serial mpi
    )

ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)

ENDIF(${PACKAGE_NAME}_ENABLE_Amesos)

IF (${PACKAGE_NAME}_ENABLE_AztecOO)
IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BKS_AztecOO_example
    SOURCES BlockKrylovSchurEpetraExGenAztecOO.cpp
    COMM serial mpi
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)
ENDIF(${PACKAGE_NAME}_ENABLE_AztecOO)

IF (${PACKAGE_NAME}_ENABLE_Belos)
IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BKS_Belos_example
    SOURCES BlockKrylovSchurEpetraExGenBelos.cpp
    COMM serial mpi
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)
ENDIF(${PACKAGE_NAME}_ENABLE_Belos)
