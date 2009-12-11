ENABLE_TESTING()

INCLUDE(PyTrilinosMakeTest)

PACKAGE_COPY_FILES_TO_BINARY_DIR(test_testutil.py
  SOURCE_FILES testutil.src
  DEST_FILES   testutil.py)

IF(PyTrilinos_ENABLE_Teuchos)
  PyTrilinos_MAKE_TEST(testTeuchos              )
  PyTrilinos_MAKE_TEST(testTeuchos_Comm         )
  PyTrilinos_MAKE_TEST(testTeuchos_ParameterList)
  PyTrilinos_MAKE_TEST(testTeuchos_ScalarTraits )
ENDIF(PyTrilinos_ENABLE_Teuchos)

IF(PyTrilinos_ENABLE_Epetra)
  PyTrilinos_MAKE_TEST(testEpetra_BlockMap      )
  PyTrilinos_MAKE_TEST(testEpetra_Comm          )
  PyTrilinos_MAKE_TEST(testEpetra_CompObject    )
  PyTrilinos_MAKE_TEST(testEpetra_CrsGraph      )
  PyTrilinos_MAKE_TEST(testEpetra_CrsMatrix     )
  PyTrilinos_MAKE_TEST(testEpetra_FEVector      )
  PyTrilinos_MAKE_TEST(testEpetra_ImportExport  )
  PyTrilinos_MAKE_TEST(testEpetra_IntSerialDense)
  PyTrilinos_MAKE_TEST(testEpetra_IntVector     )
  PyTrilinos_MAKE_TEST(testEpetra_LinearProblem )
  PyTrilinos_MAKE_TEST(testEpetra_LocalMap      )
  PyTrilinos_MAKE_TEST(testEpetra_Map           )
  PyTrilinos_MAKE_TEST(testEpetra_MapColoring   )
  PyTrilinos_MAKE_TEST(testEpetra_MultiVector   )
  PyTrilinos_MAKE_TEST(testEpetra_Object        )
  PyTrilinos_MAKE_TEST(testEpetra_SerialDense   )
  PyTrilinos_MAKE_TEST(testEpetra_Util          )
  PyTrilinos_MAKE_TEST(testEpetra_Vector        )
ENDIF(PyTrilinos_ENABLE_Epetra)

IF(PyTrilinos_ENABLE_Triutils)
  PyTrilinos_MAKE_TEST(testTriUtils)
ENDIF(PyTrilinos_ENABLE_Triutils)

IF(PyTrilinos_ENABLE_EpetraExt)
  PyTrilinos_MAKE_TEST(testEpetraExt_Coloring)
  PyTrilinos_MAKE_TEST(testEpetraExt_Property)
ENDIF(PyTrilinos_ENABLE_EpetraExt)

IF(PyTrilinos_ENABLE_AztecOO)
  PyTrilinos_MAKE_TEST(testAztecOO)
ENDIF(PyTrilinos_ENABLE_AztecOO)

IF(PyTrilinos_ENABLE_Galeri)
  PyTrilinos_MAKE_TEST(testGaleri_CrsMatrices)
  PyTrilinos_MAKE_TEST(testGaleri_Maps       )
ENDIF(PyTrilinos_ENABLE_Galeri)

IF(PyTrilinos_ENABLE_Amesos)
  PyTrilinos_MAKE_TEST(testAmesos)
ENDIF(PyTrilinos_ENABLE_Amesos)

IF(PyTrilinos_ENABLE_Ifpack)
  PyTrilinos_MAKE_TEST(testIFPACK)
ENDIF(PyTrilinos_ENABLE_Ifpack)

IF(PyTrilinos_ENABLE_Komplex)
  PyTrilinos_MAKE_TEST(testKomplex_LinearProblem)
ENDIF(PyTrilinos_ENABLE_Komplex)

IF(PyTrilinos_ENABLE_Anasazi)
  PyTrilinos_MAKE_TEST(testAnasazi)
ENDIF(PyTrilinos_ENABLE_Anasazi)

IF(PyTrilinos_ENABLE_Pliris)
ENDIF(PyTrilinos_ENABLE_Pliris)

IF(PyTrilinos_ENABLE_ML)
  PyTrilinos_MAKE_TEST(testML_Space)
  PyTrilinos_MAKE_TEST(testML_CompObject)
  PyTrilinos_MAKE_TEST(testML_TimeObject)
  PyTrilinos_MAKE_TEST(testML_MultiVector)
  PyTrilinos_MAKE_TEST(testML_PyMatrix)
  PyTrilinos_MAKE_TEST(testML_Preconditioner)
ENDIF(PyTrilinos_ENABLE_ML)

IF(PyTrilinos_ENABLE_NOX)
  PyTrilinos_MAKE_TEST(testNOX_StatusTest)
ENDIF(PyTrilinos_ENABLE_NOX)

IF(PyTrilinos_ENABLE_LOCA)
ENDIF(PyTrilinos_ENABLE_LOCA)
