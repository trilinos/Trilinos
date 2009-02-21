SET(PyTrilinos_ALL "'Epetra'")

IF(PyTrilinos_ENABLE_Teuchos)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'Teuchos'")
ENDIF(PyTrilinos_ENABLE_Teuchos)

IF(PyTrilinos_ENABLE_Triutils)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'TriUtils'")
ENDIF(PyTrilinos_ENABLE_Triutils)

IF(PyTrilinos_ENABLE_EpetraExt)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'EpetraExt'")
ENDIF(PyTrilinos_ENABLE_EpetraExt)

IF(PyTrilinos_ENABLE_AztecOO)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'AztecOO'")
ENDIF(PyTrilinos_ENABLE_AztecOO)

IF(PyTrilinos_ENABLE_Galeri)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'Galeri'")
ENDIF(PyTrilinos_ENABLE_Galeri)

IF(PyTrilinos_ENABLE_Amesos)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'Amesos'")
ENDIF(PyTrilinos_ENABLE_Amesos)

IF(PyTrilinos_ENABLE_Ifpack)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'IFPACK'")
ENDIF(PyTrilinos_ENABLE_Ifpack)

IF(PyTrilinos_ENABLE_Komplex)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'Komplex'")
ENDIF(PyTrilinos_ENABLE_Komplex)

IF(PyTrilinos_ENABLE_Anasazi)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'Anasazi'")
ENDIF(PyTrilinos_ENABLE_Anasazi)

IF(PyTrilinos_ENABLE_Pliris)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'Pliris'")
ENDIF(PyTrilinos_ENABLE_Pliris)

IF(PyTrilinos_ENABLE_ML)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'ML'")
ENDIF(PyTrilinos_ENABLE_ML)

IF(PyTrilinos_ENABLE_NOX)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'NOX'")
ENDIF(PyTrilinos_ENABLE_NOX)

IF(PyTrilinos_ENABLE_LOCA)
  SET(PyTrilinos_ALL "${PyTrilinos_ALL},\n           'LOCA'")
ENDIF(PyTrilinos_ENABLE_LOCA)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/__init__.src
  ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
