# MatrixDriver

IF (Cthulhu_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(CMake/tpetra)
ENDIF()

IF (Cthulhu_ENABLE_Epetra)
  ADD_SUBDIRECTORY(CMake/epetra)
ENDIF()
