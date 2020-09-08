
IF (Tpetra_ENABLE_CUDA)
  TRIBITS_ADD_EXECUTABLE(
    KokkosUVM
    SOURCES KokkosUVM
    COMM serial
  )
ENDIF()
