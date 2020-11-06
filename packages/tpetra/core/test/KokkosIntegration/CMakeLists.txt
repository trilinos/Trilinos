
IF (Tpetra_ENABLE_CUDA)
  TRIBITS_ADD_EXECUTABLE(
    KokkosUVM
    SOURCES KokkosUVM
    COMM serial
  )
  TRIBITS_ADD_EXECUTABLE(
    KokkosDeepcopy
    SOURCES KokkosDeepcopy
    COMM serial
  )
  TRIBITS_ADD_LIBRARY(
    ApiTest
    SOURCES ApiTest
    SHARED
    TESTONLY
  )
  TRIBITS_ADD_LIBRARY(
    Intercept
    SOURCES Intercept
    SHARED
    TESTONLY
  )
ENDIF()
