
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Zoltan2Core)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

IF (${PACKAGE_NAME}_ENABLE_Zoltan2Core)

  TRIBITS_ADD_EXECUTABLE(
    EEP_Tpetra_Isorropia_LPTrans_Ex
    SOURCES EEP_Tpetra_Isorropia_LPTrans_Ex.cpp
    IMPORTEDLIBS /Users/eeprude/install/gklib/lib/libGKlib.a -lGKlib.a
    ARGS "-L/Users/eeprude/install/gklib/lib/ -lGKlib"
    COMM mpi
  )

ENDIF()

ENDIF()
