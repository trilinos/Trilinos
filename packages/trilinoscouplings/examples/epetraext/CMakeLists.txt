
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Isorropia)
IF (${PACKAGE_NAME}_ENABLE_Isorropia)

  PACKAGE_ADD_EXECUTABLE(
    EpetraExt_Isorropia_LPTrans_Ex
    SOURCES EpetraExt_Isorropia_LPTrans_Ex.cpp
    COMM mpi
  )

ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Amesos)

  PACKAGE_ADD_EXECUTABLE(
    epetraext_amesosbtf_graphtrans_example
    SOURCES EpetraExt_AmesosBTF_GraphTrans_Ex.cpp
    COMM serial
  )

  PACKAGE_ADD_EXECUTABLE(
    epetraext_amesosbtf_matrixtrans_example
    SOURCES EpetraExt_AmesosBTF_MatrixTrans_Ex.cpp
    COMM serial
  )

ENDIF()

ASSERT_DEFINED(EpetraExt_ENABLE_PETSC)
IF(EpetraExt_ENABLE_PETSC)
PACKAGE_ADD_EXECUTABLE(
  petsc
  SOURCES EpetraExt_petsc.cpp
  )
ENDIF()
