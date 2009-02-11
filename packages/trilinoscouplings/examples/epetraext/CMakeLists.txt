
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraExt_Isorropia_LPTrans_Ex
  SOURCES EpetraExt_Isorropia_LPTrans_Ex.cpp
  ARGS
  COMM mpi
  )
