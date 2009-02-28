
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  ml_AztecMsr
  SOURCES ml_AztecMsr.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  ml_EpetraVbr
  SOURCES ml_EpetraVbr.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  ml_MatrixFree
  SOURCES ml_MatrixFree.cpp
  COMM serial mpi
  )
