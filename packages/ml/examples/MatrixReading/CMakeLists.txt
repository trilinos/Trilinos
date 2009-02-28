
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  ml_read_matrix_AIJ
  SOURCES ml_read_matrix_AIJ.cpp
  COMM serial mpi
  )
