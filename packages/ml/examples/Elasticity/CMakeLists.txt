
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  ml_example_elasticity
  SOURCES ml_example_elasticity.c place_holder.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  ml_read_elas
  SOURCES ml_read_elas.c place_holder.cpp
  COMM serial mpi
  )
