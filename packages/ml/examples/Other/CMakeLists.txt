
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  ml_example2d
  SOURCES ml_example2d.c place_holder.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  ml_example3d
  SOURCES ml_example3d.c place_holder.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  ml_readfile
  SOURCES ml_readfile.c place_holder.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  ml_read_salsa
  SOURCES ml_read_salsa.c place_holder.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  ml_recirc
  SOURCES ml_recirc.c place_holder.cpp
  COMM serial mpi
  )
