
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  ml_viz
  SOURCES ml_viz.cpp
  COMM serial mpi
  )
