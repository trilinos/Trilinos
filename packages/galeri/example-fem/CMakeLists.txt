INCLUDE(PackageAddExecutable)

PACKAGE_ADD_EXECUTABLE(
  Grid
  SOURCES Grid.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  Laplacian3D
  SOURCES Laplacian3D.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  AdvDiff2D
  SOURCES AdvDiff2D.cpp
  COMM serial mpi
  )


