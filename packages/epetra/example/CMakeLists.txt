# $Header$

# Required modules
#INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)

# Example definitions

ADD_SUBDIRECTORY(petra_power_method)
ADD_SUBDIRECTORY(verySimple)
ADD_SUBDIRECTORY(UG_Ex1)
ADD_SUBDIRECTORY(UG_ExMaps)
ADD_SUBDIRECTORY(my_example)

# Below this line are subdirectories in epetra/example that are not exposed using autotools

# Linking error
#PACKAGE_ADD_EXECUTABLE_AND_TEST(C_wrappers EXAMPLE SOURCES c_main.c)
#PACKAGE_ADD_EXECUTABLE(
  #C_wrappers
  #COMM mpi serial  
  #SOURCES c_main.c
  #DIRECTORY C_wrappers
  #)

# Requires ifpack, aztecOO
#PACKAGE_ADD_EXECUTABLE_AND_TEST(Fortran EXAMPLE SOURCES f_main.F)
#PACKAGE_ADD_EXECUTABLE(
  #Fortran
  #COMM mpi serial  
  #SOURCES f_main.F
  #DIRECTORY Fortran
  #)

# Must be run on 3 MPI processes
#PACKAGE_ADD_EXECUTABLE_AND_TEST(ImportExport EXAMPLE SOURCES cxx_main.cpp)
#PACKAGE_ADD_EXECUTABLE(
  #ImportExport
  #COMM mpi serial  
  #SOURCES cxx_main.cpp
  #DIRECTORY ImportExport
  #)

# Requires the first argument to be a Harwell-Boeing matrix file
#PACKAGE_ADD_EXECUTABLE_AND_TEST(InverseIteration EXAMPLE SOURCES cxx_main.cpp)
#PACKAGE_ADD_EXECUTABLE(
  #InverseIteration
  #COMM mpi serial  
  #SOURCES cxx_main.cpp
  #DIRECTORY InverseIteration
  #)

# Works
#PACKAGE_ADD_EXECUTABLE_AND_TEST(MapColoring EXAMPLE SOURCES cxx_main.cpp ARGS 100 200)
#PACKAGE_ADD_EXECUTABLE(
  #MapColoring
  #COMM mpi serial  
  #SOURCES cxx_main.cpp
  #DIRECTORY MapColoring
  #)

# Includes/linking errors
#PACKAGE_ADD_EXECUTABLE_AND_TEST(petra_nonlinear EXAMPLE SOURCES cc_main.cc basis.cc ARGS 51)
#PACKAGE_ADD_EXECUTABLE(
  #petra_nonlinear
  #COMM mpi serial  
  #SOURCES cxx_main.cpp
  #DIRECTORY petra_nonlinear
  #)
