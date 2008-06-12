# $Header$

# Required modules
INCLUDE(Trilinos_Add_Executable)

# Setup include directories and library dependencies
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})

# Example definitions
TRILINOS_ADD_EXECUTABLE(petra_power_method EXAMPLE SOURCES cxx_main.cpp ARGS 100)
TRILINOS_ADD_EXECUTABLE(verySimple EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(UG_Ex1 EXAMPLE SOURCES cxx_main.cpp power_method.cpp ARGS 100)
TRILINOS_ADD_EXECUTABLE(UG_ExMaps EXAMPLE SOURCES cxx_main.cpp ARGS 100)
TRILINOS_ADD_EXECUTABLE(my_example EXAMPLE SOURCES cxx_main.cpp ARGS -v)

# Below this line are subdirectories in epetra/example that are not exposed using autotools

# Linking error
#TRILINOS_ADD_EXECUTABLE(C_wrappers EXAMPLE SOURCES c_main.c)

# Requires ifpack, aztecOO
#TRILINOS_ADD_EXECUTABLE(Fortran EXAMPLE SOURCES f_main.F)

# Must be run on 3 MPI processes
#TRILINOS_ADD_EXECUTABLE(ImportExport EXAMPLE SOURCES cxx_main.cpp)

# Requires the first argument to be a Harwell-Boeing matrix file
#TRILINOS_ADD_EXECUTABLE(InverseIteration EXAMPLE SOURCES cxx_main.cpp)

# Works
#TRILINOS_ADD_EXECUTABLE(MapColoring EXAMPLE SOURCES cxx_main.cpp ARGS 100 200)

# Includes/linking errors
#TRILINOS_ADD_EXECUTABLE(petra_nonlinear EXAMPLE SOURCES cc_main.cc basis.cc ARGS 51)

