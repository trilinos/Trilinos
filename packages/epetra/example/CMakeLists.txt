# $Header$

# Required modules
INCLUDE(Trilinos_Add_Executable_And_Test)

# Setup include directories and library dependencies
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})

# Example definitions
TRILINOS_ADD_EXECUTABLE_AND_TEST(petra_power_method EXAMPLE SOURCES cxx_main.cpp ARGS 100)
TRILINOS_ADD_EXECUTABLE_AND_TEST(verySimple EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(UG_Ex1 EXAMPLE SOURCES cxx_main.cpp power_method.cpp ARGS 100)
TRILINOS_ADD_EXECUTABLE_AND_TEST(UG_ExMaps EXAMPLE SOURCES cxx_main.cpp ARGS 100)
TRILINOS_ADD_EXECUTABLE_AND_TEST(my_example EXAMPLE SOURCES cxx_main.cpp ARGS -v)

# Below this line are subdirectories in epetra/example that are not exposed using autotools

# Linking error
#TRILINOS_ADD_EXECUTABLE_AND_TEST(C_wrappers EXAMPLE SOURCES c_main.c)

# Requires ifpack, aztecOO
#TRILINOS_ADD_EXECUTABLE_AND_TEST(Fortran EXAMPLE SOURCES f_main.F)

# Must be run on 3 MPI processes
#TRILINOS_ADD_EXECUTABLE_AND_TEST(ImportExport EXAMPLE SOURCES cxx_main.cpp)

# Requires the first argument to be a Harwell-Boeing matrix file
#TRILINOS_ADD_EXECUTABLE_AND_TEST(InverseIteration EXAMPLE SOURCES cxx_main.cpp)

# Works
#TRILINOS_ADD_EXECUTABLE_AND_TEST(MapColoring EXAMPLE SOURCES cxx_main.cpp ARGS 100 200)

# Includes/linking errors
#TRILINOS_ADD_EXECUTABLE_AND_TEST(petra_nonlinear EXAMPLE SOURCES cc_main.cc basis.cc ARGS 51)

