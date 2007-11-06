INCLUDE(${CMAKE_MODULE_PATH}/Epetra_Add_Test.cmake)

# Syntax for adding a test:
#   EPETRA_ADD_TEST(<TEST_DIR_NAME> SOURCES <SOURCE1> <SOURCE2> ... ARGS <ARG1> <ARG2> ...)

# The default arguments include "-v", as specified in the macro EPETRA_ADD_TEST

EPETRA_ADD_TEST(petra_power_method SOURCES cxx_main.cpp ARGS 100)

EPETRA_ADD_TEST(verySimple SOURCES cxx_main.cpp)

EPETRA_ADD_TEST(UG_Ex1 SOURCES cxx_main.cpp power_method.cpp ARGS 100)

EPETRA_ADD_TEST(UG_ExMaps SOURCES cxx_main.cpp ARGS 100)

EPETRA_ADD_TEST(my_example SOURCES cxx_main.cpp ARGS -v)


# Below this line are subdirectories in epetra/example that are not exposed using autotools

# Linking error
#EPETRA_ADD_TEST(C_wrappers SOURCES c_main.c)

# Requires ifpack, aztecOO
#EPETRA_ADD_TEST(Fortran SOURCES f_main.F)

# Must be run on 3 MPI processes
#EPETRA_ADD_TEST(ImportExport SOURCES cxx_main.cpp)

# Requires the first argument to be a Harwell-Boeing matrix file
#EPETRA_ADD_TEST(InverseIteration SOURCES cxx_main.cpp)

# Works
#EPETRA_ADD_TEST(MapColoring SOURCES cxx_main.cpp ARGS 100 200)

# Includes/linking errors
#EPETRA_ADD_TEST(petra_nonlinear SOURCES cc_main.cc basis.cc ARGS 50)
