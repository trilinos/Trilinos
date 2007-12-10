INCLUDE(Trilinos_Add_Executable)

# Syntax for adding a test:
#   TRILINOS_ADD_EXECUTABLE(<TEST_DIR_NAME> SOURCES <SOURCE1> <SOURCE2> ... ARGS <ARG1> <ARG2> ...)

INCLUDE_DIRECTORIES(${Epetra_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_SOURCE_DIR}/src)

LINK_LIBRARIES(epetra ${TRILINOS_LAPACK_LIBRARY} ${TRILINOS_BLAS_LIBRARY})

IF(TRILINOS_ENABLE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	LINK_LIBRARIES(${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
ENDIF(TRILINOS_ENABLE_MPI)

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

