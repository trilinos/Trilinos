PROJECT(epetra C CXX Fortran)

# Flag for turning on MPI support
OPTION(ENABLE_MPI 
       "Build with MPI for parallel codes."
       OFF
)

# Flag to enable tests
OPTION(ENABLE_TESTS
       "Build tests for this package."
       OFF
)

SET(CMAKE_MODULE_PATH "${epetra_SOURCE_DIR}/cmake")

# MPI related configuration
IF(ENABLE_MPI)
  FIND_PACKAGE(MPI)
  IF(DEFINED MPI_LIBRARY AND DEFINED MPI_INCLUDE_PATH)
    SET(HAVE_MPI TRUE)
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
    ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
# The following should be set when HAVE_MPI is set
    ADD_DEFINITIONS(-DEPETRA_MPI)
    FIND_PROGRAM(MPI_EXECUTABLE
      NAMES mpiexec mpirun
      PATHS /usr/bin /usr/local/bin /usr/local/mpi/bin 
      ${MPI_LIBRARY}/../bin
      "$ENV{ProgramFiles}/MPICH/SDK/Bin"
      "$ENV{ProgramFiles}/MPICH2/Bin"
      "C:/Program Files/MPICH/SDK/Bin"
      "${MPI_LIBRARY}/../Bin"
    )
    MARK_AS_ADVANCED(MPI_EXECUTABLE)
  ELSE(DEFINED MPI_LIBRARY AND DEFINED MPI_INCLUDE_PATH)
    MESSAGE( SEND_ERROR "MPI_LIBRARY and MPI_INCLUDE_PATH must be specified")
    SET(HAVE_MPI FALSE)
  ENDIF(DEFINED MPI_LIBRARY AND DEFINED MPI_INCLUDE_PATH)
ENDIF(ENABLE_MPI)

# Use secure versions of CRT methods in Windows 
# (addresses C4996 build errors in MSVS)
IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES)
ENDIF(WIN32)

ADD_SUBDIRECTORY(src)

IF(ENABLE_TESTS)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(test)
ENDIF(ENABLE_TESTS)

# Installation
INCLUDE(CPack)
