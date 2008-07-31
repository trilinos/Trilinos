# $Header$

# 2008/07/31: rabartl: TODO: There are a bunch of cache variables like
# MPIEXEC, MPIEXEC_MAX_NUMPROCS, MPIEXEC_NUMPROC_FLAG, MPIEXEC_POSTFLAGS, and
# MPIEXEC_PREFLAGS that are getting defined by in the file MPIConfig.cmake
# (involed by FIND_PACKAGE(MPI)) that are not being used and are just sitting
# in the cache.
#
# Options:
#
# 1) Just reuse the cache variables defined by MPIConfig.cmake as needed and
# don't define new ones below.
#
# 2) Undefine all of the MPI cache variables being defined by MPIConfig.cmake


# Use CMake module to find MPI_LIBRARY and MPI_INCLUDE_PATH
FIND_PACKAGE(MPI)

# Check to see if MPI was found
IF(DEFINED MPI_LIBRARY AND DEFINED MPI_INCLUDE_PATH)

  # Found MPI, now set it up
  SET(HAVE_MPI TRUE)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
# The following should be set when HAVE_MPI is set
  ADD_DEFINITIONS(-DEPETRA_MPI)
  # 2008/06/09: rabartl: Above, EPETRA_MPI gets set up in Epetra_ConfigDefs.h.
  # It should not be set here and this top-level file should not refer to
  # a specific package like Epetra.

  # Find MPI executable (mpiexec or mpirun)
  FIND_PROGRAM(MPI_EXECUTABLE
    NAMES mpiexec mpirun
    PATHS /usr/bin /usr/local/bin /usr/local/mpi/bin 
    ${MPI_LIBRARY}/../bin
    "$ENV{ProgramFiles}/MPICH/SDK/Bin"
    "$ENV{ProgramFiles}/MPICH2/Bin"
    "C:/Program Files/MPICH/SDK/Bin"
    "${MPI_LIBRARY}/../Bin"
  )
  IF(MPI_EXECUTABLE)
    IF(${MPI_EXECUTABLE} MATCHES mpiexec)
      #The number of processors should not be hard coded
      SET(MPI_EXECUTABLE_FLAGS 
  	    -n 2
	    CACHE STRING
        "Flags for the MPI executable."
      )
      SET(MPI_NUMPROCS_FLAG 
  	    -n
	    CACHE STRING
        "Flag setting the number of processors to use."
      )
    ELSE()
      #The number of processors should not be hard coded
      SET(MPI_EXECUTABLE_FLAGS 
  	    -np 2
	    CACHE STRING
        "Flags for the MPI executable."
      )
      SET(MPI_NUMPROCS_FLAG 
  	    -np
	    CACHE STRING
        "Flag setting the number of processors to use."
      )
    ENDIF()

    # 2008/07/31: rabartl: TODO: We should consider appending the num-processor
    # flag to the MPI_EXECUABLE name (e.g. '/usr/local/mpi/bin/mpiexec -np ')
    # in order to be able to handle systems that need --num-procs=N
    # (e.g. 'mpirun --num-procs=').
    #
    # Perhaps we should use MPI_GO instead defined with all of the input
    # parameters and significant terminal whitespace (e.g. 'mpiexec -np ' and
    # 'yod -sz')

    MARK_AS_ADVANCED(MPI_EXECUTABLE_FLAGS)
  ENDIF(MPI_EXECUTABLE)
  MARK_AS_ADVANCED(MPI_EXECUTABLE)

ELSE(DEFINED MPI_LIBRARY AND DEFINED MPI_INCLUDE_PATH)

  # Did not find MPI
  MESSAGE( SEND_ERROR "MPI_LIBRARY and MPI_INCLUDE_PATH must be specified")
  SET(HAVE_MPI FALSE)

ENDIF(DEFINED MPI_LIBRARY AND DEFINED MPI_INCLUDE_PATH)
