INCLUDE(Assert_Defined)
INCLUDE(Global_Null_Set)
INCLUDE(Global_Set)
INCLUDE(Advanced_Set)


IF (TPL_MPI_INCLUDE_DIRS)

  # The MPI library has already been found or the user has specified
  # these manually.  In this case, just make sure that everything has been
  # specified correctly

  # Make sure the other variables have been defined
  ASSERT_DEFINED(TPL_MPI_LIBRARY_DIRS)
  ASSERT_DEFINED(TPL_MPI_LIBRARIES)

  # Verify that indeed we have found MPI!

  # ToDo: Implement!

ELSE()

  # Otherwise, we need to look for the MPI headers and libraries

  FIND_PACKAGE(MPI)

  IF(DEFINED MPI_LIBRARY AND DEFINED MPI_INCLUDE_PATH)

    ASSERT_DEFINED(MPI_INCLUDE_PATH)
    SET( TPL_MPI_INCLUDE_DIRS "${MPI_INCLUDE_PATH}"
      CACHE PATH "")
    MARK_AS_ADVANCED(TPL_MPI_INCLUDE_DIRS)

    ASSERT_DEFINED(MPI_LIBRARIES)
    SET( TPL_MPI_LIBRARIES "${MPI_LIBRARIES}"
      CACHE PATH "")
    MARK_AS_ADVANCED(TPL_MPI_LIBRARIES)
  
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
        SET(MPI_NUMPROCS_FLAG_DEFAULT -n)
      ELSE()
        SET(MPI_NUMPROCS_FLAG_DEFAULT -np)
      ENDIF()

      ADVANCED_SET(MPI_NUMPROCS_FLAG ${MPI_NUMPROCS_FLAG_DEFAULT}
        CACHE STRING
        "Flag setting the number of processors to use."
        )

      ADVANCED_SET(TRILINOS_MPI_GO ${MPI_EXECUTABLE} ${MPI_NUMPROCS_FLAG}
        " " CACHE STRING "The actual command uses to run MPI jobs" )
  
    ENDIF()
  
    MARK_AS_ADVANCED(MPI_EXECUTABLE)
  
  ELSE()
  
    # Did not find MPI
    MESSAGE( FATAL_ERROR
      "TPL_MPI_[INCLUDE_DIRS,LIBRARIES,LIBARY_DIRS] must be specified!")
    SET(HAVE_MPI FALSE)
  
  ENDIF()

  # Just needs to be set for consistency
  GLOBAL_NULL_SET(TPL_MPI_LIBRARY_DIRS)

ENDIF()
