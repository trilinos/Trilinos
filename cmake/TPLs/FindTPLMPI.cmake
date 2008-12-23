INCLUDE(AssertDefined)
INCLUDE(GlobalNullSet)
INCLUDE(GlobalSet)
INCLUDE(AdvancedSet)


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

  #
  # Otherwise, we need to look for the MPI headers and libraries
  #

  #
  # A) Use the CMake module to look for MPI
  #

  # Set up some default paths that FIND_[PROGRAM,LIBRARY,PATH] will look in.
  
  ADVANCED_SET(MPI_BASE_DIR "" CACHE PATH
    "Base directory for MPI implementation under which ./bin, ./include, and ./lib are found")

  SET(TMP_CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH})
  SET(TMP_CMAKE_PROGRAM_PATH ${CMAKE_PROGRAM_PATH})
  SET(TMP_CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH})

  IF (MPI_BASE_DIR)
    SET(ENV{CMAKE_INCLUDE_PATH} # Used by FIND_PATH(...)
      "${MPI_BASE_DIR}/include:${MPI_BASE_DIR}/bin:$ENV{CMAKE_INCLUDE_PATH}" )
    SET(ENV{CMAKE_PROGRAM_PATH} # Used by FIND_PROGRAM(...)
      "${MPI_BASE_DIR}/include:${MPI_BASE_DIR}/bin:$ENV{CMAKE_PROGRAM_PATH}" )
    SET(ENV{CMAKE_LIBRARY_PATH} # Used by FIND_LIBRARY(...)
      "${MPI_BASE_DIR}/include:${MPI_BASE_DIR}/bin:$ENV{CMAKE_LIBRARY_PATH}" )
  ENDIF()

  # Use the CMake module for finding MPI

  FIND_PACKAGE(MPI)
 
  #
  # B) Put together TPL_MPI_LIBRARIES, TPL_MPI_INCLUDE_DIRS
  #

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

  # Set back the path environment variables

  SET(ENV{CMAKE_INCLUDE_PATH} "${TMP_CMAKE_INCLUDE_PATH}")
  SET(ENV{CMAKE_PROGRAM_PATH} "${TMP_CMAKE_PROGRAM_PATH}")
  SET(ENV{CMAKE_LIBRARY_PATH} "${TMP_CMAKE_LIBRARY_PATH}")

ENDIF()
