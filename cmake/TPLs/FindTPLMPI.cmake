INCLUDE(TPLDeclareLibraries)

IF (MPI_USE_COMPILER_WRAPPERS)

  # The Compiler wrappers take care of all of these so just set them to null
  GLOBAL_SET(TPL_MPI_INCLUDE_DIRS)
  GLOBAL_SET(TPL_MPI_LIBRARIES)
  GLOBAL_SET(TPL_MPI_LIBRARY_DIRS)

ELSE()

  # We must get the header include dirs and libraries from the user.
  TPL_DECLARE_LIBRARIES( MPI
    REQUIRED_HEADERS mpi.h
    REQUIRED_LIBS_NAMES mpi
    )

ENDIF()
