INCLUDE(AdvancedSet)
INCLUDE(MultilineSet)
INCLUDE(PrintNonemptyVar)
INCLUDE(FindProgramPlus)
INCLUDE(TPLDeclareLibraries)


FUNCTION(PACKAGE_ARCH_EXTRACT_BASE_DIR FILE_PATH BASE_DIR)
  IF (NOT ${BASE_DIR})
    GET_FILENAME_COMPONENT( ${BASE_DIR} ${FILE_PATH} PATH )
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${BASE_DIR})
    ENDIF()
    SET(${BASE_DIR} ${${BASE_DIR}} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


MACRO(PACKAGE_ARCH_FIND_MPI_COMPILER LANG)
  IF (${PROJECT_NAME}_ENABLE_${LANG})
    IF (CMAKE_${LANG}_COMPILER)
      MESSAGE(STATUS "Leaving current CMAKE_${LANG}_COMPILER="
        "${CMAKE_${LANG}_COMPILER} alone since it was alredy set!")
    ELSE()
      FIND_PROGRAM_PLUS(
        MPI_${LANG}_COMPILER
        NAMES ${ARGN}
        PATHS ${MPI_BIN_DIR_PATHS}

        )
      PACKAGE_ARCH_EXTRACT_BASE_DIR(${MPI_${LANG}_COMPILER} MPI_BASE_DIR)
      MESSAGE(STATUS "Setting CMAKE_${LANG}_COMPILER=\${MPI_${LANG}_COMPILER}")
      SET(CMAKE_${LANG}_COMPILER "${MPI_${LANG}_COMPILER}"
        CACHE FILEPATH
        "${LANG} compiler overridden by MPI_${LANG}_COMPILER")
      PRINT_VAR(MPI_${LANG}_COMPILER)
    ENDIF()
  ENDIF()
ENDMACRO()


FUNCTION(PACKAGE_ARCH_SETUP_MPI)

  #
  # A) Get the directorie paths
  #
  
  MULTILINE_SET( DOC
    "Base directory for the MPI implementation under which"
    " the bin, include, and lib directories are found" )
  ADVANCED_SET( MPI_BASE_DIR "" CACHE PATH ${DOC} )
  PRINT_NONEMPTY_VAR(MPI_BASE_DIR)

  IF (MPI_BASE_DIR)
    SET(MPI_BIN_DIR_DEFAULT "${MPI_BASE_DIR}/bin")
  ELSE()
    SET(MPI_BIN_DIR_DEFAULT "")
  ENDIF()
  MULTILINE_SET( DOC
    "Path to the bin directory where the MPI compiler"
    " and runtime executables are found" )
  ADVANCED_SET( MPI_BIN_DIR ${MPI_BIN_DIR_DEFAULT} CACHE PATH ${DOC} )
  PRINT_NONEMPTY_VAR(MPI_BIN_DIR)

  MULTILINE_SET( DOC
    "If set to 'ON', then the MPI compiler wrappers will be used."
    "  Set MPI_[C,CXX,Fortran]_COMPILER:FILEPATH=XXX to set compilers." )
  ADVANCED_SET( MPI_USE_COMPILER_WRAPPERS ON CACHE BOOL ${DOC} )
  PRINT_VAR(MPI_USE_COMPILER_WRAPPERS)

  SET(MPI_BIN_DIR_PATHS
    ${MPI_BIN_DIR}
    /usr/local/mpi/bin
    /usr/local/bin
    /usr/bin
    "$ENV{ProgramFiles}/MPICH/SDK/Bin"
    "$ENV{ProgramFiles}/MPICH2/Bin"
    "C:/Program Files/MPICH/SDK/Bin"
    )

  #
  # B) Get the MPI compilers and/or just the raw include paths and libraries
  #

  IF (MPI_USE_COMPILER_WRAPPERS)

    # B.1) Set up to use the MPI wrappers

    PACKAGE_ARCH_FIND_MPI_COMPILER(C mpicc)

    PACKAGE_ARCH_FIND_MPI_COMPILER(CXX mpic++ mpiCC mpicxx)

    PACKAGE_ARCH_FIND_MPI_COMPILER(Fortran mpif77)

  ELSE()

    # B.2) Set up to use raw configure options

    MESSAGE(FATAL_ERROR "ToDo: Implement support for setting raw MPI options with raw compilers!")
 
    TPL_DECLARE_LIBRARIES( MPI
      REQUIRED_HEADERS mpi.h
      REQUIRED_LIBS_NAMES "mpi"
      )

  ENDIF()

  #
  # C) Get the MPI executable
  #

   FIND_PROGRAM_PLUS( MPI_EXEC
    NAMES mpiexec mpirun
    PATHS ${MPI_BIN_DIR_PATHS}
    DOC "MPI executable used to run MPI programs"
    )
  MARK_AS_ADVANCED(MPI_EXEC)

  IF(MPI_EXEC)

    GET_FILENAME_COMPONENT( MPI_EXEC_NAME "${MPI_EXEC}" PATH )

    IF(MPI_EXEC_NAME STREQUAL mpiexec)
      SET(MPI_EXEC_NUMPROCS_FLAG_DEFAULT -n)
    ELSE()
      SET(MPI_EXEC_NUMPROCS_FLAG_DEFAULT -np)
    ENDIF()
    ADVANCED_SET( MPI_EXEC_NUMPROCS_FLAG
      ${MPI_EXEC_NUMPROCS_FLAG_DEFAULT}
      CACHE STRING
      "Flag setting the number of processors to use with MPI run command." )
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_NONEMPTY_VAR(MPI_EXEC_NUMPROCS_FLAG)
    ENDIF()

    ADVANCED_SET( MPI_EXEC_MAX_NUMPROCS "4"
      CACHE STRING
      "The maximum mumber of processes to use when running MPI programs." )

    ADVANCED_SET( MPI_EXEC_PRE_NUMPROCS_FLAGS ""
      CACHE STRING
      "Extra command-line args to the MPI exec before num-procs args." )

    ADVANCED_SET( MPI_EXEC_POST_NUMPROCS_FLAGS ""
      CACHE STRING
      "Extra command-line args to the MPI exec after num-procs args." )

  ENDIF()

  #MESSAGE(FATAL_ERROR "Stopping!")

ENDFUNCTION()

# 2009/01/23: rabartl: ToDo: Above: create util FIND_PROGRAM_PATH_FIRST(...) 
# in order to implement looking in the input path first and not last
