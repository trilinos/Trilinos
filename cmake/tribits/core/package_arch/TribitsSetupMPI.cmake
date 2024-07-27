# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(AdvancedSet)
include(MultilineSet)
include(PrintNonemptyVar)
include(FindProgramPlus)


function(tribits_extract_base_dir FILE_PATH BASE_DIR)
  if (NOT ${BASE_DIR})
    get_filename_component( ${BASE_DIR} ${FILE_PATH} PATH )
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(${BASE_DIR})
    endif()
    set(${BASE_DIR} ${${BASE_DIR}} PARENT_SCOPE)
  endif()
endfunction()


macro(tribits_find_mpi_compiler LANG)
  if (${PROJECT_NAME}_ENABLE_${LANG})
    if (CMAKE_${LANG}_COMPILER)
      message(STATUS "Leaving current CMAKE_${LANG}_COMPILER="
        "${CMAKE_${LANG}_COMPILER} since it is already set!")
    else()
      find_program_plus(
        MPI_${LANG}_COMPILER
        NAMES ${ARGN}
        PATHS ${MPI_BIN_DIR_PATHS}
        )
      tribits_extract_base_dir(${MPI_${LANG}_COMPILER} MPI_BASE_DIR)
      message(STATUS "Setting CMAKE_${LANG}_COMPILER=\${MPI_${LANG}_COMPILER}")
      set(CMAKE_${LANG}_COMPILER "${MPI_${LANG}_COMPILER}"
        CACHE FILEPATH
        "${LANG} compiler overridden by MPI_${LANG}_COMPILER")
      print_var(CMAKE_${LANG}_COMPILER)
    endif()
  endif()
endmacro()


function(tribits_setup_mpi)

  #
  # A) Get the directory paths
  #

  multiline_set( DOC
    "Base directory for the MPI implementation under which"
    " the bin, include, and lib directories are found" )
  advanced_set( MPI_BASE_DIR "" CACHE PATH ${DOC} )
  print_nonempty_var(MPI_BASE_DIR)

  if (MPI_BASE_DIR)
    set(MPI_BIN_DIR_DEFAULT "${MPI_BASE_DIR}/bin")
  else()
    set(MPI_BIN_DIR_DEFAULT "")
  endif()
  multiline_set( DOC
    "Path to the bin directory where the MPI compiler"
    " and runtime executables are found" )
  advanced_set( MPI_BIN_DIR ${MPI_BIN_DIR_DEFAULT} CACHE PATH ${DOC} )
  print_nonempty_var(MPI_BIN_DIR)

  multiline_set( DOC
    "If set to 'ON', then the MPI compiler wrappers will be used."
    "  Set MPI_[C,CXX,Fortran]_COMPILER:FILEPATH=XXX to set compilers." )
  advanced_set( MPI_USE_COMPILER_WRAPPERS ON CACHE BOOL ${DOC} )
  print_var(MPI_USE_COMPILER_WRAPPERS)

  file(TO_CMAKE_PATH "$ENV{ProgramFiles}" PROGRAM_FILES)
  if(MPI_BIN_DIR)
    set(MPI_BIN_DIR_PATHS ${MPI_BIN_DIR})
  else()
    set(MPI_BIN_DIR_PATHS
      /usr/local/mpi/bin
      /usr/local/bin
      /usr/bin
      "${PROGRAM_FILES}/Microsoft HPC Pack 2008 SDK/Bin"
      "C:/Program Files/Microsoft HPC Pack 2008 SDK/Bin"
      "${PROGRAM_FILES}/MPICH/SDK/Bin"
      "${PROGRAM_FILES}/MPICH2/Bin"
      "C:/Program Files/MPICH/SDK/Bin"
      "C:/Program Files/MPICH2/Bin"
      )
  endif()

  #
  # B) Get the MPI compilers and/or just the raw include paths and libraries
  #

  if (MPI_USE_COMPILER_WRAPPERS)

    # B.1) Set up to use the MPI wrappers

    tribits_find_mpi_compiler(C mpicc)

    tribits_find_mpi_compiler(CXX  mpicxx mpic++ mpiCC)

    tribits_find_mpi_compiler(Fortran mpif90 mpif77)

  else()

    # B.2) Set up to use raw configure options

    advanced_set( MPI_COMPILE_FLAGS ""
      CACHE STRING
      "List of general compiler flags (excluding include directories)." )

    advanced_set( MPI_LINK_FLAGS ""
      CACHE STRING
      "Link Flags for MPI executables." )

    # NOTE: Test rest of the flags will be set up by the
    # FindTPLMPI.cmake module!

  endif()

  #
  # C) Get the MPI executable
  #

   find_program_plus( MPI_EXEC
    NAMES mpiexec mpirun
    PATHS ${MPI_BIN_DIR_PATHS}
    DOC "MPI executable used to run MPI programs"
    )
  mark_as_advanced(MPI_EXEC)

  if(MPI_EXEC)

    get_filename_component( MPI_EXEC_NAME "${MPI_EXEC}" PATH )

    if(MPI_EXEC_NAME STREQUAL mpiexec)
      set(MPI_EXEC_NUMPROCS_FLAG_DEFAULT -n)
    else()
      set(MPI_EXEC_NUMPROCS_FLAG_DEFAULT -np)
    endif()
    advanced_set( MPI_EXEC_NUMPROCS_FLAG
      ${MPI_EXEC_NUMPROCS_FLAG_DEFAULT}
      CACHE STRING
      "Flag setting the number of processors to use with MPI run command." )
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_nonempty_var(MPI_EXEC_NUMPROCS_FLAG)
    endif()

    advanced_set( MPI_EXEC_DEFAULT_NUMPROCS "4"
      CACHE STRING
      "The default number of processes to use when running MPI programs." )

    if ( "${MPI_EXEC_MAX_NUMPROCS_DEFAULT}" STREQUAL "" )
      set(MPI_EXEC_MAX_NUMPROCS_DEFAULT 4)
    endif()
    advanced_set( MPI_EXEC_MAX_NUMPROCS ${MPI_EXEC_MAX_NUMPROCS_DEFAULT}
      CACHE STRING
      "The maximum number of processes to use when running MPI programs.  Tests with more procs are excluded." )

    advanced_set( MPI_EXEC_PRE_NUMPROCS_FLAGS ""
      CACHE STRING
      "Extra command-line args to the MPI exec before num-procs args." )

    advanced_set( MPI_EXEC_POST_NUMPROCS_FLAGS ""
      CACHE STRING
      "Extra command-line args to the MPI exec after num-procs args." )

  endif()

  #message(FATAL_ERROR "Stopping!")

endfunction()

# 2009/01/23: rabartl: ToDo: Above: create util find_program_path_first(...)
# in order to implement looking in the input path first and not last
