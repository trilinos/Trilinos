

SET(CMAKE_EXECUTABLE_SUFFIX ".exe")


#
# Process the COMM arguments
#
# NOTE: The COMM array arguments is passed as ${ARGN}
#

FUNCTION( PACKAGE_PROCESS_COMM_ARGS  ADD_SERIAL_FEATURE  ADD_MPI_FEATURE )

  SET(COMM_ARRAY ${ARGN})

  IF (COMM_ARRAY)
    SET(${ADD_SERIAL_FEATURE} OFF PARENT_SCOPE)
    SET(${ADD_MPI_FEATURE} OFF PARENT_SCOPE)
    FOREACH(COMM ${COMM_ARRAY})
      IF(COMM STREQUAL "serial")
        SET(${ADD_SERIAL_FEATURE} ON PARENT_SCOPE)
      ELSEIF (COMM STREQUAL "mpi")
        SET(${ADD_MPI_FEATURE} ON PARENT_SCOPE)
      ELSE()
        MESSAGE(SEND_ERROR "Error, the COMM value '${COMM}' is not valid."
          "  Only 'mpi' and 'serial' are allowed.")
      ENDIF()
    ENDFOREACH()
  ELSE()
    SET(${ADD_MPI_FEATURE} ON PARENT_SCOPE)
    SET(${ADD_SERIAL_FEATURE} ON PARENT_SCOPE)
  ENDIF()

  IF (TPL_ENABLE_MPI)
    SET(${ADD_SERIAL_FEATURE} OFF PARENT_SCOPE)
  ELSE()
    SET(${ADD_MPI_FEATURE} OFF PARENT_SCOPE)
  ENDIF()

ENDFUNCTION()


FUNCTION( PACKAGE_CREATE_NAME_FROM_CURRENT_SOURCE_DIRECTORY DIRECTORY_NAME )
    SET(DIRECTORY_NAME "")
    #Get the unique part of the path for this test directory
    STRING(REGEX REPLACE ${PACKAGE_SOURCE_DIR} "" unique_dir_path ${CMAKE_CURRENT_SOURCE_DIR})
    
    #strip off the preceeding "/"
    STRING(LENGTH ${unique_dir_path} udp_length)
    MATH(EXPR last_index "${udp_length}-1")
    STRING(SUBSTRING ${unique_dir_path} 1 ${last_index} unique_dir_path)

    #Make the name acceptable for filesystems. This may need to be made compatible with windows
    #since they use a "\" instead of a "/" for directory delimiters. I'm not sure how this will
    #react if we encounter a directory name with a space in it.
    STRING(REGEX REPLACE "/" "_" DIRECTORY_NAME ${unique_dir_path})

    #PRINT_VAR(DIRECTORY_NAME)
    SET(DIRECTORY_NAME ${DIRECTORY_NAME} PARENT_SCOPE)
ENDFUNCTION()
