# ROL CMake Utilities
# Provides compatibility between TriBITS and standalone CMake builds

# Function to add executable and test, compatible with both TriBITS and standalone CMake
function(ROL_ADD_EXECUTABLE_AND_TEST TARGET_NAME)
  # Parse arguments
  set(options ADD_DIR_TO_NAME)
  set(oneValueArgs NUM_MPI_PROCS PASS_REGULAR_EXPRESSION COMM)
  set(multiValueArgs SOURCES ARGS)
  
  cmake_parse_arguments(ROLET "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  
  # Check if we're in TriBITS environment
  if(COMMAND TRIBITS_ADD_EXECUTABLE_AND_TEST)
    # Forward to TriBITS function with all original arguments
    tribits_add_executable_and_test(${TARGET_NAME} ${ARGN})
  else()
    # Standalone CMake implementation
    
    # Handle ADD_DIR_TO_NAME by prepending current directory name
    set(ACTUAL_TARGET_NAME ${TARGET_NAME})
    if(ROLET_ADD_DIR_TO_NAME)
      get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
      set(ACTUAL_TARGET_NAME "${DIR_NAME}_${TARGET_NAME}")
    endif()
    
    # Create the executable
    add_executable(${ACTUAL_TARGET_NAME} ${ROLET_SOURCES})
    
    # Link against ROL if available
    if(TARGET rol)
      target_link_libraries(${ACTUAL_TARGET_NAME} rol)
    endif()
    
    # Add the test
    add_test(NAME ${ACTUAL_TARGET_NAME} COMMAND ${ACTUAL_TARGET_NAME} ${ROLET_ARGS})
    
    # Set test properties
    if(ROLET_PASS_REGULAR_EXPRESSION)
      set_tests_properties(${ACTUAL_TARGET_NAME} PROPERTIES 
        PASS_REGULAR_EXPRESSION "${ROLET_PASS_REGULAR_EXPRESSION}")
    endif()
    
    # Handle MPI if specified (basic implementation)
    if(ROLET_NUM_MPI_PROCS AND ROLET_NUM_MPI_PROCS GREATER 1)
      find_package(MPI QUIET)
      if(MPI_FOUND)
        set_tests_properties(${ACTUAL_TARGET_NAME} PROPERTIES 
          PROCESSORS ${ROLET_NUM_MPI_PROCS})
      endif()
    endif()
  endif()
endfunction()

# Function to add executable only (no test)
function(ROL_ADD_EXECUTABLE TARGET_NAME)
  # Parse arguments
  set(options ADD_DIR_TO_NAME)
  set(oneValueArgs)
  set(multiValueArgs SOURCES)
  
  cmake_parse_arguments(ROLE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  
  # Check if we're in TriBITS environment
  if(COMMAND TRIBITS_ADD_EXECUTABLE)
    # Forward to TriBITS function with all original arguments
    tribits_add_executable(${TARGET_NAME} ${ARGN})
  else()
    # Standalone CMake implementation
    
    # Handle ADD_DIR_TO_NAME by prepending current directory name
    set(ACTUAL_TARGET_NAME ${TARGET_NAME})
    if(ROLE_ADD_DIR_TO_NAME)
      get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
      set(ACTUAL_TARGET_NAME "${DIR_NAME}_${TARGET_NAME}")
    endif()
    
    add_executable(${ACTUAL_TARGET_NAME} ${ROLE_SOURCES})
    
    # Link against ROL if available
    if(TARGET rol)
      target_link_libraries(${ACTUAL_TARGET_NAME} rol)
    endif()
  endif()
endfunction()

# Function to copy files to binary directory
function(ROL_COPY_FILES_TO_BINARY_DIR TARGET_NAME)
  # Parse arguments
  set(options)
  set(oneValueArgs SOURCE_DIR DEST_DIR)
  set(multiValueArgs SOURCE_FILES)
  
  cmake_parse_arguments(ROLCOPY "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  
  # Check if we're in TriBITS environment
  if(COMMAND TRIBITS_COPY_FILES_TO_BINARY_DIR)
    # Forward to TriBITS function with all original arguments
    tribits_copy_files_to_binary_dir(${TARGET_NAME} ${ARGN})
  else()
    # Standalone CMake implementation
    if(NOT ROLCOPY_SOURCE_DIR)
      set(ROLCOPY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    
    if(NOT ROLCOPY_DEST_DIR)
      set(ROLCOPY_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()
    
    foreach(FILE ${ROLCOPY_SOURCE_FILES})
      configure_file(
        ${ROLCOPY_SOURCE_DIR}/${FILE}
        ${ROLCOPY_DEST_DIR}/${FILE}
        COPYONLY
      )
    endforeach()
  endif()
endfunction()

# Function to include directories, compatible with both TriBITS and standalone CMake
function(ROL_INCLUDE_DIRECTORIES)
  # Check if we're in TriBITS environment
  if(COMMAND TRIBITS_INCLUDE_DIRECTORIES)
    # Forward to TriBITS function with all arguments
    tribits_include_directories(${ARGN})
  else()
    # Standalone CMake implementation
    cmake_parse_arguments(
      INC
      "REQUIRED_DURING_INSTALLATION_TESTING"
      ""
      ""
      ${ARGN}
    )
    include_directories(${INC_UNPARSED_ARGUMENTS})
  endif()
endfunction()