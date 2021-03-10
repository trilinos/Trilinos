# Run tests in parallel iff MPI enabled 
IF (MYAPP_MPI)
  set(NUM_MPI_PROCS 4)
  set(MPIGO mpiexec -np ${NUM_MPI_PROCS})
ELSE()
  set(MPIGO "")
ENDIF()

function(add_app_test)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.xml
    ${CMAKE_CURRENT_BINARY_DIR}/input.xml COPYONLY)
  get_filename_component(testName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  add_test(${testName} ${MPIGO} ${MyApp_BINARY_DIR}/src/MyApp.exe)
  if (MYAPP_MPI)
    set_tests_properties(${testName} PROPERTIES PROCESSORS ${NUM_MPI_PROCS})
  endif()
endfunction()

# Add test problem subdirectories. Can be conditional on build
add_subdirectory(ex1)
IF(MYAPP_EPETRA)
  add_subdirectory(ex2)
ENDIF()

