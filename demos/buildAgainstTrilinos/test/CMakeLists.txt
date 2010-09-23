# Run tests in parallel iff MPI enabled 
IF (MYAPP_MPI)
  set(MPIGO mpiexec -np 4)
ELSE()
  set(MPIGO  )
ENDIF()

set(MyTest.exe ${MPIGO} ${MyApp_BINARY_DIR}/src/MyApp.exe)

# Add test problem subdirectories. Can be conditional on build
add_subdirectory(ex1)
IF(MYAPP_EPETRA)
  add_subdirectory(ex2)
ENDIF()

