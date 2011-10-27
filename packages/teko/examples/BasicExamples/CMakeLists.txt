INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
   InverseFactoryOperator-Epetra
   SOURCES 
      inverse-factory-operator.cpp
   COMM mpi
   )

PACKAGE_ADD_EXECUTABLE(
   StratimikosSolve
   SOURCES 
      stratimikos-solve.cpp
   COMM mpi
   )

PACKAGE_COPY_FILES_TO_BINARY_DIR(BasicExamples
  SOURCE_FILES 
    simple_example.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  EXEDEPS InverseFactoryOperator-Epetra
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(BasicExamples-Strat
  SOURCE_FILES 
    strat_example.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  EXEDEPS StratimikosSolve
  )
