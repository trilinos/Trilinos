INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
   BuilderPreconditioner_Step1_Example
   SOURCES 
      example-test.cpp
   COMM mpi
   )
