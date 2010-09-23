ADD_SUBDIRECTORY(step1)

INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
   BuildPreconditioner_Driver_Example
   SOURCES 
      example-driver.cpp
   COMM mpi
   )
