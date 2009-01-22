INCLUDE(PackageCopyFilesToBinaryDir)
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../base/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../base/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../support-Trilinos)



PACKAGE_ADD_EXECUTABLE_AND_TEST(
  beam
  SOURCES beam.cpp
  ARGS "-i beam_crs_aztecoo.i"
  COMM serial mpi
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(copyexamplebeamfeifiles
  SOURCE_FILES beam_crs_aztecoo.i
  )
