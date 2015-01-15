
SET(HEADERS "")
SET(SOURCES "")


SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

# The performance test currently depends on KokkosCore for its timers.
# Otherwise, it does not use anything else in KokkosCore directly.
#
# FIXME (mfh 15 Jan 2015) The test currently does not build for me if
# the Kokkos refactor version of Tpetra is disabled.
IF (Trilinos_ENABLE_KokkosCore AND Tpetra_ENABLE_Kokkos_Refactor)
  TRIBITS_ADD_EXECUTABLE(
    Performance-CGSolve
    SOURCES
      ${SOURCES}
    COMM mpi
  )

#  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyXmlForTpetraPerfCgSolve
#    SOURCE_DIR ${Tpetra_SOURCE_DIR}/test/
#    SOURCE_FILES Tpetra_PerformanceTests.xml
#    EXEDEPS Performance-CGSolve
#  )
ENDIF ()


