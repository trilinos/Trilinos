
SET(HEADERS "")
SET(SOURCES "")


SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

# The performance test currently depends on KokkosCore for its timers.
# Otherwise, it does not use anything else in KokkosCore directly.
#
IF (Tpetra_INST_DOUBLE)
  # This test makes sense only for Scalar == double
  TRIBITS_ADD_EXECUTABLE(
    Performance-CGSolve
    SOURCES
      ${SOURCES}
    COMM mpi
  )
ENDIF()

#TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyXmlForTpetraPerfCgSolve
#  SOURCE_DIR ${Tpetra_SOURCE_DIR}/test/
#  SOURCE_FILES Tpetra_PerformanceTests.xml
#  EXEDEPS Performance-CGSolve
#)


