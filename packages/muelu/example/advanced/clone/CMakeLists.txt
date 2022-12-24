INCLUDE(XpetraAddTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Clone example
#
# FIXME (mfh 11 Feb 2015) I am removing KokkosClassic::ThrustGPUNode,
# so this example won't build.  We need to test this example, though!
# We should figure out how to test it in a sensible way without
# needing ThrustGPUNode.
#
# IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND TpetraClassic_ENABLE_ThrustGPUNode)
#   TRIBITS_ADD_EXECUTABLE_AND_TEST(
#     Clone-Tpetra
#     SOURCES Clone-Tpetra.cpp
#     NUM_MPI_PROCS 1
#     COMM serial mpi
#     )
# ENDIF()
