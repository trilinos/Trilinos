INCLUDE(CheckCXXSourceCompiles)

FUNCTION(TSQR_CHECK_MPI_COMM THING VARNAME)
  # Test whether mpi.h defines MPI_COMM_${THING} as an MPI
  # communicator.  We do this by attempting to compile C++ code that
  # uses this symbol as a communicator input of MPI_Comm_rank().
  set (_SYMBOL "MPI_COMM_${THING}")
  IF (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
    message (STATUS "Checking whether MPI defines the symbol ${_SYMBOL}")
  ENDIF()
  SET(SOURCE
  "
#include <mpi.h>

int main(int argc, char* argv[]) {
  int worldRank, otherRank;
  MPI_Init (&argc, &argv);
  // Sanity check; MPI_COMM_WORLD had better be there.
  MPI_Comm_rank (MPI_COMM_WORLD, &worldRank);
  MPI_Comm_rank (${_SYMBOL}, &otherRank);
  MPI_Finalize ();
  return 0;
}
"
  )
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()


FUNCTION(TSQR_CHECK_MPI_COMM_SUBSETS)
  # Test whether mpi.h defines MPI_COMM_NETWORK, MPI_COMM_NODE,
  # MPI_COMM_SOCKET, and / or MPI_COMM_CACHE.  These are communicators
  # which are subsets of MPI_COMM_WORLD.  MPI_COMM_NETWORK includes
  # exactly one MPI rank per node, MPI_COMM_NODE includes only the MPI
  # ranks on a node, MPI_COMM_SOCKET includes only the MPI ranks that
  # inhabit the same socket on a node, and MPI_COMM_CACHE includes
  # only the MPI ranks that share a common cache on a node.
  foreach (_THING "NETWORK" "NODE" "SOCKET" "CACHE")
    TSQR_CHECK_MPI_COMM ("${_THING}" "HAVE_MPI_COMM_${_THING}")
  endforeach ()
ENDFUNCTION()
