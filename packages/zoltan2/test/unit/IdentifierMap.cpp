// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER

// TODO: doxygen comments

#include <vector>
#include <ostream>
#include <Zoltan2_IdentifierMap.hpp>
#include <Teuchos_GlobalMPISession.hpp>   // So we don't have to #ifdef HAVE_MPI
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#define nobjects 10000

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int nprocs = session.getNProc();
  int rank = session.getRank();

  int numLocalObjects = nobjects /  nprocs;
  int leftOver = nobjects % nprocs;

  if (rank < leftOver) numLocalObjects++;

  vector<int> localIDs(numLocalObjects);
  vector<long> globalIDs(numLocalObjects);

  long base = nobjects * rank;

  for (i=0; i < numLocalObjects; i++){
    globalIDs[i] = base + i;
    localIDs[i] = i;
  }

  Teuchos::RCP<const Comm<long>> comm = Teuchos::DefaultComm<long>::getComm();

  IdentifierMap<long, int> map(*comm, globalIDs, localIDs);

  std::cout << "PASS" << std::endl;
}
