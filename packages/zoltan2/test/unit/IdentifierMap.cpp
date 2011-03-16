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

#define APPGID_IS_NOT_GNO

#include <vector>
#include <ostream>
#include <Teuchos_GlobalMPISession.hpp>   // So we don't have to #ifdef HAVE_MPI
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Zoltan2_Partitioner.hpp>

#define nobjects 10000

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int nprocs = session.getNProc();
  int rank = session.getRank();

  int numLocalObjects = nobjects /  nprocs;
  int leftOver = nobjects % nprocs;

  if (rank < leftOver) numLocalObjects++;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  Teuchos::RCP< std::vector<long> > gids = Teuchos::rcp(new std::vector<long>(numLocalObjects));
  Teuchos::RCP< std::vector<int> > lids = Teuchos::rcp(new std::vector<int>(numLocalObjects));

  long base = nobjects * rank;

  for (int i=0; i < numLocalObjects; i++){
    (*gids)[i] = base + i;
    (*lids)[i] = i;
  }

  Z2::IdentifierMap<long, int> idmap(comm, gids, lids);

  std::cout << "PASS" << std::endl;
}
