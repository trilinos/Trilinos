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

// Test the conversion of MPI communicators to Teuchos::MPIComm objects.
// Test creation of sub-communicators.

#include <iostream>
#include <vector>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_Environment.hpp>

#ifdef HAVE_MPI
#include "mpi.h"
#include <Teuchos_DefaultMpiComm.hpp>
static int call_zoltan2(MPI_Comm callerComm, int rank, int size)
{
  // The zoltan2 library will work with a copy of the caller's communicator.
  MPI_Comm comm;
  MPI_Comm_dup(callerComm,&comm);
  MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
  Zoltan2::Environment env;

  // env not really implemented yet - just set it up
  env._errorOStream = &std::cerr;
  env._errorCheckLevel = 0;

  // Get the Teuchos communicator that Zoltan will use.

  Teuchos::RCP<Teuchos::MpiComm<int> > tcomm = Zoltan2::getTeuchosMpiComm<int>(comm);

std::cout << tcomm.strong_count() << " " << tcomm.weak_count() << std::endl;
std::cout << tcomm.is_valid_ptr() << std::endl;

  if ( (tcomm->getRank() != rank) ||
       (tcomm->getSize() != size) ){
    std::cerr << "Bad tcomm" << std::endl;
    return 1;
  }

  tcomm->barrier();

  int origRank = rank;

  std::cout << "Proc " << origRank << " of " << size << std::endl;

  int subGroup = rank%2;
  std::vector<int> ranks(size);
  for (int i=subGroup, j=0; i < size; i+=2){
    ranks[j++] = i;
  }

  Teuchos::RCP<Teuchos::MpiComm<int> > subComm = 
    Zoltan2::getTeuchosMpiSubComm(tcomm, ranks, env);

  subComm->barrier();

  rank = subComm->getRank();
  size = subComm->getSize();

  std::cout << "Proc " << rank << " of " << size << std::endl;

  return 0;
}
#endif

using namespace std;

int main(int argc, char *argv[])
{
  int fail=0;
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  fail = call_zoltan2(MPI_COMM_WORLD, rank, size);
  
#else
  std::cout << "PASS" << std::endl; 
#endif
  return fail;
}



