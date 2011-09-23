// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

// Test the conversion of MPI communicators to Teuchos::MPIComm objects.
// Test creation of sub-communicators.

#include <iostream>
#include <vector>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_Environment.hpp>

#ifdef HAVE_MPI
#include "mpi.h"
#include <Teuchos_DefaultMpiComm.hpp>
static int call_zoltan2(MPI_Comm callerComm, int rank, int size)
{
  int fail=0, myFail=0;
  // We're mimicing the zoltan2 library here, and it will
  // work with a copy of the caller's communicator.
  MPI_Comm comm;
  MPI_Comm_dup(callerComm,&comm);
  MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
  Zoltan2::Environment env;

  // Just take whatever the default environment is.
  env.commitParameters();

  // Get the Teuchos communicator that Zoltan will use.

  Teuchos::RCP<Teuchos::MpiComm<int> > tcomm; 

  try {
    tcomm = Zoltan2::getTeuchosMpiComm<int>(comm);
  }
  catch (std::exception &e){
    myFail = 1;
  }

  if ( !myFail && ((tcomm->getRank()!=rank) || (tcomm->getSize()!=size))){
    myFail = 1;
  }

  MPI_Allreduce(&myFail, &fail, 1, MPI_INT, MPI_MAX, callerComm);
  if (fail){
    if (!rank)
      std::cerr << "Failure in Zoltan2::getTeuchosMpiComm" << std::endl;
    return fail;
  }

  // Create a sub communicator by giving a member list
  int numSubGroups = 2;
  int subGroup = rank%numSubGroups;
  int newRank = 0;
  int newSize = 0;
  std::vector<int> ranks;

  for (int i=subGroup; i < size; i+=numSubGroups, newSize++){
    if (i == rank) newRank = newSize;
    ranks.push_back(i);
  }

  Teuchos::RCP<Teuchos::MpiComm<int> > subComm;

  try{
    subComm = Zoltan2::getTeuchosMpiSubComm(tcomm, ranks, env);
  }
  catch (std::exception &e){
    myFail = 1;
  }

  if (!myFail && 
     ((subComm->getRank()!=newRank) || (subComm->getSize()!=newSize))){
    myFail = 1;
  }

  Teuchos::reduceAll<int, int>(*tcomm, Teuchos::REDUCE_MAX, 1, &myFail, &fail);
  if (fail){
    if (!rank){
      std::cerr << "Failure in Zoltan2::getTeuchosMpiSubComm with list";
      std::cerr << std::endl;
    }
    return fail;
  }

  // Create a sub communicator by giving a color

  try{
    subComm = Zoltan2::getTeuchosMpiSubComm(tcomm, subGroup, env);
  }
  catch (std::exception &e){
    myFail = 1;
  }

  if (!myFail && 
     ((subComm->getRank()!=newRank) || (subComm->getSize()!=newSize))){
    myFail = 1;
  }

  Teuchos::reduceAll<int, int>(*tcomm, Teuchos::REDUCE_MAX, 1, &myFail, &fail);
  if (fail){
    if (!rank){
      std::cerr << "Failure in Zoltan2::getTeuchosMpiSubComm with color";
      std::cerr << std::endl;
    }
    return fail;
  }

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

  if (!rank){
    if (!fail)
      std::cout << "PASS" << std::endl; 
    else
      std::cout << "FAIL" << std::endl; 
  }

  MPI_Finalize();
  
#else
  std::cout << "NOTRUN" << std::endl; 
#endif

  return fail;
}



