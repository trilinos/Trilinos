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

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Zoltan2_Environment.hpp>

#include <iostream>
#include <vector>

static int call_zoltan2(const Teuchos::Comm<int>& callerComm)
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_MAX;

  int rank = callerComm.getRank();
  int size = callerComm.getSize();

  int fail=0, myFail=0;

  // We're mimicing the zoltan2 library here, and it will
  // work with a copy of the caller's communicator.
  Teuchos::RCP< Teuchos::Comm<int> > comm = callerComm.duplicate();
  Zoltan2::Environment env;

  // Just take whatever the default environment is.
  env.commitParameters();

  if ( !myFail && ((comm->getRank()!=callerComm.getRank()) ||
                   (comm->getSize()!=callerComm.getSize()))) {
    myFail = 1;
  }

  Teuchos::reduceAll<int, int>(callerComm, REDUCE_MAX, myFail, outArg(fail));
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

  Teuchos::RCP<Teuchos::Comm<int> > subComm;

  try{
    subComm = comm->createSubcommunicator(ranks);
  }
  catch (std::exception &e){
    myFail = 1;
  }

  if (!myFail && 
     ((subComm->getRank()!=newRank) || (subComm->getSize()!=newSize))){
    myFail = 1;
  }

  Teuchos::reduceAll<int, int>(*comm, REDUCE_MAX, myFail, outArg(fail));
  if (fail){
    if (!rank){
      std::cerr << "Failure in Zoltan2::getTeuchosMpiSubComm with list";
      std::cerr << std::endl;
    }
    return fail;
  }

  // Create a sub communicator by giving a color

  try{
    subComm = comm->split(subGroup, 0);
  }
  catch (std::exception &e){
    myFail = 1;
  }

  if (!myFail && 
     ((subComm->getRank()!=newRank) || (subComm->getSize()!=newSize))){
    myFail = 1;
  }

  Teuchos::reduceAll<int, int>(*comm, REDUCE_MAX, myFail, outArg(fail));
  if (fail){
    if (!rank){
      std::cerr << "Failure in Zoltan2::getTeuchosMpiSubComm with color";
      std::cerr << std::endl;
    }
    return fail;
  }

  return 0;
}

using namespace std;

int main(int argc, char *argv[])
{
  // Initialize MPI if necessary.
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int fail=0;
  // Get a communicator, it will be the right kind regardless of whether we
  // are using MPI or not.
  Teuchos::RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  fail = call_zoltan2(*comm);

  if (!comm->getRank()){
    if (!fail)
      std::cout << "PASS" << std::endl; 
    else
      std::cout << "FAIL" << std::endl; 
  }
  
  return fail;
}



