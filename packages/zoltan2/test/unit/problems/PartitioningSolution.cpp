// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER
//
// Test the PartitioningSolution class.
//
// Normally a Solution is created by a Problem.  
// We create a few Solutions in this unit test.

#include <Zoltan2_PartitioningSolution.hpp>
#include <ErrorHandlingForTests.hpp>

typedef long myid_t;
typedef int lno_t;
typedef float scalar_t;

typedef Zoltan2::BasicUserTypes<scalar_t, myid_t, lno_t, myid_t> user_t;

using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::arcp;

double epsilon = 10e-6;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  int fail=0, gfail=0;

  /////////////
  // A default environment
  RCP<const Zoltan2::Environment> env = Zoltan2::getDefaultEnvironment();

  /////////////
  // A simple identifier map.

  myid_t *myGids = new myid_t [10];
  for (int i=0, x=rank*10; i < 10; i++){
    myGids[i] = x++;
  }

  ArrayRCP<myid_t> gidArray(myGids, 0, 10, true);

  RCP<const Zoltan2::IdentifierMap<user_t> > idMap = 
    rcp(new Zoltan2::IdentifierMap<user_t>(env, comm, gidArray)); 

  /////////////
  // Part sizes: some parts are double the size of others.

  int weightDim = 1;     // the default as far as the Solution is concerned
  int numLocalParts = 1; // the default in the Environment's parameter list
  int numGlobalParts = nprocs; // also the default in parameter list
  float partSize = 1.0;
  if (rank % 2) partSize = 2.0;

  Array<size_t> partIdArray(numLocalParts);
  Array<float> partSizeArray(numLocalParts);

  partIdArray[0] = rank;           // my part
  partSizeArray[0] = partSize;     // size of my part

  // Normalized part size for every part, for checking later on

  float *normalizedPartSizes = new float [numGlobalParts];
  float sumSizes=0;
  for (int i=0; i < numGlobalParts; i++){
    normalizedPartSizes[i] = 1.0;
    if (i % 2) normalizedPartSizes[i] = 2.0;
    sumSizes += normalizedPartSizes[i];
  }
  for (int i=0; i < numGlobalParts; i++)
    normalizedPartSizes[i] /= sumSizes;

  // We need one partIdArray and one partSizeArray for each weight dimension.

  Array<Array<size_t> > partArrayList;
  partArrayList.push_back(partIdArray);

  Array<Array<float> > partSizeList;
  partSizeList.push_back(partSizeArray);

  /////////////
  // Create a solution object with part size information, and check it.

  RCP<Zoltan2::PartitioningSolution<user_t> > solution;

  try{
    solution = rcp(new Zoltan2::PartitioningSolution<user_t>(
      env,                // application environment info
      comm,               // problem communicator
      idMap,              // problem identifiers (global Ids, local Ids)
      weightDim,
      partArrayList,      // For each weight dim, a list of part Ids
      partSizeList));      // For each weight dim, a list of part sizes
  }
  catch (std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "constructor call 1", 1);

  // Test the Solution queries that are used by algorithms

  if (solution->getGlobalNumberOfParts() != numGlobalParts)
    fail=2;

  if (!fail && solution->getLocalNumberOfParts() != 1)
    fail=3;

  if (!fail && !solution->oneToOnePartDistribution())
    fail=4;

  if (!fail && solution->getPartDistribution() != NULL)
    fail=5;

  if (!fail && solution->getProcDistribution() != NULL)
    fail=6;
      
  if (!fail && 
        (nprocs>1 && solution->criteriaHasUniformPartSizes(0) ||
         nprocs==1 && !solution->criteriaHasUniformPartSizes(0)) )
    fail=8;

  if (!fail){
    for (int partId=0; !fail && partId < numGlobalParts; partId++){
      float psize = solution->getCriteriaPartSize(0, partId);

      if ( psize < normalizedPartSizes[partId] - epsilon ||
           psize > normalizedPartSizes[partId] + epsilon )
        fail=9;
    }
  }

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }

  // Test the Solution set method that is called by algorithms

  size_t *partAssignments = new size_t [10];
  for (int i=0; i < 10; i++){
    partAssignments[i] = myGids[i] % numGlobalParts;  // round robin
  }
  ArrayRCP<size_t> partList = arcp(partAssignments, 0, 10);

  float *imbalances = new float [weightDim];
  for (int i=0; i < weightDim; i++){
    imbalances[i] = 1.0;    // for now just set it to perfect balance
  }
  ArrayRCP<float> metrics = arcp(imbalances, 0, weightDim);

  try{
    solution->setParts(gidArray.view(0, 10), partList, metrics); 
  }
  catch (std::exception &e){
    fail=10;
  }

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }

  // Test the Solution get methods that may be called by users 
  // or migration functions.

  if (solution->getNumberOfIds() != 10)
    fail = 11;

  if (!fail){
    const myid_t *gids = solution->getGlobalIdList();
    for (int i=0; !fail && i < 10; i++){
      if (gids[i] != myGids[i])
        fail = 12;
    }
  }

  if (!fail){
    const size_t *parts = solution->getPartList();
    for (int i=0; !fail && i < 10; i++){
      if (parts[i] != myGids[i] % numGlobalParts)
        fail = 13;
    }
  }

  if (!fail){
    const float *val = solution->getImbalance();
    for (int i=0; !fail && i < weightDim; i++){
      if (val[i] != 1.0)
        fail = 14;
    }
  }

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }

  if (rank==0)
    std::cout << "PASS" << std::endl;
  
  ///////////////////////////////////////////////////////////////////
  //  TODO:  
  /////////////
  // Create a solution object without part size information, and check it.
  /////////////
  // Test multiple weights.
  /////////////
  // Test multiple parts per process.
  /////////////
  // Specify a list of parts of size 0.  (The rest should be uniform.)

}
