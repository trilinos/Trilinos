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

#define NUM_LOCAL_IDS 10

void testSolution( RCP<const Environment> &env,
      RCP<const Comm<int> > &comm, RCP<const IdentifierMap<User_t> > &idMap,
      int wdim, int numGlobalParts, int numLocalParts,
      ArrayRCP<myid_t> &gidArray, 
      ArrayView<Array<size_t> > reqPartIds,
      ArrayView<Array<float> > reqPartSizes,
      float **correctPartSizes )    // null if uniform part sizes
{
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  int fail=0, gfail=0;
  double epsilon = 10e-6;

  bool specifyPartSizes = false;
  for (int i=0; i < wdim; i++){
    if (reqPartIds[i].size() > 0){
      specifyPartSizes = true;
      break;
    }
  }

  ////////////////////////////////////////////////////////////
  // Instantiate a partitioning solution.  
  ////////////////////////////////////////////////////////////

  RCP<Zoltan2::PartitioningSolution<user_t> > solution;

  if (specifyPartSizes){
    try{
      solution = rcp(new Zoltan2::PartitioningSolution<user_t>(
        env,                // application environment info
        comm,               // problem communicator
        idMap,              // problem identifiers (global Ids, local Ids)
        wdim,
        partArrayList,      // For each weight dim, a list of part Ids
        partSizeList));      // For each weight dim, a list of part sizes
    }
    catch (std::exception &e){
      fail=1;
    }
  else
    try{
      solution = rcp(new Zoltan2::PartitioningSolution<user_t>(
        env,                // application environment info
        comm,               // problem communicator
        idMap,              // problem identifiers (global Ids, local Ids)
        wdim));
    }
    catch (std::exception &e){
      fail=100;
    }
  }

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "constructor call 1", 1);

  ////////////////////////////////////////////////////////////
  // Test the Solution queries that are used by algorithms
  ////////////////////////////////////////////////////////////

  if (solution->getGlobalNumberOfParts() != numGlobalParts)
    fail=2;

  if (!fail && solution->getLocalNumberOfParts() != numLocalParts)
    fail=3;

  bool oneToOneDist = (numGlobalParts == nprocs);

  if (oneToOneDist){
    if (!fail && !solution->oneToOnePartDistribution())
      fail=4;

    if (!fail && solution->getPartDistribution() != NULL)
      fail=5;

    if (!fail && solution->getProcDistribution() != NULL)
      fail=6;
  }
  else{
    if (!fail && solution->oneToOnePartDistribution())
      fail=7;

    if (!fail){
      int *partDist = solution->getPartDistribution();
      if (!partDist)
        fail = 8;
      else{
        if (rank==0){
          std::cout << "Owners of parts: ";
          for (int i=0; i <= numGlobalParts; i++)
            std::cout << partDist[i] << " ";
          std::cout << std::endl;
        }
      }
      size_t *procDist = solution->getProcDistribution();
      if (!procDist)
        fail = 9;
      else{
        if (rank==0){
          std::cout << "First part owned by each process: ";
          for (int i=0; i <= nprocs; i++)
            std::cout << procDist[i] << " ";
          std::cout << std::endl;
        }
      }
    }
  }

  for (int w=0; !fail && w < wdim; w++){
    bool uniformSizes = ((numGlobalParts == 1) || (correctPartSizes[w]==NULL));

    if (uniformSizes == solution->criteriaHasUniformPartSizes(w))
      fail=w*100 + 8;

    if (!fail){
      for (int partId=0; !fail && partId < numGlobalParts; partId++){
        float psize = solution->getCriteriaPartSize(w, partId);
        if ( psize < correctPartSizes[w][partId] - epsilon ||
             psize > correctPartSizes[w][partId] + epsilon )
          fail=w*100 + 9;
      }
    }
  }

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }

  ////////////////////////////////////////////////////////////
  // Test the Solution set method that is called by algorithms
  ////////////////////////////////////////////////////////////

  size_t *partAssignments = new size_t [NUM_LOCAL_IDS];
  for (int i=0; i < NUM_LOCAL_IDS; i++){
    partAssignments[i] = gidArray[i] % numGlobalParts;  // round robin
  }
  ArrayRCP<size_t> partList = arcp(partAssignments, 0, NUM_LOCAL_IDS);

  float *imbalances = new float [wdim];
  for (int i=0; i < wdim; i++){
    imbalances[i] = 1.0;    // for now just set it to perfect balance
  }
  ArrayRCP<float> metrics = arcp(imbalances, 0, wdim);

  try{
    solution->setParts(gidArray.view(0, NUM_LOCAL_IDS), partList, metrics); 
  }
  catch (std::exception &e){
    fail=10;
  }

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }

  ////////////////////////////////////////////////////////////
  // Finally, test the Solution get methods that may be called 
  // by users or migration functions.
  ////////////////////////////////////////////////////////////

  if (solution->getNumberOfIds() != NUM_LOCAL_IDS)
    fail = 11;

  if (!fail){
    const myid_t *gids = solution->getGlobalIdList();
    for (int i=0; !fail && i < NUM_LOCAL_IDS; i++){
      if (gids[i] != gidArray[i])
        fail = 12;
    }
  }

  if (!fail){
    const size_t *parts = solution->getPartList();
    for (int i=0; !fail && i < NUM_LOCAL_IDS; i++){
      if (parts[i] != gidArray[i] % numGlobalParts)
        fail = 13;
    }
  }

  if (!fail){
    const float *val = solution->getImbalance();
    for (int i=0; !fail && i < wdim; i++){
      if (val[i] != 1.0)
        fail = 14;
    }
  }

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }
}


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  int fail=0, gfail=0;

  /////////////
  // A problem parameter list
  Teuchos::ParameterList params;
  Teuchos::ParameterList &partitioningParams = params.sublist("partitioning");

  partitioningParams.set("num_global_parts", nprocs);
  partitioningParams.set("num_local_parts", 1);

  /////////////
  // An environment
  RCP<const Zoltan2::Environment> env = 
    rcp(new Zoltan2::Environment(params, comm);

  env->commitParameters();

  /////////////
  // A simple identifier map.

  myid_t *myGids = new myid_t [NUM_LOCAL_IDS];
  for (int i=0, x=rank*NUM_LOCAL_IDS; i < NUM_LOCAL_IDS; i++){
    myGids[i] = x++;
  }

  ArrayRCP<myid_t> gidArray(myGids, 0, NUM_LOCAL_IDS, true);

  RCP<const Zoltan2::IdentifierMap<user_t> > idMap = 
    rcp(new Zoltan2::IdentifierMap<user_t>(env, comm, gidArray)); 

  /////////////////////////////////////////////////////////////
  // First test: 
  //    Weight dimension one.
  //    One part per process.
  //    Some part sizes are double the others.

  int weightDim = 1;
  int numLocalParts = 1;
  int numGlobalParts = nprocs;

  Array<size_t> partIdArray(numLocalParts);
  Array<float> partSizeArray(numLocalParts);

  partIdArray[0] = rank;           // my part

  if (rank%2 ==0)
    partSizeArray[0] = 1.0;     // size of my part
  else
    partSizeArray[0] = 2.0;     // size of my part

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

  testSolution(env, comm, idMap, weightDim, numGlobalParts, numLocalParts,
   gidArray, partArrayList, partSizeList, &normalizedPartSizes);

  delete [] normalizedPartSizes;

  /////////////////////////////////////////////////////////////
  // First test: 
  //    Weight dimension one.
  //    Several parts per process.
  //    All part sizes are the same.

  int weightDim = 1;
  int numGlobalParts = nprocs * 5;
  int numLocalParts = 0;

  partitioningParams.set("num_global_parts", numGlobalParts);
  partitioningParams.set("num_local_parts", 0);

  env = rcp(new Zoltan2::Environment(params, comm);
  env->commitParameters();

  Array<size_t> emptyPartIdArray;
  Array<float> emptyPartSizeArray;

  float *dummyCorrectPartSizes = NULL;

  partArrayList[0] = emptyPartIdArray;
  partSizeList[0] = emptyPartSizeArray;

  testSolution(env, comm, idMap, weightDim, numGlobalParts, numLocalParts,
   gidAddary, partArrayList, partSizeList, &normalizedPartSizes);


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
