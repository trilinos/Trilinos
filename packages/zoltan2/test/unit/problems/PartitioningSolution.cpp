// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Test the PartitioningSolution class.
//
// Normally a Solution is created by a Problem.  
// We create a few Solutions in this unit test.

// TODO test ::convertSolutionToImportList

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_TestHelpers.hpp>

typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> user_t;
typedef Zoltan2::BasicIdentifierInput<user_t> idInput_t;

using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::arcp;

double epsilon = 10e-6;

typedef zoltan2_partId_t partId_t;

void makeArrays(int wdim, int *lens, partId_t **ids, scalar_t **sizes,
  ArrayRCP<ArrayRCP<partId_t> > &idList, ArrayRCP<ArrayRCP<scalar_t> > &sizeList)
{
  ArrayRCP<partId_t> *idArrays = new ArrayRCP<partId_t> [wdim];
  ArrayRCP<scalar_t> *sizeArrays = new ArrayRCP<scalar_t> [wdim];

  for (int w=0; w < wdim; w++){
    idArrays[w] = arcp(ids[w], 0, lens[w], false);
    sizeArrays[w] = arcp(sizes[w], 0, lens[w], false);
  }

  idList = arcp(idArrays, 0, wdim, true);
  sizeList = arcp(sizeArrays, 0, wdim, true);
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  int fail=0, gfail=0;

  ////////////////
  // Arrays to hold part Ids and part Sizes for each weight dimension.

  int numIdsPerProc = 10;
  int maxWeightDim = 3;
  int maxNumPartSizes = nprocs;
  int *lengths = new int [maxWeightDim];
  partId_t **idLists = new partId_t * [maxWeightDim];
  scalar_t **sizeLists = new scalar_t * [maxWeightDim];

  for (int w=0; w < maxWeightDim; w++){
    idLists[w] = new partId_t [maxNumPartSizes];
    sizeLists[w] = new scalar_t [maxNumPartSizes];
  }

  /////////////
  // A default environment
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  /////////////
  // A simple identifier map.

  gno_t *myGids = new gno_t [numIdsPerProc];
  for (int i=0, x=rank*numIdsPerProc; i < numIdsPerProc; i++){
    myGids[i] = x++;
  }

  ArrayRCP<const gno_t> gidArray(myGids, 0, numIdsPerProc, true);

  RCP<const Zoltan2::IdentifierMap<user_t> > idMap = 
    rcp(new Zoltan2::IdentifierMap<user_t>(env, comm, gidArray)); 

  /////////////
  // TEST:
  // One weight dimension, one part per proc.
  // Some part sizes are 2 and some are 1.

  int numGlobalParts = nprocs;
  int weightDim = 1;

  ArrayRCP<ArrayRCP<partId_t> > ids;
  ArrayRCP<ArrayRCP<scalar_t> > sizes;

  memset(lengths, 0, sizeof(int) * maxWeightDim);

  lengths[0] = 1;                    // We give a size for 1 part.
  idLists[0][0] = rank;              // The part is my part.
  sizeLists[0][0] = rank%2 + 1.0;    // The size is 1.0 or 2.0

  makeArrays(1, lengths, idLists, sizeLists, ids, sizes);

  // Normalized part size for every part, for checking later on

  scalar_t *normalizedPartSizes = new scalar_t [numGlobalParts];
  scalar_t sumSizes=0;
  for (int i=0; i < numGlobalParts; i++){
    normalizedPartSizes[i] = 1.0;
    if (i % 2) normalizedPartSizes[i] = 2.0;
    sumSizes += normalizedPartSizes[i];
  }
  for (int i=0; i < numGlobalParts; i++)
    normalizedPartSizes[i] /= sumSizes;

  /////////////
  // Create a solution object with part size information, and check it.

  RCP<Zoltan2::PartitioningSolution<idInput_t> > solution;

  try{
    solution = rcp(new Zoltan2::PartitioningSolution<idInput_t>(
      env,                // application environment info
      comm,               // problem communicator
      idMap,              // problem identifiers (global Ids, local Ids)
      weightDim,                  // weight dimension
      ids.view(0,weightDim),      // part ids
      sizes.view(0,weightDim))); // part sizes
  }
  catch (std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "constructor call 1", 1);

  // Test the Solution queries that are used by algorithms

  if (solution->getTargetGlobalNumberOfParts() != size_t(numGlobalParts))
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
        ((nprocs>1 && solution->criteriaHasUniformPartSizes(0)) ||
         (nprocs==1 && !solution->criteriaHasUniformPartSizes(0))) )
    fail=8;

  if (!fail){
    for (int partId=0; !fail && partId < numGlobalParts; partId++){
      scalar_t psize = solution->getCriteriaPartSize(0, partId);

      if ( psize < normalizedPartSizes[partId] - epsilon ||
           psize > normalizedPartSizes[partId] + epsilon )
        fail=9;
    }
  }

  delete [] normalizedPartSizes;

  gfail = globalFail(comm, fail);
  if (gfail){
    printFailureCode(comm, fail);   // exits after printing "FAIL"
  }

  // Test the Solution set method that is called by algorithms

  partId_t *partAssignments = new partId_t [numIdsPerProc];
  for (int i=0; i < numIdsPerProc; i++){
    partAssignments[i] = myGids[i] % numGlobalParts;  // round robin
  }
  ArrayRCP<partId_t> partList = arcp(partAssignments, 0, numIdsPerProc);

  try{
    solution->setParts(gidArray, partList);
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

  if (solution->getLocalNumberOfIds() != size_t(numIdsPerProc))
    fail = 11;

  if (!fail){
    const gno_t *gids = solution->getIdList();
    for (int i=0; !fail && i < numIdsPerProc; i++){
      if (gids[i] != myGids[i])
        fail = 12;
    }
  }

  if (!fail){
    const partId_t *parts = solution->getPartList();
    for (int i=0; !fail && i < numIdsPerProc; i++){
      if (parts[i] != myGids[i] % numGlobalParts)
        fail = 13;
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

  delete [] lengths;
  for (int w=0; w < maxWeightDim; w++){
    delete [] idLists[w];
    delete [] sizeLists[w];
  }
  delete [] idLists;
  delete [] sizeLists;
}
