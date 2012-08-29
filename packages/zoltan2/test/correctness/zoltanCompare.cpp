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

/*! \file zoltanCompare.cpp
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

#include <Tpetra_MultiVector.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;

//
// A few of the RCB tests done by Zoltan in nightly testing.
//

#if 0
#define NUMTESTS 19

static int testNumProcs[NUMTESTS] = {
2,2,
3,3,3,3,
4,4,4,4,4,4,4,
5,
6,6,6,6,
8
};

static string testArgs[NUMTESTS*3] = {
"simple", "no", "no",
"vwgt2", "no", "no",

"vwgt", "no", "no",
"bug", "no", "no",
"drake", "no", "no",
"onedbug", "no", "no",

"ewgt", "no", "no", 
"grid20x19", "no", "no", 
"grid20x19", "yes", "no",
"grid20x19", "no", "yes",
"nograph", "no", "no", 
"simple", "no", "no", 
"simple", "yes", "no",

"brack2_3", "no", "no",

"hammond2", "no", "no",
"degenerateAA", "no", "no",
"degenerate", "no", "no",
"degenerate", "no", "yes",

"hammond", "no", "no"
};
#else
#define NUMTESTS 2
static int testNumProcs[NUMTESTS] = {4,4};

static string testArgs[NUMTESTS*3] = {
"grid20x19", "no", "no",
"simple", "no", "no"};
#endif

typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tMatrix_t;
typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Zoltan2::XpetraCrsMatrixInput<tMatrix_t> inputAdapter_t;

int runRCB(const RCP<const Comm<int> > &comm,
  string fname, bool average_cuts, bool rectilinear_blocks,
  int numGlobalParts)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  

  // Read this test data from the Zoltan(1) test directory.

  RCP<UserInputForTests> uinput;
  try{
    uinput = rcp(new UserInputForTests(zoltanTestDirectory, fname, comm, true));
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: UserInputForTests" << std::endl;
    return 1;
  }

  RCP<tMatrix_t> matrix;
  try{
    matrix = uinput->getTpetraCrsMatrix();
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: get matrix" << std::endl;
    return 1;
  }

  RCP<const tMatrix_t> matrixConst = rcp_const_cast<const tMatrix_t>(matrix);

  RCP<tMVector_t> coords;
  try{
   coords = uinput->getCoordinates();
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: get coordinates" << std::endl;
    return 1;
  }

  int coordDim = (coords.is_null() ? 0 : coords->getNumVectors());

  RCP<tMVector_t> weights;
  try{
   weights = uinput->getWeights();
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: get weights" << std::endl;
    return 1;
  }

  int weightDim = (weights.is_null() ? 0 : weights->getNumVectors());

  // Create an input adapter.

  RCP<inputAdapter_t> ia;

  try{
    ia = rcp(new inputAdapter_t(matrixConst, coordDim, weightDim));
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: input adapter" << std::endl;
    return 1;
  }

  for (int dim=0; dim < coordDim; dim++){
    ia->setRowCoordinates(dim, coords->getData(dim).getRawPtr(), 1);
  }

  for (int dim=0; dim < weightDim; dim++)
    ia->setRowWeights(dim, weights->getData(dim).getRawPtr(), 1);

 // Parameters

  Teuchos::ParameterList params;
  params.set("timer_output_stream" , "std::cout");
  //params.set("debug_level" , "verbose_detailed_status");

  params.set("algorithm", "rcb");
  params.set("partitioning_objective", "multicriteria_balance_total_maximum");
  if (rank == 0)
    std::cout << "algorithm = rcb" << std::endl;

  double tolerance = 1.1;
  params.set("imbalance_tolerance", tolerance );
  if (rank == 0)
    std::cout << "imbalance_tolerance = " << tolerance << std::endl;

  if (nprocs == 1){
    params.set("num_global_parts", numGlobalParts);
    std::cout << "num_global_parts = " << numGlobalParts << std::endl;
  }

  params.set("bisection_num_test_cuts", 1);
  if (rectilinear_blocks){
    params.set("rectilinear_blocks", "yes");
    if (rank == 0)
      std::cout << "rectilinear_blocks = yes" << std::endl;
  }
  if (average_cuts){
    params.set("average_cuts", "yes");
    if (rank == 0)
      std::cout << "average_cuts = yes" << std::endl;
  }

  if (rank == 0){
    std::cout << "coordinate dimension: " << coordDim << std::endl;
    std::cout << "weight dimension: " << weightDim << std::endl;
    if (weightDim > 1)
      std::cout << 
        "objective: multicriteria_balance_total_maximum (2-norm)" << std::endl;
  }

  // Create the problem.

  RCP<Zoltan2::PartitioningProblem<inputAdapter_t> > problem;
  try{
    problem = rcp(new Zoltan2::PartitioningProblem<inputAdapter_t>(
      ia.getRawPtr(), &params));
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: problem" << std::endl;
    return 1;
  }

  try{
    problem->solve();
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: solve" << std::endl;
    return 1;
  }

  if (rank == 0){
    problem->printMetrics(cout);
  }

  problem->printTimers();

  return 0;
}
  
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  bool ac = false, rb = false;

  if (getenv("DEBUGME")){
    std::cout << getpid() << std::endl;
    sleep(10);
  }

  int fail=0;

  if (argc > 1){
   
    Teuchos::CommandLineProcessor cmdp (false, false);
  
    string inputFile("none"), average_cuts("no"), rectilinear_blocks("no");
  
    cmdp.setOption("inputFile", &inputFile, 
      "root of file name: \"grid20x19\" for \"grid20x19_coord.mtx\"");
    cmdp.setOption("average_cuts", &average_cuts, 
      "yes or no");
    cmdp.setOption("rectilinear_blocks", &rectilinear_blocks, 
      "yes or no");
  
    try{
      cmdp.parse(argc, argv);
    }
    catch(...){
      if (rank == 0)
        std::cout << "FAIL: arguments" << std::endl;
      return 1;
    }
  
    if (inputFile == string("none"))
      return 0;
  
    if (average_cuts == string("yes"))
      ac = true;
    if (rectilinear_blocks == string("yes"))
      rb = true;

    fail = runRCB(comm, inputFile, ac, rb, nprocs);
  }
  else{         // do all the Zoltan tests
    int numRan = 0;
    for (int i=0,ii=0; i < NUMTESTS; i++, ii+=3){
      int numProcs = testNumProcs[i];
      if ((nprocs == 1) || (nprocs == numProcs)){
        numRan++;
        ac = (testArgs[ii+1] == string("yes"));
        rb = (testArgs[ii+2] == string("yes"));
        fail = runRCB(comm, testArgs[ii], ac, rb, numProcs);

        // AlltoAll hangs second time around on 3 or 5 procs.
        // On s861036 and on octopi.
        // Tried many re-writes of AlltoAll using both Teuchos
        // and MPI.
        // TODO
 
        if ((nprocs == 3) || (nprocs == 5))
          break;

      }
    }
    if (numRan == 0){
      fail = runRCB(comm, "grid20x19", "yes", "yes", nprocs);
    }
  }
  
  if (rank == 0 && !fail)
    std::cout << "PASS" << std::endl;
  
  return 0;
}

