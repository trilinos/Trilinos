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
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <zoltan_cpp.h>

#include <Tpetra_MultiVector.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;

//
// A few of the RCB tests done by Zoltan in nightly testing.
//

#define NUMTESTS 22

static int testNumProcs[] = {
2,2,
3,3,3,3,3,3,
4,4,4,4,4,4,4,4,
5,
6,6,6,6,
8
};

static string testArgs[] = {
// Filename  AverageCuts  RectilinearBlocks
"simple",       "no",          "no",
"vwgt2",        "no",          "no",

"bug",          "no",          "no",
"drake",        "no",          "no",
"onedbug",      "no",          "no",
"simple",       "no",          "no",
"vwgt",         "no",          "no",
"vwgt2",        "no",          "no",

"ewgt",         "no",          "no", 
"grid20x19",    "no",          "no", 
"grid20x19",    "yes",         "no",
"grid20x19",    "no",          "yes",
"nograph",      "no",          "no", 
"simple",       "no",          "no", 
"simple",       "yes",         "no",
"vwgt2",        "no",          "no",

"brack2_3",     "no",          "no",

"hammond2",     "no",          "no",
"degenerateAA", "no",          "no",
"degenerate",   "no",          "no",
"degenerate",   "no",          "yes",

"hammond",      "no",          "no"
};

static string objectives[] = {
  "balance_object_count",
  "balance_object_weight",
  "multicriteria_balance_total_maximum"
};


typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> tMatrix_t;
typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> vectorAdapter_t;
typedef Zoltan2::XpetraCrsMatrixAdapter<tMatrix_t,tMVector_t> matrixAdapter_t;

////////////////////////////////////////////////////////////////////////////////
// Zoltan callbacks

template <typename MV>
int znumobj(void *data, int *ierr) 
{
  *ierr = ZOLTAN_OK;
  MV *vec = (MV *) data;
  return vec->getLocalLength();
}

template <typename MV, typename scalar_t>
void zobjlist(void *data, int ngid, int nlid, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
              int nwgts, float *wgts, int *ierr)
{
  *ierr = ZOLTAN_OK;
  MV *vec = (MV *) data;
  int n = vec->getLocalLength();
  for (int i = 0; i < n; i++) {
    gids[i] = vec->getMap()->getGlobalElement(i);
    lids[i] = i;
  }
  for (int w = 0; w < nwgts; w++) {
    ArrayRCP<const scalar_t> wvec = vec->getData(w);
    for (int i = 0; i < n; i++)
      wgts[i*nwgts+w] = wvec[i];
  }
}

template <typename MV>
int znumgeom(void *data, int *ierr) 
{
  *ierr = ZOLTAN_OK;
  MV *cvec = (MV *) data;
  return cvec->getNumVectors();
}

template <typename MV, typename scalar_t>
void zgeom(void *data, int ngid, int nlid, int nobj, 
          ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
          int ndim, double *coords, int *ierr)
{
  *ierr = ZOLTAN_OK;
  MV *vec = (MV *) data;
  for (int d = 0; d < ndim; d++) {
    ArrayRCP<const scalar_t> cvec = vec->getData(d);
    for (int i = 0; i < nobj; i++) {
      coords[lids[i]*ndim+d] = cvec[lids[i]];
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// Function to compute both Zoltan2 and Zoltan partitions and print metrics

int runRCB(
  const RCP<const Comm<int> > &comm,
  int numGlobalParts,
  int testCnt
)
{
  int me = comm->getRank();
  int np = comm->getSize();
  

  // Read this test data from the Zoltan(1) test directory.

  UserInputForTests *uinput;
  try{
    uinput = new UserInputForTests(zoltanTestDirectory, testArgs[testCnt*3],
                                   comm, true);
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: UserInputForTests " << e.what() << std::endl;
    return 1;
  }

  RCP<tMatrix_t> matrix;
  try{
    matrix = uinput->getUITpetraCrsMatrix();
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: get matrix " << e.what() << std::endl;
    return 1;
  }

  RCP<const tMatrix_t> matrixConst = rcp_const_cast<const tMatrix_t>(matrix);

  RCP<tMVector_t> coords;
  try{
   coords = uinput->getUICoordinates();
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: get coordinates " << e.what() << std::endl;
    return 1;
  }

  int coordDim = (coords.is_null() ? 0 : coords->getNumVectors());

  RCP<tMVector_t> weights;
  try{
   weights = uinput->getUIWeights();
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: get weights " << e.what() << std::endl;
    return 1;
  }

  int nWeights = (weights.is_null() ? 0 : weights->getNumVectors());

  // Create input adapters for the matrix and its coordinates

  matrixAdapter_t *ia;

  try{
    ia = new matrixAdapter_t(matrixConst, nWeights);
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: matrix adapter " << e.what() << std::endl;
    return 1;
  }

  for (int idx=0; idx < nWeights; idx++)
    ia->setRowWeights(weights->getData(idx).getRawPtr(), 1, idx);

  vectorAdapter_t *ca = NULL;

  try{
    ca = new vectorAdapter_t(coords);
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: vector adapter " << e.what() << std::endl;
    return 1;
  }

  ia->setCoordinateInput(ca);
  
 // Parameters

  Teuchos::ParameterList params;
  params.set("timer_output_stream" , "std::cout");
  params.set("compute_metrics", "true");
  // params.set("debug_level" , "verbose_detailed_status");

  params.set("algorithm", "rcb");
  params.set("partitioning_objective", objectives[(nWeights > 2 ? 2 : nWeights)]);

  double tolerance = 1.1;
  params.set("imbalance_tolerance", tolerance );
  params.set("num_global_parts", numGlobalParts);
  params.set("bisection_num_test_cuts", 1);
  params.set("average_cuts", testArgs[testCnt*3+1]);
  params.set("rectilinear", testArgs[testCnt*3+2]);

  if (me == 0) {
    std::cout << "Test " << testCnt << " filename            = "
              << testArgs[testCnt*3] << std::endl;
    std::cout << "Test " << testCnt << " num processors      = "
              << np << std::endl;
    std::cout << "Test " << testCnt << " algorithm           = rcb"
              << std::endl;
    std::cout << "Test " << testCnt << " num_global_parts    = "
              << numGlobalParts << std::endl;
    std::cout << "Test " << testCnt << " imbalance_tolerance = "
              << tolerance << std::endl;
    std::cout << "Test " << testCnt << " coordinate dim      = "
              << coordDim << std::endl;
    std::cout << "Test " << testCnt << " num weights per ID  = "
              << nWeights << std::endl;
    std::cout << "Test " << testCnt << " partition objective = "
              << objectives[(nWeights > 2 ? 2 : nWeights)] << std::endl;
    std::cout << "Test " << testCnt << " average_cuts        = "
              << testArgs[testCnt*3+1] << std::endl;
    std::cout << "Test " << testCnt << " rectilinear_blocks  = "
              << testArgs[testCnt*3+2] << std::endl;
  }

  // Create the problem.

  Zoltan2::PartitioningProblem<matrixAdapter_t> *problem;
#ifdef HAVE_ZOLTAN2_MPI
  // TPLs may want an MPI communicator

  const Teuchos::MpiComm<int> *tmpicomm =
                 dynamic_cast<const Teuchos::MpiComm<int> *>(comm.getRawPtr());
  MPI_Comm mpiComm = *(tmpicomm->getRawMpiComm());

  try{
    problem = new Zoltan2::PartitioningProblem<matrixAdapter_t>(ia, &params,mpiComm);
  }
#else
  try{
    problem = new Zoltan2::PartitioningProblem<matrixAdapter_t>(ia, &params);
  }
#endif
  catch(std::exception &e){
    if (me == 0)
      std::cout << "FAIL: problem " << e.what() << std::endl;
    return 1;
  }

  try{
    problem->solve();
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "FAIL: solve " << e.what() << std::endl;
    return 1;
  }

  // Now run the same partitioning using Zoltan
#ifdef HAVE_ZOLTAN2_MPI
  Zoltan zz(mpiComm);
#else
  Zoltan zz;
#endif
  char tmp[56];
  zz.Set_Param("LB_METHOD", "RCB");
  
  sprintf(tmp, "%d", numGlobalParts);
  zz.Set_Param("NUM_GLOBAL_PARTS", tmp);
  sprintf(tmp, "%d", nWeights);
  zz.Set_Param("OBJ_WEIGHT_DIM", tmp);
  sprintf(tmp, "%f", tolerance);
  zz.Set_Param("IMBALANCE_TOL", tmp);
  zz.Set_Param("RETURN_LISTS", "PART");
  zz.Set_Param("FINAL_OUTPUT", "1");
  zz.Set_Param("CHECK_GEOM", "0");
  if (testArgs[testCnt*3+1] == "yes") zz.Set_Param("AVERAGE_CUTS", "1");
  if (testArgs[testCnt*3+2] == "yes") zz.Set_Param("RCB_RECTILINEAR_BLOCKS", "1");

  zz.Set_Num_Obj_Fn(znumobj<tMVector_t>, (void *) coords.getRawPtr());
  if (nWeights)
    zz.Set_Obj_List_Fn(zobjlist<tMVector_t,zscalar_t>, (void *) weights.getRawPtr());
  else
    zz.Set_Obj_List_Fn(zobjlist<tMVector_t,zscalar_t>, (void *) coords.getRawPtr());
  zz.Set_Num_Geom_Fn(znumgeom<tMVector_t>, (void *) coords.getRawPtr());
  zz.Set_Geom_Multi_Fn(zgeom<tMVector_t,zscalar_t>, (void *) coords.getRawPtr());

  int changes, ngid, nlid;
  int numd, nump;
  ZOLTAN_ID_PTR dgid = NULL, dlid = NULL, pgid = NULL, plid = NULL;
  int *dproc = NULL, *dpart = NULL, *pproc = NULL, *ppart = NULL;
  zz.LB_Partition(changes, ngid, nlid, numd, dgid, dlid, dproc, dpart,
                                       nump, pgid, plid, pproc, ppart);
  zz.LB_Free_Part(&pgid, &plid, &pproc, &ppart);

  if (me == 0){
    problem->printMetrics(cout);
  }

  problem->printTimers();

  delete ia;
  delete ca;
  delete problem;
  delete uinput;

  return 0;
}
  
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int me = comm->getRank();
  int np = comm->getSize();

  int fail=0;

  Array<int> ranks(np);
  for (int i = 0; i < np; i++) ranks[i] = i;

  for (int i=0; i < NUMTESTS; i++) {
    int nTestProcs = testNumProcs[i];
    if (nTestProcs > np) {
      if (me == 0) {
        std::cout << "Skipping test " << i << " on " << testArgs[i*3]
                  << "; required number of procs " << nTestProcs 
                  << " is greater than available procs " << np << std::endl;
      }
      continue;
    }
    RCP<const Comm<int> > testcomm;
    if (nTestProcs == np)
      testcomm = comm;
    else
      testcomm = comm->createSubcommunicator(ranks.view(0,nTestProcs));

    if (me < nTestProcs) {
      fail = runRCB(testcomm, nTestProcs, i);

      // AlltoAll hangs second time around on 3 or 5 procs.
      // On s861036 and on octopi.
      // Tried many re-writes of AlltoAll using both Teuchos
      // and MPI.
      // TODO
      // if ((np == 3) || (np == 5))
      //   break;
    }
  }
  
  if (me == 0 && !fail)
    std::cout << "PASS" << std::endl;
  
  return 0;
}

