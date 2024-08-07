// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file zoltanCompare.cpp
 *  Compares zoltan execution through Zoltan2 with direct zoltan execution
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

#include <Tpetra_MultiVector.hpp>
#include <zoltan.h>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;

//
// A few of the tests done by Zoltan in nightly testing.
//

enum testFields {
  TESTNAMEOFFSET = 0,
  TESTMETHODOFFSET,
  TESTOBJWGTOFFSET,
  TESTNUMPROCS,
  TESTNUMARGS
};


#define NUMTESTS 30
static string testArgs[] = {
// Filename  LB_Method   ObjWeightDim   NumProcs
"simple",       "phg",          "0",      "2",
"simple",       "rcb",          "0",      "2",
"vwgt2",        "rcb",          "2",      "2",

"bug",          "rcb",          "1",      "3",
"drake",        "rcb",          "0",      "3",
"onedbug",      "rcb",          "0",      "3",
"simple",       "rcb",          "0",      "3",
"vwgt",         "rcb",          "1",      "3",
"vwgt",         "phg",          "1",      "3",
"vwgt2",        "rcb",          "2",      "3",

"simple",       "default",      "0",      "4",
"ewgt",         "hsfc",         "0",      "4",
"grid20x19",    "hsfc",         "0",      "4",
"grid20x19",    "hsfc",         "0",      "4",
"grid20x19",    "hsfc",         "0",      "4",
"nograph",      "rib",          "0",      "4",
"nograph",      "phg",          "0",      "4",
"simple",       "rib",          "0",      "4",
"simple",       "rib",          "0",      "4",
"vwgt2",        "rib",          "2",      "4",
"simple",       "rib",          "0",      "4",
"simple",       "phg",          "0",      "4",

"brack2_3",     "rcb",          "3",      "5",

"hammond2",     "rcb",          "2",      "6",
"degenerateAA", "rcb",          "0",      "6",
"degenerate",   "rcb",          "0",      "6",
"degenerate",   "rcb",          "0",      "6",

"hammond",      "rcb",          "0",      "8",
"hammond",      "phg",          "0",      "8",
"vwgt2",        "rcb",          "2",      "8"
};

typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tGraph_t;
typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> tMatrix_t;
typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> vectorAdapter_t;
typedef Zoltan2::XpetraCrsMatrixAdapter<tMatrix_t,tMVector_t> matrixAdapter_t;
typedef Zoltan2::EvaluatePartition<matrixAdapter_t> quality_t;

// Number of ZOLTAN_ID in a zgno_t (for NUM_GID_ENTRIES)
static constexpr int znGidEnt = sizeof(zgno_t) / sizeof(ZOLTAN_ID_TYPE);
#define SET_ZID(n,a,b)                                            \
   {int ZOLTAN_ID_LOOP;                                                 \
    for (ZOLTAN_ID_LOOP = 0; ZOLTAN_ID_LOOP < (n); ZOLTAN_ID_LOOP++)    \
      (a)[ZOLTAN_ID_LOOP] = (b)[ZOLTAN_ID_LOOP];                        \
   }


////////////////////////////////////////////////////////////////////////////////
// Zoltan callbacks

static int znumobj(void *data, int *ierr) 
{
  *ierr = ZOLTAN_OK;
  tMVector_t *vec = (tMVector_t *) data;
  return vec->getLocalLength();
}

static void zobjlist(void *data, int ngid, int nlid, 
                     ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                     int nwgts, float *wgts, int *ierr)
{
  *ierr = ZOLTAN_OK;
  tMVector_t *vec = (tMVector_t *) data;
  size_t n = vec->getLocalLength();

  for (size_t i = 0; i < n; i++) {
    zgno_t vgid = vec->getMap()->getGlobalElement(i);
    ZOLTAN_ID_PTR vgidptr = (ZOLTAN_ID_PTR) &vgid;
    SET_ZID(znGidEnt, &(gids[i*znGidEnt]), vgidptr);
    lids[i] = i;
  }

  for (int w = 0; w < nwgts; w++) {
    ArrayRCP<const zscalar_t> wvec = vec->getData(w);
    for (size_t i = 0; i < n; i++)
      wgts[i*nwgts+w] = wvec[i];
  }
}

static int znumgeom(void *data, int *ierr) 
{
  *ierr = ZOLTAN_OK;
  tMVector_t *cvec = (tMVector_t *) data;
  return cvec->getNumVectors();
}

static void zgeom(void *data, int ngid, int nlid, int nobj, 
                  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                  int ndim, double *coords, int *ierr)
{
  *ierr = ZOLTAN_OK;
  tMVector_t *vec = (tMVector_t *) data;
  for (int d = 0; d < ndim; d++) {
    ArrayRCP<const zscalar_t> cvec = vec->getData(d);
    for (int i = 0; i < nobj; i++) {
      coords[lids[i]*ndim+d] = cvec[lids[i]];
    }
  }
}

static void zhgsize(void *data, int *nLists, int *nPins, int *format, int *ierr)
{
  tMatrix_t *matrix = (tMatrix_t *) data;
  *nLists = matrix->getLocalNumRows();
  *nPins = matrix->getLocalNumEntries();
  *format = ZOLTAN_COMPRESSED_VERTEX;
  *ierr = ZOLTAN_OK;
}

static void zhg(void *data, int ngid, int nLists, int nPins, int format,
                ZOLTAN_ID_PTR listGids, int *offsets, ZOLTAN_ID_PTR pinGids, 
                int *ierr)
{
  tMatrix_t *matrix = (tMatrix_t *) data;
  RCP<const tGraph_t> graph = matrix->getCrsGraph();
  zlno_t nrows = graph->getLocalNumRows();
  
  offsets[0] = 0;
  for (zlno_t i = 0; i < nrows; i++) {

    zgno_t tmp = graph->getRowMap()->getGlobalElement(i);
    ZOLTAN_ID_PTR ztmp = (ZOLTAN_ID_PTR) &tmp;
    SET_ZID(znGidEnt, &(listGids[i*znGidEnt]), ztmp);

    size_t nEntries = graph->getNumEntriesInLocalRow(i);
    offsets[i+1] = offsets[i] + nEntries;

    
    typename tMatrix_t::local_inds_host_view_type colind;
    graph->getLocalRowView(i, colind);

    for (size_t j = 0; j < nEntries; j++) {
      tmp = graph->getColMap()->getGlobalElement(colind[j]);
      ztmp = (ZOLTAN_ID_PTR) &tmp;
      SET_ZID(znGidEnt, &(pinGids[(offsets[i]+j)*znGidEnt]), ztmp);
    }
  }

  *ierr = ZOLTAN_OK;
}

////////////////////////////////////////////////////////////////////////////////
// Function to compute both Zoltan2 and Zoltan partitions and print metrics

int run(
  const RCP<const Comm<int> > &comm,
  int numGlobalParts,
  int testCnt,
  std::string *thisTest
)
{
#ifdef HAVE_ZOLTAN2_MPI
  // Zoltan needs an MPI comm
  const Teuchos::MpiComm<int> *tmpicomm =
               dynamic_cast<const Teuchos::MpiComm<int> *>(comm.getRawPtr());
  MPI_Comm mpiComm = *(tmpicomm->getRawMpiComm());
#endif

  int me = comm->getRank();
  int np = comm->getSize();
  double tolerance = 1.05;


  //////////////////////////////////////////////
  // Read test data from Zoltan's test directory
  //////////////////////////////////////////////

  UserInputForTests *uinput;
  try{
    uinput = new UserInputForTests(zoltanTestDirectory,
                                   thisTest[TESTNAMEOFFSET],
                                   comm, true);
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: UserInputForTests "
           << e.what() << std::endl;
    return 1;
  }

  RCP<tMatrix_t> matrix;
  try{
    matrix = uinput->getUITpetraCrsMatrix();
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: get matrix "
           << e.what() << std::endl;
    return 1;
  }

  RCP<const tMatrix_t> matrixConst = rcp_const_cast<const tMatrix_t>(matrix);

  RCP<tMVector_t> coords;
  try{
   coords = uinput->getUICoordinates();
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: get coordinates "
           << e.what() << std::endl;
    return 1;
  }

  RCP<tMVector_t> weights;
  try{
   weights = uinput->getUIWeights();
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: get weights "
           << e.what() << std::endl;
    return 1;
  }
  int nWeights = atoi(thisTest[TESTOBJWGTOFFSET].c_str());

  if (me == 0) {
    std::cout << "Test " << testCnt << " filename            = "
         << thisTest[TESTNAMEOFFSET] << std::endl;
    std::cout << "Test " << testCnt << " num processors      = "
         << np << std::endl;
    std::cout << "Test " << testCnt << " zoltan method       = "
         << thisTest[TESTMETHODOFFSET] << std::endl;
    std::cout << "Test " << testCnt << " num_global_parts    = "
         << numGlobalParts << std::endl;
    std::cout << "Test " << testCnt << " imbalance_tolerance = "
         << tolerance << std::endl;
    std::cout << "Test " << testCnt << " num weights per ID  = "
         << nWeights << std::endl;
  }

  /////////////////////////////////////////
  // PARTITION USING ZOLTAN DIRECTLY
  /////////////////////////////////////////

  if (me == 0) std::cout << "Calling Zoltan directly" << std::endl;

# ifdef HAVE_ZOLTAN2_MPI
    Zoltan zz(mpiComm);
# else
    Zoltan zz;
# endif

  char tmp[56];
  zz.Set_Param("LB_METHOD", thisTest[TESTMETHODOFFSET]);
  
  sprintf(tmp, "%d", znGidEnt);
  zz.Set_Param("NUM_GID_ENTRIES", tmp);
  sprintf(tmp, "%d", numGlobalParts);
  zz.Set_Param("NUM_GLOBAL_PARTS", tmp);
  sprintf(tmp, "%d", nWeights);
  zz.Set_Param("OBJ_WEIGHT_DIM", tmp);
  sprintf(tmp, "%f", tolerance);
  zz.Set_Param("IMBALANCE_TOL", tmp);
  zz.Set_Param("RETURN_LISTS", "PART");
  zz.Set_Param("FINAL_OUTPUT", "1");
  zz.Set_Param("SEED", "1111");
  zz.Set_Param("LB_APPROACH", "PARTITION");

  zz.Set_Num_Obj_Fn(znumobj, (void *) coords.getRawPtr());
  if (nWeights)
    zz.Set_Obj_List_Fn(zobjlist, (void *) weights.getRawPtr());
  else
    zz.Set_Obj_List_Fn(zobjlist, (void *) coords.getRawPtr());
  zz.Set_Num_Geom_Fn(znumgeom, (void *) coords.getRawPtr());
  zz.Set_Geom_Multi_Fn(zgeom, (void *) coords.getRawPtr());
  zz.Set_HG_Size_CS_Fn(zhgsize, (void *) &(*matrix));
  zz.Set_HG_CS_Fn(zhg, (void *) &(*matrix));

  int changes, ngid, nlid;
  int numd, nump;
  ZOLTAN_ID_PTR dgid = NULL, dlid = NULL, pgid = NULL, plid = NULL;
  int *dproc = NULL, *dpart = NULL, *pproc = NULL, *ppart = NULL;

  int ierr = zz.LB_Partition(changes, ngid, nlid,
                             numd, dgid, dlid, dproc, dpart,
                             nump, pgid, plid, pproc, ppart);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: direct Zoltan call" << std::endl;
    zz.LB_Free_Part(&pgid, &plid, &pproc, &ppart);
    return 1;
  }

for(int i = 0; i < nump; i++) {
  std::cout << me << " KDD Z1 " << pgid[i] << " " << plid[i] << " " << ppart[i] << std::endl;
}

  /////////////////////////////////////////
  // PARTITION USING ZOLTAN THROUGH ZOLTAN2
  /////////////////////////////////////////

  if (me == 0) std::cout << "Calling Zoltan through Zoltan2" << std::endl;

  matrixAdapter_t *ia;
  try{
    ia = new matrixAdapter_t(matrixConst, nWeights);
  }
  catch(std::exception &e){
    if (me == 0)
      std::cout << "Test " << testCnt << ":  FAIL: matrix adapter "
           << e.what() << std::endl;
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
      std::cout << "Test " << testCnt << ":  FAIL: vector adapter "
           << e.what() << std::endl;
    return 1;
  }
  ia->setCoordinateInput(ca);
  
  Teuchos::ParameterList params;
  params.set("timer_output_stream" , "std::cout");
  // params.set("debug_level" , "verbose_detailed_status");

  params.set("algorithm", "zoltan");
  params.set("imbalance_tolerance", tolerance );
  params.set("num_global_parts", numGlobalParts);

  if (thisTest[TESTMETHODOFFSET] != "default") {
    // "default" tests case of no Zoltan parameter sublist
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD",thisTest[TESTMETHODOFFSET]);
    zparams.set("LB_APPROACH", "PARTITION");
    zparams.set("SEED", "1111");
  }

  Zoltan2::PartitioningProblem<matrixAdapter_t> *problem;
# ifdef HAVE_ZOLTAN2_MPI
    try{
      problem = new Zoltan2::PartitioningProblem<matrixAdapter_t>(ia, &params,
                                                                  mpiComm);
    }
# else
    try{
      problem = new Zoltan2::PartitioningProblem<matrixAdapter_t>(ia, &params);
    }
# endif
  catch(std::exception &e){
    std::cout << "Test " << testCnt << " FAIL: problem " << e.what() << std::endl;
    return 1;
  }

  try {
    problem->solve();
  }
  catch(std::exception &e){
    std::cout << "Test " << testCnt << " FAIL: solve " << e.what() << std::endl;
    return 1;
  }

for(int i = 0; i < nump; i++) {
  std::cout << me << " KDD Z2 " << coords->getMap()->getGlobalElement(i)  << " " << i << " " << problem->getSolution().getPartListView()[i] << std::endl;
}

  // create metric object

 RCP<quality_t>metricObject = rcp(new quality_t(ia,&params,problem->getComm(),
						&problem->getSolution()));
  if (me == 0){
    metricObject->printMetrics(std::cout);
  }
  problem->printTimers();

  /////////////////////////////////////////
  // COMPARE RESULTS
  /////////////////////////////////////////
  size_t nObj = coords->getLocalLength();
  const int *z2parts = problem->getSolution().getPartListView();
  int diffcnt = 0, gdiffcnt = 0;
  for (size_t i = 0; i < nObj; i++) {
    if (z2parts[plid[i]] != ppart[i]) {
      diffcnt++;
      std::cout << me << " DIFF for " << i << " (" 
           << coords->getMap()->getGlobalElement(i) << "):  "
           << "Z2 = " << z2parts[i] << "; Z1 = " << ppart[plid[i]] << std::endl;
    }
  }

  /////////////////////////////////////////
  // CLEAN UP
  /////////////////////////////////////////
  zz.LB_Free_Part(&pgid, &plid, &pproc, &ppart);
  delete ia;
  delete ca;
  delete problem;
  delete uinput;

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &diffcnt, &gdiffcnt);
  if (gdiffcnt > 0) {
    if (me == 0) 
      std::cout << "Test " << testCnt << " "
           << thisTest[TESTNAMEOFFSET] << " "
           << thisTest[TESTMETHODOFFSET] << " "
           << thisTest[TESTOBJWGTOFFSET] << " "
           << " FAIL: comparison " << std::endl;
    return 1;
  }

  return 0;
}
  
////////////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int me = comm->getRank();
  int np = comm->getSize();

  int fail=0;

  Array<int> ranks(np);
  for (int i = 0; i < np; i++) ranks[i] = i;

  for (int i=0; i < NUMTESTS; i++) {
    std::string *thisTest = &(testArgs[i*TESTNUMARGS]);
    int nTestProcs = atoi(thisTest[TESTNUMPROCS].c_str());
    if (nTestProcs > np) {
      if (me == 0) {
        std::cout << "Skipping test " << i << " on "
             << thisTest[TESTNAMEOFFSET]
             << "; required number of procs " << nTestProcs 
             << " is greater than available procs " << np << std::endl;
      }
      continue;
    }

    // Make a communicator of appropriate size for the test
    RCP<const Comm<int> > testcomm;
    if (nTestProcs == np)
      testcomm = comm;
    else
      testcomm = comm->createSubcommunicator(ranks.view(0,nTestProcs));

    // Run the test if in the communicator
    if (me < nTestProcs) {
      fail += run(testcomm, nTestProcs, i, thisTest);
    }
  }
  
  if (me == 0 && !fail)
    std::cout << "PASS" << std::endl;
  
  return 0;
}

