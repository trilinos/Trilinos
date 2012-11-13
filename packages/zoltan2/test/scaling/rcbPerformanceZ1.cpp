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

/*! \file rcbPerformanceZ1.cpp
    \brief rcbPerformance with Zoltan1

    Geometry is a uniform mesh.
    \todo  get the imbalance when done
*/

#include <Zoltan2_TestHelpers.hpp>

#ifdef HAVE_ZOLTAN2_ZOLTAN
#include <zoltan.h>
#include <Teuchos_CommandLineProcessor.hpp>

#include <vector>
#include <ostream>
#include <sstream>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bad_alloc;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::CommandLineProcessor;

typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Tpetra::Map<lno_t, gno_t, node_t> tMap_t;

static ArrayRCP<ArrayRCP<scalar_t> > weights;
static RCP<tMVector_t> coordinates;

// Zoltan1 query functions

int getNumObj(void *data, int *ierr)
{
  *ierr = 0;
  return coordinates->getLocalLength();
}

int getDim(void *data, int *ierr)
{
  *ierr = 0;
  return 3;
}

void getObjList(void *data, int numGid, int numLid,
  gid_t * gids, gid_t * lids, 
  int wgt_dim, float *obj_wgts, int *ierr)
{
  *ierr = 0;
  size_t localLen = coordinates->getLocalLength();
  const gno_t *ids = coordinates->getMap()->getNodeElementList().getRawPtr();
  gno_t *idsNonConst = const_cast<gno_t *>(ids);

  if (sizeof(gid_t) == sizeof(gno_t)){
    memcpy(gids, idsNonConst, sizeof(gid_t) * localLen);
  }
  else{
    for (size_t i=0; i < localLen; i++)
      gids[i] = static_cast<gid_t>(idsNonConst[i]);
  }

  if (wgt_dim > 0){
    float *wgts = obj_wgts;
    for (size_t i=0; i < localLen; i++)
      for (int w=0; w < wgt_dim; w++)
        *wgts++ = static_cast<float>(weights[w][i]);
  }
}

void getCoordinates(void *data, int numGid, int numLid,
  int numObj, gid_t * gids, gid_t * lids,
  int dim, double *coords, int *ierr)
{
  // I know that Zoltan asks for coordinates in gid order.
  *ierr = 0;
  double *val = coords;
  const scalar_t *x = coordinates->getData(0).getRawPtr();
  const scalar_t *y = coordinates->getData(1).getRawPtr();
  const scalar_t *z = coordinates->getData(2).getRawPtr();
  for (lno_t i=0; i < numObj; i++){
    *val++ = static_cast<double>(x[i]);
    *val++ = static_cast<double>(y[i]);
    *val++ = static_cast<double>(z[i]);
  }
}



enum weightTypes{
  upDown,
  roundRobin,
  increasing,
  numWeightTypes
};

ArrayRCP<scalar_t> makeWeights(
  const RCP<const Teuchos::Comm<int> > & comm,
  lno_t len, weightTypes how, scalar_t scale, int rank)
{
  scalar_t *wgts = new scalar_t [len];
  if (!wgts)
    throw bad_alloc();

  ArrayRCP<scalar_t> wgtArray(wgts, 0, len, true);

  if (how == upDown){
    scalar_t val = scale + rank%2;
    for (lno_t i=0; i < len; i++)
      wgts[i] = val;
  }
  else if (how == roundRobin){
    for (int i=0; i < 10; i++){
      scalar_t val = (i + 10)*scale;
      for (int j=i; j < len; j += 10)
         wgts[j] = val;
    }
  }
  else if (how == increasing){
    scalar_t val = scale + rank;
    for (lno_t i=0; i < len; i++)
      wgts[i] = val;
  }

  return wgtArray;
}

/*! \brief Create a mesh of approximately the desired size.
 *
 *  We want 3 dimensions close to equal in length.
 */
const RCP<tMVector_t> getMeshCoordinates(
    const RCP<const Teuchos::Comm<int> > & comm,
    gno_t numGlobalCoords)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  gno_t xdim = static_cast<int>(floor(xdimf));
  gno_t ydim = xdim;
  gno_t zdim = numGlobalCoords / (xdim*ydim);
  gno_t num=xdim*ydim*zdim;
  gno_t diff = numGlobalCoords - num;
  gno_t newdiff = 0;

  while (diff > 0){
    if (zdim > xdim && zdim > ydim){
      zdim++;
      newdiff = diff - (xdim*ydim);
      if (newdiff < 0)
        if (diff < -newdiff)
          zdim--;
    }
    else if (ydim > xdim && ydim > zdim){
      ydim++;
      newdiff = diff - (xdim*zdim);
      if (newdiff < 0)
        if (diff < -newdiff)
          ydim--;
    }
    else{
      xdim++;
      newdiff = diff - (ydim*zdim);
      if (newdiff < 0)
        if (diff < -newdiff)
          xdim--;
    }

    diff = newdiff;
  }

  num=xdim*ydim*zdim;
  diff = numGlobalCoords - num;
  if (diff < 0)
    diff /= -numGlobalCoords;
  else
    diff /= numGlobalCoords;

  if (rank == 0){
    if (diff > .01)
      cout << "Warning: Difference " << diff*100 << " percent" << endl;
    cout << "Mesh size: " << xdim << "x" << ydim << "x" <<
      zdim << ", " << num << " vertices." << endl;
  }

  // Divide coordinates.

  gno_t numLocalCoords = num / nprocs;
  gno_t leftOver = num % nprocs;
  gno_t gid0 = 0;

  if (rank <= leftOver)
    gid0 = gno_t(rank) * (numLocalCoords+1);
  else
    gid0 = (leftOver * (numLocalCoords+1)) + 
           ((gno_t(rank) - leftOver) * numLocalCoords);

  if (rank < leftOver)
    numLocalCoords++;

  gno_t gid1 = gid0 + numLocalCoords;

  gno_t *ids = new gno_t [numLocalCoords];
  if (!ids)
    throw bad_alloc();
  ArrayRCP<gno_t> idArray(ids, 0, numLocalCoords, true);

  for (gno_t i=gid0; i < gid1; i++)
    *ids++ = i;   

  RCP<const tMap_t> idMap = rcp(
    new tMap_t(num, idArray.view(0, numLocalCoords), 0, comm));

  // Create a Tpetra::MultiVector of coordinates.

  scalar_t *x = new scalar_t [numLocalCoords*3]; 
  if (!x)
    throw bad_alloc();
  ArrayRCP<scalar_t> coordArray(x, 0, numLocalCoords*3, true);

  scalar_t *y = x + numLocalCoords;
  scalar_t *z = y + numLocalCoords;

  gno_t xStart = 0;
  gno_t yStart = 0;
  gno_t xyPlane = xdim*ydim;
  gno_t zStart = gid0 / xyPlane;
  gno_t rem = gid0 % xyPlane;
  if (rem > 0){
    yStart = rem / xdim;
    xStart = rem % xdim;
  }

  lno_t next = 0;
  for (scalar_t zval=zStart; next < numLocalCoords && zval < zdim; zval++){
    for (scalar_t yval=yStart; next < numLocalCoords && yval < ydim; yval++){
      for (scalar_t xval=xStart; next < numLocalCoords && xval < xdim; xval++){
        x[next] = xval;
        y[next] = yval;
        z[next] = zval;
        next++;
      }
      xStart = 0;
    }
    yStart = 0;
  }

  // TODO:  KDD ArrayRCP<ArrayView<>> has to go!

  ArrayView<const scalar_t> xArray(x, numLocalCoords);
  ArrayView<const scalar_t> yArray(y, numLocalCoords);
  ArrayView<const scalar_t> zArray(z, numLocalCoords);
  ArrayRCP<ArrayView<const scalar_t> > coordinateArrays =
    arcp(new ArrayView<const scalar_t> [3], 0, 3);
  coordinateArrays[0] = xArray;
  coordinateArrays[1] = yArray;
  coordinateArrays[2] = zArray;

  ArrayRCP<const ArrayView<const scalar_t> > constCoords =
   coordinateArrays.getConst();

  RCP<tMVector_t> meshCoords = rcp(new tMVector_t(
    idMap, constCoords.view(0,3), 3));

  return meshCoords;
}


int main(int argc, char *argv[])
{
  // MEMORY_CHECK(true, "Before initializing MPI");

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  MEMORY_CHECK(rank==0 || rank==nprocs-1, "After initializing MPI");

  if (rank==0)
    cout << "Number of processes: " << nprocs << endl;

  // Default values
  double numGlobalCoords = 1000;
  int weightDim = 0;
  int debugLevel=2;          // for timing
  string memoryOn("memoryOn");
  string memoryOff("memoryOff");
  bool doMemory=false;
  int numGlobalParts = nprocs;
  int dummyTimer=0;
  bool remap=0;

  CommandLineProcessor commandLine(false, true);
  commandLine.setOption("size", &numGlobalCoords, 
    "Approximate number of global coordinates.");
  commandLine.setOption("numParts", &numGlobalParts, 
    "Number of parts (default is one per proc).");
  commandLine.setOption("weightDim", &weightDim, 
    "Number of weights per coordinate, zero implies uniform weights.");
  commandLine.setOption("debug", &debugLevel, "Zoltan1 debug level");
  commandLine.setOption("remap", "no-remap", &remap,
    "Zoltan1 REMAP parameter; disabled by default for scalability testing");
  commandLine.setOption("timers", &dummyTimer, "ignored");
  commandLine.setOption(memoryOn.c_str(), memoryOff.c_str(), &doMemory,
    "do memory profiling");

  string balanceCount("balance_object_count");
  string balanceWeight("balance_object_weight");
  string mcnorm1("multicriteria_minimize_total_weight");
  string mcnorm2("multicriteria_balance_total_maximum");
  string mcnorm3("multicriteria_minimize_maximum_weight");

  string objective(balanceWeight);   // default

  string doc(balanceCount);
  doc.append(": ignore weights\n");

  doc.append(balanceWeight);
  doc.append(": balance on first weight\n");

  doc.append(mcnorm1);
  doc.append(": given multiple weights, balance their total.\n");

  doc.append(mcnorm3);
  doc.append(": given multiple weights, balance the maximum for each coordinate.\n");

  doc.append(mcnorm2);
  doc.append(": given multiple weights, balance the L2 norm of the weights.\n");

  commandLine.setOption("objective", &objective,  doc.c_str());

  CommandLineProcessor::EParseCommandLineReturn rc = 
    commandLine.parse(argc, argv);

  if (rc != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL){
    if (rc == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED){
      if (rank==0)
        cout << "PASS" << endl;
      return 1;
    }
    else{
      if (rank==0)
        cout << "FAIL" << endl;
      return 0;
    }
  }

  //MEMORY_CHECK(doMemory && rank==0, "After processing parameters");

  gno_t globalSize = static_cast<gno_t>(numGlobalCoords);

  coordinates = getMeshCoordinates(comm, globalSize);
  size_t numLocalCoords = coordinates->getLocalLength();

#if 0
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p==rank){
      cout << "Rank " << rank << ", " << numLocalCoords << "coords" << endl;
      const scalar_t *x = coordinates->getData(0).getRawPtr();
      const scalar_t *y = coordinates->getData(1).getRawPtr();
      const scalar_t *z = coordinates->getData(2).getRawPtr();
      for (lno_t i=0; i < numLocalCoords; i++)
        cout << " " << x[i] << " " << y[i] << " " << z[i] << endl;
    }
    cout.flush();
    comm->barrier();
  }
#endif

  if (weightDim > 0){

    weights = arcp(new ArrayRCP<scalar_t> [weightDim],
      0, weightDim, true);

    int wt = 0;
    scalar_t scale = 1.0;
    for (int i=0; i < weightDim; i++){
      weights[i] = 
        makeWeights(comm, numLocalCoords, weightTypes(wt++), scale, rank);

      if (wt == numWeightTypes){
        wt = 0;
        scale++;
      }
    }
  }

  MEMORY_CHECK(doMemory && rank==0, "After creating input");

  // Now call Zoltan to partition the problem.

  float ver;
  int aok = Zoltan_Initialize(argc, argv, &ver);

  if (aok != 0){
    printf("sorry...\n");
    exit(0);
  }

  struct Zoltan_Struct *zz;
  zz = Zoltan_Create(MPI_COMM_WORLD);
  
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "CHECK_GEOM", "0");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); // compiled with ULONG option
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "PART");
  Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.1");
  std::ostringstream oss;
  oss << numGlobalParts;
  Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", oss.str().c_str());
  oss.str("");
  oss << debugLevel;
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", oss.str().c_str());

  if (remap)
    Zoltan_Set_Param(zz, "REMAP", "1");
  else
    Zoltan_Set_Param(zz, "REMAP", "0");

  if (objective != balanceCount){
    oss.str("");
    oss << weightDim;
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", oss.str().c_str());

    if (objective == mcnorm1)
      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "1");
    else if (objective == mcnorm2)
      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "2");
    else if (objective == mcnorm3)
      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "3");
  }
  else{
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
  }

  Zoltan_Set_Num_Obj_Fn(zz, getNumObj, NULL);
  Zoltan_Set_Obj_List_Fn(zz, getObjList,NULL);
  Zoltan_Set_Num_Geom_Fn(zz, getDim, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, getCoordinates,NULL);

  int changes, numGidEntries, numLidEntries, numImport, numExport;
  gid_t * importGlobalGids, * importLocalGids;
  gid_t * exportGlobalGids, * exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  MEMORY_CHECK(doMemory && rank==0, "Before Zoltan_LB_Partition");

  if (rank == 0) std::cout << "Calling Zoltan_LB_Partition" << std::endl;
  aok = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the vertices I must send */
        &exportLocalGids,   /* Local IDs of the vertices I must send */
        &exportProcs,    /* Process to which I send each of the vertices */
        &exportToPart);  /* Partition to which each vertex will belong */

  if (rank == 0) std::cout << "Returned from Zoltan_LB_Partition" << std::endl;
  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_LB_Partition");

  /* Print the load-balance stats here */

  scalar_t *sumWgtPerPart = new scalar_t[numGlobalParts];
  scalar_t *gsumWgtPerPart = new scalar_t[numGlobalParts];
  for (int i = 0; i < numGlobalParts; i++) sumWgtPerPart[i] = 0.;

  for (size_t i = 0; i < numLocalCoords; i++)
    sumWgtPerPart[exportToPart[i]] += (weightDim ? weights[0][i]: 1.);

  reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, numGlobalParts,
                           sumWgtPerPart, gsumWgtPerPart);

  scalar_t maxSumWgtPerPart = 0.;
  scalar_t minSumWgtPerPart = std::numeric_limits<scalar_t>::max();
  scalar_t totWgt = 0.;
  int maxSumWgtPart=0, minSumWgtPart=0;
  for (int i = 0; i < numGlobalParts; i++) {
    if (gsumWgtPerPart[i] > maxSumWgtPerPart) {
      maxSumWgtPerPart = gsumWgtPerPart[i];
      maxSumWgtPart = i;
    }
    if (gsumWgtPerPart[i] < minSumWgtPerPart) {
      minSumWgtPerPart = gsumWgtPerPart[i];
      minSumWgtPart = i;
    }
    totWgt += gsumWgtPerPart[i];
  }

  if (rank == 0)
    std::cout << std::endl << std::endl
              << "Part loads (per part for " << numGlobalParts << " parts):"
              << std::endl
              << "   min = " << minSumWgtPerPart
                             << " in part " << minSumWgtPart << std::endl
              << "   max = " << maxSumWgtPerPart
                             << " in part " << maxSumWgtPart << std::endl
              << "   tot = " << totWgt << std::endl
              << "   avg = " << totWgt / numGlobalParts 
              << std::endl << std::endl << std::endl;

  delete [] sumWgtPerPart;
  delete [] gsumWgtPerPart;

  Zoltan_Destroy(&zz);
  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_Destroy");

  if (rank==0){
    if (aok != 0)
      std::cout << "FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return 0;
}

#else
#include <iostream>
int main(int argc, char *argv[])
{
  std::cout << "Test did not run due to faulty configuration." << std::endl;
  std::cout << "PASS" << std::endl;
}
#endif
