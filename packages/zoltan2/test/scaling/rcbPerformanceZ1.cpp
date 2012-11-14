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
*/

#include "Zoltan2_config.h"

#ifdef HAVE_ZOLTAN2_ZOLTAN
#include <zoltan.h>

#include <Zoltan2_Util.hpp>
#define MEMORY_CHECK(iPrint, msg) \
  if (iPrint){ \
    long kb = Zoltan2::getProcessKilobytes(); \
    std::cout.width(10); \
    std::cout.fill('*'); \
    std::cout << kb << " KB, " << msg << std::endl; \
    std::cout.width(0); \
    std::cout.fill(' '); \
  }

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Kokkos_DefaultNode.hpp>

#include <vector>
#include <string>
#include <ostream>
#include <sstream>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bad_alloc;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
using Teuchos::CommandLineProcessor;

typedef int lno_t;
typedef long gno_t;
typedef double scalar_t;
typedef Kokkos::DefaultNode::DefaultNodeType node_t;

typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Tpetra::Map<lno_t, gno_t, node_t> tMap_t;

//////////////////////////////////////////////////////////////////////////////
// Data structure for data
typedef struct dots {
  vector<vector<float> > weights;
  tMVector_t *coordinates;
} DOTS;

//////////////////////////////////////////////////////////////////////////////
// Zoltan1 query functions

int getNumObj(void *data, int *ierr)
{
  *ierr = 0;
  DOTS *dots = (DOTS *) data;
  return dots->coordinates->getLocalLength();
}

//////////////////////////
int getDim(void *data, int *ierr)
{
  *ierr = 0;
  DOTS *dots = (DOTS *) data;
  return dots->coordinates->getNumVectors();
}

//////////////////////////
void getObjList(void *data, int numGid, int numLid,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
  int wgt_dim, float *obj_wgts, int *ierr)
{
  *ierr = 0;
  DOTS *dots = (DOTS *) data;

  size_t localLen = dots->coordinates->getLocalLength();
  const gno_t *ids = 
               dots->coordinates->getMap()->getNodeElementList().getRawPtr();

  if (sizeof(ZOLTAN_ID_TYPE) == sizeof(gno_t)) 
    memcpy(gids, ids, sizeof(ZOLTAN_ID_TYPE) * localLen);
  else 
    for (size_t i=0; i < localLen; i++)
      gids[i] = static_cast<ZOLTAN_ID_TYPE>(ids[i]);

  if (wgt_dim > 0){
    float *wgts = obj_wgts;
    for (size_t i=0; i < localLen; i++)
      for (int w=0; w < wgt_dim; w++)
        *wgts++ = dots->weights[w][i];
  }
}

//////////////////////////
void getCoordinates(void *data, int numGid, int numLid,
  int numObj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
  int dim, double *coords, int *ierr)
{
  // I know that Zoltan asks for coordinates in gid order.
  *ierr = 0;
  DOTS *dots = (DOTS *) data;
  double *val = coords;
  const scalar_t *x = dots->coordinates->getData(0).getRawPtr();
  const scalar_t *y = dots->coordinates->getData(1).getRawPtr();
  const scalar_t *z = dots->coordinates->getData(2).getRawPtr();
  for (int i=0; i < numObj; i++){
    *val++ = static_cast<double>(x[i]);
    *val++ = static_cast<double>(y[i]);
    *val++ = static_cast<double>(z[i]);
  }
}


//////////////////////////////////////////////////////////////////////////////

enum weightTypes{
  upDown,
  roundRobin,
  increasing,
  numWeightTypes
};

void makeWeights(
  const RCP<const Teuchos::Comm<int> > & comm,
  vector<float> &wgts, weightTypes how, float scale, int rank)
{
  lno_t len = wgts.size();
  if (how == upDown){
    float val = scale + rank%2;
    for (lno_t i=0; i < len; i++)
      wgts[i] = val;
  }
  else if (how == roundRobin){
    for (int i=0; i < 10; i++){
      float val = (i + 10)*scale;
      for (int j=i; j < len; j += 10)
         wgts[j] = val;
    }
  }
  else if (how == increasing){
    float val = scale + rank;
    for (lno_t i=0; i < len; i++)
      wgts[i] = val;
  }
}

//////////////////////////////////////////////////////////////////////////////
/* Create a mesh of approximately the desired size.
 *
 *  We want 3 dimensions close to equal in length.
 */
tMVector_t* makeMeshCoordinates(
    const RCP<const Teuchos::Comm<int> > & comm,
    gno_t numGlobalCoords)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  gno_t xdim = static_cast<gno_t>(floor(xdimf));
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
  ArrayView<gno_t> idArray(ids, numLocalCoords);

  for (gno_t i=gid0; i < gid1; i++)
    *ids++ = i;   

  RCP<const tMap_t> idMap = rcp(new tMap_t(num, idArray, 0, comm));

  delete [] ids;

  // Create a Tpetra::MultiVector of coordinates.

  scalar_t *x = new scalar_t [numLocalCoords*3]; 
  if (!x) throw bad_alloc();

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
  for (scalar_t zval=zStart; next < numLocalCoords && zval < zdim; zval+=1.){
    for (scalar_t yval=yStart; next < numLocalCoords && yval < ydim; yval+=1.){
      for (scalar_t xval=xStart; next < numLocalCoords && xval < xdim;xval+=1.){
        x[next] = xval;
        y[next] = yval;
        z[next] = zval;
        next++;
      }
      xStart = 0;
    }
    yStart = 0;
  }

  ArrayView<const scalar_t> xArray(x, numLocalCoords*3);
  tMVector_t *dots = new tMVector_t(idMap, xArray, numLocalCoords, 3);

  return dots;
}


int main(int argc, char *argv[])
{
  // MEMORY_CHECK(true, "Before initializing MPI");

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  DOTS dots;

  MEMORY_CHECK(rank==0 || rank==nprocs-1, "After initializing MPI");

  if (rank==0)
    cout << "Number of processes: " << nprocs << endl;

  // Default values
  gno_t numGlobalCoords = 1000;
  int weightDim = 0;
  int debugLevel=0;
  string memoryOn("memoryOn");
  string memoryOff("memoryOff");
  bool doMemory=false;
  int numGlobalParts = nprocs;
  int dummyTimer=0;
  bool remap=0;

  string balanceCount("balance_object_count");
  string balanceWeight("balance_object_weight");
  string mcnorm1("multicriteria_minimize_total_weight");
  string mcnorm2("multicriteria_balance_total_maximum");
  string mcnorm3("multicriteria_minimize_maximum_weight");
  string objective(balanceWeight);   // default

  // Process command line input
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

  string doc(balanceCount);
  doc.append(": ignore weights\n");
  doc.append(balanceWeight);
  doc.append(": balance on first weight\n");
  doc.append(mcnorm1);
  doc.append(": given multiple weights, balance their total.\n");
  doc.append(mcnorm3);
  doc.append(": given multiple weights, "
             "balance the maximum for each coordinate.\n");
  doc.append(mcnorm2);
  doc.append(": given multiple weights, balance the L2 norm of the weights.\n");
  commandLine.setOption("objective", &objective,  doc.c_str());

  CommandLineProcessor::EParseCommandLineReturn rc = 
    commandLine.parse(argc, argv);
  if (rc != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    if (rc == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
      if (rank==0) cout << "PASS" << endl;
      return 1;
    }
    else {
      if (rank==0) cout << "FAIL" << endl;
      return 0;
    }
  }

  //MEMORY_CHECK(doMemory && rank==0, "After processing parameters");

  // Create the data structure
  dots.coordinates = makeMeshCoordinates(comm, numGlobalCoords);
  size_t numLocalCoords = dots.coordinates->getLocalLength();

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
 
    dots.weights.resize(weightDim);

    int wt = 0;
    float scale = 1.0;
    for (int i=0; i < weightDim; i++){
      dots.weights[i].resize(numLocalCoords);
      makeWeights(comm, dots.weights[i], weightTypes(wt++), scale, rank);

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
    printf("Zoltan_Initialize failed\n");
    exit(0);
  }

  struct Zoltan_Struct *zz;
  zz = Zoltan_Create(MPI_COMM_WORLD);
  
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "CHECK_GEOM", "0");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "PART");
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

  Zoltan_Set_Num_Obj_Fn(zz, getNumObj, &dots);
  Zoltan_Set_Obj_List_Fn(zz, getObjList, &dots);
  Zoltan_Set_Num_Geom_Fn(zz, getDim, &dots);
  Zoltan_Set_Geom_Multi_Fn(zz, getCoordinates, &dots);

  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids;
  ZOLTAN_ID_PTR exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  MEMORY_CHECK(doMemory && rank==0, "Before Zoltan_LB_Partition");

  if (rank == 0) std::cout << "Calling Zoltan_LB_Partition" << std::endl;
  aok = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
                            &numImport, &importGlobalGids, &importLocalGids,
                            &importProcs, &importToPart,
                            &numExport, &exportGlobalGids, &exportLocalGids,
                            &exportProcs, &exportToPart);
  if (rank == 0) std::cout << "Returned from Zoltan_LB_Partition" << std::endl;

  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_LB_Partition");

  /* Print the load-balance stats here */

  float *sumWgtPerPart = new float[numGlobalParts];
  float *gsumWgtPerPart = new float[numGlobalParts];
  for (int i = 0; i < numGlobalParts; i++) sumWgtPerPart[i] = 0.;

  for (size_t i = 0; i < numLocalCoords; i++)
    sumWgtPerPart[exportToPart[i]] += (weightDim ? dots.weights[0][i]: 1.);

  Teuchos::reduceAll<int, float>(*comm, Teuchos::REDUCE_SUM, numGlobalParts,
                                  sumWgtPerPart, gsumWgtPerPart);

  float maxSumWgtPerPart = 0.;
  float minSumWgtPerPart = std::numeric_limits<float>::max();
  float totWgt = 0.;
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

  delete dots.coordinates;
  for (int i = 0; i < weightDim; i++)
    dots.weights[i].clear();
  dots.weights.clear();

  MEMORY_CHECK(doMemory && rank==0, "After destroying input");

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
