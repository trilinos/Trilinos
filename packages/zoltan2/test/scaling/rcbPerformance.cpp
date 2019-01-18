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

/*! \file rcbPerformance.cpp
    \brief A test that can do large scale problems and time them.

    Geometry is a uniform mesh.
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

#include <Teuchos_CommandLineProcessor.hpp>

#include <vector>

using std::vector;
using std::bad_alloc;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::CommandLineProcessor;

typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Tpetra::Map<zlno_t, zgno_t, znode_t> tMap_t;

enum weightTypes{
  upDown,
  roundRobin,
  increasing,
  numWeightTypes
};

ArrayRCP<zscalar_t> makeWeights(
  const RCP<const Teuchos::Comm<int> > & comm,
  zlno_t len, weightTypes how, zscalar_t scale, int rank)
{
  zscalar_t *wgts = new zscalar_t [len];
  if (!wgts)
    throw bad_alloc();

  ArrayRCP<zscalar_t> weights(wgts, 0, len, true);

  if (how == upDown){
    zscalar_t val = scale + rank%2;
    for (zlno_t i=0; i < len; i++)
      wgts[i] = val;
  }
  else if (how == roundRobin){
    for (int i=0; i < 10; i++){
      zscalar_t val = (i + 10)*scale;
      for (zlno_t j=i; j < len; j += 10)
         weights[j] = val;
    }
  }
  else if (how == increasing){
    zscalar_t val = scale + rank;
    for (zlno_t i=0; i < len; i++)
      wgts[i] = val;
  }

  return weights;
}

/*! \brief Create a mesh of approximately the desired size.
 *
 *  We want 3 dimensions close to equal in length.
 */
const RCP<tMVector_t> getMeshCoordinates(
    const RCP<const Teuchos::Comm<int> > & comm,
    zgno_t numGlobalCoords)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  ssize_t xdim = static_cast<ssize_t>(floor(xdimf));
  ssize_t ydim = xdim;
  ssize_t zdim = numGlobalCoords / (xdim*ydim);
  ssize_t num=xdim*ydim*zdim;
  ssize_t diff = numGlobalCoords - num;
  ssize_t newdiff = 0;

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
      std::cout << "Warning: Difference " << diff*100 << " percent" << std::endl;
    std::cout << "Mesh size: " << xdim << "x" << ydim << "x" <<
      zdim << ", " << num << " vertices." << std::endl;
  }

  // Divide coordinates.

  ssize_t numLocalCoords = num / nprocs;
  ssize_t leftOver = num % nprocs;
  ssize_t gid0 = 0;

  if (rank <= leftOver)
    gid0 = zgno_t(rank) * (numLocalCoords+1);
  else
    gid0 = (leftOver * (numLocalCoords+1)) + 
           ((zgno_t(rank) - leftOver) * numLocalCoords);

  if (rank < leftOver)
    numLocalCoords++;

  ssize_t gid1 = gid0 + numLocalCoords;

  zgno_t *ids = new zgno_t [numLocalCoords];
  if (!ids)
    throw bad_alloc();
  ArrayRCP<zgno_t> idArray(ids, 0, numLocalCoords, true);

  for (ssize_t i=gid0; i < gid1; i++)
    *ids++ = zgno_t(i);   

  RCP<const tMap_t> idMap = rcp(
    new tMap_t(num, idArray.view(0, numLocalCoords), 0, comm));

  // Create a Tpetra::MultiVector of coordinates.

  zscalar_t *x = new zscalar_t [numLocalCoords*3]; 
  if (!x)
    throw bad_alloc();
  ArrayRCP<zscalar_t> coordArray(x, 0, numLocalCoords*3, true);

  zscalar_t *y = x + numLocalCoords;
  zscalar_t *z = y + numLocalCoords;

  zgno_t xStart = 0;
  zgno_t yStart = 0;
  zgno_t xyPlane = xdim*ydim;
  zgno_t zStart = gid0 / xyPlane;
  zgno_t rem = gid0 % xyPlane;
  if (rem > 0){
    yStart = rem / xdim;
    xStart = rem % xdim;
  }

  zlno_t next = 0;
  for (zscalar_t zval=zStart; next < numLocalCoords && zval < zdim; zval++){
    for (zscalar_t yval=yStart; next < numLocalCoords && yval < ydim; yval++){
      for (zscalar_t xval=xStart; next < numLocalCoords && xval < xdim; xval++){
        x[next] = xval;
        y[next] = yval;
        z[next] = zval;
        next++;
      }
      xStart = 0;
    }
    yStart = 0;
  }

  ArrayView<const zscalar_t> xArray(x, numLocalCoords);
  ArrayView<const zscalar_t> yArray(y, numLocalCoords);
  ArrayView<const zscalar_t> zArray(z, numLocalCoords);
  ArrayRCP<ArrayView<const zscalar_t> > coordinates =
    arcp(new ArrayView<const zscalar_t> [3], 0, 3);
  coordinates[0] = xArray;
  coordinates[1] = yArray;
  coordinates[2] = zArray;

  ArrayRCP<const ArrayView<const zscalar_t> > constCoords =
   coordinates.getConst();

  RCP<tMVector_t> meshCoords = rcp(new tMVector_t(
    idMap, constCoords.view(0,3), 3));

  return meshCoords;
}


int main(int narg, char *arg[])
{
  // MEMORY_CHECK(true, "Before initializing MPI");

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  MEMORY_CHECK(rank==0, "After initializing MPI");

  if (rank==0)
    std::cout << "Number of processes: " << nprocs << std::endl;

  // Default values
  double numGlobalCoords = 1000;
  int numTestCuts = 1;
  int nWeights = 0;
  string timingType("no_timers");
  string debugLevel("basic_status");
  string memoryOn("memoryOn");
  string memoryOff("memoryOff");
  string memoryProcs("0");
  bool doMemory=false;
  int numGlobalParts = nprocs;

  CommandLineProcessor commandLine(false, true);
  commandLine.setOption("size", &numGlobalCoords, 
    "Approximate number of global coordinates.");
  commandLine.setOption("testCuts", &numTestCuts, 
    "Number of test cuts to make when looking for bisector.");
  commandLine.setOption("numParts", &numGlobalParts, 
    "Number of parts (default is one per proc).");
  commandLine.setOption("nWeights", &nWeights, 
    "Number of weights per coordinate, zero implies uniform weights.");

  string balanceCount("balance_object_count");
  string balanceWeight("balance_object_weight");
  string mcnorm1("multicriteria_minimize_total_weight");
  string mcnorm2("multicriteria_balance_total_maximum");
  string mcnorm3("multicriteria_minimize_maximum_weight");

  string objective(balanceWeight);   // default

  string doc(balanceCount); doc.append(": ignore weights\n");

  doc.append(balanceWeight); doc.append(": balance on first weight\n");

  doc.append(mcnorm1);
  doc.append(": given multiple weights, balance their total.\n");

  doc.append(mcnorm3);
  doc.append(": given multiple weights, balance the maximum for each coordinate.\n");

  doc.append(mcnorm2);
  doc.append(": given multiple weights, balance the L2 norm of the weights.\n");

  commandLine.setOption("objective", &objective,  doc.c_str());

  commandLine.setOption("timers", &timingType,
    "no_timers, micro_timers, macro_timers, both_timers, test_timers");

  commandLine.setOption("debug", &debugLevel,
   "no_status, basic_status, detailed_status, verbose_detailed_status");

  commandLine.setOption(memoryOn.c_str(), memoryOff.c_str(), &doMemory,
    "do memory profiling");

  commandLine.setOption("memoryProcs", &memoryProcs,
   "list of processes that output memory usage");

  CommandLineProcessor::EParseCommandLineReturn rc = 
    commandLine.parse(narg, arg);

  if (rc != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL){
    if (rc == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED){
      if (rank==0)
        std::cout << "PASS" << std::endl;
      return 1;
    }
    else{
      if (rank==0)
        std::cout << "FAIL" << std::endl;
      return 0;
    }
  }

  //MEMORY_CHECK(doMemory && rank==0, "After processing parameters");

  zgno_t globalSize = static_cast<zgno_t>(numGlobalCoords);

  RCP<tMVector_t> coordinates = getMeshCoordinates(comm, globalSize);
  size_t numLocalCoords = coordinates->getLocalLength();

#if 0
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p==rank){
      std::cout << "Rank " << rank << ", " << numLocalCoords << "coords" << std::endl;
      const zscalar_t *x = coordinates->getData(0).getRawPtr();
      const zscalar_t *y = coordinates->getData(1).getRawPtr();
      const zscalar_t *z = coordinates->getData(2).getRawPtr();
      for (zlno_t i=0; i < numLocalCoords; i++)
        std::cout << " " << x[i] << " " << y[i] << " " << z[i] << std::endl;
    }
    std::cout.flush();
    comm->barrier();
  }
#endif

  Array<ArrayRCP<zscalar_t> > weights(nWeights);

  if (nWeights > 0){
    int wt = 0;
    zscalar_t scale = 1.0;
    for (int i=0; i < nWeights; i++){
      weights[i] = 
        makeWeights(comm, numLocalCoords, weightTypes(wt++), scale, rank);
      if (wt == numWeightTypes){
        wt = 0;
        scale++;
      }
    }
  }

  MEMORY_CHECK(doMemory && rank==0, "After creating input");

  // Create an input adapter.
  const RCP<const tMap_t> &coordmap = coordinates->getMap();
  ArrayView<const zgno_t> ids = coordmap->getNodeElementList();
  const zgno_t *globalIds = ids.getRawPtr();
  
  size_t localCount = coordinates->getLocalLength();
  typedef Zoltan2::BasicVectorAdapter<tMVector_t> inputAdapter_t;
  RCP<inputAdapter_t> ia;
  
  if (nWeights == 0){
    ia = rcp(new inputAdapter_t (localCount, globalIds, 
      coordinates->getData(0).getRawPtr(), coordinates->getData(1).getRawPtr(),
      coordinates->getData(2).getRawPtr(), 1,1,1));
  }
  else{
    vector<const zscalar_t *> values(3);
    for (int i=0; i < 3; i++)
      values[i] = coordinates->getData(i).getRawPtr();
    vector<int> valueStrides(0);  // implies stride is one
    vector<const zscalar_t *> weightPtrs(nWeights);
    for (int i=0; i < nWeights; i++)
      weightPtrs[i] = weights[i].getRawPtr();
    vector<int> weightStrides(0); // implies stride is one

    ia = rcp(new inputAdapter_t (localCount, globalIds, 
      values, valueStrides, weightPtrs, weightStrides));
  }

  MEMORY_CHECK(doMemory && rank==0, "After creating input adapter");

  // Parameters

  Teuchos::ParameterList params;

  if (timingType != "no_timers"){
    params.set("timer_output_stream" , "std::cout");
    params.set("timer_type" , timingType);
  }

  if (doMemory){
    params.set("memory_output_stream" , "std::cout");
    params.set("memory_procs" , memoryProcs);
  }

  params.set("debug_output_stream" , "std::cerr");
  params.set("debug_procs" , "0");

  if (debugLevel != "basic_status"){
    params.set("debug_level" , debugLevel);
  }

  params.set("algorithm", "rcb");
  params.set("partitioning_objective", objective);
  double tolerance = 1.1;
  params.set("imbalance_tolerance", tolerance );

  if (numGlobalParts != nprocs)
    params.set("num_global_parts" , numGlobalParts);

  if (rank==0){
    std::cout << "Number of parts: " << numGlobalParts << std::endl;
  }

  // Create a problem, solve it, and display the quality.

  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&(*ia), &params);

  problem.solve();

  comm->barrier();

  problem.printTimers();

  comm->barrier();

  if (rank == 0){
    std::cout << "PASS" << std::endl;
  }

  return 0;
}

