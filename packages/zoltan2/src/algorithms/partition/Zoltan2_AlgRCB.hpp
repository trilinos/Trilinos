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

/*! \file Zoltan2_AlgRCB.hpp
    \brief Contains the recursive coordinate bisection algorthm.
*/

#ifndef _ZOLTAN2_ALGRCB_HPP_
#define _ZOLTAN2_ALGRCB_HPP_

#include <Zoltan2_AlgRCB_methods.hpp>
#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_Exceptions.hpp>

#include <Teuchos_ParameterList.hpp>

namespace Zoltan2{

/*! \brief Recursive coordinate bisection algorithm.
 *
 *  \param env   library configuration and problem parameters
 *  \param comm the communicator for the problem
 *  \param coords    a CoordinateModel with user data
 *  \param solution  a PartitioningSolution, on input it 
 *      contains part information, on return it also contains 
 *      the solution and quality metrics.
 *                    
 *   \todo timing and memory usage profiling and debug messages
 *   \todo  catch errors and pass back
 *   \todo write the rcb tree back to the solution
 *   \todo for "repartition", start with the tree in the solution
 *   \todo  work on performance issues.  Some of the global
 *               communication can probably be consolidated
 *                into fewer messages.
 *
 * The algorithm is documented in \ref rcbPage.  Please document
 * changes at this page.
 */

template <typename Adapter>
void AlgRCB(
  const RCP<const Environment> &env,
  RCP<Comm<int> > &problemComm,
  const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &coords, 
  RCP<PartitioningSolution<Adapter> > &solution
) 
{
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL

  Z2_THROW_EXPERIMENTAL("Zoltan2 RCB is strictly experimental software "
                        "due to performance problems in its use of Tpetra.")

#else  // INCLUDE_ZOLTAN2_EXPERIMENTAL

  typedef typename Adapter::node_t node_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  // Make a copy of communicator because
  // we subdivide the communicator during the algorithm.

  RCP<Comm<int> > comm = problemComm->duplicate();

  std::bitset<NUM_RCB_PARAMS> params;

  env->debug(DETAILED_STATUS, string("Entering AlgPartRCB"));

  const Teuchos::ParameterList &pl = env->getParameters();

  env->timerStart(BOTH_TIMERS, "RCB set up");

  ////////////////////////////////////////////////////////
  // Partitioning problem parameters of interest:
  //    objective
  //    imbalance_tolerance

  multiCriteriaNorm mcnorm = normBalanceTotalMaximum;
  string obj;

  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("partitioning_objective");
  if (pe)
    obj = pe->getValue(&obj);

  if (!pe){
    params.set(rcb_balanceWeight);
    mcnorm = normBalanceTotalMaximum;
  }
  else if (obj == string("balance_object_count")){
    params.set(rcb_balanceCount);
  }
  else if (obj == string("multicriteria_minimize_total_weight")){
    params.set(rcb_minTotalWeight);
    mcnorm = normMinimizeTotalWeight;
  }
  else if (obj == string("multicriteria_minimize_maximum_weight")){
    params.set(rcb_minMaximumWeight);
    mcnorm = normMinimizeMaximumWeight;
  }
  else if (obj == string("multicriteria_balance_total_maximum")){
    params.set(rcb_balanceTotalMaximum);
    mcnorm = normBalanceTotalMaximum;
  }
  else{
    params.set(rcb_balanceWeight);
    mcnorm = normBalanceTotalMaximum;
  }

  scalar_t imbalanceTolerance = .1;
  pe = pl.getEntryPtr("imbalance_tolerance");
  if (pe){
    double tol;
    tol = pe->getValue(&tol);
    imbalanceTolerance = tol - 1.0;
  }

  if (imbalanceTolerance <= 0)
    imbalanceTolerance = 10e-4;  // TODO - what's a good choice

  ////////////////////////////////////////////////////////
  // Geometric partitioning problem parameters of interest:
  //    average_cuts
  //    rectilinear_blocks
  //    bisection_num_test_cuts (experimental)

  int val = 0;
  pe = pl.getEntryPtr("average_cuts");
  if (pe)
    val = pe->getValue(&val);

  if (val == 1)
    params.set(rcb_averageCuts);

  val = 0;
  pe = pl.getEntryPtr("rectilinear_blocks");
  if (pe)
    val = pe->getValue(&val);

  if (val == 1)
    params.set(rcb_rectilinearBlocks);

  int numTestCuts = 1;
  pe = pl.getEntryPtr("bisection_num_test_cuts");
  if (pe)
    numTestCuts = pe->getValue(&numTestCuts);

  ////////////////////////////////////////////////////////
  // From the CoordinateModel we need:
  //    coordinate values
  //    coordinate weights, if any
  //    coordinate global Ids

  typedef StridedData<lno_t, scalar_t> input_t;

  int coordDim = coords->getCoordinateDim();
  int weightDim = coords->getCoordinateWeightDim();
  size_t numLocalCoords = coords->getLocalNumCoordinates();
  global_size_t numGlobalCoords = coords->getGlobalNumCoordinates();

  ArrayView<const gno_t> gnos;
  ArrayView<input_t>     xyz;
  ArrayView<input_t>     wgts;

  coords->getCoordinates(gnos, xyz, wgts);

  Array<ArrayRCP<const scalar_t> > values(coordDim);
  for (int dim=0; dim < coordDim; dim++){
    ArrayRCP<const scalar_t> ar;
    xyz[dim].getInputArray(ar);
    values[dim] = ar;
  }

  int criteriaDim = (weightDim ? weightDim : 1);
  bool ignoreWeights = params.test(rcb_balanceCount);

  if (criteriaDim > 1 && ignoreWeights)
    criteriaDim = 1;

  ArrayRCP<bool> uniformWeights(new bool [criteriaDim], 0, criteriaDim, true);
  Array<ArrayRCP<const scalar_t> > weights(criteriaDim);

  if (weightDim == 0 || ignoreWeights)
    uniformWeights[0] = true;
  else{
    for (int wdim = 0; wdim < weightDim; wdim++){
      if (wgts[wdim].size() == 0){
        uniformWeights[wdim] = true;
      }
      else{
        uniformWeights[wdim] = false;
        ArrayRCP<const scalar_t> ar;
        wgts[wdim].getInputArray(ar);
        weights[wdim] = ar;
      }
    }
  }

  if (env->doStatus() && (numGlobalCoords < 500)){
    ostringstream oss;
    oss << "Problem: ";
    for (size_t i=0; i < numLocalCoords; i++){
      oss << gnos[i] << " (";
      for (int dim=0; dim < coordDim; dim++)
        oss << (xyz[dim])[i] << " ";
      oss << ") ";
    }

    env->debug(VERBOSE_DETAILED_STATUS, oss.str());
  }

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.
  // If the part sizes for a given criteria are not uniform,
  // then they are values that sum to 1.0.

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

  Array<bool> uniformParts(criteriaDim);
  Array<ArrayRCP<scalar_t> > partSizes(criteriaDim);

  for (int wdim = 0; wdim < criteriaDim; wdim++){
    if (solution->criteriaHasUniformPartSizes(wdim)){
      uniformParts[wdim] = true;
    }
    else{
      scalar_t *tmp = new scalar_t [numGlobalParts];
      env->localMemoryAssertion(__FILE__, __LINE__, numGlobalParts, tmp) ;
    
      for (size_t i=0; i < numGlobalParts; i++){
        tmp[i] = solution->getCriteriaPartSize(wdim, i);
      }

      partSizes[wdim] = arcp(tmp, 0, numGlobalParts);
    }
  }

  // It may not be possible to solve the partitioning problem
  // if we have multiple weight dimensions with part size
  // arrays that differ. So let's be aware of this possibility.

  bool multiplePartSizeSpecs = false;

  if (criteriaDim > 1){
    for (int wdim1 = 0; wdim1 < criteriaDim; wdim1++)
      for (int wdim2 = wdim1+1; wdim2 < criteriaDim; wdim2++)
        if (!solution->criteriaHaveSamePartSizes(wdim1, wdim2)){
          multiplePartSizeSpecs = true;
          break;
        }
  }
  
  if (multiplePartSizeSpecs)
    params.set(rcb_multiplePartSizeSpecs);

  ////////////////////////////////////////////////////////
  // Create the distributed data for the algorithm.
  //
  // It is a multivector containing one vector for each coordinate
  // dimension, plus a vector for each weight dimension that is not
  // uniform.

  typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> mvector_t;

  int multiVectorDim = coordDim;
  for (int wdim = 0; wdim < criteriaDim; wdim++)
    if (!uniformWeights[wdim]) multiVectorDim++;

  gno_t gnoMin, gnoMax;
  coords->getIdentifierMap()->getGnoRange(gnoMin, gnoMax);

  RCP<map_t> map;
  try{
    map = rcp(new map_t(numGlobalCoords, gnos, gnoMin, comm));
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  typedef ArrayView<const scalar_t> coordList_t;

  coordList_t *avList = new coordList_t [multiVectorDim];

  for (int dim=0; dim < coordDim; dim++)
    avList[dim] = values[dim].view(0, numLocalCoords);

  for (int wdim=0, idx=coordDim; wdim < criteriaDim; wdim++)
    if (!uniformWeights[wdim])
      avList[idx++] = weights[wdim].view(0, numLocalCoords);

  ArrayRCP<const ArrayView<const scalar_t> > vectors =
    arcp(avList, 0, multiVectorDim);

  RCP<mvector_t> mvector;

  try{
    mvector = rcp(new mvector_t(
      map, vectors.view(0, multiVectorDim), multiVectorDim));
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->timerStop(BOTH_TIMERS, "RCB set up");

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  // The algorithm
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  partId_t part0 = 0;
  partId_t part1 = numGlobalParts-1;
  int sanityCheck = numGlobalParts;
  int groupSize = comm->getSize();
  int rank = comm->getRank();

  long imbalanceReductionFactor(1);
  long nparts = numGlobalParts;
  while ((nparts >>= 1) != 0) imbalanceReductionFactor++;

  imbalanceTolerance /= imbalanceReductionFactor;

  int iteration = 1;

  env->memory("RCB algorithm set up");

  env->timerStart(MACRO_TIMERS, "Parallel RCB");

  while (part1>part0 && groupSize>1 && numGlobalCoords>0 && sanityCheck--){

    ////////////////////////////////////////////////////////
    // Which coordinates are left and which are right?

    Array<unsigned char> lrflags(numLocalCoords);
    scalar_t cutValue=0;  // TODO eventually save this for user
    int cutDimension=0;
    scalar_t imbalance=0, weightLeft=0, weightRight=0;
    partId_t leftHalfNumParts=0;

    env->timerStart(MICRO_TIMERS, "Find cut", iteration, 2);

    try{
      determineCut<mvector_t, Adapter>(env, comm, 
        params, numTestCuts, imbalanceTolerance,
        coordDim, mvector, 
        uniformWeights.view(0,criteriaDim), mcnorm, solution,
        part0, part1,
        lrflags.view(0, numLocalCoords), 
        cutDimension, cutValue, imbalance, leftHalfNumParts,
        weightLeft, weightRight);
    }
    Z2_FORWARD_EXCEPTIONS

    env->timerStop(MICRO_TIMERS, "Find cut", iteration, 2);

    // Do we have empty left or right halves?

    bool skipLeft = (weightLeft == 0);
    bool skipRight = (weightRight == 0);

    ////////////////////////////////////////////////////////
    // Migrate the multivector of data.

    int leftHalfNumProcs=0;

    env->timerStart(MICRO_TIMERS, "Migrate", iteration, 2);

    try{ // on return mvector has my new data

      migrateData<mvector_t>( env, comm, lrflags.view(0,numLocalCoords), 
        mvector, leftHalfNumProcs);
    }
    Z2_FORWARD_EXCEPTIONS

    env->timerStop(MICRO_TIMERS, "Migrate", iteration, 2);

    env->localBugAssertion(__FILE__, __LINE__, "num procs in half",
      leftHalfNumProcs > 0 && leftHalfNumProcs < groupSize,
      BASIC_ASSERTION);

    bool inLeftHalf = (rank < leftHalfNumProcs);

    if ((inLeftHalf && skipLeft) || (!inLeftHalf && skipRight)){
      groupSize = 1;
      numLocalCoords = 0;
      continue;
    }

    ////////////////////////////////////////////////////////
    // Divide into two subgroups.

    env->timerStart(MICRO_TIMERS, "Create sub group, sub data", iteration, 2);

    int *ids = NULL;

    if (rank < leftHalfNumProcs){
      groupSize = leftHalfNumProcs;
      ids = new int [groupSize];
      env->localMemoryAssertion(__FILE__, __LINE__, groupSize, ids);
      for (int i=0; i < groupSize; i++)
        ids[i] = i;
      part1 = part0 + leftHalfNumParts - 1;
    }
    else {
      groupSize = comm->getSize() - leftHalfNumProcs;
      rank -= leftHalfNumProcs;
      ids = new int [groupSize];
      env->localMemoryAssertion(__FILE__, __LINE__, groupSize, ids);
      for (int i=0; i < groupSize; i++)
        ids[i] = i + leftHalfNumProcs;
      part0 += leftHalfNumParts;
    }

    ArrayView<const int> idView(ids, groupSize);
    comm = comm->createSubcommunicator(idView);

    delete [] ids;

    ////////////////////////////////////////////////////////
    // Create a new multivector for my smaller group.

    ArrayView<const gno_t> gnoList = mvector->getMap()->getNodeElementList();
    size_t localSize = mvector->getLocalLength();
  
    // Tpetra will calculate the globalSize.
    size_t globalSize = Teuchos::OrdinalTraits<size_t>::invalid();
  
    RCP<map_t> subMap;
    try{
      subMap= rcp(new map_t(globalSize, gnoList, 0, comm));
    }
    Z2_THROW_OUTSIDE_ERROR(*env)

    coordList_t *avSubList = new coordList_t [multiVectorDim];
  
    for (int dim=0; dim < multiVectorDim; dim++)
      avSubList[dim] = mvector->getData(dim).view(0, localSize);
  
    ArrayRCP<const ArrayView<const scalar_t> > subVectors =
      arcp(avSubList, 0, multiVectorDim);
  
    RCP<mvector_t> subMvector;
  
    try{
      subMvector = rcp(new mvector_t(
        subMap, subVectors.view(0, multiVectorDim), multiVectorDim));
    }
    Z2_THROW_OUTSIDE_ERROR(*env)

    env->timerStop(MICRO_TIMERS, "Create sub group, sub data", iteration, 2);
  
    mvector = subMvector;

    numLocalCoords = mvector->getLocalLength();
    numGlobalCoords = mvector->getGlobalLength();

    iteration++;

    env->memory("New subgroup data created");
  } 

  env->timerStop(MACRO_TIMERS, "Parallel RCB");

  env->localBugAssertion(__FILE__, __LINE__, "partitioning failure", 
    sanityCheck, BASIC_ASSERTION);

  ArrayRCP<partId_t> partId;

  if (numLocalCoords > 0){
    partId_t *tmp = new partId_t [numLocalCoords];
    env->localMemoryAssertion(__FILE__, __LINE__, numLocalCoords, tmp);
    partId = arcp(tmp, 0, numLocalCoords, true);
  }

  env->memory("Solution array created");

  if ((part1 > part0) && (numLocalCoords > 0)){ // Serial partitioning

    // scalar_t cutValue;   TODO
    // int cutDimension;
    // scalar_t imbalance;

    env->timerStart(MACRO_TIMERS, "Serial RCB");

    try{
      ArrayView<lno_t> emptyIndex;

      serialRCB<mvector_t, Adapter>(env, 1, params,
        numTestCuts, imbalanceTolerance,
        coordDim, mvector, emptyIndex, 
        uniformWeights.view(0,criteriaDim), solution,
        part0, part1, partId.view(0,numLocalCoords));
    }
    Z2_FORWARD_EXCEPTIONS

    env->timerStop(MACRO_TIMERS, "Serial RCB");

  }
  else{
    for (lno_t i=0; i < partId.size(); i++)
      partId[i] = part0;
  }

  ////////////////////////////////////////////////////////
  // Done: update the solution

  ArrayRCP<const gno_t> gnoList = 
    arcpFromArrayView(mvector->getMap()->getNodeElementList());

  if (env->getDebugLevel() >= VERBOSE_DETAILED_STATUS && 
     (numGlobalCoords < 500)){
    ostringstream oss;
    oss << "Solution: ";
    for (gno_t i=0; i < gnoList.size(); i++)
      oss << gnoList[i] << " (" << partId[i] << ") ";
    
    env->debug(VERBOSE_DETAILED_STATUS, oss.str());
  }

  solution->setParts(gnoList, partId, false);
#endif // INCLUDE_ZOLTAN2_EXPERIMENTAL
}

} // namespace Zoltan2

#endif
