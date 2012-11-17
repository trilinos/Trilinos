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

/*! \file Zoltan2_AlgRCB_methods.hpp
    \brief Methods written for the RCB algorithm
*/

#ifndef _ZOLTAN2_ALGRCB_METHODS_HPP_
#define _ZOLTAN2_ALGRCB_METHODS_HPP_

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Metric.hpp>

#include <Tpetra_Vector.hpp>

#include <sstream>
#include <string>
#include <cmath>
#include <bitset>

namespace Teuchos{

/*! \brief A Teuchos::MPIComm reduction operation.
 *
 *  Given a list of values T, some number are added, the next group
 *      is min'ed, the final group is max'ed.
 */

template <typename Ordinal, typename T>
  class Zoltan2_RCBOperation : public ValueTypeReductionOp<Ordinal,T>
{
private:
  Ordinal numSum_, numMin_, numMax_;

public:
  /*! \brief Default Constructor
   */
  Zoltan2_RCBOperation():numSum_(0), numMin_(0), numMax_(0){}

  /*! \brief Constructor
   *   \param nsum  the count of how many sums will be computed at the
   *             start of the list.
   *   \param nmin  following the sums, this many minimums will be computed.
   *   \param nmax  following the minimums, this many maximums will be computed.
   */
  Zoltan2_RCBOperation(Ordinal nsum, Ordinal nmin, Ordinal nmax):
    numSum_(nsum), numMin_(nmin), numMax_(nmax){}
  
  /*! \brief Implement Teuchos::ValueTypeReductionOp interface
   */
  void reduce( const Ordinal count, const T inBuffer[], T inoutBuffer[]) const
  {
    Ordinal next=0;
    for (Ordinal i=0; i < numSum_; i++, next++)
      inoutBuffer[next] += inBuffer[next];

    for (Ordinal i=0; i < numMin_; i++, next++)
      if (inoutBuffer[next] > inBuffer[next])
        inoutBuffer[next] = inBuffer[next];

    for (Ordinal i=0; i < numMax_; i++, next++)
      if (inoutBuffer[next] < inBuffer[next])
        inoutBuffer[next] = inBuffer[next];
  }

  /*! \brief Convenience method.
   */
  ArrayRCP<T> getInitializedBuffer(T minExtreme, T maxExtreme) const
  {
    Ordinal len = numSum_+numMin_+numMax_;
    T *buf = new T [len];
    if (!buf)
      throw std::bad_alloc();

    Ordinal next=0;
    for (Ordinal i=0; i < numSum_; i++)
      buf[next++] = 0;
    for (Ordinal i=0; i < numMin_; i++)
      buf[next++] = minExtreme;
    for (Ordinal i=0; i < numMax_; i++)
      buf[next++] = maxExtreme;

    ArrayRCP<T> returnBuf = arcp(buf, 0, len, true);
    return returnBuf;
  }
};
} // namespace Teuchos

namespace Zoltan2{

/*! \brief The boolean parameters of interest to the RCB algorithm.
 */
enum rcbParams {
  rcb_balanceCount,            /*!< objective = balance_object_count */
  rcb_balanceWeight,          /*!< objective = balance_object_weight */
  rcb_minTotalWeight,      /*!< objective = mc_minimize_total_weight */
  rcb_minMaximumWeight,  /*!< objective = mc_minimize_maximum_weight */
  rcb_balanceTotalMaximum, /*!< objective = mc_balance_total_maximum */
  rcb_averageCuts,          /*!< averageCuts = yes */
  rcb_rectilinearBlocks,    /*!< rectilinearBlocks = yes */
  rcb_multiplePartSizeSpecs,  /*!< multicriteria w/differing part sizes */
  NUM_RCB_PARAMS
};

/*! \brief During partitioning, flags are stored in unsigned char arrays.
 *
 *  Flag is also used to store a region number, but there are at most
 *  251 regions.  Therefore region number will not conflict with leftFlag
 *  and rightFlag. (Number of regions is number of test cuts plus one.)
 */

enum leftRightFlag{
  unsetFlag = 0xfd,       /*!< 253 */
  leftFlag = 0xfe,        /*!< 254 */
  rightFlag = 0xff        /*!< 255 */
};

// TODO: RCB has several helper methods.
// Do we want a convention for naming algorithm methods?  Do
// we want sub-namespaces for algorithms? 

/*! \brief Determine the fraction of work in the left half.
 *
 *   \param env the Environment for the application.
 *   \param part0  is the first part of the parts to bisected.
 *   \param part1  is the last part of the parts to bisected.
 *   \param solution which contains the part sizes.
 *   \param fractionLeft on return <tt>fractionLeft[w]</tt> is the
 *             fraction of work wanted in the left half for weight
 *             dimension \c w.
 *   \param numPartsLeftHalf on return is the number of parts in the
 *               left half.
 */

template <typename Adapter>
  void getFractionLeft(
    const RCP<const Environment> &env,
    partId_t part0,
    partId_t part1,
    const RCP<PartitioningSolution<Adapter> > &solution,
    ArrayRCP<double> &fractionLeft,
    partId_t &numPartsLeftHalf)
{
  env->debug(DETAILED_STATUS, string("Entering getFractionLeft"));
  partId_t numParts = part1 - part0 + 1;
  // TODO In LDRD, substitute call to machine model for
  // TODO computation of numPartsLeftHalf below.
  numPartsLeftHalf = numParts / 2;
  partId_t left0 = part0;   // First part in left half
  partId_t left1 = left0 + numPartsLeftHalf - 1;  // Last part in left half
  partId_t right0 = left1 + 1;  // First part in right half
  partId_t right1 = part1;  // Last part in right half

  int weightDim = solution->getNumberOfCriteria();
  fractionLeft = arcp(new double [weightDim], 0, weightDim);

  for (int wdim=0; wdim<weightDim; wdim++){
    if (solution->criteriaHasUniformPartSizes(wdim)){
      fractionLeft[wdim] = double(numPartsLeftHalf) / double(numParts);
    }
    else{
      fractionLeft[wdim] = 0;
      for(int partId=left0; partId <= left1; partId++){
        fractionLeft[wdim] += solution->getCriteriaPartSize(wdim, partId);
      }
      double total = fractionLeft[wdim];
      for(int partId=right0; partId <= right1; partId++){
        total += solution->getCriteriaPartSize(wdim, partId);
      }
      fractionLeft[wdim] /= total;
    }
  }
  env->debug(DETAILED_STATUS, string("Exiting getFractionLeft"));
}

/*! \brief Choose a coordinate dimension for the next cut.
 *
 *   \param env the Environment for the application.
 *   \param comm the communicator
 *   \param coordDim the first \c coordDim vectors in the \c vectors
 *              list are coordinates, the rest are weights.
 *   \param vectors lists of coordinates and non-uniform weights
 *   \param index is the index into the \c vectors arrays for the
 *              coordinates to be included in the partitioning.
 *              If <tt>index.size()</tt> is zero, then include
 *              all coordinates.
 *   \param dimension on return is the dimension to cut
 *   \param minCoord on the return is the minimum coordinate value in
 *            this dimension.
 *   \param maxCoord on the return is the maximum coordinate value in
 *            this dimension.
 */

template <typename mvector_t>
  void getCutDimension(
    const RCP<const Environment> &env,
    const RCP<Comm<int> > &comm,
    int coordDim,
    const RCP<mvector_t> &vectors,
    ArrayView<typename mvector_t::local_ordinal_type> index,
    int &dimension,                        // output
    typename mvector_t::scalar_type &minCoord,      // output
    typename mvector_t::scalar_type &maxCoord)      // output
{
  env->debug(DETAILED_STATUS, string("Entering getCutDimension"));
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;

  int nprocs = comm->getSize();
  bool useIndices = index.size() > 0;
  lno_t numLocalCoords = 0;

  if (useIndices)
    numLocalCoords = index.size();
  else
    numLocalCoords = vectors->getLocalLength();

  Array<scalar_t> spans(coordDim*2);
  int next = 0;

  if (useIndices){
    if (numLocalCoords > 0){
      for (int dim=0; dim < coordDim; dim++){
        const scalar_t *val = vectors->getData(dim).getRawPtr();
        scalar_t v = val[index[0]];
        scalar_t min = v;
        scalar_t max = v;
        for (lno_t i=1; i < numLocalCoords; i++){
          v = val[index[i]];
          if (v > max)
            max = v;
          else if (v < min)
            min = v;
        }
  
        spans[next++] = min;
        spans[next++] = max * -1.0;
      }
    }
    else{
      for (int dim=0; dim < coordDim; dim++){
        spans[next++] = numeric_limits<scalar_t>::max();
        spans[next++] = numeric_limits<scalar_t>::min() * -1.0;
      }
    }
  }
  else{
    for (int dim=0; dim < coordDim; dim++){
      const scalar_t *val = vectors->getData(dim).getRawPtr();
      pair<scalar_t, scalar_t> minMax = 
        z2LocalMinMax<scalar_t>(val, numLocalCoords);
      spans[next++] = minMax.first;
      spans[next++] = minMax.second * -1.0;
    }
  }

  if (nprocs > 1){
    Array<scalar_t> globalValues(coordDim*2);
    try{
      reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_MIN, coordDim*2,
        spans.getRawPtr(), globalValues.getRawPtr());
    }
    Z2_THROW_OUTSIDE_ERROR(*env)

    for (int i=0; i < coordDim*2; i++)
      spans[i] = globalValues[i];
  }

  scalar_t maxSpan = 0;
  dimension = -1;
  next = 0;

  for (int dim=0; dim < coordDim; dim++){
    scalar_t min = spans[next++];
    scalar_t max = spans[next++] * -1.0;
    scalar_t newSpan = max - min;
    if (newSpan > maxSpan){
      maxSpan = newSpan;
      dimension = dim;
      minCoord = min;
      maxCoord = max;
    }
  }
  env->debug(DETAILED_STATUS, string("Exiting getCutDimension"));
}

/*! \brief Migrate coordinates and weights to new processes.
 *
 *   \param env the Environment for the application.
 *   \param comm the communicator for this step.
 *   \param lrflags each element in this array must have the value leftFlag 
 *         or rightFlag to
 *        indicate whether the corresponding coordinate is on the left of
 *        on the right of the cut. 
 *   \param vectors is a list of the data to be migrated.  On return,
 *        \c vectors contains the new data belonging to this process.
 *   param leftNumProcs on return is the number of processes receiving 
 *         the left half.
 *
 *   \return the number of processes in the left half is returned.
 */

template <typename mvector_t>
  void migrateData(
    const RCP<const Environment> &env, const RCP<Comm<int> > &comm,
    ArrayView<unsigned char> lrflags,
    RCP<mvector_t> &vectors,    // on return is the new data
    int &leftNumProcs)          // on return is num procs with left data
{
  env->debug(DETAILED_STATUS, string("Entering migrateData"));
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;

  int nprocs = comm->getSize();
  size_t nobj = vectors->getLocalLength();
  size_t nGlobalObj = vectors->getGlobalLength();

  env->localBugAssertion(__FILE__, __LINE__, "migrateData input", 
    nprocs>1 && size_t(lrflags.size())==nobj, DEBUG_MODE_ASSERTION);

  gno_t myNumLeft= 0, numLeft;
  for (size_t i=0; i < nobj; i++)
    if (lrflags[i] == leftFlag)
      myNumLeft++;

  try{
    reduceAll<int, gno_t>(
      *comm, Teuchos::REDUCE_SUM, 1, &myNumLeft, &numLeft);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)
    
  scalar_t leftFraction = scalar_t(numLeft)/scalar_t(nGlobalObj);
  leftNumProcs = static_cast<int>(floor((nprocs*leftFraction)+.5));

  if (leftNumProcs == 0)
    leftNumProcs++;
  else if (leftNumProcs == nprocs)
    leftNumProcs--;

  ///////////////////////////////////////////////////////
  // Get a list of my new global numbers.

  int *sendCount = new int [nprocs];
  env->localMemoryAssertion(__FILE__, __LINE__, nprocs, sendCount) ;
  memset(sendCount, 0, sizeof(int) * nprocs);
  ArrayView<int> sendCountView(sendCount, nprocs);
  ArrayView<gno_t> sendBufView;

  if (nobj > 0){
    int *procId = new int [nobj];
    env->localMemoryAssertion(__FILE__, __LINE__, nobj, procId) ;
    int leftProc0 = 0;
    int rightProc0 = leftProc0 + leftNumProcs;
  
    int nextLeftProc = leftProc0;
    int nextRightProc = rightProc0;
    int *p = procId;
  
    for (size_t i=0; i < nobj; i++){
      if (lrflags[i] == leftFlag){
        if (nextLeftProc == rightProc0)
          nextLeftProc = leftProc0;
  
        sendCount[nextLeftProc]++;
        *p++ = nextLeftProc++;
      }
      else{
        if (nextRightProc == nprocs)
          nextRightProc = rightProc0;
        sendCount[nextRightProc]++;
  
        *p++ = nextRightProc++;
      }
    }
  
    gno_t *sendOffset = new gno_t [nprocs];
    env->localMemoryAssertion(__FILE__, __LINE__, nprocs, sendOffset) ;
    sendOffset[0] = 0;
    for (int i=0; i < nprocs-1; i++)
      sendOffset[i+1] = sendOffset[i] + sendCount[i];

    gno_t *sendBuf = new gno_t [nobj];
    env->localMemoryAssertion(__FILE__, __LINE__, nobj, sendBuf) ;
    sendBufView = ArrayView<gno_t>(sendBuf, nobj);

    ArrayView<const gno_t> gnoList = vectors->getMap()->getNodeElementList();

    for (size_t i=0; i < nobj; i++){
      int proc = procId[i];
      lno_t offset = sendOffset[proc]++;

      sendBuf[offset] = gnoList[i];
    }

    delete [] sendOffset;
    delete [] procId;
  }

  ArrayRCP<gno_t> recvBuf;
  ArrayRCP<int> recvCount;

  try{
    AlltoAllv<gno_t>(*comm, *env,
      sendBufView, sendCountView,
      recvBuf, recvCount);
  }
  Z2_FORWARD_EXCEPTIONS

  if (nobj > 0){
    delete [] sendBufView.getRawPtr();
    delete [] sendCountView.getRawPtr();
  }

  ///////////////////////////////////////////////////////
  // Migrate the multivector of data.

  gno_t numMyNewGnos = 0;
  for (int i=0; i < nprocs; i++)
    numMyNewGnos += recvCount[i];

  RCP<const mvector_t> newMultiVector;
  RCP<const mvector_t> constInput = rcp_const_cast<const mvector_t>(vectors);

  try{
    newMultiVector = XpetraTraits<mvector_t>::doMigration(
      constInput, numMyNewGnos, recvBuf.getRawPtr());
  }
  Z2_FORWARD_EXCEPTIONS

  vectors = rcp_const_cast<mvector_t>(newMultiVector);
  env->memory("Former problem data replaced with new data");
  env->debug(DETAILED_STATUS, string("Exiting migrateData"));
}

template <typename lno_t, typename scalar_t>
  scalar_t getCoordWeight(lno_t id, 
    multiCriteriaNorm mcnorm, 
    ArrayView<StridedData<lno_t, scalar_t> > weights)
{
  scalar_t wgtValue = 1.0;
  size_t weightDim = weights.size();
    
  if (weightDim > 1){
    Array<scalar_t> coordWeights(weightDim, 1.0);
    for (size_t wdim=0; wdim < weightDim; wdim++){
      if (weights[wdim].size() > 0)
        coordWeights[wdim] = weights[wdim][id];
    }

    wgtValue = normedWeight<scalar_t>(coordWeights.view(0,weightDim), mcnorm);
  }
  else if (weights[0].size() > 0){
    wgtValue = weights[0][id];
  }

  return wgtValue;
}

/*! \brief Solve partitioning if there are empty parts.
 *  \param fractionLeft  amount of work in left half
 *  \param minCoord       minimum coordinate value
 *  \param maxCoord       maximum coordinate value
 *  \param lrf   on return, if one side is empty,
 *            \c lrf will be leftFlag or rightFlag
 *            to indicate which side is empty.
 *  \param imbalance on return, if one half is empty, the
 *                     imbalance will be set to 0.0 (perfect balance)
 *  \param cutValue on return, if one half is empty, the
 *                    cutValue will be set to minCoord or
 *                    maxCoord depending on which half is empty.
 *
 *  \return true if all coords were moved to one half because
 *            the other half is empty, false otherwise.
 */

template <typename scalar_t>
  bool emptyPartsCheck( const RCP<const Environment> &env,
   const ArrayView<double> fractionLeft, 
   scalar_t minCoord, scalar_t maxCoord,
   leftRightFlag &lrf,
   scalar_t &cutValue)
{
  env->debug(DETAILED_STATUS, string("Entering emptyPartsCheck"));
  // initialize return values
  lrf = leftFlag;
  cutValue = 0.0;

  size_t weightDim = fractionLeft.size();
  int numEmptyRight = 0, numEmptyLeft = 0;

  for (size_t i=0; i < weightDim; i++){
    if (fractionLeft[i] == 0.0)
      numEmptyLeft++;
    else if (fractionLeft[i] == 1.0)
      numEmptyRight++;
  }

  if (numEmptyRight + numEmptyLeft == 0)
    return false;

  env->localInputAssertion(__FILE__, __LINE__,
    "partitioning not solvable - conflicting part size criteria",
    (numEmptyRight==0) || (numEmptyLeft==0), BASIC_ASSERTION);

  if (numEmptyLeft > 0){
    lrf = leftFlag;
    cutValue = minCoord;
  }
  else{
    lrf = rightFlag;
    cutValue = maxCoord;
  }

  env->debug(DETAILED_STATUS, string("Exiting emptyPartsCheck"));
  return true;
}

/*! \brief Move boundary coordinates to the right.
 *
 *  \param env   the environment
 *  \param comm   the communicator
 *  \param totalWeightLeft  the total weight in the left region
 *             at the end of the last iteration.
 *  \param targetWeightLeft  the ideal weight for the left part.
 *  \param cutLocation   index into sums of boundary sum.
 *                         each boundary.
 *  \todo document parameters
 */

template <typename lno_t, typename gno_t, typename scalar_t>
  void testCoordinatesOnRightBoundary(
    const RCP<const Environment> &env,
    const RCP<Comm<int> > &comm, 
    scalar_t totalWeightLeft,
    scalar_t targetLeftScalar,
    int cutLocation,
    ArrayView<scalar_t> localSums,
    ArrayView<scalar_t> globalSums,
    ArrayView<lno_t> index,
    ArrayView<StridedData<lno_t, scalar_t> > weights,
    multiCriteriaNorm mcnorm,
    ArrayView<unsigned char> lrFlags,
    scalar_t &globalWeightMovedRight)   // output
{
  env->debug(DETAILED_STATUS, string("Entering testCoordinatesOnRightBoundary"));
  int nprocs = comm->getSize();
  int rank = comm->getRank();

  globalWeightMovedRight = 0.0;

  scalar_t localBoundarySum = localSums[cutLocation];

  scalar_t total = totalWeightLeft;
  for (int i=0; i <= cutLocation; i++)
    total += globalSums[i];

  scalar_t totalMoveRight = total - targetLeftScalar;
  scalar_t localMoveRight = localBoundarySum;
  scalar_t actualWeightMovedRight = 0.0;
  
  Array<scalar_t> scansum(nprocs+1, 0.0);
  Teuchos::gatherAll<int, scalar_t>(*comm, 1, 
    &localBoundarySum, nprocs, scansum.getRawPtr()+1);
  for (int i=2; i<=nprocs; i++)
    scansum[i] += scansum[i-1];

  if (localBoundarySum > 0.0){

    scalar_t sumMine = scansum[rank]; // sum of ranks preceding me
    scalar_t diff = scansum[nprocs] - sumMine;
    localMoveRight = 0;

    if (diff <= totalMoveRight)
      localMoveRight = localBoundarySum;  // all
    else{
      scalar_t leftPart = diff - totalMoveRight;
      if (leftPart < localBoundarySum)
        localMoveRight = localBoundarySum - leftPart;
    }
  }

  if (localMoveRight > 0.0){

    bool moveAll =  (localMoveRight >= localBoundarySum);

    if (moveAll){
      actualWeightMovedRight = localBoundarySum;
      for (lno_t i=0; i < lrFlags.size(); i++){
        if (lrFlags[i] == cutLocation){
          lrFlags[i] = rightFlag;
        }
      }
    }
    else{
      int weightDim = weights.size();
      int useIndices = index.size() > 0;

      for (lno_t i=0; i < lrFlags.size(); i++){
        if (lrFlags[i] == cutLocation){
          lrFlags[i] = rightFlag;

          lno_t idx = (useIndices ? index[i] : i);

          actualWeightMovedRight += getCoordWeight<lno_t, scalar_t>(
            idx, mcnorm, weights.view(0, weightDim));

          if (actualWeightMovedRight >= localMoveRight)
            break;
        }
      } // next coordinate
    }
  }

  try{
    reduceAll<int, scalar_t>(
      *comm, Teuchos::REDUCE_SUM, 1, &actualWeightMovedRight,
      &globalWeightMovedRight);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->debug(DETAILED_STATUS, string("Exiting testCoordinatesOnRightBoundary"));
  return;
}

/*! \brief Find the point in space that divides the data evenly with
 *     respect to the weights, part sizes, and the user's objective.
 *
 *   \param env the Environment for the application.
 *   \param comm the communicator for this step.
 *   \param params a bit map of boolean parameters.
 *   \param numTestCuts the number of test cuts to make in one round.
 *   \param tolerance the maximum acceptable imbalance (0,1).
 *   \param cutDim  the dimension of the coordinates to cut.
 *   \param coordDim the first \c coordDim vectors in the \c vectors
 *              list are coordinates, the rest are weights.
 *   \param vectors lists of coordinates and non-uniform weights
 *   \param index is the index into the \c vectors arrays for the
 *              coordinates to be included in the partitioning.
 *              If <tt>index.size()</tt> is zero, then all coordinates
 *              are included.
 *   \param fractionLeft  the size of the left part for each weight,
 *                  the right part should measure <tt>1.0 - fractionLeft</tt>.
 *   \param uniformWeights element \c w is true if weights for weight
 *                 dimension \c w are all 1.0.
 *   \param coordGlobalMin the global minimum of coordinates in dimension
 *                                \c cutDim
 *   \param coordGlobalMax the global maximum of coordinates in dimension
 *                                \c cutDim
 *   \param cutValue  on return this is the computed cut location.
 *   \param lrflags on return lists the values leftFlag or rightFlag to
 *        indicate whether the corresponding coordinate is on the left of
 *        on the right of the cut.  Allocated by caller.  In particular,
 *        if <tt>index.size()</tt> is non-zero, then
 *        <tt>lrflags[i]</tt> is the flag for 
 *        <tt>coordinate vectors[index[i]]</tt>.  Otherwise it is the
 *        flag for <tt>coordinate vectors[i]</tt>.
 *        
 *   \param totalWeightLeft on return is the global total weight 
 *                  left of the cut.
 *   \param totalWeightRight on return is the global total weight 
 *                  right of the cut.
 *   \param localCountLeft on return is the local number of objects
 *                  left of the cut.
 *   \param imbalance on return is the imbalance for the 
 *              computed partitioning (0, 1).
 * 
 *  \todo a separate simpler function when weightDim <= 1
 *  \todo During the first global comm, ensure all procs have same values
 *           for cutDim and fractionLeft.
 */

template <typename mvector_t>
  void BSPfindCut(
    const RCP<const Environment> &env,
    const RCP<Comm<int> > &comm,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts,
    typename mvector_t::scalar_type tolerance,
    int cutDim,
    int coordDim,
    const RCP<mvector_t> &vectors,
    ArrayView<typename mvector_t::local_ordinal_type> index,
    ArrayView<double> fractionLeft,
    ArrayView<bool> uniformWeights,
    typename mvector_t::scalar_type coordGlobalMin,
    typename mvector_t::scalar_type coordGlobalMax,
    typename mvector_t::scalar_type &cutValue,         // output
    ArrayView<unsigned char> lrFlags,                  // output
    typename mvector_t::scalar_type &totalWeightLeft,  // output
    typename mvector_t::scalar_type &totalWeightRight, // output
    typename mvector_t::local_ordinal_type &localCountLeft, // output
    typename mvector_t::scalar_type &imbalance)        // output
{
  env->debug(DETAILED_STATUS, string("Entering BSPfindCut"));

  // initialize output
  bool useIndices = index.size() > 0;
  int numAllCoords = vectors->getLocalLength();
  int numCoords = (useIndices ? index.size() : numAllCoords);
  cutValue = totalWeightLeft = totalWeightRight = imbalance = 0.0;
  localCountLeft = 0;

  for (int i=0; i < numCoords; i++)
    lrFlags[i] = unsetFlag;

  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  bool multiplePartSizeSpecs = params.test(rcb_multiplePartSizeSpecs);
  bool rectilinearBlocks = params.test(rcb_rectilinearBlocks);
  bool averageCuts = params.test(rcb_averageCuts);

  // A coordinate is considered to be on a cut if it is within
  // this distance of the cut.

  double epsilon = (coordGlobalMax - coordGlobalMin) * 10e-9;

  // Find the coordinate values and weights.

  int weightDim = uniformWeights.size();   // at least one
  int numNonUniformWeights = 0;
  for (int i=0; i < weightDim; i++){
    if (!uniformWeights[i])
      numNonUniformWeights++;
  }

  if (env->getDebugLevel() >= DETAILED_STATUS){
    ostringstream info;
    info << "Weight dim " << weightDim << ", Fraction left:";
    for (int i=0; i < weightDim; i++)
      info << " " << fractionLeft[i];
    info << endl << "Dimension " << cutDim << " [";
    info << coordGlobalMin << ", " << coordGlobalMax << "]";
    info << endl << "# test cuts " << numTestCuts;
    info << ", tolerance " << tolerance << endl;
    env->debug(DETAILED_STATUS, info.str());
  }

  const scalar_t *coordValue = vectors->getData(cutDim).getRawPtr();

  // An empty input_t object implies uniform weights.

  input_t *wgtinfo = new input_t [weightDim];
  env->localMemoryAssertion(__FILE__, __LINE__, weightDim, wgtinfo);
  ArrayRCP<input_t> weight(wgtinfo, 0, weightDim, true);

  if (numNonUniformWeights > 0){
    for (int wdim = 0, widx=coordDim; wdim < weightDim; wdim++){
      if (!uniformWeights[wdim]){
        weight[wdim] = input_t(vectors->getData(widx++), 1);
      }
    }
  }

  // Multicriteria norm

  multiCriteriaNorm mcnorm = normBalanceTotalMaximum;

  if (params.test(rcb_minMaximumWeight))
    mcnorm = normMinimizeMaximumWeight;
  else if (params.test(rcb_minTotalWeight))
    mcnorm = normMinimizeTotalWeight;
  
  // Goal is globally find one cut that comes close to leaving
  // partSizeLeft*totalWeight on the left side.

  Epetra_SerialDenseVector partSizeLeft( 
    View, fractionLeft.getRawPtr(), weightDim);
  
  // Where do we make the first test cuts?
  //
  //   min     1     2     3     max
  //
  // Space is [min, max].  
  // If there are three test cuts: 1, 2 and 3.
  //   4 regions: [min,1] [1,2] [2,3] [3,max].
  //   5 boundaries: min, 1, 2, 3, and max.

  int numRegions = numTestCuts + 1;
  int numBoundaries = numTestCuts + 2;
  int endBoundary = numBoundaries - 1;
  vector<scalar_t> boundaries(numBoundaries);
  vector<scalar_t> searchBoundaries(numBoundaries);

  int numSums = numBoundaries+numRegions;

  Teuchos::Zoltan2_RCBOperation<int, scalar_t> reductionOp(
    numSums,      // number of sums
    numRegions,   // number of mins
    numRegions);  // number of maxes

  typename std::vector<scalar_t>::iterator foundCut;

  bool done=false;
  bool fail=false;
  scalar_t min = coordGlobalMin;
  scalar_t max = coordGlobalMax;
  lno_t numRemaining = numCoords;
  lno_t prevNumRemaining = numCoords;
  size_t numGlobalPoints = vectors->getGlobalLength();
  size_t sanityCheck = numGlobalPoints;

  double totalWeight = 0;
  double targetLeftScalar = 0;
  double targetLeftNorm = 0;
  Epetra_SerialDenseVector targetLeftVector(weightDim);

  while (!done && !fail && sanityCheck--){

    // Create regions into which coordinates will be placed.
    scalar_t diff = (max - min) / numRegions;
    boundaries[0] = min;
    boundaries[endBoundary] = max;
    for (int i=0; i < endBoundary; i++){
      searchBoundaries[i+1] = boundaries[i+1] = boundaries[i] + diff;
    }

    // Move ends slightly so we catch points on boundary.
    searchBoundaries[0] = min - epsilon;
    searchBoundaries[endBoundary] = max + epsilon;

    // Save region and boundary sums, and region min and max.
 
    ArrayRCP<scalar_t> localSums = reductionOp.getInitializedBuffer(
      searchBoundaries[endBoundary], searchBoundaries[0]);
 
    ArrayRCP<scalar_t> globalSums = reductionOp.getInitializedBuffer(
      searchBoundaries[endBoundary], searchBoundaries[0]);

    scalar_t *sums = localSums.getRawPtr();
    scalar_t *regionMin = sums + numSums;
    scalar_t *regionMax = regionMin + numRegions;

    if (numRemaining > 0){

      // Assign each of my points to a region.
      // lower_bound() finds the first cut f, such that f >= coordValue[i].
      // So for now, objects that are on the cut boundary go into the
      // region on the "left" side.

      for (lno_t i=0; i < numCoords; i++){
  
        if (lrFlags[i] != unsetFlag)
          continue;

        int inRegion = 0;
        int idx = (useIndices ? index[i] : i);
        scalar_t value = coordValue[idx];

        if (numRegions > 2){
       
          foundCut = std::lower_bound(
            searchBoundaries.begin(), searchBoundaries.end(), value);
          
          env->localBugAssertion(__FILE__, __LINE__, "search cuts", 
            foundCut != searchBoundaries.end(), BASIC_ASSERTION);

          inRegion = foundCut - searchBoundaries.begin() - 1;        
        }
        else{
          if (value <= boundaries[1])
            inRegion = 0;
          else
            inRegion = 1;
        }

        int sumIdx = 1 + inRegion*2;

        if (value >= boundaries[inRegion+1]-epsilon){  
          // "on" right boundary of this region
          sumIdx++;
        }
        else if (inRegion==0 && (value < min + epsilon)){
          // in region 0 but far left boundary
          sumIdx--;
        }

        lrFlags[i] = (unsigned char)sumIdx;

        if (value < regionMin[inRegion])
          regionMin[inRegion] = value;
        if (value > regionMax[inRegion])
          regionMax[inRegion] = value;

        sums[sumIdx] += getCoordWeight<lno_t, scalar_t>(idx, 
              mcnorm, weight.view(0, weightDim));

      } // next coord
    }

    try{
      reduceAll<int, scalar_t>(*comm, reductionOp,
        localSums.size(), localSums.getRawPtr(), globalSums.getRawPtr());
    }
    Z2_THROW_OUTSIDE_ERROR(*env)

    sums = globalSums.getRawPtr();
    regionMin = sums + numSums;
    regionMax = regionMin + numRegions;

    if (env->getDebugLevel() >= DETAILED_STATUS){
      ostringstream info;
      info << "  Region " << min << " - " << max << endl;
      info << "  Remaining to classify: " << numRemaining << endl;
      info << "  Boundaries: ";
      for (int i=0; i < numBoundaries; i++)
        info << boundaries[i] << " ";
      info << endl << "  For searching: ";
      for (int i=0; i < numBoundaries; i++)
        info << searchBoundaries[i] << " ";
      info << endl << "  Global sums: ";
      double sumTotal=0;
      for (int i=0; i < numSums; i++){
        sumTotal += sums[i];
        info << sums[i] << " ";
      }
      info << " total: " << sumTotal << endl;
      env->debug(DETAILED_STATUS, info.str());
    }

    if (totalWeight == 0){   // first time through only

      for (int i=0; i < numSums; i++)
        totalWeight += sums[i];

      partSizeLeft.Scale(totalWeight);
      targetLeftVector = partSizeLeft;

      targetLeftScalar = targetLeftVector[0];
      targetLeftNorm = targetLeftVector.Norm2();
      totalWeightLeft = 0;
      totalWeightRight = 0;
    }

    int cutLocation=0;
    scalar_t testDiff=0, prevTestDiff=0, target=0;

    if (multiplePartSizeSpecs){
      // more complex: if we have multiple weight dimensions, the
      //   weights are non-uniform, and the part sizes requested
      //   for each each weight dimension differ, then we may not
      //   be able to reach the imbalance tolerance.
      //
      // TODO: discuss how to (or whether to) handle this case.
      // 
      // For now we look at this imbalance:
      // 
      //    |target - actual|^2 / |target|^2
  
      target = targetLeftNorm;
      Epetra_SerialDenseVector testVec(weightDim);
      for (int i=0; i < weightDim; i++)
        testVec[i] = totalWeightLeft;
      Epetra_SerialDenseVector diffVec = testVec;
      diffVec.Scale(-1.0);
      diffVec += targetLeftVector;

      testDiff = diffVec.Norm2(); // imbalance numerator
      prevTestDiff = testDiff;
      cutLocation= -1;

      while (++cutLocation< numSums){

        for (int i=0; i < weightDim; i++)
          testVec[i] += sums[cutLocation];
  
        diffVec = testVec;
        diffVec.Scale(-1.0);
        diffVec += targetLeftVector;
  
        testDiff = diffVec.Norm2();
        
        if (testDiff >= target)
          break;

        prevTestDiff = testDiff;
      }
    }
    else{    // the part sizes for each weight dimension are the same
  
      target = targetLeftScalar;
      testDiff = totalWeightLeft; 
      prevTestDiff = testDiff;
      cutLocation = -1;
  
      while (++cutLocation < numSums){
     
        testDiff += sums[cutLocation];
        if (testDiff >= target)
          break;
  
        prevTestDiff = testDiff;
      }
    }

    scalar_t diffLeftCut = target - prevTestDiff;
    scalar_t diffRightCut = testDiff - target;

    if (diffLeftCut < diffRightCut){
      imbalance = diffLeftCut / target;
      if (imbalance <= tolerance){
        env->debug(DETAILED_STATUS, "  Done, tolerance met");
        done = true;
        cutLocation--;
      }
    } 
    else{
      imbalance = diffRightCut / target;
      if (imbalance <= tolerance){
        env->debug(DETAILED_STATUS, "  Done, tolerance met");
        done = true;
     }
    }

    bool cutLocIsRegion = (cutLocation % 2 == 1);
    bool cutLocIsBoundary = !cutLocIsRegion;

    if (env->getDebugLevel() >= DETAILED_STATUS){
      ostringstream info;
      info << "  Best cut location: " << cutLocation;
      if (cutLocIsRegion) info << " just after a region." << endl;
      else info << " just after a boundary." << endl;
      env->debug(DETAILED_STATUS, info.str());
    }

    if (!done && cutLocIsBoundary){

      done = true;    // can not divide space any more

      env->debug(DETAILED_STATUS, "  Done, cutting at a region boundary");

      if (rectilinearBlocks){
        // Can not divide boundary points into two
        // different regions to achieve balance.
        fail = true;
      }
      else {
        // Maybe by moving some of the points right, we can
        // obtain a better balance.  If so, the lrFlag for
        // the points moved right will be updated here.

        scalar_t globalWeightMovedRight(0);

        testCoordinatesOnRightBoundary<lno_t, gno_t, scalar_t>(
          env, comm, totalWeightLeft, targetLeftScalar, cutLocation, 
          localSums.view(0, numSums), globalSums.view(0, numSums),
          index, weight.view(0, weightDim), mcnorm, lrFlags,
          globalWeightMovedRight);

        scalar_t newSum = testDiff - globalWeightMovedRight;
        globalSums[cutLocation] -= globalWeightMovedRight;

        if (newSum > target)
          imbalance = (target - newSum) / target;
        else
          imbalance = (newSum - target) / target;

        if (imbalance > tolerance)
          fail = true;
      }
    }

    int rightmostLeftNum=0, leftmostRightNum=0;

    if (!done){
      // Best cut is following a region.  Narrow down the boundaries.

      int regionNumber = (cutLocation - 1) / 2;

      min = regionMin[regionNumber];
      max = regionMax[regionNumber];

      rightmostLeftNum = cutLocation - 1;
      leftmostRightNum = cutLocation + 1;
    }
    else{
      rightmostLeftNum = cutLocation;
      leftmostRightNum = cutLocation + 1;

      if (cutLocIsRegion && averageCuts){
        int regionNumber = (cutLocation - 1) / 2;
        cutValue = (
          boundaries[regionNumber+1] + // boundary to right
          regionMax[regionNumber])     // rightmost point in region
          / 2.0;
      }
    }

    for (int i=0; i <= rightmostLeftNum; i++){
      totalWeightLeft += globalSums[i];
    }

    for (lno_t i=0; i < numCoords; i++){
      if (lrFlags[i] != leftFlag && lrFlags[i] != rightFlag){
        if (lrFlags[i] <= rightmostLeftNum){
          lrFlags[i] = leftFlag;
          numRemaining--;
          localCountLeft++;
        }
        else if (lrFlags[i] >= leftmostRightNum){
          lrFlags[i] = rightFlag;
          numRemaining--;
        }
        else{
          lrFlags[i] = unsetFlag;   // still to be determined
        }
      }
    }

    if (env->getDebugLevel() >= VERBOSE_DETAILED_STATUS && numCoords < 100){
      // For large numCoords, building this message
      // takes an extraordinarily long time.
      ostringstream ossLeft;
      ostringstream ossRight;
      ossLeft << "left: ";
      ossRight << "right: ";
      for (lno_t i=0; i < numCoords; i++){
        if (lrFlags[i] == unsetFlag)
          continue;
        lno_t idx = (useIndices ? index[i] : i);
        scalar_t val = coordValue[idx];
        if (lrFlags[i] == leftFlag)
          ossLeft << val << " ";
        else if (lrFlags[i] == rightFlag)
          ossRight << val << " ";
        else 
          env->localBugAssertion(__FILE__, __LINE__, 
            "left/right flags", false, BASIC_ASSERTION);
      }
      ostringstream msg;
      msg << ossLeft.str() << endl << ossRight.str() << endl;
      env->debug(VERBOSE_DETAILED_STATUS, msg.str());
    }

    prevNumRemaining = numRemaining;
  }  // while !done

  totalWeightRight = totalWeight - totalWeightLeft;

  if (fail)
    env->debug(BASIC_STATUS, "Warning: tolerance not achieved in sub-step");

  // TODO: If fail the warn that tolerance was not met.

  env->globalInputAssertion(__FILE__, __LINE__, "partitioning not solvable",
    done, DEBUG_MODE_ASSERTION, comm);

  env->memory("End of bisection");

  if (env->getDebugLevel() >= DETAILED_STATUS){
    ostringstream oss;
    oss << "Exiting BSPfindCut, ";
    oss << "# iterations: " << numGlobalPoints - sanityCheck;    
    env->debug(DETAILED_STATUS, oss.str());
  }

  return;
}

/*! \brief Divide the coordinates into a "left" half and "right" half.
 *
 *   \param env the Environment for the application.
 *   \param comm the communicator for this step.
 *   \param params a bit map of boolean parameters.
 *   \param numTestCuts the number of test cuts to make in one round.
 *   \param tolerance the maximum acceptable imbalance (0,1).
 *   \param coordDim the first \c coordDim vectors in the \c vectors
 *              list are coordinates, the rest are weights.
 *   \param vectors lists of coordinates and non-uniform weights
 *   \param uniformWeights element \c w is true if weights for weight
 *                 dimension \c w are all 1.0.
 *   \param solution for obtaining the part sizes
 *   \param part0  is the first part of the parts to bisected.
 *   \param part1  is the last part of the parts to bisected.
 *   \param lrflags on return has the value leftFlag or rightFlag to
 *        indicate whether the corresponding coordinate is on the left of
 *        on the right of the cut.  Allocated by caller.
 *   \param cutDimension on return coordinate dimension that was cut.
 *   \param cutValue  on return this is the computed cut location.
 *   \param imbalance on return is the imbalance for the computed partitioning
 *           for cutDim and fractionLeft (0,1).
 *   \param numPartsLeftHalf on return the number of parts in the left half.
 *   \param weightLeftHalf on return the total weight in the left half.
 *   \param weightRightHalf on return the total weight in the right half.
 */

template <typename mvector_t, typename Adapter>
  void determineCut(
    const RCP<const Environment> &env,
    const RCP<Comm<int> > &comm,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts,
    typename mvector_t::scalar_type tolerance,
    int coordDim, 
    const RCP<mvector_t> &vectors,
    const ArrayView<bool> uniformWeights,
    multiCriteriaNorm mcnorm,
    const RCP<PartitioningSolution<Adapter> > &solution,
    partId_t part0, 
    partId_t part1,
    ArrayView<unsigned char> lrflags,  // output
    int &cutDimension,                  // output
    typename mvector_t::scalar_type &cutValue,   // output
    typename mvector_t::scalar_type &imbalance,  // output
    partId_t &numPartsLeftHalf,                  // output
    typename mvector_t::scalar_type &weightLeftHalf,  // output
    typename mvector_t::scalar_type &weightRightHalf  // output
    )
{
  env->debug(DETAILED_STATUS, string("Entering determineCut"));
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  lno_t numLocalCoords = vectors->getLocalLength();

  // initialize return values
  cutDimension = 0;
  cutValue = imbalance = weightLeftHalf = weightRightHalf = 0.0;
  numPartsLeftHalf = 0;

  ///////////////////////////////////////////////////////
  // Pick a cut direction.

  scalar_t globalMinCoord, globalMaxCoord;
  ArrayView<lno_t> emptyIndex;

  getCutDimension<mvector_t>(env, comm, coordDim, vectors, emptyIndex,
    cutDimension, globalMinCoord, globalMaxCoord);

  ///////////////////////////////////////////////////////
  // Compute part sizes for the two parts.

  ArrayRCP<double> fractionLeft;
  size_t weightDim = uniformWeights.size();

  getFractionLeft<Adapter>(env, part0, part1, solution,
    fractionLeft, numPartsLeftHalf);

  // Special case of empty left or empty right.

  leftRightFlag lrf;

  bool emptyParts = emptyPartsCheck(env,
    fractionLeft.view(0, weightDim), // input
    globalMinCoord, globalMaxCoord,  // input
    lrf, cutValue);                  // output

  if (emptyParts){
    
    for (lno_t i=0; i < lrflags.size(); i++)
      lrflags[i] = lrf;

    imbalance = 0.0;                // perfect balance
    scalar_t totalWeight = 0.0;
    int numNonUniform = 0;

    for (size_t i=0; i < weightDim; i++)
      if (!uniformWeights[i])
        numNonUniform++;

    int wgt1 = vectors->getNumVectors() - numNonUniform;

    if (weightDim == 1){
      if (numNonUniform == 0)
        totalWeight = numLocalCoords;
      else{
        const scalar_t *val = vectors->getData(wgt1).getRawPtr();
        for (lno_t i=0; i < numLocalCoords; i++)
          totalWeight += val[i];
      }
    }
    else{  // need to add up total normed weight
      Array<input_t> wgts(weightDim);
      for (size_t i=0; i < weightDim; i++){
        if (!uniformWeights[i]){
          wgts[i] = input_t(vectors->getData(wgt1++), 1);
        }
      }

      partId_t numParts, numNonemptyParts;
      ArrayRCP<MetricValues<scalar_t> > metrics;
      ArrayRCP<scalar_t> weightSums;
    
      globalSumsByPart<scalar_t, unsigned char, lno_t>(
        env, comm, lrflags, 
        wgts.view(0, weightDim), mcnorm,
        numParts, numNonemptyParts, metrics, weightSums);

      totalWeight = weightSums[1];  // sum normed weight
    }

    if (lrf == leftFlag){
      numPartsLeftHalf = part1 - part0 + 1;
      weightLeftHalf = totalWeight;
      weightRightHalf = 0;
    }
    else{
      numPartsLeftHalf = 0;
      weightRightHalf = totalWeight;
      weightLeftHalf = 0;
    }

    env->debug(DETAILED_STATUS, string("Exiting determineCut"));
    return;
  }

  ///////////////////////////////////////////////////////
  // Divide the coordinates into balanced left and right
  // halves.

  ArrayView<lno_t> emptyIndices;
  lno_t localCountLeft;

  try{
    BSPfindCut<mvector_t>( env, comm,
      params, numTestCuts, tolerance,
      cutDimension, coordDim, vectors, emptyIndices,
      fractionLeft.view(0, weightDim), uniformWeights.view(0, weightDim),
      globalMinCoord, globalMaxCoord,
      cutValue, lrflags.view(0, numLocalCoords),
      weightLeftHalf, weightRightHalf, localCountLeft, imbalance);
  }
  Z2_FORWARD_EXCEPTIONS
  env->debug(DETAILED_STATUS, string("Exiting determineCut"));
}


/*! \brief Perform RCB on the local process only.
 *
 *   \param env the Environment for the application.
 *   \param depth  depth of recursion, for debugging, 
 *           call with "1" first time.
 *   \param params a bit map of boolean parameters.
 *   \param numTestCuts the number of test cuts to make in one round.
 *   \param tolerance the maximum acceptable imbalance (0,1).
 *   \param coordDim the first \c coordDim vectors in the \c vectors
 *              list are coordinates, the rest are weights.
 *   \param vectors lists of coordinates and non-uniform weights
 *   \param index is the index into the \c vectors arrays for the
 *              coordinates to be included in the partitioning.
 *       If index.size() is zero then indexing will not be used.
 *   \param uniformWeights element \c w is true if weights for weight
 *                 dimension \c w are all 1.0.
 *   \param solution for obtaining part sizes.
 *   \param part0  is the first part of the parts to bisected.
 *   \param part1  is the last part of the parts to bisected.
 *   \param partNum on return <tt>partNum[i]</tt> is the new
 *                part number for coordinate \c i.
 */

template <typename mvector_t, typename Adapter>
 void serialRCB(
    const RCP<const Environment> &env,
    int depth,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, 
    typename mvector_t::scalar_type tolerance, 
    int coordDim,
    const RCP<mvector_t> &vectors, 
    ArrayView<typename mvector_t::local_ordinal_type> index,
    const ArrayView<bool> uniformWeights,
    const RCP<PartitioningSolution<Adapter> > &solution,
    partId_t part0, 
    partId_t part1,
    ArrayView<partId_t> partNum)   // output
{
  env->debug(DETAILED_STATUS, string("Entering serialRCB"));
  env->timerStart(MICRO_TIMERS, "serialRCB", depth, 2);

  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;

  RCP<Comm<int> > comm(new Teuchos::SerialComm<int>);  

  lno_t numLocalCoords=0;
  bool useIndices;
  bool firstCall=false;

  if (index.size() == 0){
    // First time through there are no indices.
    useIndices = false;
    numLocalCoords = vectors->getLocalLength();
    firstCall = true;
  }
  else{
    useIndices = true;
    numLocalCoords = index.size();
  }

  if (env->getDebugLevel() >= DETAILED_STATUS){
    ostringstream info;
    info << "  Number of coordinates: " << numLocalCoords << endl;
    info << "  Use index: " << useIndices << endl;
    env->debug(DETAILED_STATUS, info.str());
  }

  ///////////////////////////////////////////////////////
  // Are we done?

  if (part1 == part0){
    if (useIndices)
      for (lno_t i=0; i < numLocalCoords; i++)
        partNum[index[i]] = part0;
    else
      for (lno_t i=0; i < numLocalCoords; i++)
        partNum[i] = part0;

    env->memory("serial RCB end");
    env->timerStop(MICRO_TIMERS, "serialRCB", depth, 2);
    env->debug(DETAILED_STATUS, string("Exiting serialRCB"));

    return;
  }

  ///////////////////////////////////////////////////////
  // Pick a cut direction

  int cutDimension;
  scalar_t minCoord, maxCoord;

  try{
    getCutDimension(env, comm, coordDim, vectors, index,
      cutDimension, minCoord, maxCoord);
  }
  Z2_FORWARD_EXCEPTIONS

  ///////////////////////////////////////////////////////
  // Compute relative sizes of the two halves.

  ArrayRCP<double> fractionLeft;
  partId_t numPartsLeftHalf;

  try{
    getFractionLeft<Adapter>(env, part0, part1, solution,
      fractionLeft, numPartsLeftHalf);
  }
  Z2_FORWARD_EXCEPTIONS

  // Check for special case of empty left or empty right.

  int weightDim = uniformWeights.size();
  scalar_t imbalance, cutValue;  //unused for now
  leftRightFlag lrf;

  bool emptyPart = emptyPartsCheck(env, fractionLeft.view(0, weightDim), 
    minCoord, maxCoord, lrf, cutValue);

  if (emptyPart){  // continue only on the side that is not empty

    if (lrf == rightFlag)  // right is empty
      part1 = part0 + numPartsLeftHalf - 1;
    else
      part0 = part0 + numPartsLeftHalf;

    if (part0 == part1){
      if (useIndices){
        for (lno_t i=0; i < numLocalCoords; i++){
          partNum[index[i]] = part0;
        }
      }
      else{
        for (lno_t i=0; i < numLocalCoords; i++){
          partNum[i] = part0;
        }
      }
      
      imbalance = 0.0;       // perfect
      env->timerStop(MICRO_TIMERS, "serialRCB", depth, 2);
      env->debug(DETAILED_STATUS, string("Exiting serialRCB"));
      return;
    }
  }

  ///////////////////////////////////////////////////////
  // Divide into balanced left and right halves.

  scalar_t totalLeft=0, totalRight=0;
  lno_t localCountLeft=0, localCountRight=0;
  unsigned char *flags = new unsigned char [numLocalCoords];
  env->localMemoryAssertion(__FILE__, __LINE__, numLocalCoords, flags) ;
  ArrayRCP<unsigned char> lrflags(flags, 0, numLocalCoords, true);

  try{
    BSPfindCut<mvector_t>( env, comm,
      params, numTestCuts, tolerance,
      cutDimension, coordDim, vectors, index,
      fractionLeft.view(0, weightDim), uniformWeights.view(0, weightDim),
      minCoord, maxCoord,
      cutValue, lrflags.view(0, numLocalCoords),
      totalLeft, totalRight, localCountLeft, imbalance);
  }
  Z2_FORWARD_EXCEPTIONS

  if (firstCall)
    env->memory("serial RCB start");

  ///////////////////////////////////////////////////////
  // Adjust indices for left half and right half

 localCountRight = numLocalCoords - localCountLeft;

  if (localCountLeft){

    lno_t *newIndex = new lno_t [localCountLeft];
    env->localMemoryAssertion(__FILE__, __LINE__, localCountLeft, newIndex);
    ArrayView<lno_t> leftIndices(newIndex, localCountLeft);

    if (useIndices){
      for (lno_t i=0, ii=0; i < numLocalCoords; i++)
        if (lrflags[i] == leftFlag)
          newIndex[ii++] = index[i];
    }
    else{
      for (lno_t i=0, ii=0; i < numLocalCoords; i++)
        if (lrflags[i] == leftFlag)
          newIndex[ii++] = i;
    }

    partId_t newPart1 = part0 + numPartsLeftHalf - 1;

    env->timerStop(MICRO_TIMERS, "serialRCB", depth, 2);

    serialRCB<mvector_t, Adapter>(env, depth+1, params, numTestCuts, tolerance, 
      coordDim, vectors, leftIndices,
      uniformWeights.view(0, weightDim), solution,
      part0, newPart1, partNum);

    env->timerStart(MICRO_TIMERS, "serialRCB", depth, 2);

    delete [] newIndex;
  }


  if (localCountRight){

    lno_t *newIndex = new lno_t [localCountRight];
    env->localMemoryAssertion(__FILE__, __LINE__, localCountRight, newIndex);
    ArrayView<lno_t> rightIndices(newIndex, localCountRight);

    if (useIndices){
      for (lno_t i=0, ii=0; i < numLocalCoords; i++)
        if (lrflags[i] == rightFlag){
          newIndex[ii++] = index[i];
        }
    }
    else{
      for (lno_t i=0, ii=0; i < numLocalCoords; i++)
        if (lrflags[i] == rightFlag)
          newIndex[ii++] = i;
    }

    partId_t newPart0 = part0 + numPartsLeftHalf;

    env->timerStop(MICRO_TIMERS, "serialRCB", depth, 2);

    serialRCB<mvector_t, Adapter>(env, depth+1, params, numTestCuts, tolerance, 
      coordDim, vectors, rightIndices,
      uniformWeights.view(0, weightDim), solution,
      newPart0, part1, partNum);

    env->timerStart(MICRO_TIMERS, "serialRCB", depth, 2);

    delete [] newIndex;
  }

  env->timerStop(MICRO_TIMERS, "serialRCB", depth, 2);
  env->debug(DETAILED_STATUS, string("Exiting serialRCB"));
}

}// namespace Zoltan2

#endif

