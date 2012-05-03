// @HEADER
//***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_AlgRCB_methods.hpp
    \brief Methods written for the RCB algorithm
*/

#ifndef _ZOLTAN2_ALGRCB_METHODS_HPP_
#define _ZOLTAN2_ALGRCB_METHODS_HPP_

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_XpetraTraits.hpp>

#include <Tpetra_Vector.hpp>

#include <sstream>
#include <string>
#include <cmath>
#include <bitset>

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
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;

  env->timerStart("getCutDimension");

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
  env->timerStop("getCutDimension");
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
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;

  env->timerStart("migrateData");

  int nprocs = comm->getSize();
  size_t nobj = vectors->getLocalLength();
  size_t nGlobalObj = vectors->getGlobalLength();

  env->localBugAssertion(__FILE__, __LINE__, "migrateData input", 
    nprocs>1 && lrflags.size()==nobj, DEBUG_MODE_ASSERTION);

  gno_t myNumLeft= 0, numLeft;
  for (lno_t i=0; i < nobj; i++)
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

  lno_t *sendCount = new lno_t [nprocs];
  env->localMemoryAssertion(__FILE__, __LINE__, nprocs, sendCount) ;
  memset(sendCount, 0, sizeof(int) * nprocs);
  ArrayView<lno_t> sendCountView(sendCount, nprocs);
  ArrayView<gno_t> sendBufView;

  if (nobj > 0){
    int *procId = new int [nobj];
    env->localMemoryAssertion(__FILE__, __LINE__, nobj, procId) ;
    int leftProc0 = 0;
    int rightProc0 = leftProc0 + leftNumProcs;
  
    int nextLeftProc = leftProc0;
    int nextRightProc = rightProc0;
    int *p = procId;
  
    for (lno_t i=0; i < nobj; i++){
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
  
    lno_t *sendOffset = new lno_t [nprocs];
    env->localMemoryAssertion(__FILE__, __LINE__, nprocs, sendOffset) ;
    sendOffset[0] = 0;
    for (int i=0; i < nprocs-1; i++)
      sendOffset[i+1] = sendOffset[i] + sendCount[i];

    gno_t *sendBuf = new gno_t [nobj];
    env->localMemoryAssertion(__FILE__, __LINE__, nobj, sendBuf) ;
    sendBufView = ArrayView<gno_t>(sendBuf, nobj);

    ArrayView<const gno_t> gnoList = vectors->getMap()->getNodeElementList();

    for (lno_t i=0; i < nobj; i++){
      int proc = procId[i];
      lno_t offset = sendOffset[proc]++;

      sendBuf[offset] = gnoList[i];
    }

    delete [] sendOffset;
    delete [] procId;
  }

  ArrayRCP<gno_t> recvBuf;
  ArrayRCP<lno_t> recvCount;

  try{
    AlltoAllv<gno_t, lno_t>(*comm, *env,
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

  lno_t numMyNewGnos = 0;
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

  env->timerStop("migrateData");
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
    for (int wdim=0; wdim < weightDim; wdim++){
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
  // initialize return values
  lrf = leftFlag;
  cutValue = 0.0;

  size_t weightDim = fractionLeft.size();
  int numEmptyRight = 0, numEmptyLeft = 0;

  for (int i=0; i < weightDim; i++){
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

  return true;
}

/*! \brief Move boundary coordinates to the right if necessary.
 *
 *  \param env   the environment
 *  \param comm   the communicator
 *  \param rectilinearBlocks should be true if points on a boundary
 *             can be split over two regions, false otherwise.
 *  \param totalWeightLeft  the total weight in the left region
 *             at the end of the last iteration.
 *  \param targetWeightLeft  the ideal weight for the left part.
 *  \param rightBoundary the index into the boundaries array of
 *                        bounary in question.
 *  \param boundaries the scalar value of each boundary.
 *  \param boundarySum the array of local sums of coordinates on
 *                         each boundary.
 *  \todo document parameters
 *
 *  All coordinates that are on a boundary have been placed in 
 *  the region to the left of the boundary.  The total weight
 *  on the left of this boundary exceeds the target weight for
 *  the left part by more than the tolerance. If the imbalance
 *  can be improved by moving some or all of them to the region on
 *  the right, then do so.
 */

template <typename lno_t, typename gno_t, typename scalar_t>
  void testCoordinatesOnRightBoundary(
    const RCP<const Environment> &env,
    const RCP<Comm<int> > &comm, 
    bool rectilinearBlocks,
    scalar_t totalWeightLeft,
    scalar_t targetLeftScalar,
    int rightBoundary,
    vector<scalar_t> &boundaries,
    ArrayView<scalar_t> boundarySum,
    ArrayView<StridedData<lno_t, scalar_t> > weights,
    multiCriteriaNorm mcnorm,
    ArrayView<const scalar_t> coords,
    scalar_t epsilon,
    ArrayView<unsigned char> lrFlags,
    ArrayView<lno_t> index,
    ArrayView<scalar_t> regionSums,
    scalar_t &globalWeightMovedRight,
    gno_t &localCountMovedRight)
{
  globalWeightMovedRight = 0.0;
  localCountMovedRight = 0;

  env->timerStart("testCoordinatesOnRightBoundary");

  scalar_t localBoundarySum = boundarySum[rightBoundary];
  scalar_t globalBoundarySum = 0;

  try{
    reduceAll<int, scalar_t>( *comm, Teuchos::REDUCE_SUM, 1,
      &localBoundarySum, &globalBoundarySum); 
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  if (globalBoundarySum == 0.0){  // no boundary coordinates to move
    env->timerStop("testCoordinatesOnRightBoundary");
    return;
  }

  int weightDim = weights.size();
  int useIndices = index.size() > 0;
  int numLocalCoords = (useIndices ? index.size() : coords.size());
  int currentRegion = rightBoundary - 1;

  // Check the total weight on the left if we move some
  // or all boundary coordinates to the right.

  scalar_t testLeft = totalWeightLeft;
  for (int i=0; i <= currentRegion; i++)
    testLeft += regionSums[i];

  scalar_t testSum = testLeft - globalBoundarySum;
  scalar_t totalMoveRight = 0.0;

  if (testSum < targetLeftScalar){
    // With coords in the left region, totalWeightLeft
    // exceeded the targetLeftScalar.  So balance may
    // be achievable if boundary coords are split.
    if (!rectilinearBlocks)
      totalMoveRight = targetLeftScalar - testSum;
  }
  else{
    // Move boundary coords to right to get to balance sooner.
    totalMoveRight = globalBoundarySum;
  }

  if (totalMoveRight == 0.0)
    env->timerStop("testCoordinatesOnRightBoundary");
    return;

  // We are moving some or all coordinates which are on the right
  // boundary of the new region over to the right.

  scalar_t localMoveRight = localBoundarySum;
  scalar_t actualWeightMovedRight = 0.0;

  if (totalMoveRight < globalBoundarySum){
 
    int nprocs = comm->getSize();
    int rank = comm->getRank();
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
  }

  if (localMoveRight > 0.0){

    // Check coordinates in the region to the left of right boundary.
    // If they are on the right boundary, move them to the right until 
    // the weight adds up to localMoveRight.

    bool moveAll =  (localMoveRight >= localBoundarySum);

    for (size_t i=0; i < numLocalCoords; i++){
      if (lrFlags[i] == currentRegion){
        lno_t idx = (useIndices ? index[i] : i);
        if (coords[idx] >= boundaries[rightBoundary]-epsilon){ 

          lrFlags[i] = rightFlag;
          localCountMovedRight++;

          scalar_t w = getCoordWeight<lno_t, scalar_t>(idx, 
            mcnorm, weights.view(0, weightDim));

          actualWeightMovedRight += w;

          if (!moveAll && (actualWeightMovedRight >= localMoveRight))
            break;
        }
      }
    } // next coordinate
  }

  try{
    reduceAll<int, scalar_t>(
      *comm, Teuchos::REDUCE_SUM, 1, &actualWeightMovedRight,
      &globalWeightMovedRight);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->timerStop("testCoordinatesOnRightBoundary");

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
  if (env->doStatus())
    env->debug(DETAILED_STATUS, string("Entering BSPfindCut"));

  int nprocs = comm->getSize();
  if (nprocs > 1)
    env->timerStart("BSPfindCut");
  else
    env->timerStart("BSPfindCut - serial");

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

  if (env->doStatus()){
    ostringstream info;
    info << "Weight dim " << weightDim << ", Fraction left:";
    for (int i=0; i < weightDim; i++)
      info << " " << fractionLeft[i];
    info << endl << "Dimension " << cutDim << " [";
    info << coordGlobalMin << ", " << coordGlobalMax << "]";
    info << endl << "# test cuts " << numTestCuts;
    info << ", tolerance " << tolerance << endl;
    env->debug(VERBOSE_DETAILED_STATUS, info.str());
  }

  const scalar_t *coordValue = vectors->getData(cutDim).getRawPtr();

  // An empty input_t object implies uniform weights.

  input_t *info = new input_t [weightDim];
  env->localMemoryAssertion(__FILE__, __LINE__, weightDim, info);
  ArrayRCP<input_t> weight(info, 0, weightDim, true);

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
  typename std::vector<scalar_t>::iterator foundCut;

  bool done=false;
  bool fail=false;
  scalar_t min = coordGlobalMin;
  scalar_t max = coordGlobalMax;
  lno_t numRemaining = numCoords;
  size_t sanityCheck = vectors->getGlobalLength();

  double totalWeight = 0;
  double targetLeftScalar = 0;
  double targetLeftNorm = 0;
  Epetra_SerialDenseVector targetLeftVector(weightDim);

  while (!done && !fail && sanityCheck--){

    // Create regions into which coordinates will be placed.

    scalar_t diff = (max - min) / numRegions;
    boundaries[0] = min;
    for (int i=0; i < numRegions; i++)
      boundaries[i+1] = boundaries[i] + diff;

    boundaries[endBoundary] += epsilon;

    Array<scalar_t> boundarySum(numBoundaries, 0.0);// total weight on boundary
    Array<scalar_t> regionMin(numRegions, max+1.0); // lowest value in region
    Array<scalar_t> regionMax(numRegions, min-1.0); // highest value in region

    if (numRemaining > 0){

      // Assign each of my points to a region.
      // lower_bound() finds the first cut f, such that f >= coordValue[i].
      // So for now, objects that are on the cut boundary go into the
      // region on the "left" side.

      for (size_t i=0; i < numCoords; i++){
  
        if (lrFlags[i] != unsetFlag)
          continue;

        int inRegion = 0;
        int idx = (useIndices ? index[i] : i);
        scalar_t value = coordValue[idx];
  
        if (numRegions > 2){
       
          foundCut = std::lower_bound(boundaries.begin(), boundaries.end(), 
              value);
          
          env->localBugAssertion(__FILE__, __LINE__, "search cuts", 
            foundCut != boundaries.end(), BASIC_ASSERTION);

          inRegion = foundCut - boundaries.begin() - 1;        
        }
        else{
          if (value <= boundaries[1])
            inRegion = 0;
          else
            inRegion = 1;
        }

        lrFlags[i] = (unsigned char)inRegion;

        if (averageCuts){
          if (value < regionMin[inRegion])
            regionMin[inRegion] = value;
          if (value > regionMax[inRegion])
            regionMax[inRegion] = value;
        }

        if (value >= boundaries[inRegion+1]-epsilon){  // "on" boundary
          scalar_t w = getCoordWeight<lno_t, scalar_t>(idx, 
            mcnorm, weight.view(0, weightDim));
          boundarySum[inRegion+1] += w;
        }
      }
    }

    partId_t numParts, numNonemptyParts;
    ArrayRCP<MetricValues<scalar_t> > metrics;
    ArrayRCP<scalar_t> weightSums;

    // Get the global sums in each region

    globalSumsByPart<scalar_t, unsigned char, lno_t>(
       env, comm,                 // environment
       lrFlags,                   // part assignments
       0, numRegions - 1,         // parts to be included in count
       weight.view(0, weightDim), // weights to include
       mcnorm,                    // multicriteria norm
       numParts,         // output: number of parts
       numNonemptyParts, // output: number that are non-empty
       metrics,          // output: MetricValues objects
       weightSums);      // output: lists of sums per part

    scalar_t *regionSums = NULL;
    if (numNonUniformWeights > 0)
      regionSums = weightSums.getRawPtr() + numParts; // normed weight sum
    else 
      regionSums = weightSums.getRawPtr();            // object count sums

    if (totalWeight == 0){   // first time through only

      for (int i=0; i < numParts; i++)
        totalWeight += regionSums[i];

      partSizeLeft.Scale(totalWeight);
      targetLeftVector = partSizeLeft;

      targetLeftScalar = targetLeftVector[0];
      targetLeftNorm = targetLeftVector.Norm2();
      totalWeightLeft = 0;
      totalWeightRight = 0;
    }

    int regionNum=0;
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

      scalar_t testDiff = diffVec.Norm2(); // imbalance numerator
      scalar_t prevTestDiff = testDiff;
      regionNum = 0;

      while (regionNum < numParts){

        for (int i=0; i < weightDim; i++)
          testVec[i] += regionSums[regionNum];
  
        diffVec = testVec;
        diffVec.Scale(-1.0);
        diffVec += targetLeftVector;
  
        scalar_t testDiff = diffVec.Norm2();
        
        if (testDiff >= target)
          break;

        prevTestDiff = testDiff;
        regionNum++;
      }
  
      if (regionNum == numParts)
        fail = true;
    }
    else{    // the part sizes for each weight dimension are the same
  
      target = targetLeftScalar;
      testDiff = totalWeightLeft; 
      prevTestDiff = testDiff;
      regionNum = 0;
  
      while (regionNum < numParts){
     
        testDiff += regionSums[regionNum];
        if (testDiff >= target)
          break;
  
        prevTestDiff = testDiff;
        regionNum++;
      }
    }
  
    int leftBoundary = regionNum;
    int rightBoundary = regionNum+1;
    int rightmostLeftRegion=0, leftmostRightRegion=0;
  
    if (!fail){
      scalar_t diffLeftCut = target - prevTestDiff;
      scalar_t diffRightCut = testDiff - target;
  
      if (diffLeftCut < diffRightCut){
        imbalance = diffLeftCut / target;
        if (imbalance <= tolerance){
          done = true;
          min = max = cutValue = boundaries[leftBoundary];
          rightmostLeftRegion = regionNum-1;
          leftmostRightRegion = regionNum;
        }
      }
      else{
        imbalance = diffRightCut / target;
        if (imbalance <= tolerance){
          done = true;
          min = max = cutValue = boundaries[rightBoundary];
          rightmostLeftRegion = regionNum;
          leftmostRightRegion = regionNum+1;
        }
      }
    
      if (!done) {
        if (leftBoundary != 0)
          min = boundaries[leftBoundary];
        if (rightBoundary != endBoundary)
          max = boundaries[rightBoundary];
        rightmostLeftRegion = regionNum-1;
        leftmostRightRegion = regionNum+1;
      }
    }    // if !fail
  

    if (!fail && done && averageCuts){
      scalar_t localExtrema[2], globalExtrema[2];
      localExtrema[0] = regionMax[rightmostLeftRegion];
      localExtrema[1] = regionMin[leftmostRightRegion] * -1.0;

      try{
        reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_MAX, 2,
          localExtrema, globalExtrema);
      }
      Z2_THROW_OUTSIDE_ERROR(*env)

      scalar_t globalLeftOfCut = globalExtrema[0];
      scalar_t globalRightOfCut = globalExtrema[1] * -1.0;

      cutValue = (globalRightOfCut + globalLeftOfCut) * 0.5;
    }

    if (!done && !multiplePartSizeSpecs){

      scalar_t globalWeightMovedRight = 0.0;
      gno_t localCountMovedRight = 0;

      // See if moving coordinates on the boundary will help.
      // TODO: This can be fixed to work for multiplePartSizeSpecs

      ArrayView<lno_t> indexArray;
      if (useIndices)
        indexArray = ArrayView<lno_t>(index.getRawPtr(), numCoords);

      ArrayView<const scalar_t> coords(coordValue, numAllCoords);
      ArrayView<scalar_t> sums(regionSums, numRegions);

      try{
        testCoordinatesOnRightBoundary<lno_t, gno_t, scalar_t>( 
          env, comm, rectilinearBlocks,
          totalWeightLeft, target, 
          rightBoundary, boundaries, boundarySum.view(0, numBoundaries),
          weight.view(0, weightDim), mcnorm, 
          coords, epsilon, lrFlags.view(0, numCoords), indexArray,
          sums,
          globalWeightMovedRight,           // output
          localCountMovedRight);            // output
      }
      Z2_FORWARD_EXCEPTIONS

      if (globalWeightMovedRight > 0.0){

        regionSums[rightmostLeftRegion+1] -= globalWeightMovedRight;
        numRemaining -= localCountMovedRight;

        // Check to see if moving the coordinates balanced them.

        scalar_t testWeightLeft = totalWeightLeft;
        for (int i=0; i <= rightmostLeftRegion+1; i++)
          testWeightLeft += regionSums[i];
  
        scalar_t newDiff = 0.0;
        if (testWeightLeft < targetLeftScalar)
          newDiff = targetLeftScalar - testWeightLeft;
        else
          newDiff = testWeightLeft - targetLeftScalar;
  
        imbalance = newDiff / targetLeftScalar;
        if (imbalance <= tolerance){
          done = true;
          min = max = cutValue = boundaries[rightBoundary];
          rightmostLeftRegion++;
        }
      }
    }

    if (!fail){

      for (int i=0; i <= rightmostLeftRegion; i++){
        totalWeightLeft += regionSums[i];
      }

      for (size_t i=0; i < numCoords; i++){
        if (lrFlags[i] != leftFlag && lrFlags[i] != rightFlag){
          if (lrFlags[i] <= rightmostLeftRegion){
            lrFlags[i] = leftFlag;
            numRemaining--;
            localCountLeft++;
          }
          else if (lrFlags[i] >= leftmostRightRegion){
            lrFlags[i] = rightFlag;
            numRemaining--;
          }
          else{
            lrFlags[i] = unsetFlag;   // still to be determined
          }
        }
      }

      if (env->doStatus()){
        ostringstream ossLeft;
        ostringstream ossRight;
        ossLeft << "left: ";
        ossRight << "right: ";
        for (size_t i=0; i < numCoords; i++){
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
    }   // if !fail

  }  // while !done

  totalWeightRight = totalWeight - totalWeightLeft;

  env->globalInputAssertion(__FILE__, __LINE__, "partitioning not solvable",
    done && !fail, DEBUG_MODE_ASSERTION, comm);

  if (env->doStatus())
    env->debug(DETAILED_STATUS, string("Exiting BSPfindCut"));

  if (nprocs > 1)
    env->timerStop("BSPfindCut");
  else
    env->timerStop("BSPfindCut - serial");

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
  // initialize return values
  size_t numLocalCoords = vectors->getLocalLength();
  cutDimension = 0;
  cutValue = imbalance = weightLeftHalf = weightRightHalf = 0.0;
  numPartsLeftHalf = 0;

  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;
  typedef StridedData<lno_t, scalar_t> input_t;

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

    for (int i=0; i < weightDim; i++)
      if (!uniformWeights[i])
        numNonUniform++;

    int wgt1 = vectors->getNumVectors() - numNonUniform;

    if (weightDim == 1){
      if (numNonUniform == 0)
        totalWeight = numLocalCoords;
      else{
        const scalar_t *val = vectors->getData(wgt1).getRawPtr();
        for (size_t i=0; i < numLocalCoords; i++)
          totalWeight += val[i];
      }
    }
    else{  // need to add up total normed weight
      Array<input_t> wgts(weightDim);
      for (int i=0; i < weightDim; i++){
        if (!uniformWeights[i]){
          wgts[i] = input_t(vectors->getData(wgt1++), 1);
        }
      }

      partId_t numParts, numNonemptyParts;
      ArrayRCP<MetricValues<scalar_t> > metrics;
      ArrayRCP<scalar_t> weightSums;
    
      globalSumsByPart<scalar_t, unsigned char, lno_t>(
        env, comm, lrflags, 
        0, 0,                    // says to ignore flags
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
}


/*! \brief Perform RCB on the local process only.
 *
 *   \param env the Environment for the application.
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
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;

  RCP<Comm<int> > comm(new Teuchos::SerialComm<int>);  

  int numLocalCoords=0;
  bool useIndices;

  if (index.size() == 0){
    // First time through there are no indices.
    useIndices = false;
    numLocalCoords = vectors->getLocalLength();
  }
  else{
    useIndices = true;
    numLocalCoords = index.size();
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

    serialRCB<mvector_t, Adapter>(env, params, numTestCuts, tolerance, 
      coordDim, vectors, leftIndices,
      uniformWeights.view(0, weightDim), solution,
      part0, newPart1, partNum);

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

    serialRCB<mvector_t, Adapter>(env, params, numTestCuts, tolerance, 
      coordDim, vectors, rightIndices,
      uniformWeights.view(0, weightDim), solution,
      newPart0, part1, partNum);

    delete [] newIndex;
  }
}

}// namespace Zoltan2

#endif

