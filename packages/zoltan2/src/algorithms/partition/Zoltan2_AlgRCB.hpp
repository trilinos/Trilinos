// @HEADER
//***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_AlgRCB.hpp
    \brief Contains the recursive coordinate bisection algorthm.
*/

#ifndef _ZOLTAN2_ALGRCB_HPP_
#define _ZOLTAN2_ALGRCB_HPP_

#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Metric.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_GetParameter.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Tpetra_Vector.hpp>

#include <sstream>
#include <string>
#include <cmath>
#include <bitset>

namespace Zoltan2{

/*! \brief The boolean parameters of interest to the RCB algorithm.
 */
enum rcbParams {
  rcb_fastSolution,         /*!< speed_versus_quality = speed */
  rcb_goodSolution,         /*!< speed_versus_quality = quality */
  rcb_balanceSolution,     /*!< speed_versus_quality = balance */
  rcb_lowMemory,            /*!< memory_versus_speed = memory */
  rcb_lowRunTime,           /*!< memory_versus_speed = speed */
  rcb_balanceMemoryRunTime, /*!< memory_versus_speed = balance */
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
  leftFlag = 0xfe,     /*!< 254 */
  rightFlag = 0xff     /*!< 255 */
};

template <typename mvector_t>
  void getCutDimension( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm, 
    int coordDim, const RCP<mvector_t> &vectors, 
    ArrayView<typename mvector_t::local_ordinal_type> index,
    int &dimension, 
    typename mvector_t::scalar_type &minCoord, 
    typename mvector_t::scalar_type &maxCoord);

template <typename mvector_t>
 void serialRCB( const RCP<const Environment> &env,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, typename mvector_t::scalar_type tolerance, 
    int coordDim, const RCP<mvector_t> &vectors, 
    ArrayView<typename mvector_t::local_ordinal_type> index,
    const ArrayView<bool> uniformWeights,
    const ArrayView<ArrayRCP<typename mvector_t::scalar_type > > partSizes,
    partId_t part0, partId_t part1, ArrayView<partId_t> partNum);

template <typename mvector_t>
  void BSPfindCut( const RCP<const Environment> &env,
    const RCP<const Teuchos::Comm<int> > &comm, 
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, 
    typename mvector_t::scalar_type tolerance, int cutDim,
    int coordDim, const RCP<mvector_t> &vectors,
    ArrayView<typename mvector_t::local_ordinal_type> index,
    ArrayView<double> &fractionLeft, 
    ArrayView<bool> uniformWeights,
    typename mvector_t::scalar_type coordGlobalMin, 
    typename mvector_t::scalar_type coordGlobalMax, 
    typename mvector_t::scalar_type &cutValue, 
    ArrayView<unsigned char> lrflags, 
    typename mvector_t::scalar_type &imbalance);

template <typename mvector_t>
  void determineCut( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, typename mvector_t::scalar_type tolerance,
    int coordDim, const RCP<mvector_t> &vectors, 
    const ArrayView<bool> uniformWeights,
    const ArrayView<ArrayRCP<typename mvector_t::scalar_type > > partSizes,
    partId_t part0, partId_t part1,
    ArrayView<unsigned char> lrflags,
    int &cutDimension, 
    typename mvector_t::scalar_type &cutValue, 
    typename mvector_t::scalar_type &imbalance,
    partId_t &numPartsLeftHalf);

template <typename mvector_t>
  void migrateData( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    ArrayView<unsigned char> lrflags,
    RCP<mvector_t> &vectors, int &numProcsLeftHalf);

template <typename scalar_t>
  void getFractionLeft( const RCP<const Environment> &env,
    partId_t part0, partId_t part1,
    const ArrayView<ArrayRCP<scalar_t> > partSizes,
    ArrayRCP<double> &fractionLeft, partId_t &numPartsLeftHalf);

/*! \brief Recursive coordinate bisection algorithm.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm  the communicator for the problem
 *  \param coords    a CoordinateModel with user data
 *  \param solution  a PartitioningSolution, on input it 
 *      contains part information, on return it also contains 
 *      the solution and quality metrics.
 *                    
 *   \todo timing and memory usage profiling and debug messages
 *   \todo  catch errors and pass back
 *   \todo write the rcb tree back to the solution
 *   \todo for "repartition", start with the tree in the solution
 *   \todo implement rectilinear_blocks and average_cuts
 *   \todo  work on performance issues.  Some of the global
 *               communication can probably be consolidated
 *                into fewer messages.
 */

template <typename Adapter>
void AlgRCB(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
  const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &coords, 
  RCP<PartitioningSolution<Adapter> > &solution
) 
{
  typedef typename Adapter::node_t node_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  std::bitset<NUM_RCB_PARAMS> params;

  bool isSet;
  string strChoice;
  int intChoice;

  ////////////////////////////////////////////////////////
  // Library parameters of interest:
  //
  //    speed_versus_quality
  //    memory_versus_speed

  if (env->doStatus())
    env->debug(DETAILED_STATUS, string("Entering AlgPartRCB"));

  const Teuchos::ParameterList &pl = env->getParameters();

  getParameterValue<string>(pl, "memory_versus_speed", 
    isSet, strChoice);

  if (isSet && strChoice==string("memory"))
    params.set(rcb_lowMemory);
  else if (isSet && strChoice==string("speed"))
    params.set(rcb_lowRunTime);
  else
    params.set(rcb_balanceMemoryRunTime);

  getParameterValue<string>(pl, "speed_versus_quality",
    isSet, strChoice);

  if (isSet || strChoice==string("speed"))
    params.set(rcb_fastSolution);
  else if (isSet || strChoice==string("quality"))
    params.set(rcb_goodSolution);
  else 
    params.set(rcb_balanceSolution);

  ////////////////////////////////////////////////////////
  // Partitioning problem parameters of interest:
  //    objective
  //    imbalance_tolerance

  scalar_t imbalanceTolerance;

  getParameterValue<string>(pl, "partitioning",
    "objective", isSet, strChoice);

  if (isSet && strChoice == string("balance_object_count"))
    params.set(rcb_balanceCount);
  else if (isSet && strChoice == 
    string("multicriteria_minimize_total_weight"))
    params.set(rcb_minTotalWeight);
  else if (isSet && strChoice == 
    string("multicriteria_minimize_maximum_weight"))
    params.set(rcb_minMaximumWeight);
  else if (isSet && strChoice == 
    string("multicriteria_balance_total_maximum"))
    params.set(rcb_balanceTotalMaximum);
  else
    params.set(rcb_balanceWeight);

  double tol;

  getParameterValue<double>(pl, "partitioning",
    "imbalance_tolerance", isSet, tol);

  if (!isSet)
    imbalanceTolerance = .1;
  else
    imbalanceTolerance = tol - 1.0;

  if (imbalanceTolerance <= 0)
    imbalanceTolerance = 10e-4;  // TODO - what's a good choice

  ////////////////////////////////////////////////////////
  // Geometric partitioning problem parameters of interest:
  //    average_cuts
  //    rectilinear_blocks
  //    bisection_num_test_cuts (experimental)

  getParameterValue<int>(pl, "partitioning", "geometric",
    "average_cuts", isSet, intChoice);

  if (isSet && intChoice==1)
    params.set(rcb_averageCuts);

  getParameterValue<int>(pl, "partitioning", "geometric",
    "rectilinear_blocks", isSet, intChoice);

  if (isSet && intChoice==1)
    params.set(rcb_rectilinearBlocks);

  getParameterValue<int>(pl, "partitioning", "geometric",
    "bisection_num_test_cuts", isSet, intChoice);

  int numTestCuts = 5;
  if (isSet)
    numTestCuts = intChoice;

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

  ArrayRCP<bool> uniformWeights(new bool [criteriaDim], 0, criteriaDim, true);
  Array<ArrayRCP<const scalar_t> > weights(criteriaDim);

  if (weightDim == 0)             // uniform weights are implied
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

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.
  // If the part sizes for a given criteria are not uniform,
  // then they are values that sum to 1.0.

  size_t numGlobalParts = solution->getGlobalNumberOfParts();

  Array<bool> uniformParts(criteriaDim);
  Array<ArrayRCP<scalar_t> > partSizes(criteriaDim);

  for (int wdim = 0; wdim < criteriaDim; wdim++){
    if (solution->criteriaHasUniformPartSizes(wdim)){
      uniformParts[wdim] = true;
    }
    else{
      scalar_t *tmp = new scalar_t [numGlobalParts];
      env->localMemoryAssertion(__FILE__, __LINE__, numGlobalParts, tmp) ;
    
      for (partId_t i=0; i < numGlobalParts; i++){
        tmp[i] = solution->getCriteriaPartSize(i, wdim);
      }

      partSizes[wdim] = arcp(tmp, 0, numGlobalParts);
    }
  }

  // It may not be possible to solve the partitioning problem
  // if we have multiple weight dimensions with part size
  // arrays that differ. So let's be aware of this possibility.

  bool multiplePartSizeSpecs = false;

  if (weightDim > 1){
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
    map = rcp(new map_t(numGlobalCoords, gnos, gnoMin, problemComm));
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

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  // The algorithm
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  RCP<const Comm<int> > comm = rcp_const_cast<Comm<int> >(problemComm);
  partId_t part0 = 0;
  partId_t part1 = numGlobalParts-1;
  int sanityCheck = numGlobalParts;
  int groupSize = comm->getSize();
  int rank = comm->getRank();

  scalar_t imbalanceReductionFactor = 1.0;
  if (numGlobalParts > 3)
    imbalanceReductionFactor = log(scalar_t(numGlobalParts));

  imbalanceTolerance /= imbalanceReductionFactor;

  while (part1>part0 && groupSize>1 && sanityCheck--){

    ////////////////////////////////////////////////////////
    // Which coordinates are left and which are right?

    Array<unsigned char> lrflags(numLocalCoords);
    scalar_t cutValue;  // TODO eventually save this for user
    int cutDimension;
    scalar_t imbalance;
    partId_t leftHalfNumParts;

    try{
      determineCut<mvector_t>(env, comm, 
        params, numTestCuts, imbalanceTolerance,
        coordDim, mvector, 
        uniformWeights.view(0,criteriaDim), partSizes.view(0, criteriaDim), 
        part0, part1,
        lrflags.view(0, numLocalCoords), 
        cutDimension, cutValue, imbalance, leftHalfNumParts);
    }
    Z2_FORWARD_EXCEPTIONS

    ////////////////////////////////////////////////////////
    // Migrate the multivector of data.

    int leftHalfNumProcs;

    try{
      // on return mvector has my new data
      migrateData<mvector_t>( env, comm, lrflags.view(0,numLocalCoords), 
        mvector, leftHalfNumProcs);
    }
    Z2_FORWARD_EXCEPTIONS

    env->localBugAssertion(__FILE__, __LINE__, "num procs in half",
      leftHalfNumProcs > 0 && leftHalfNumProcs < groupSize,
      BASIC_ASSERTION);

    ////////////////////////////////////////////////////////
    // Divide into two subgroups.

    int *ids = NULL;

    if (rank < leftHalfNumProcs){
      groupSize = leftHalfNumProcs;
      ids = new int [groupSize];
      env->localMemoryAssertion(__FILE__, __LINE__, groupSize, ids);
      for (int i=0; i < groupSize; i++)
        ids[i] = i;
      part1 = part0 + leftHalfNumParts - 1;
    }
    else{
      groupSize = comm->getSize() - leftHalfNumProcs;
      rank -= leftHalfNumProcs;
      ids = new int [groupSize];
      env->localMemoryAssertion(__FILE__, __LINE__, groupSize, ids);
      for (int i=0; i < groupSize; i++)
        ids[i] = i + leftHalfNumProcs;
      part0 += leftHalfNumParts;
    }

    ArrayView<const int> idView(ids, groupSize);
    RCP<Comm<int> > subComm(comm->createSubcommunicator(idView));
    comm = rcp_const_cast<const Comm<int> >(subComm);

    ////////////////////////////////////////////////////////
    // Create a new multivector for my smaller group.

    ArrayView<const gno_t> gnoList = mvector->getMap()->getNodeElementList();
    size_t localSize = mvector->getLocalLength();
    size_t globalSize = localSize;
  
    pair<gno_t, gno_t> minMax = 
      z2LocalMinMax<gno_t>(gnoList.getRawPtr(), localSize);
    gno_t localMin = minMax.first;
    gno_t globalMin = localMin;

    if (groupSize > 1){
      try{
        reduceAll<int, size_t>(
          *comm, Teuchos::REDUCE_SUM, 1, &localSize, &globalSize);
      }
      Z2_THROW_OUTSIDE_ERROR(*env)
  
      try{
        reduceAll<int, gno_t>(
          *comm, Teuchos::REDUCE_MIN, 1, &localMin, &globalMin);
      }
      Z2_THROW_OUTSIDE_ERROR(*env)
    }
  
    RCP<map_t> subMap;
    try{
      subMap= rcp(new map_t(globalSize, gnoList, globalMin, comm));
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
  
    mvector = subMvector;
    numLocalCoords = mvector->getLocalLength();
  } 

  env->localBugAssertion(__FILE__, __LINE__, "partitioning failure", 
    sanityCheck, BASIC_ASSERTION);

  partId_t *tmp = new partId_t [numLocalCoords];
  env->localMemoryAssertion(__FILE__, __LINE__, numLocalCoords, tmp);
  ArrayRCP<partId_t> partId(tmp, 0, numLocalCoords, true);

  if (part1 > part0){    // Serial partitioning

    // scalar_t cutValue;   TODO
    // int cutDimension;
    // scalar_t imbalance;

    if (multiplePartSizeSpecs)
      numTestCuts = 3;
    else
      numTestCuts = 1;

    try{
      ArrayView<lno_t> emptyIndex;

      serialRCB<mvector_t>(env, params,
        numTestCuts, imbalanceTolerance,
        coordDim, mvector, emptyIndex, 
        uniformWeights.view(0,criteriaDim), partSizes.view(0, criteriaDim), 
        part0, part1, partId.view(0,numLocalCoords));
    }
    Z2_FORWARD_EXCEPTIONS
  }
  else{
    for (lno_t i=0; i < partId.size(); i++)
      partId[i] = part0;
  }

  ////////////////////////////////////////////////////////
  // Done: Compute quality metrics and update the solution
  
  ArrayRCP<MetricValues<scalar_t> > metrics;
  partId_t numParts, numNonemptyParts;

  multiCriteriaNorm mcnorm = normBalanceTotalMaximum;

  if (params.test(rcb_minMaximumWeight))
    mcnorm = normMinimizeMaximumWeight;
  else if (params.test(rcb_minTotalWeight))
    mcnorm = normMinimizeTotalWeight;

  ArrayRCP<input_t> objWgt(new input_t [criteriaDim], 0, criteriaDim, true);
  Array<ArrayView<scalar_t> > partSizeArrays(criteriaDim);

  for (int wdim = 0, widx=coordDim; wdim < criteriaDim; wdim++){
    if (!uniformWeights[wdim]){
      objWgt[wdim] = input_t(mvector->getData(widx++), 1);
    }
    if (partSizes[wdim].size() > 0)
      partSizeArrays[wdim] = 
        ArrayView<scalar_t>(partSizes[wdim].getRawPtr(), numGlobalParts);
  }

  objectMetrics<scalar_t, lno_t>(
    env, problemComm, numGlobalParts,                  // input
    partSizeArrays.view(0, criteriaDim),               // input
    partId.view(0, numLocalCoords),                    // input
    objWgt.view(0, criteriaDim), mcnorm,               // input
    numParts, numNonemptyParts, metrics);              // output

  ArrayRCP<const gno_t> gnoList = 
    arcpFromArrayView(mvector->getMap()->getNodeElementList());

  solution->setParts(gnoList, partId, metrics);
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
 *   \param fractionLeft  the size of the left part for each weight
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
 *   \param imbalance on return is the imbalance for the 
 *              computed partitioning (0, 1).
 * 
 *  \todo need a serial method that creates n parts
 *  \todo a separate simpler function when weightDim <= 1
 *  \todo During the first global comm, ensure all procs have same values
 *           for cutDim and fractionLeft.
 */

template <typename mvector_t>
  void BSPfindCut(
    const RCP<const Environment> &env,
    const RCP<const Teuchos::Comm<int> > &comm,
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
    ArrayView<unsigned char> lrFlags,         // output
    typename mvector_t::scalar_type &imbalance)        // output
{
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;
  typedef typename mvector_t::node_type node_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  bool useIndices = index.size() > 0;

  // Find the coordinate values and weights.

  int numAllCoords = vectors->getLocalLength();
  int numCoords = 0;

  if (useIndices)
    numCoords = index.size();
  else 
    numCoords = numAllCoords;

  int weightDim = uniformWeights.size();
  int numNonUniformWeights = 0;
  for (int i=0; i < weightDim; i++){
    if (!uniformWeights[i])
      numNonUniformWeights++;
  }

  set<double> differentSizes;
  for (int i=0; i < weightDim; i++)
    differentSizes.insert(fractionLeft[i]);

  // TODO: within some epsilon
  int numDifferentPartSizes = differentSizes.size();

  if (numDifferentPartSizes > 1 && numTestCuts < 3)
    numTestCuts = 3;

  const scalar_t *coordValue = vectors->getData(cutDim).getRawPtr();

  ArrayRCP<input_t> weight;
  if (weightDim > 0){
    input_t *info = new input_t [weightDim];
    env->localMemoryAssertion(__FILE__, __LINE__, weightDim, info);
    weight = arcp(info, 0, weightDim);
  }

  for (int wdim = 0, widx=coordDim; wdim < weightDim; wdim++){
    if (!uniformWeights[wdim]){
      weight[wdim] = input_t(vectors->getData(widx++), 1);
    }
  }

  // Multicriteria norm

  multiCriteriaNorm mcnorm = normBalanceTotalMaximum;

  if (params.test(rcb_minMaximumWeight))
    mcnorm = normMinimizeMaximumWeight;
  else if (params.test(rcb_minTotalWeight))
    mcnorm = normMinimizeTotalWeight;
  
  // Goal is globally find one cut that comes close to partSizeLeft.

  Epetra_SerialDenseVector partSizeLeft( 
    View, fractionLeft.getRawPtr(), weightDim);
  
  // Where do we make the first test cuts?
  //
  //   min     1     2     3     max
  //
  // Space is [min, max].  
  // Three test cuts: 1, 2 and 3.
  // 4 regions: [min,1] [1,2] [2,3] [3,max].
  // 5 boundaries: min, 1, 2, 3, and max.

  int numRegions = numTestCuts + 1;
  int numBoundaries = numTestCuts + 2;
  int endBoundary = numBoundaries - 1;
  std::vector<scalar_t> boundaries(numBoundaries);
  typename std::vector<scalar_t>::iterator foundCut;

  memset(lrFlags.getRawPtr(), 0, numCoords); // 0 at top of loop means unset
  bool done=false;
  bool fail=false;
  scalar_t min = coordGlobalMin;
  scalar_t max = coordGlobalMax;
  lno_t numRemaining = numCoords;
  lno_t sanityCheck = numCoords;

  double totalWeight = 0, totalWeightLeft=0;
  double targetLeftScalar = 0;
  double targetLeftNorm = 0;
  Epetra_SerialDenseVector targetLeftVector(weightDim);

  int firstNonLeftRegion = 0;

  while (!done && !fail && sanityCheck--){

    scalar_t diff = (max - min) / numRegions;
    boundaries[0] = min;
    for (int i=0; i < numRegions; i++)
      boundaries[i+1] = boundaries[i] + diff;

    // Be sure to catch all remaining within the boundaries.
    boundaries[endBoundary] += 1.0;

    if (numRemaining > 0){

      // Assign each of my points to a region.
      // lower_bound finds the first cut f, such that f >= coordValue[i].
      // So for now, objects that are on the cut boundary go into the
      // region on the "left" side.

      if (useIndices){
        for (size_t i=0; i < numCoords; i++){
    
          if (lrFlags[i] != 0) 
            continue;
    
          if (numRegions > 2){
            foundCut = std::lower_bound(boundaries.begin(), boundaries.end(), 
                coordValue[index[i]]);
            
            env->localBugAssertion(__FILE__, __LINE__, "search cuts", 
              foundCut != boundaries.end(), BASIC_ASSERTION);
            
            lrFlags[i] = (unsigned char)(foundCut - boundaries.begin() - 1);
          }
          else{
            lrFlags[i] = (coordValue[index[i]] <= boundaries[1] ? 0 : 1);
          }
        }
      }
      else {
        for (size_t i=0; i < numCoords; i++){
    
          if (lrFlags[i] != 0)
            continue;
    
          if (numRegions > 2){
      
            foundCut = std::lower_bound(boundaries.begin(), boundaries.end(), 
                coordValue[i]);

            env->localBugAssertion(__FILE__, __LINE__, "search cuts", 
              foundCut != boundaries.end(), BASIC_ASSERTION);
            
            lrFlags[i] = (unsigned char)(foundCut - boundaries.begin() - 1);
          }
          else{
            lrFlags[i] = (coordValue[i] <= boundaries[1] ? 0 : 1);
          }
        }
      }
    }

    partId_t numParts, numNonemptyParts;
    ArrayRCP<MetricValues<scalar_t> > metrics;
    ArrayRCP<scalar_t> weightSums;

    // Get the global sums in each region

    if (params.test(rcb_balanceCount) || numNonUniformWeights == 0){

      // Ignore weights and sum the objects in each part.

      Array<StridedData<lno_t, scalar_t> > noWeights(weightDim);

      globalSumsByPart<scalar_t, unsigned char, lno_t>(
        env, comm, 
        lrFlags,                      // part assignments
        0, numRegions - 1,            // parts to be included in count
        noWeights.view(0, weightDim),
        mcnorm,                       // ignored
        numParts,
        numNonemptyParts,
        metrics,
        weightSums);    // numParts sums of objects in each part
    }
    else{
      globalSumsByPart<scalar_t, unsigned char, lno_t>(
        env, comm, 
        lrFlags,                      // part assignments
        0, numRegions - 1,            // parts to be included in count
        weight.view(0, weightDim),
        mcnorm,
        numParts,
        numNonemptyParts,
        metrics,
        weightSums); // object sums, normed weight sums, and if
                     // weightDim>1 individual weight sums
    }

    // If the objective is to balance object count, get the
    // object count imbalance.  Otherwise get the normed 
    // weight imbalance.

    scalar_t *values = weightSums.getRawPtr();   // object count sums

    if (!params.test(rcb_balanceCount) && numNonUniformWeights > 0){
      values += numParts;   // normed weight sums
    }

    if (totalWeight == 0){   // first time through only

      for (int i=0; i < numParts; i++)
        totalWeight += values[i];

      partSizeLeft.Scale(totalWeight);
      targetLeftVector = partSizeLeft;

      targetLeftScalar = targetLeftVector[0];
      targetLeftNorm = targetLeftVector.Norm2();
      totalWeightLeft = 0;
    }

    if (numDifferentPartSizes > 1){   

      // more complex: if we have multiple weight dimensions, the
      //   weights are non-uniform, and the part sizes requested
      //   for each each weight dimension differ, then we may not
      //   be able to reach the imbalance tolerance.
      //
      // TODO: discuss how to handle this case.
  
      Epetra_SerialDenseVector testVec(weightDim);
      for (int i=0; i < weightDim; i++)
        testVec[i] = totalWeightLeft;
  
      // If we cut before the first region, what would be 
      //  difference |target - actual|^2
  
      Epetra_SerialDenseVector diffVec = testVec;
      diffVec.Scale(-1.0);
      diffVec += targetLeftVector;
      scalar_t minDiff = diffVec.Norm2();
  
      int num = 0;
      int bestNum = -1;
  
      // Find the first cut where the difference does not improve.

      while (num < numParts){   // region "num"

        // |target-actual| if cut is after region num

        for (int i=0; i < weightDim; i++)
          testVec[i] += values[num];
  
        testVec.Scale(-1.0);
        diffVec = testVec;
        diffVec += targetLeftVector;
       
        scalar_t diff = diffVec.Norm2(); 

        if (diff < minDiff){
          minDiff = diff;
          bestNum = num;
        }
        else
          break;

        num++;
      }

      //
      // |target-actual| / |target|
      //
      imbalance = minDiff / targetLeftNorm;

      if (imbalance <= tolerance){
        cutValue = boundaries[bestNum + 1];
        min = max = cutValue;
        firstNonLeftRegion = bestNum + 1;
        done = true;
      }
      else if ((bestNum >= 0) && (bestNum <= numParts-2)){
        // this only works if number of regions is at least three
        min = boundaries[bestNum];    // left of this region
        max = boundaries[bestNum+2];  // right of next region
        if (bestNum+2 == endBoundary)
          max -= 1.0;
        firstNonLeftRegion = bestNum;
      }
      else{
        fail = 1;  // We can't solve this.  Norms don't have a min.
      }
    }
    else{    // the part sizes for each weight dimension are the same

      scalar_t testLeft = totalWeightLeft; 
      scalar_t prevLeft = testLeft;
      int num = 0;

      while (num < numParts){
   
        testLeft += values[num];
        if (testLeft > targetLeftScalar)
          break;

        prevLeft = testLeft;
        num++;
      }

      scalar_t diffLeftCut = targetLeftScalar - prevLeft;
      scalar_t diffRightCut = testLeft - targetLeftScalar;

      if (diffLeftCut < diffRightCut){
        imbalance = diffLeftCut / targetLeftScalar;
        if (imbalance <= tolerance){
          cutValue = boundaries[num];
          min = max = cutValue;
          done = true;
          firstNonLeftRegion = num;
        }
      }
      else{
        imbalance = diffRightCut / targetLeftScalar;
        if (imbalance <= tolerance){
          cutValue = boundaries[num+1];
          min = max = cutValue;
          done = true;
          firstNonLeftRegion = num+1;
        }
      }

      if (!done){
        min = boundaries[num];
        max = boundaries[num+1];
        if (num+1 == endBoundary)
          max -= 1.0;
        firstNonLeftRegion = num;
      }
    }

    if (!fail){

      for (int i=0; i < firstNonLeftRegion; i++){
        totalWeightLeft += values[i];
      }
  
      if (useIndices){
        for (size_t i=0; i < numCoords; i++){
          if (lrFlags[i] != leftFlag && lrFlags[i] != rightFlag){
            if (coordValue[index[i]] <= min){
              lrFlags[i] = leftFlag;
              numRemaining--;
            }
            else if (coordValue[index[i]] > max){
              lrFlags[i] = rightFlag;
              numRemaining--;
            }
            else{
              lrFlags[i] = 0;   // still to be determined
            }
          }
        }
      }
      else{
        for (size_t i=0; i < numCoords; i++){
          if (lrFlags[i] != leftFlag && lrFlags[i] != rightFlag){
            if (coordValue[i] <= min){
              lrFlags[i] = leftFlag;
              numRemaining--;
            }
            else if (coordValue[i] > max){
              lrFlags[i] = rightFlag;
              numRemaining--;
            }
            else{
              lrFlags[i] = 0;   // still to be determined
            }
          }
        }
      }
    }

  }    // next set of test cuts

  env->globalInputAssertion(__FILE__, __LINE__, "partitioning not solvable",
    done && !fail, DEBUG_MODE_ASSERTION, comm);
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
 *   \param partSizes is a list of part sizes for each weight dimension.
 *          Sizes should sum to 1.0 for each dimension.  If
 *          <tt>partSizes[w].size()</tt> is zero it means that
 *          part sizes for weight dimension \c w are all the same.
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
 */

template <typename mvector_t>
  void determineCut(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, typename mvector_t::scalar_type tolerance,
    int coordDim, const RCP<mvector_t> &vectors,
    const ArrayView<bool> uniformWeights,
    const ArrayView<ArrayRCP<typename mvector_t::scalar_type> > partSizes,
    partId_t part0, 
    partId_t part1,
    ArrayView<unsigned char> lrflags,   // output
    int &cutDimension,                  // output
    typename mvector_t::scalar_type &cutValue,   // output
    typename mvector_t::scalar_type &imbalance,  // output
    partId_t &numPartsLeftHalf)                  // output
{
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;

  ///////////////////////////////////////////////////////
  // Pick a cut direction.

  scalar_t globalMinCoord, globalMaxCoord;
  ArrayView<lno_t> emptyIndex;

  getCutDimension<mvector_t>(env, comm, coordDim, vectors, emptyIndex,
    cutDimension, globalMinCoord, globalMaxCoord);

  ///////////////////////////////////////////////////////
  // Compute part sizes for the two parts.

  ArrayRCP<double> fractionLeft;
  size_t weightDim = partSizes.size();

  getFractionLeft<scalar_t>(env, part0, part1, 
    partSizes.view(0, weightDim), fractionLeft, numPartsLeftHalf);

  ///////////////////////////////////////////////////////
  // Divide the coordinates into balanced left and right
  // halves.

  size_t numLocalCoords = vectors->getLocalLength();
  ArrayView<lno_t> emptyIndices;

  try{
    BSPfindCut<mvector_t>( env, comm,
      params, numTestCuts, tolerance,
      cutDimension, coordDim, vectors, emptyIndices,
      fractionLeft.view(0, weightDim), uniformWeights.view(0, weightDim),
      globalMinCoord, globalMaxCoord,
      cutValue, lrflags.view(0, numLocalCoords),
      imbalance);
  }
  Z2_FORWARD_EXCEPTIONS
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
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    ArrayView<unsigned char> lrflags,
    RCP<mvector_t> &vectors,    // on return is the new data
    int &leftNumProcs)          // on return is num procs with left data
{
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;
  typedef typename mvector_t::global_ordinal_type gno_t;

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
  leftNumProcs = static_cast<int>(nprocs * leftFraction);

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
 *   \param partSizes is a list of part sizes for each weight dimension.
 *          Sizes should sum to 1.0 for each dimension.  If
 *          <tt>partSizes[w].size()</tt> is zero it means that
 *          part sizes for weight dimension \c w are all the same.
 *   \param part0  is the first part of the parts to bisected.
 *   \param part1  is the last part of the parts to bisected.
 *   \param partNum on return <tt>partNum[i]</tt> is the new
 *                part number for coordinate \c i.
 */

template <typename mvector_t>
 void serialRCB(
    const RCP<const Environment> &env,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, 
    typename mvector_t::scalar_type tolerance, 
    int coordDim,
    const RCP<mvector_t> &vectors, 
    ArrayView<typename mvector_t::local_ordinal_type> index,
    const ArrayView<bool> uniformWeights,
    const ArrayView<ArrayRCP<typename mvector_t::scalar_type > > partSizes,
    partId_t part0, 
    partId_t part1,
    ArrayView<partId_t> partNum)   // output
{
  typedef typename mvector_t::scalar_type scalar_t;
  typedef typename mvector_t::local_ordinal_type lno_t;

  RCP<const Comm<int> > comm(new Teuchos::SerialComm<int>);  

  int numLocalCoords=0;
  bool useIndices;

  if (index.size() == 0){
    useIndices = false;
    numLocalCoords = vectors->getLocalLength();
  }
  else{
    useIndices = true;
    numLocalCoords = index.size();
  }

  ///////////////////////////////////////////////////////
  // Are we done?

  if (numLocalCoords == 0)
    return;

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
  // Compute part sizes for the two parts.

  ArrayRCP<double> fractionLeft;
  partId_t numPartsLeftHalf;
  int weightDim = uniformWeights.size();

  try{
    getFractionLeft<scalar_t>(env, part0, part1, 
      partSizes.view(0, weightDim), 
      fractionLeft, numPartsLeftHalf);
  }
  Z2_FORWARD_EXCEPTIONS

  ///////////////////////////////////////////////////////
  // Divide into balanced left and right halves.

  scalar_t imbalance, cutValue;  //unused for now

  unsigned char *newFlags = new unsigned char [numLocalCoords];
  env->localMemoryAssertion(__FILE__, __LINE__, numLocalCoords, newFlags);
  ArrayRCP<unsigned char> lrflags(newFlags, 0, numLocalCoords, true);

  try{
    BSPfindCut<mvector_t>( env, comm,
      params, numTestCuts, tolerance,
      cutDimension, coordDim, vectors, index,
      fractionLeft.view(0, weightDim), uniformWeights.view(0, weightDim),
      minCoord, maxCoord,
      cutValue, lrflags.view(0, numLocalCoords),
      imbalance);
  }
  Z2_FORWARD_EXCEPTIONS

  ///////////////////////////////////////////////////////
  // Adjust indices for left half and right half

  int leftCount = 0;
  for (int i=0; i < numLocalCoords; i++)
    if (lrflags[i] == leftFlag)
     leftCount++;

  int rightCount = numLocalCoords - leftCount;

  ArrayRCP<lno_t> leftIndices;
  ArrayRCP<lno_t> rightIndices;

  if (leftCount){
    lno_t *tmpLeft = new lno_t [leftCount];
    env->localMemoryAssertion(__FILE__, __LINE__, leftCount, tmpLeft);
    leftIndices = arcp(tmpLeft, 0, leftCount);
  }

  if (rightCount){
    lno_t *tmpRight = new lno_t [rightCount];
    env->localMemoryAssertion(__FILE__, __LINE__, rightCount, tmpRight);
    rightIndices = arcp(tmpRight, 0, rightCount);
  }

  int nextLeft = 0;
  int nextRight = 0;

  if (useIndices)
    for (int i=0; i < numLocalCoords; i++){
      if (lrflags[i] == leftFlag)
        leftIndices[nextLeft++] = index[i];
      else
        rightIndices[nextRight++] = index[i];
    }
  else{
    for (int i=0; i < numLocalCoords; i++){
      if (lrflags[i] == leftFlag)
        leftIndices[nextLeft++] = i;
      else
        rightIndices[nextRight++] = i;
    }
  }

  int leftPart0 = part0;
  int leftPart1 = part0 + numPartsLeftHalf - 1;
  int rightPart0 = leftPart1 + 1;
  int rightPart1 = part1;

  serialRCB(env, params, numTestCuts, tolerance, 
    coordDim, vectors, leftIndices.view(0, leftCount),
    uniformWeights.view(0, weightDim), partSizes.view(0, weightDim),
    leftPart0, leftPart1, partNum);

  serialRCB(env, params, numTestCuts, tolerance, 
    coordDim, vectors, rightIndices.view(0, rightCount),
    uniformWeights.view(0, weightDim), partSizes.view(0, weightDim),
    rightPart0, rightPart1, partNum);
}

/*! \brief Determine the fraction of work in the left half.
 *
 *   \param env the Environment for the application.
 *   \param part0  is the first part of the parts to bisected.
 *   \param part1  is the last part of the parts to bisected.
 *   \param partSizes is a list of part sizes for each weight dimension.
 *          Sizes should sum to 1.0 for each dimension.  If
 *          <tt>partSizes[w].size()</tt> is zero it means that
 *          part sizes for weight dimension \c w are all the same.
 *   \param fractionLeft on return <tt>fractionLeft[w]</tt> is the
 *             fraction of work wanted in the left half for weight
 *             dimension \c w.
 *   \param numPartsLeftHalf on return is the number of parts in the
 *               left half.
 */


template <typename scalar_t>
  void getFractionLeft(
    const RCP<const Environment> &env,
    partId_t part0,
    partId_t part1,
    const ArrayView<ArrayRCP<scalar_t> > partSizes,
    ArrayRCP<double> &fractionLeft,
    partId_t &numPartsLeftHalf)
{
  partId_t numParts = part1 - part0 + 1;
  numPartsLeftHalf = numParts / 2;
  partId_t left0 = part0;
  partId_t left1 = left0 + numPartsLeftHalf - 1;
  partId_t right0 = left1 + 1;
  partId_t right1 = part1;

  int weightDim = partSizes.size();
  fractionLeft = arcp(new double [weightDim], 0, weightDim);

  for (int wdim=0; wdim<weightDim; wdim++){
    if (partSizes[wdim].size() == 0){
      fractionLeft[wdim] = scalar_t(numPartsLeftHalf) / scalar_t(numParts);
    }
    else{
      fractionLeft[wdim] = 0;
      for(int partId=left0; partId <= left1; partId++){
        fractionLeft[wdim] += partSizes[wdim][partId];
      }
      double total = fractionLeft[wdim];
      for(int partId=right0; partId <= right1; partId++){
        total += partSizes[wdim][partId];
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
    const RCP<const Comm<int> > &comm,
    int coordDim,
    const RCP<mvector_t> &vectors,
    ArrayView<typename mvector_t::local_ordinal_type> index,
    int &dimension,                        // output
    typename mvector_t::scalar_type &minCoord,      // output
    typename mvector_t::scalar_type &maxCoord)      // output
{
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
        for (lno_t i=0; i < numLocalCoords; i++){
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
}

}  // namespace Zoltan2
#endif
