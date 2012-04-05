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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Vector.hpp>

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
  NUM_RCB_PARAMS
};

/*! \brief During partitioning flags are stored in unsigned char arrays.
 *  Flag is also used to store a region number, but there are at most
 *  251 regions.  Therefore region number will not conflict with leftFlag
 *  and rightFlag. (Number of regions is number of test cuts plus one.)
 */

enum leftRightFlag{
  leftFlag = 0xfe,     /*!< 254 */
  rightFlag = 0xff     /*!< 255 */
};

typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> mvector_t;

template <typename mvector_t>
  void getCutDimension( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm, 
    int coordDim, const RCP<mvector_t> &vectors, 
    ArrayView<mvector_t::local_ordinal_type> index,
    int &dimension, 
    mvector_t::scalar_type &minCoord, mvector_t::scalar_type &maxCoord);

template <typename mvector_t>
 void serialRCB<mvector_t(
    const RCP<const Environment> &env,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, mvector_t::scalar_type imbalanceTolerance, 
    int coordDim, const RCP<mvector_t> &vectors, 
    ArrayView<mvector_t::local_ordinal_type> index,
    const ArrayView<bool> uniformWeights,
    const typename ArrayView<ArrayRCP<mvector_t::scalar_type> > partSizes,
    partId_t part0, partId_t part1, ArrayView<partId_t> partNum);

template <typename mvector_t>
  void BSPfindCut( const RCP<const Environment> &env,
    const RCP<const Teuchos::Comm<int> > &comm, 
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, mvector_t::scalar_type imbalanceTolerance, int cutDim,
    int coordDim, const RCP<mvector_t> &vectors,
    ArrayView<mvector_t::local_ordinal_type> index,
    ArrayView<mvector_t::scalar_type> &fractionLeft, 
    ArrayView<bool> uniformWeights,
    mvector_t::scalar_type coordGlobalMin, 
    mvector::scalar_type coordGlobalMax, 
    mvector_t::scalar_type &cutValue, ArrayView<unsigned char> lrflags, 
    mvector_t::scalar_type &imbalance);

template <typename mvector_t>
  void determineCut( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, mvector_t::scalar_type imbalanceTolerance,
    int coordDim, const RCP<mvector_t> &vectors, 
    const ArrayView<bool> uniformWeights,
    const typename ArrayView<ArrayRCP<mvector_t::scalar_type> > partSizes,
    partId_t part0, partId_t part1,
    ArrayView<unsigned char> lrflags,
    int &cutDimension, 
    mvector_t::scalar_type &cutValue, mvector_t::scalar_type &imbalance);

template <typename mvector_t>
  void migrateData( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    int numParts, ArrayView<unsigned char> lrflags,
    const RCP<mvector_t> &vectors);

template <typename scalar_t>
  getFractionLeft( const RCP<const Environment> &env,
    partId_t part0, partId_t part1,
    const ArrayView<ArrayRCP<scalar_t> > partSizes,
    ArrayRCP<scalar_t> &fractionLeft);

/*! \brief Recursive coordinate bisection algorithm.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm  the communicator for the problem
 *  \param coords    a CoordinateModel with user data
 *  \param solution  a PartitioningSolution, on input it 
 *      contains part information, on return it also contains 
 *      the solution and quality metrics.
 *                    
 *   \todo timing and memory usage profiling
 *   \todo  catch errors and pass back
 *   \todo write the rcb tree back to the solution
 *   \todo for "repartition", start with the tree in the solution
 *   \todo for now we balance the first weight, so we need to add
 *             the multicriteria options as Zoltan1 does it.
 *   \todo incorporate part sizes as Zoltan1 does it.
 * \todo implement rectilinear_blocks and average_cuts
 *   \todo  work on performance issues.  Some of the global
 *               communication can probably be consolidated
 *                into fewer messages.
 */

template <typename Adapter>
void AlgRCB(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
  const RCP<const CoordinateModel<Adapter> > &coords, 
  RCP<PartitioningSolution<typename Adapter::user_t> > &solution
) 
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  int rank = env->myRank_;
  int nprocs = env->numProcs_;

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

  env->getValue<string>(env->getParameters(), "memory_versus_speed", 
    isSet, strChoice);

  if (isSet && strChoice==string("memory"))
    params.set(rcb_lowMemory);
  else if (isSet && strChoice==string("speed"))
    params.set(rcb_lowRunTime);
  else
    params.set(rcb_balanceMemoryRunTime);

  env->getValue<string>(env->getParameters(), "speed_versus_quality",
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

  env->getValue<string>(
     env->getList(env->getParameters(), "partitioning"), 
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

  env->getValue<double>(
     env->getList(env->getParameters(), "partitioning"), 
    "imbalance_tolerance", isSet, tol);

  if (!isSet)
    imbalanceTolerance = 1.1;
  else
    imbalanceTolerance = tol;

  ////////////////////////////////////////////////////////
  // Geometric partitioning problem parameters of interest:
  //    average_cuts
  //    rectilinear_blocks
  //    bisection_num_test_cuts (experimental)

  env->getValue<int>(
    env->getList(
      env->getList(env->getParameters(), "partitioning"), "geometric"),
    "average_cuts", isSet, intChoice);

  if (isSet && intChoice==1)
    params.set(rcb_averageCuts);

  env->getValue<int>(
    env->getList(
      env->getList(env->getParameters(), "partitioning"), "geometric"),
    "rectilinear_blocks", isSet, intChoice);

  if (isSet && intChoice==1)
    params.set(rcb_rectilinearBlocks);

  env->getValue<int>(
    env->getList(
      env->getList(env->getParameters(), "partitioning"), "geometric"),
    "bisection_num_cuts", isSet, intChoice);

  int numTestCuts = 3;
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

  Array<bool> uniformWeights(weightDim);
  Array<ArrayRCP<const scalar_t> > weights(weightDim);

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

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.
  // If the part sizes for a given criteria are not uniform, then
  // they are values that sum to 1.0.

  size_t numGlobalParts = solution->getGlobalNumberOfParts();

  Array<bool> uniformParts(weightDim);
  Array<ArrayRCP<scalar_t> > partSizes(weightDim);

  for (int wdim = 0; wdim < weightDim; wdim++){
    if (solution->criteriaHasUniformPartSizes(wdim)){
      uniformParts[wdim] = true;
    }
    else{
      uniformParts[wdim] = false;
      scalar_t *tmp = new scalar_t [numGlobalParts];
      env->localMemoryAssertion(__FILE__, __LINE__, numGlobalParts, tmp) ;
      for (partId_t i=0; i < numGlobalParts; i++){
        tmp[i] = solution->getCriteriaPartSize(i, wdim);
      }
      partSizes[wdim] = arcp(tmp, 0, numGlobalParts);
    }
  }

  bool multiplePartSizeSpecs = false;

  // TODO if weightDim > 1 figure out if user asked for
  //       different part size specifications for different
  //      weight dimensions. 

  ////////////////////////////////////////////////////////
  // Create the distributed data for the algorithm.
  //
  // It is a multivector containing one vector for each coordinate
  // dimension, plus a vector for each weight dimension that is not
  // uniform.

  typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> mvector_t;

  int multiVectorDim = coordDim;
  for (int wdim = 0; wdim < weightDim; wdim++)
    if (!uniformWeights[wdim]) multiVectorDim++;

  gno_t gnoMin, gnoMax;
  coords->getIdentifierMap()->getGnoRange(gnoMin, gnoMax);

  RCP<map_t> map;
  try{
    map = rcp(new map_t(numGlobalCoords, gnos, gnoMin, problemComm));
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  ArrayRCP<const ArrayView<const scalar_t> > vectors =
    arcp(new ArrayView<const scalar_t> [multiVectorDim], 0, multiVectorDim);

  for (int dim=0; dim < coordDim; dim++)
    vectors[dim] = values[dim].view(0, numLocalCoords);

  for (int wdim=0, idx=coordDim; wdim < weightDim; wdim++)
    if (!uniformWeights[wdim])
      vectors[idx++] = weights[wdim].view(0, numLocalCoords);

  RCP<mvector_t> mvector;

  try{
    mvector = rcp(new mvector_t(
      map, vectors.view(0, multiVectorDim), multiVectorDim));
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  ////////////////////////////////////////////////////////
  // The algorithm

  bool done=false;
  const RCP<const Comm<int> > &comm = problemComm;
  partId_t part0 = 0;
  partId_t part1 = numGlobalParts-1;
  int sanityCheck = numGlobalParts;

  if (comm->size() == 1)
    done = true;

  while (!done && sanityCheck--){

    // Determine which coordinates are left and which
    // are right.

    Array<unsigned char> lrflags(mvector->getLocalLength());
    scalar_t cutValue;  // TODO eventually save this for user
    int cutDimension;
    scalar_t imbalance;

    try{
      determineCut<mvector_t>(env, comm, 
        params, numTestCuts, imbalanceTolerance,
        coordDim, mvector, uniformWeights, partSizes, part0, part1,
        lrflags, cutDimension, cutValue, imbalance);
    }
    Z2_FORWARD_EXCEPTIONS

    
    // Migrate the multivector of data.

    int numParts = part1 - part0 + 1;
    int numLeftHalf = numParts / 2;
    int leftNumProcs;

    try{
      migrateData<mvector__t>( env, comm, 
        numParts, numLeftHalf, lrflags, 
        mvector,      // on return, mvector is the new data
        leftNumProcs) // on return, number of procs in left half
    }
    Z2_FORWARD_EXCEPTIONS

    // Divide into two subgroups

    const RCP<const Comm<int> > subComm;
    lno_t groupSize = 0;
    int *ids = NULL;

    if (comm->getRank()< leftNumProcs){
      groupSize = leftNumProcs;
      if (groupSize > 1) {
        ids = new int [groupSize];
        env->localMemoryAssertion(__FILE__, __LINE__, ids, groupSize);
        for (int i=0; i < groupSize; i++)
          ids[i] = i;
      }
      part1 = part0 + numLeftHalf - 1;
    }
    else{
      groupSize = comm->getSize() - leftNumProcs;
      if (groupSize > 1) {
        ids = new int [groupSize];
        env->localMemoryAssertion(__FILE__, __LINE__, ids, groupSize);
        for (int i=0; i < groupSize; i++)
          ids[i] = i + leftNumProcs;
      }

      part0 += numLeftHalf;
    }

    if ((part0 == part1) || (groupSize == 1)){
      done = true;
    }
    else {
      ArrayView<int> idView(ids, groupSize);
      subComm = comm->createSubcommunicator(idView);
      comm = subComm;
    }
  }

  env->localBugAssertion(__FILE__, __LINE__, "partitioning failure", 
    done, BASIC_ASSERTION);

  Array<partId_t> partId(mvector->getLocalLength());

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
        coordDim, mvector, emptyIndex, uniformWeights, partSizes,
        part0, part1, partId);
    }
    Z2_FORWARD_EXCEPTIONS
  }
  else{
    for (lno_t i=0; i < partId.size(); i++)
      partId[i] = part0;
  }

  // Now we have a part assignment.  Compute the imbalance and update the
  // solution.

}

/*! \brief Find the point in space that divides the data evenly with
 *     respect to the weights, part sizes, and the user's objective.
 *
 *   \param env the Environment for the application.
 *   \param comm the communicator for this step.
 *   \param params a bit map of boolean parameters.
 *   \param numTestCuts the number of test cuts to make in one round.
 *   \param imbalanceTolerance the maximum acceptable imbalance.
 *   \param cutDim  the dimension of the coordinates to cut.
 *   \param coordDim the first \c coordDim vectors in the \c vectors
 *              list are coordinates, the rest are weights.
 *   \param vectors lists of coordinates and non-uniform weights
 *   \param index is the index into the \c vectors arrays for the
 *              coordinates to be included in the partitioning.
 *   \param fractionLeft  the size of the left part for each weight
 *   \param uniformWeights element \c w is true if weights for weight
 *                 dimension \c w are all 1.0.
 *   \param coordGlobalMin the global minimum of coordinates in dimension
 *                                \c cutDim
 *   \param coordGlobalMax the global maximum of coordinates in dimension
 *                                \c cutDim
 *   \param cutValue  on return this is the computed cut location.
 *   \param lrflags on return has the value leftFlag or rightFlag to
 *        indicate whether the corresponding coordinate is on the left of
 *        on the right of the cut.  Allocated by caller.
 *   \param imbalance on return is the imbalance for the computed partitioning
 * 
 *  \todo need a serial method that creates n parts
 *  \todo the imbalance tolerance needs to be tighter because we are
 *          only making one cut in the series.
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
    mvector_t::scalar_type imbalanceTolerance,
    int cutDim,
    int coordDim,
    const RCP<mvector_t> &vectors,
    ArrayView<mvector_t::local_ordinal_type> index,
    ArrayView<mvector_t::scalar_type> fractionLeft,
    ArrayView<bool> uniformWeights,
    mvector_t::scalar_type coordGlobalMin,
    mvector_t::scalar_type coordGlobalMax,
    mvector_t::scalar_type &cutValue,         // output
    ArrayView<unsigned char> lrFlags,         // output
    mvector_t::scalar_type &imbalance)        // output
{
  typename typedef mvector_t::local_ordinal_type lno_t;
  typename typedef mvector_t::global_ordinal_type gno_t;
  typename typedef mvector_t::node_type node_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  bool useIndices = index.size() > 0;

  // Find the coordinate values and weights.

  int numCoords = 0;

  if (useIndices)
    numCoords = index.size();
  else 
    numCoords = vectors->getLocalLength();

  int weightDim = uniformWeights.size();
  int numNonUniformWeights = 0;
  for (int i=0; i < weightDim; i++){
    if (!uniformWeights[i])
      numNonUniformWeights++;
  }

  set<scalar_t> differentSizes;
  for (int i=0; i < weightDim; i++)
    set.insert(fractionLeft[i]);

  int numDifferentPartSizes = set.size(); // TODO: within some epsilon

  if (numDifferentPartSizes > 1 && numTestCuts < 3)
    numTestCuts = 3;

  const scalar_t *coordValue = vectors->getData(cutDim).getRawPtr();

  ArrayRCP<input_t> weight(new input_t [weightDim], 0, weightDim);

  for (int wdim = 0, widx=coordDim; wdim < weightDim; wdim++){
    if (!uniformWeights[wdim]){
      ArrayView<scalar_t> v(vectors->getData(widx++).getRawPtr(), numCoords);
      weight[wdim] = input_t(v, 1);
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
  
  // TODO adjust tolerance due to fact we are making cuts incrementally
  scalar_t tolerance = imbalanceTolerance - 1.0;
  if (tolerance < 0)
    tolerance = 1e-2;

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
  std::vector<scalar_t> testCuts(numBoundaries);
  std::vector<scalar_t>::iterator foundCut;

  memset(lrFlags.getRawPtr(), 0, numCoords); // 0 at top of loop means unset

  bool done=false;
  bool fail=false;
  scalar_t min = coordGlobalMin;
  scalar_t max = coordGlobalMax;
  lno_t numRemaining = numCoords;

  scalar_t totalWeight = 0, totalWeightLeft=0, totalWeightRight=0;
  scalar_t targetLeftScalar = 0;
  Epetra_SerialDenseVector targetLeftVector(weightDim);

  unsigned char partNumMin = 0; 
  unsigned char partNumMax = 0; 
  int firstNonLeftRegion = 0;

  while (!done && !fail){

    scalar_t diff = (max - min) / numRegions;
    testCuts[0] = min;
    for (int i=1; i <= numTestCuts; i++)
      testCuts[i+1] = testCuts[i] + diff;

    // Catch objects on the boundary
    testCuts[0] -= 1.0;
    testCuts[numBoundaries-1] = max + 1.0;

    if (numRemaining > 0){

      // Assign each of my points to a region.
      // lower_bound find sthe first cut f, such that f >= coordValue[i].
      // So for now, objects that are on the cut boundary go into the
      // region on the "left" side.

      if (useIndices){
        for (size_t i=0; i < numCoords; i++){
    
          if (lrFlags[i] != 0) 
            continue;
    
          if (numRegions > 2){
            foundCut = std::lower_bound(testCuts.begin(), testCuts.end(), 
                coordValue[index[i]]);
            
            env->localBugAssertion(__FILE__, __LINE__, "search cuts", 
              foundCut != testCuts.end(), BASIC_ASSERTION);
            
            lrFlags[i] = (unsigned char)(foundCut - testCuts.begin() - 1);
          }
          else{
            lrFlags[i] = (coordValue[index[i]] <= testCuts[1] ? 0 : 1);
          }
        }
      }
      else {
        for (size_t i=0; i < numCoords; i++){
    
          if (lrFlags[i] != 0)
            continue;
    
          if (numRegions > 2){
      
            foundCut = std::lower_bound(testCuts.begin(), testCuts.end(), 
                coordValue[i]);
            
            env->localBugAssertion(__FILE__, __LINE__, "search cuts", 
              foundCut != testCuts.end(), BASIC_ASSERTION);
            
            lrFlags[i] = (unsigned char)(foundCut - testCuts.begin() - 1);
          }
          else{
            lrFlags[i] = (coordValue[i] <= testCuts[1] ? 0 : 1);
          }
        }
      }
    }

    unsigned char numParts, numNonemptyParts;
    ArrayRCP<MetricValues<scalar_t> > metrics;
    ArrayRCP<scalar_t> weightSums;

    // Get the global sums in each region

    if (params.test(rcb_balanceCount) || numNonUniformWeights == 0){

      // Ignore weights and sum the objects in each part.

      Array<StridedInput<lno_t, scalar_t> > noWeights;

      globalSumsByPart<scalar_t, unsigned char, lno_t>(
        env, comm, 
        lrFlags,                      // part assignments
        partNumMin,
        partNumMax,
        noWeights.view(0, 1),
        mcNorm,                       // ignored
        numParts,
        numNonemptyParts,
        metrics,
        weightSums);    // numParts sums of objects in each part
    }
    else{
      globalSumsByPart<scalar_t, unsigned char, lno_t>(
        env, comm, 
        lrFlags,                      // part assignments
        partNumMin,
        partNumMax,
        weight.view(0, weightDim),
        mcNorm,
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
      for (int i=0; i < numParts, i++)
        totalWeight += values[i];

      totalWeightLeft = 0;
      targetLeftVector = partSizeLeft.Scale(totalWeight);
      targetLeftScalar = targetLeftVector[0];

      partNumMax = numParts-1;
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

        // |target-actual|^2 if cut is after region num

        for (int i=0; i < weightDim; i++)
          testVec[i] += values[num];
  
        diffVec = testVec.Scale(-1.0);
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

      if (minDiff <= tolerance){
        cutValue = testCuts[bestNum + 1];
        imbalance = 1.0 + minDiff;
        min = max = cutValue;
        firstNonLeftRegion = bestNum + 1;
        done = true;
      }
      else if ((bestNum >= 0) && (bestNum <= numParts-2)){
        // this only works if number of regions is at least three
        min = testCuts[bestNum];    // left of this region
        max = testCuts[bestNum+2];  // right of next region
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
        scalar_t imbalance = diffLeftCut / targetLeftScalar;
        if (imbalance <= tolerance){
          cutValue = testCuts[num];
          imbalance = 1.0 + diffLeftCut;
          min = max = cutValue;
          done = true;
          firstNonLeftRegion = num;
        }
      }
      else{
        scalar_t imbalance = diffRightCut / targetLeftScalar;
        if (imbalance <= tolerance){
          cutValue = testCuts[num+1];
          imbalance = 1.0 + diffRightCut;
          min = max = cutValue;
          done = true;
          firstNonLeftRegion = num+1;
        }
      }

      if (!done){
        min = testCuts[num];
        max = testCuts[num+1];
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
    fail==false, DEBUG_MODE_ASSERTION, comm);
}

/*! \brief Divide the coordinates into a "left" half and "right" half.
 *
 *   \param env the Environment for the application.
 *   \param comm the communicator for this step.
 *   \param params a bit map of boolean parameters.
 *   \param numTestCuts the number of test cuts to make in one round.
 *   \param imbalanceTolerance the maximum acceptable imbalance.
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
 *           for cutDim and fractionLeft.
 */

template <typename mvector_t>
  void determineCut(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const std::bitset<NUM_RCB_PARAMS> &params,
    int numTestCuts, mvector_t::scalar_type imbalanceTolerance,
    int coordDim, const RCP<mvector_t> &vectors,
    const ArrayView<bool> uniformWeights,
    const typename ArrayView<ArrayRCP<mvector_t::scalar_type> > partSizes,
    partId_t part0, 
    partId_t part1,
    ArrayView<unsigned char> lrflags,   // output
    int &cutDimension,                  // output
    mvector_t::scalar_type &cutValue,   // output
    mvector_t::scalar_type &imbalance)  // output
{
  typename typedef mvector_t::scalar_type scalar_t;
  typename typedef mvector_t::local_ordinal_type lno_t;
  typename typedef mvector_t::global_ordinal_type gno_t;
  typename typedef mvector_t::node_type node_t;

  ///////////////////////////////////////////////////////
  // Pick a cut direction.

  int cutDimension;
  scalar_t globalMinCoord, globalMaxCoord;
  ArrayView<lno_t> emptyIndex;

  getCutDimension<mvector_t>(env, comm, coordDim, mvector, emptyIndex,
    cutDimension, globalMinCoord, globalMaxCoord);

  ///////////////////////////////////////////////////////
  // Compute part sizes for the two parts.

  ArrayRCP<scalar_t> fractionLeft;

  getFractionLeft<scalar_t>(env, part0, part1, partSizes, fractionLeft);

  ///////////////////////////////////////////////////////
  // Divide the coordinates into balanced left and right
  // halves.

  int weightDim = uniformWeights.size();
  ArrayView<lno_t> emptyIndices;

  try{
    BSPfindCut<scalar_t, mvector_t>( env, comm,
      params, numTestCuts, imbalanceTolerance,
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
 *   \param numParts  the total number of parts represented.
 *   \param numLeftHalf  then number of parts in the new left half.
 *   \param lrflags on return has the value leftFlag or rightFlag to
 *        indicate whether the corresponding coordinate is on the left of
 *        on the right of the cut.  Allocated by caller.
 *   \param vectors is a list of the data to be migrated.  On return,
 *        \c vectors is new data belonging to this process.
 */
template <typename mvector_t>
  void migrateData(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    int numParts, int numLeftHalf,
    ArrayView<unsigned char> lrflags,
    const RCP<mvector_t> &vectors)    // on return is the new data
{
  typename typedef mvector_t::scalar_type scalar_t;
  typename typedef mvector_t::local_ordinal_type lno_t;
  typename typedef mvector_t::global_ordinal_type gno_t;
  typename typedef mvector_t::node_type node_t;

  int nprocs = comm->getSize();
  int rank = comm->getRank();
  int numVectors = vectors->getNumVectors();
  int numLocalCoords = vectors->getLocalLength();

  scalar_t leftFraction = scalar_t(numLeftHalf)/scalar_t(numParts);
  int leftNumProcs = nprocs * leftFraction;

  ///////////////////////////////////////////////////////
  // Get a list of my new global numbers.

  lno_t *sendCount = new lno_t [nprocs];
  env->localMemoryAssertion(__FILE__, __LINE__, nprocs, sendCount) ;
  memset(sendCount, 0, sizeof(int) * nprocs);
  ArrayView<lno_t> sendCountView(sendCount, nprocs);
  ArrayView<gno_t> sendBufView;

  if (numLocalCoords > 0){
    int *procId = new int [numLocalCoords];
    env->localMemoryAssertion(__FILE__, __LINE__, numLocalCoords, procId) ;
    int leftProc0 = 0;
    int rightProc0 = leftProc0 + leftNumProcs;
  
    int nextLeftProc = leftProc0;
    int nextRightProc = rightProc0;
    int *p = procId;
  
    for (lno_t i=0; i < numLocalCoords; i++){
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

    gno_t *sendBuf = new gno_t [numLocalCoords];
    env->localMemoryAssertion(__FILE__, __LINE__, numLocalCoords, sendBuf) ;
    sendBufView = ArrayView<gno_t>(sendBuf, numLocalCoords);

    ArrayView<const gno_t> gnoList = vectors->getMap()->getNodeElementList();

    for (lno_t i=0; i < numLocalCoords; i++){
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

  if (numLocalCoords > 0){
    delete [] sendBufView.getRawPtr();
    delete [] sendCountView.getRawPtr();
  }

  ///////////////////////////////////////////////////////
  // Migrate the multivector of data.

  lno_t numMyNewGnos = 0;
  for (int i=0; i < nprocs; i++)
    numMyNewGnos += recvCount;

  RCP<const mvector_t> newMultiVector;
  RCP<const mvector_t> constInput = rcp_const_cast<const mvector_t>(vectors);

  try{
    newMultiVector = XpetraTraits<mvector_t>::doMigration(
      constInput, numMyNewGnos, recvBuf.getRawPtr());
  }
  Z2_FORWARD_EXCEPTIONS

  vectors = newMultiVector;
}

/*! \brief Perform RCB on the local process only.
 *
 *   \param env the Environment for the application.
 *   \param params a bit map of boolean parameters.
 *   \param numTestCuts the number of test cuts to make in one round.
 *   \param imbalanceTolerance the maximum acceptable imbalance.
 *   \param coordDim the first \c coordDim vectors in the \c vectors
 *              list are coordinates, the rest are weights.
 *   \param vectors lists of coordinates and non-uniform weights
 *   \param index is the index into the \c vectors arrays for the
 *              coordinates to be included in the partitioning.
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
    mvector_t::scalar_type imbalanceTolerance, 
    int coordDim,
    const RCP<mvector_t> &vectors, 
    ArrayView<mvector_t::local_ordinal_type> index,
    const ArrayView<bool> uniformWeights,
    const typename ArrayView<ArrayRCP<mvector_t::scalar_type > > partSizes,
    partId_t part0, 
    partId_t part1,
    ArrayView<partId_t> partNum)   // output
{
  typename typedef mvector_t::scalar_type scalar_t;
  typename typedef mvector_t::local_ordinal_type lno_t;

  int numLocalCoords;
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

  ArrayRCP<scalar_t> fractionLeft;

  try{
    getFractionLeft<scalar_t>(env, part0, part1, partSizes, fractionLeft);
  }
  Z2_FORWARD_EXCEPTIONS

  ///////////////////////////////////////////////////////
  // Divide into balanced left and right halves.

  int weightDim = uniformWeights.size();
  int numLocalCoords = vectors->getLocalLength();
  RCP<const Comm<int> > comm = Teuchos::SerialComm();  

  scalar_t imbalance, cutValue;  //unused for now

  unsigned char *newFlags = new unsigned char [numLocalCoords];
  env->localMemoryAssertion(__FILE__, __LINE__, numLocalCoords, newFlags);
  ArrayRCP<unsigned char> lrflags(newFlags, 0, numLocalCoords);

  try{
    BSPfindCut<mvector_t>( env, comm,
      params, numTestCuts, imbalanceTolerance,
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

  lno_t *tmpLeft = new lno_t [leftCount];
  env->localMemoryAssertion(__FILE__, __LINE__, leftCount, tmpLeft);

  lno_t *tmpRight = new lno_t [rightCount];
  env->localMemoryAssertion(__FILE__, __LINE__, rightCount, tmpRight);

  ArrayRCP<lno_t> leftIndices(tmpLeft, 0, leftCount);
  ArrayRCP<lno_t> rightIndices(tmpRight, 0, rightCount);
  int nextLeft = nextRight = 0;

  for (int i=0; i < numLocalCoords; i++){
    if (lrflags[i] == leftFlag)
      leftIndices[nextLeft++] = index[i];
    else
      rightIndices[nextRight++] = index[i];
  }

  int numParts = part1 - part0 + 1;
  int numLeftHalf = numParts / 2;

  int leftPart0 = part0;
  int leftPart1 = part0 + numLeftHalf - 1;
  int rightPart0 = leftPart1 + 1;
  int rightPart1 = part1;

  serialRCB(env, params, numTestCuts, imbalanceTolerance, 
    coordDim, vectors, leftIndices.view(0, leftCount),
    uniformWeights, partSizes,
    leftPart0, leftPart1, partNum);

  serialRCB(env, params, numTestCuts, imbalanceTolerance, 
    coordDim, vectors, rightIndices.view(0, rightCount),
    uniformWeights, partSizes,
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
 */


template <typename scalar_t>
  getFractionLeft(
    const RCP<const Environment> &env,
    partId_t part0,
    partId_t part1,
    const ArrayView<ArrayRCP<scalar_t> > partSizes,
    ArrayRCP<scalar_t> &fractionLeft)
{
  partId_t numParts = part1 - part0 + 1;
  partId_t numLeftHalf = numParts / 2;
  partId_t numRightHalf = numParts - numLeftHalf;
  partId_t left0 = part0;
  partId_t left1 = left0 + numLeftHalf - 1;
  partId_t right0 = left1 + 1;
  partId_t right1 = part1;

  int weightDim = partSizes.size();
  fractionLeft = arcp(new scalar_t [weightDim], 0, weightDim);

  for (int wdim=0; wdim<weightDim; wdim++){
    if (partSizes[wdim].size() == 0){
      fractionLeft[wdim] = scalar_t(numLeftHalf) / scalar_t(numParts);
    }
    else{
      fractionLeft[wdim] = 0;
      for(int partId=left0; partId <= left1; partId++){
        fractionLeft[wdim] += partSizes[wdim][partId];
      }
      scalar_t total = fractionLeft[wdim];
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
    ArrayView<mvector_t::local_ordinal_type> index,
    int &dimension,                        // output
    mvector_t::scalar_type &minCoord,      // output
    mvector_t::scalar_type &maxCoord)      // output
{
  typename typedef mvector_t::scalar_type scalar_t;
  typename typedef mvector_t::local_ordinal_type lno_t;
  typename typedef mvector_t::global_ordinal_type gno_t;
  typename typedef mvector_t::node_type node_t;

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
      pair<scalar_t, scalar_t> minMax =
        z2LocalMinMax<scalar_t>(values[dim].getRawPtr(), numLocalCoords);
      spans[next++] = minMax.first();
      spans[next++] = minMax.second() * -1.0;
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
  int cutDimension = 0;
  int next = 0;

  for (int dim=0; dim < coordDim; dim++){
    min = span[next++];
    max = span[next++] * -1.0;
    scalar_t newSpan = max - min;
    if (newSpan > maxSpan){
      maxSpan = newSpan;
      cutDimension = dim;
      minCoord = min;
      maxCoord = max;
    }
  }
}

}  // namespace Zoltan2
#endif
