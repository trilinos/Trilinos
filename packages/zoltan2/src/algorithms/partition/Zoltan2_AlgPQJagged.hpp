// @HEADER
//***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_AlgRCB.hpp
\brief Contains the recursive coordinate bisection algorthm.
 */

#ifndef _ZOLTAN2_ALGPQJagged_HPP_
#define _ZOLTAN2_ALGPQJagged_HPP_

#include <Zoltan2_AlgRCB_methods.hpp>
#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_Metric.hpp>             // won't need this
#include <Zoltan2_GetParameter.hpp>

#include <Teuchos_ParameterList.hpp>

#define EPS_SCALE 5
#ifdef HAVE_ZOLTAN2_OMP
#include <omp.h>
#endif
#include <mpi.h>


//#define RCBCODE
#define mpi
#define LEAF_IMBALANCE_FACTOR 0.1

//imbalance calculation. Wreal / Wexpected - 1
#define imbalanceOf(Wachieved, totalW, expectedRatio) \
    (Wachieved) / ((totalW) * (expectedRatio)) - 1




namespace Zoltan2{

//diffclock for temporary timing experiments.
double diffclock(clock_t clock1,clock_t clock2)
{
  double diffticks=clock1-clock2;
  double diffms=(diffticks*10)/CLOCKS_PER_SEC;
  return diffms;
}



/*! \brief A helper class containing array representation of
 *  coordinate linked lists.
 */

template <typename lno_t, typename size_tt>
class pqJagged_PartVertices{
private:
  lno_t *linkedList; //initially filled with -1's.
  lno_t *partBegins; //initially filled with -1's.
  lno_t *partEnds; //initially filled with -1's.
public:

  //default constructor
  pqJagged_PartVertices(){};

  /*! \brief  The memory is provided to class via set function.
   *  \param linkedList_ is the array with size as the number of coordinates. Assumes all array is filled -1's.
   *  Each element of array points to the next element of the array in the linked list.
   *  \param partBegins_ is the array with size as the number of parts to be divided in current coordinate dimension. Assumes that array is filled with -1's.
   *  Holds the beginning of each part.
   *  \param partEnds_ is the array with size as the number of parts to be divided in current coordinate dimension. Assumes that array is filled with -1's.
   *  Holds the end coordinate of each part.
   */

  void set(lno_t *linkedList_, lno_t *partBegins_, lno_t *partEnds_){
    linkedList = linkedList_;
    partBegins = partBegins_;
    partEnds = partEnds_;
  }

  //user is responsible from providing the correct number of part counts
  /*! \brief Inserting a coordinate to a particular part.
   * Since, class does not have the size information,
   * it is user's responsibility to provide indices for partNo and coordinateIndex that are in the range.
   * \param partNo is the part number that the coordinate is inserted.
   * \param coordinateIndex is index of coordinate to be inserted.
   */
  void inserToPart (int partNo, lno_t coordinateIndex){

    switch (partEnds[partNo]){
    case -1: // this means partBegins[partNo] is also -1.
      partBegins[partNo] = coordinateIndex;
      partEnds[partNo] = coordinateIndex;
      break;
    default:
      linkedList[coordinateIndex] = partBegins[partNo];
      partBegins[partNo] = coordinateIndex;
      break;
    }


  }

  /*! \brief
   * linkedList getter function.
   */
  lno_t *getLinkedList(){ return linkedList;}

  /*! \brief
   * partBegins getter function.
   */
  lno_t *getPartBegins(){ return partBegins;}
  /*! \brief
   * partEnds getter function.
   */
  lno_t *getPartEnds(){ return partEnds;}

};





/*! \brief
 * Function that calculates the next pivot position,
 * according to given coordinates of upper bound and lower bound, the weights at upper and lower bounds, and the expected weight.
 * \param cutUpperBounds is the pointer to the array holding the upper bounds coordinates of the cuts.
 * \param cutLowerBounds is the pointer to the array holding the lower bound coordinates of the cuts.
 * \param currentCutIndex is the index of the current cut, the indices used for lower and upper bound arrays.
 * \param cutUpperWeight is the pointer to the array holding the weights at the upper bounds of each cut.
 * \param cutLowerWeight is the pointer to the array holding the weights at the lower bounds of each cut.
 * \param ew is the expected weight that should be placed on the left of the cut line.
 */
template <typename scalar_t>
inline scalar_t pivotPos (scalar_t * cutUpperBounds, scalar_t *cutLowerBounds,size_t currentCutIndex, scalar_t *cutUpperWeight, scalar_t *cutLowerWeight, scalar_t ew){

  return ((cutUpperBounds[currentCutIndex] - cutLowerBounds[currentCutIndex]) /
      (cutUpperWeight[currentCutIndex] - cutLowerWeight[currentCutIndex]))  * (ew - cutLowerWeight[currentCutIndex]) + cutLowerBounds[currentCutIndex];
}


/*! \brief  Returns the parameters such as:
 * Partitioning problem parameters of interest:
 *	Partitioning objective
 *	imbalance_tolerance
 *
 *Geometric partitioning problem parameters of interest:
 *	average_cuts
 *	rectilinear_blocks
 *	bisection_num_test_cuts (experimental)
 *	\param pl is the ParameterList object read from the environment.
 *	\param float-like value representing imbalance tolerance. An output of the function.
 *	(for imbalance 0.03, the user currently inputs 1.03. However, This function will return 0.03 for that case.)
 *	\param mcnorm multiCriteriaNorm. An output of the function. //not used in pqJagged algorithm currently.
 *	\param  params holding the bits for objective problem description. An output of the function. //not used in pqJagged algorithm currently..
 *	\param numTestCuts. An output of the function. //not used in pqJagged algorithm currently.
 *	\param ignoreWeights is the boolean value to treat the each coordinate as uniform regardless of the given input weights. Output of the function.
 */
template <typename T>
void pqJagged_getParameters(const Teuchos::ParameterList &pl, T &imbalanceTolerance,
    multiCriteriaNorm &mcnorm, std::bitset<NUM_RCB_PARAMS> &params,  int &numTestCuts, bool &ignoreWeights){

  bool isSet;
  string strChoice;
  int intChoice;


  getParameterValue<string>(pl, "partitioning",
      "objective", isSet, strChoice);

  if (isSet && strChoice == string("balance_object_count"))
    params.set(rcb_balanceCount);
  else if (isSet && strChoice ==
      string("multicriteria_minimize_total_weight")){
    params.set(rcb_minTotalWeight);
    mcnorm = normMinimizeTotalWeight;
  }
  else if (isSet && strChoice ==
      string("multicriteria_minimize_maximum_weight")){
    params.set(rcb_minMaximumWeight);
    mcnorm = normMinimizeMaximumWeight;
  }
  else if (isSet && strChoice ==
      string("multicriteria_balance_total_maximum")){
    params.set(rcb_balanceTotalMaximum);
    mcnorm = normBalanceTotalMaximum;
  }
  else{
    params.set(rcb_balanceWeight);
    mcnorm = normBalanceTotalMaximum;
  }

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
  if (isSet)
    numTestCuts = intChoice;
  ignoreWeights = params.test(rcb_balanceCount);


}


/*! \brief  Returns the input coordinate value parameters.
 * \param coords is the coordinate model representing the input.
 * \param coordDim is the output of the function representing the dimension count of the input.
 * \param weightDim is the output of the function representing the dimension of the weights.
 * \param numLocalCoords is the output representing the count of the local coordinates.
 * \param numGlobalCoords is the output representing the count of the global coordinates.
 * \param criteriaDim is the output representing the multi-objective count.
 * \param ignoreWeights is the boolean input of the function, to treat the each coordinate
 * as uniform regardless of the given input weights. Output of the function.
 *
 */
template <typename Adapter>
void pqJagged_getCoordinateValues( const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &coords, int &coordDim,
    int &weightDim, size_t &numLocalCoords, global_size_t &numGlobalCoords, int &criteriaDim, const bool &ignoreWeights){

  coordDim = coords->getCoordinateDim();
  weightDim = coords->getCoordinateWeightDim();
  numLocalCoords = coords->getLocalNumCoordinates();
  numGlobalCoords = coords->getGlobalNumCoordinates();
  criteriaDim = (weightDim ? weightDim : 1);
  if (criteriaDim > 1 && ignoreWeights)
    criteriaDim = 1;
}


/*! \brief  Function returning the input values for the problem such as the coordinates,
 * weights and desiredPartSizes.
 * \param env environment.
 * \param coords is the coordinate model representing the input.
 * \param solution is the partitioning solution object.
 * \param params is the  bitset to represent multiple parameters. //not used currently by pqJagged.
 * \param coordDim is an integer value to represent the count of coordinate dimensions.
 * \param weightDim is an integer value to represent the count of weight dimensions.
 * \param numLocalCoords is a size_t value to represent the count of local coordinates.
 * \param numGlobalParts is a size_t value to represent the total part count. //not used currently inside the pqJagged algorithm.
 * \param pqJagged_multiVectorDim  ...//not used by pqJagged algorithm.
 * \param pqJagged_values is the output representing the coordinates of local points.
 *  Its size is coordDim x numLocalCoords and allocated before the function.
 * \param criteriaDim ...//not used currently inside the pqJagged algorithm.
 * \param pqJagged_weights is the two dimensional array output, representing the weights of
 * coordinates in each weight dimension. Sized weightDim x numLocalCoords, and allocated before the function.
 * \param pqJagged_gnos is the ArrayView output representing the global indices of each vertex. No allocation is needed.
 * \param ignoreWeights is the boolean input of the function, to treat the each coordinate
 * \param pqJagged_uniformWeights is the boolean array representing whether or not coordinates have uniform weight in each weight dimension.
 * \param pqJagged_uniformParts is the boolean array representing whether or not the desired part weights are uniform in each weight dimension.
 * \param pqJagged_partSizes is the two dimensional float-like array output that represents the ratio of each part.
 *
 */
template <typename Adapter, typename scalar_t, typename gno_t>
void pqJagged_getInputValues(
    const RCP<const Environment> &env, const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &coords,
    RCP<PartitioningSolution<Adapter> > &solution,
    std::bitset<NUM_RCB_PARAMS> &params,
    const int &coordDim,
    const int &weightDim,
    const size_t &numLocalCoords, size_t &numGlobalParts, int &pqJagged_multiVectorDim,
    scalar_t **pqJagged_values, const int &criteriaDim, scalar_t **pqJagged_weights, ArrayView<const gno_t> &pqJagged_gnos, bool &ignoreWeights,
    bool *pqJagged_uniformWeights, bool *pqJagged_uniformParts, scalar_t **pqJagged_partSizes
){
  typedef typename Adapter::node_t node_t;
  typedef typename Adapter::lno_t lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  ArrayView<const gno_t> gnos;
  ArrayView<input_t>     xyz;
  ArrayView<input_t>     wgts;

  coords->getCoordinates(gnos, xyz, wgts);
  pqJagged_gnos = gnos;


  //std::cout << std::endl;

  for (int dim=0; dim < coordDim; dim++){
    ArrayRCP<const scalar_t> ar;
    xyz[dim].getInputArray(ar);
    //pqJagged coordinate values assignment
    pqJagged_values[dim] =  (scalar_t *)ar.getRawPtr();
  }

  if (weightDim == 0 || ignoreWeights){

    pqJagged_uniformWeights[0] = true;
  }
  else{
    for (int wdim = 0; wdim < weightDim; wdim++){
      if (wgts[wdim].size() == 0){
        pqJagged_uniformWeights[wdim] = true;
      }
      else{
        ArrayRCP<const scalar_t> ar;
        wgts[wdim].getInputArray(ar);
        pqJagged_uniformWeights[wdim] = false;
        pqJagged_weights[wdim] = (scalar_t *) ar.getRawPtr();
      }
    }
  }

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.
  // If the part sizes for a given criteria are not uniform,
  // then they are values that sum to 1.0.

  numGlobalParts = solution->getTargetGlobalNumberOfParts();

  for (int wdim = 0; wdim < criteriaDim; wdim++){
    if (solution->criteriaHasUniformPartSizes(wdim)){
      pqJagged_uniformParts[wdim] = true;
    }
    else{
      scalar_t *tmp = new scalar_t [numGlobalParts];
      env->localMemoryAssertion(__FILE__, __LINE__, numGlobalParts, tmp) ;
      for (size_t i=0; i < numGlobalParts; i++){
        tmp[i] = solution->getCriteriaPartSize(wdim, i);
      }
      pqJagged_partSizes[wdim] = tmp;

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

  pqJagged_multiVectorDim = coordDim;
  for (int wdim = 0; wdim < criteriaDim; wdim++){
    if (!pqJagged_uniformWeights[wdim]) pqJagged_multiVectorDim++;
  }

}


/*! \brief Printing the input values, to check configuration.
 *
 */
template <typename scalar_t, typename gno_t>
void pqJagged_printInput(int coordDim, int weightDim, size_t numLocalCoords, global_size_t numGlobalCoords,
    int criteriaDim, scalar_t **pqJagged_values, scalar_t **pqJagged_weights,
    bool *pqJagged_uniformParts, bool *pqJagged_uniformWeights, gno_t *pqJagged_gnos,
    bool ignoreWeights,size_t numGlobalParts, scalar_t **pqJagged_partSizes){

  std::cout << "numLocalCoords:" << numLocalCoords << std::endl;
  std::cout << "coordDim:" << coordDim << std::endl;
  for(int i = 0; i < numLocalCoords; ++i){
    for (int ii = 0; ii < coordDim; ++ii){
      std::cout <<  pqJagged_values[ii][i] << " ";
    }
    std::cout << std::endl;
  }


  std::cout << "criteriaDim:" << criteriaDim << std::endl;
  std::cout << "weightDim:" << weightDim << std::endl;
  if(weightDim){
    for(int i = 0; i < numLocalCoords; ++i){
      for (int ii = 0; ii < weightDim; ++ii){
        std::cout <<  pqJagged_weights[ii][i] << " ";
      }
      std::cout << std::endl;
    }
  }

  std::cout << "pqJagged_uniformWeights:" << pqJagged_uniformWeights[0] << std::endl;
  for(int i = 0; i < criteriaDim; ++i){
    std::cout << pqJagged_uniformWeights[i] << " ";
  }
  std::cout << std::endl;


  std::cout << "gnos" << std::endl;
  for(int i = 0; i < numLocalCoords; ++i){
    std::cout <<  pqJagged_gnos[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "ignoreWeights:" << ignoreWeights << std::endl;

  std::cout << "pqJagged_uniformParts:" << pqJagged_uniformParts[0] << std::endl;
  for(int i = 0; i < criteriaDim; ++i){
    std::cout << pqJagged_uniformParts[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "pqJagged_partSizes:" << std::endl;
  std::cout << "numGlobalParts:" << numGlobalParts << std::endl;
  for(int i = 0; i < criteriaDim; ++i){
    if(!pqJagged_uniformParts[i])
      for(int ii = 0; ii < numGlobalParts; ++ii){
        std::cout << pqJagged_partSizes[i][ii] << " ";
      }
    std::cout << std::endl;
  }
}


/*! \brief Function returning the available thread number by the processor.
 *
 */
int pqJagged_getNumThreads(){
  int numThreads = 1;


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(numThreads)
  {
    numThreads = omp_get_num_threads();
  }

#endif


  return numThreads;

}


/*! \brief Function to determine the minimum and maximum coordinate
 * in the given set of points.
 * \param comm is the Comm object. //Pure mpi all reduce is used, therefore, it is not used currently.
 * \param pqJagged_coordinates float-like array representing the coordinates in a single dimension. Sized as numLocalCoords.
 * \param minCoordinate is the output to represent the minimumCoordinate in  given range of coordinates.
 * \param maxCoordinate is the output to represent the maximum coordinate in the given range of coordinates.
 * \param env environment object. //not used currently.
 * \param numThreads is the integer value to represent the number of threads available for each processor.
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param coordinateBegin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinateEnd is the end index of the given partition on partitionedPointPermutations.
 * \param max_min_array provided work array sized numThreads * 2.
 */
template <typename scalar_t, typename lno_t>
void pqJagged_getMinMaxCoord(RCP<Comm<int> > &comm, scalar_t *pqJagged_coordinates, scalar_t &minCoordinate, scalar_t &maxCoordinate,
    const RCP<const Environment> &env, int numThreads,
    lno_t *partitionedPointPermutations, lno_t coordinateBegin, lno_t coordinateEnd, scalar_t *max_min_array /*sized nothreads * 2*/, scalar_t maxScalar, scalar_t minScalar){


  if(coordinateBegin >= coordinateEnd)
  {
    minCoordinate = maxScalar;
    maxCoordinate = minScalar;
  }
  else {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
    {

      int myId = 0;
#ifdef HAVE_ZOLTAN2_OMP
      myId = omp_get_thread_num();
#endif
      scalar_t myMin, myMax;

      //size_t j = firstIteration ? 0 : partitionedPointCoordinates[coordinateBegin];
      size_t j = partitionedPointPermutations[coordinateBegin];
      myMin = myMax = pqJagged_coordinates[j];
      //cout << "coordinateBegin:" << coordinateBegin << " end:" << coordinateEnd << endl;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(lno_t j = coordinateBegin + 1; j < coordinateEnd; ++j){
        //int i = firstIteration ? j:partitionedPointCoordinates[coordinateBegin];
        int i = partitionedPointPermutations[j];
        //cout << pqJagged_coordinates[i] << endl;
        if(pqJagged_coordinates[i] > myMax) myMax = pqJagged_coordinates[i];
        if(pqJagged_coordinates[i] < myMin) myMin = pqJagged_coordinates[i];
      }
      max_min_array[myId] = myMin;
      max_min_array[myId + numThreads] = myMax;
      //cout << "myMin:" << myMin << " myMax:" << myMax << endl;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#endif



#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single nowait
#endif
      {
        minCoordinate = max_min_array[0];
        for(int i = 1; i < numThreads; ++i){
          if(max_min_array[i] < minCoordinate) minCoordinate = max_min_array[i];
        }
      }
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single nowait
#endif
      {
        maxCoordinate = max_min_array[numThreads];
        for(int i = numThreads + 1; i < numThreads * 2; ++i){
          if(max_min_array[i] > maxCoordinate) maxCoordinate = max_min_array[i];
        }
      }

    }
  }
  //cout <<" mx:" << maxCoordinate << " mn:" << minCoordinate << endl;

#ifdef mpi
  scalar_t minm = 0;scalar_t maxm = 0;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
  {

    //TODO:Merge these communications.
    try{
      reduceAll<int,scalar_t>(
          *comm, Teuchos::REDUCE_MIN,
          1, &minCoordinate, &minm
      );

      reduceAll<int,scalar_t>(
          *comm, Teuchos::REDUCE_MAX,
          1, &maxCoordinate, &maxm
      );
    }
    Z2_THROW_OUTSIDE_ERROR(*env)

    //MPI_Allreduce ( &minCoordinate, &minm, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    //MPI_Allreduce ( &maxCoordinate, &maxm, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  }
  //cout << "minm:" << minm << " maxm:" << maxm << endl;

  //cout << "sizeof scalar_t:" << sizeof(scalar_t)  << endl;
  minCoordinate = minm;
  maxCoordinate = maxm;
#endif
}


template <typename scalar_t>
void pqJagged_getCutCoord_Weights(
    scalar_t minCoordinate, scalar_t maxCoordinate,
    bool pqJagged_uniformParts, scalar_t *pqJagged_partSizes /*p sized, weight ratios of each part*/,
    size_t noCuts/*p-1*/ ,
    scalar_t *cutCoordinates /*p - 1 sized, coordinate of each cut line*/,
    scalar_t *cutPartRatios /*cumulative weight ratios, at left side of each cut line. p-1 sized*/){

  scalar_t coordinateRange = maxCoordinate - minCoordinate;
  if(pqJagged_uniformParts){
    scalar_t uniform = 1. / (noCuts + 1);
    scalar_t slice = uniform * coordinateRange;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < noCuts; ++i){
      cutPartRatios[i] =  uniform * (i + 1);
      cutCoordinates[i] = minCoordinate + slice * (i + 1);
    }
    cutPartRatios[noCuts] = 1;
  }
  else {
    cutPartRatios[0] = pqJagged_partSizes[0];
    cutCoordinates[0] = coordinateRange * cutPartRatios[0];
    for(size_t i = 1; i < noCuts; ++i){
      cutPartRatios[i] = pqJagged_partSizes[i] + cutPartRatios[i - 1];
      cutCoordinates[i] = coordinateRange * cutPartRatios[i];
    }
  }
}


/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 * \param total_part_count is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
 * \param totalPartWeights is the array holding the weight of parts. Assumes there are 2*P - 1 parts.
 * \param isDone is the boolean array to determine if the correct position for a cut line is found.
 * \param cutPartRatios are the desired cumulative part ratios, sized P.
 * \param totalWeight is the global total weight in the current range of coordinates.
 * \param imbalanceTolerance is the maximum allowed imbalance ratio.
 * \param allDone is the number of cut lines that are not completed.
 * \param cutUpperBounds is the array holding the upper bound coordinate for each cut line. Sized P - 1.
 * \param cutLowerBounds is the array holding the lower bound coordinate for each cut line. Sized P - 1.
 * \param cutCoordinates is the array holding the coordinates of each cut line. Sized P - 1.
 * \param noCuts is the number of cut lines. P - 1.
 * \param maxCoordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
 * \param minCoordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
 * \param leftClosestDistance is the array holding the distances to the closest points to the cut lines from left.
 * \param rightClosestDistance is the array holding the distances to the closest points to the cut lines from right.
 * \param cutLowerWeight is the array holding the weight of the parts at the left of lower bound coordinates.
 * \param cutUpperWeight is the array holding the weight of the parts at the left of upper bound coordinates.
 * \param cutCoordinatesWork is the work array, sized P - 1.
 */
template <typename scalar_t>
void getNewCoordinates_simple(
    const size_t &total_part_count, const scalar_t * totalPartWeights,
    bool *isDone, const scalar_t *cutPartRatios,
    const scalar_t &totalWeight, const scalar_t &imbalanceTolerance,
    size_t &allDone, scalar_t *cutUpperBounds, scalar_t *cutLowerBounds,
    scalar_t *&cutCoordinates, const size_t &noCuts,
    const scalar_t &maxCoordinate, const scalar_t &minCoordinate,
    scalar_t *leftClosestDistance, scalar_t *rightClosestDistance,
    scalar_t * cutLowerWeight,scalar_t * cutUpperWeight, scalar_t *&cutCoordinatesWork,
    float *nonRectelinearPart,
    bool allowNonRectelinearPart){

  scalar_t seenW = 0;
  float expected = 0;
  scalar_t leftImbalance = 0, rightImbalance = 0;

  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();

  /*
for (size_t i = 0; i < noCuts; i++){
cout << "i:" << i << " pw:" << totalPartWeights[i * 2] << endl;
}
   */

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
  for (size_t i = 0; i < noCuts; i++){


    //if a left and right closes point is not found, set the distance to 0.
    if(leftClosestDistance[i] < 0) leftClosestDistance[i] = 0;
    if(rightClosestDistance[i] < 0) rightClosestDistance[i] = 0;


    //current weight of the part at the left of the cut line.
    seenW = totalPartWeights[i * 2];

    //cout << "i:" << i << " w: "<< seenW << endl;

    //if already determined at previous iterations, do nothing.
    if(isDone[i]) {
      cutCoordinatesWork[i] = cutCoordinates[i];
      continue;
    }
    expected = cutPartRatios[i];
    /*
if(i == 2){
cout << "seen: " << seenW << " expected:" << totalWeight * expected << " cutLine:" << totalPartWeights[i * 2 + 1] << endl;
cout << "upper:" << cutUpperBounds[i] << " lower:" <<  cutLowerBounds[i] << endl;
}
     */
    leftImbalance = imbalanceOf(seenW, totalWeight, expected);
    rightImbalance = imbalanceOf(totalWeight - seenW, totalWeight, 1 - expected);

    bool isLeftValid = leftImbalance <= imbalanceTolerance && leftImbalance >= -imbalanceTolerance;
    bool isRightValid = rightImbalance <= imbalanceTolerance && rightImbalance >= -imbalanceTolerance;


    //cout << "leftImbalance: " << leftImbalance << imbalanceTolerance << endl;
    //cout << "rightImbalance:" << rightImbalance << imbalanceTolerance <<  endl;
    //if the cut line reaches to desired imbalance.
    if(isLeftValid && isRightValid){
      isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
      allDone -=1;

      cutCoordinatesWork [ i] = cutCoordinates[i];
    }
    //if upper bound and lower bound reaches to the same point.
    //the imbalance cannot be satisfied. Accept as it is.
    /*else if(cutUpperBounds[i] != -1 && (cutUpperBounds[i] - cutLowerBounds[i]) < _EPSILON ){
isDone[i] = true;
__sync_fetch_and_sub(&allDone, 1);
cutCoordinatesWork [ i] = cutCoordinates[i];
}
     */
    //if left imbalance is lower than 0, then cut line should be moved to right.

    else if(leftImbalance < 0){

      scalar_t ew = totalWeight * expected;

      //cout << "allow:" << allowNonRectelinearPart << " " << totalPartWeights[i * 2 + 1] << " " << ew << endl;
      if(allowNonRectelinearPart){
        if (totalPartWeights[i * 2 + 1] >= ew){
          isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          allDone -=1;
          cutCoordinatesWork [ i] = cutCoordinates[i];
          //cout << "left:" << totalPartWeights[i * 2] << " with line:" << totalPartWeights[i * 2 + 1 ]<< " expected:" << ew << endl;

          nonRectelinearPart[i] = 1 - (totalPartWeights[i * 2 + 1] - ew) / (totalPartWeights[i * 2 + 1] - totalPartWeights[i * 2]);
          //cout << "nonRectelinearPart[i]:" << nonRectelinearPart[i] << endl;
          continue;
        }
      }
      //TODO if boundary is enough binary search on processors -- only if different part assignment in the same line is allowed.

      //TODO for now cutUpperBounds and cutLowerBounds are initialized with -1, which is not okay. Think about a way without using another assignment array.

      //if it is the first iteration.
      if (cutUpperBounds[i] == -1){
        bool done = false;
        //compare it with other cut lines.
        for (size_t ii = i + 1; ii < noCuts ; ++ii){
          scalar_t pw = totalPartWeights[ii * 2];

          if(pw >= ew){
            //if pw are equal, that coordinate is the desired solution.
            if(pw == ew){
              isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
              allDone -=1;
              cutCoordinatesWork [i] = cutCoordinates[ii];
            }
            else {

              //if the part weight at the left of a cut line is bigger than the expected weight,
              //choose it as the upper bound.
              cutUpperBounds[i] = cutCoordinates[ii] - leftClosestDistance[ii]; //epsilon to shift the cut line. currently no shift.
              //cutUpperWeight[i] = totalPartWeights [2 * ii + 1];
              cutUpperWeight[i] = totalPartWeights [2 * ii];

              //lower bound becomes the line before the upper bound.
              cutLowerBounds[i] = cutCoordinates [ ii -1] + rightClosestDistance[ii - 1];
              //cutLowerWeight[i] = totalPartWeights [2 * ii - 1];
              cutLowerWeight[i] = totalPartWeights [2 * ii - 2];


              //calculate new position.
              scalar_t newPivot = pivotPos<scalar_t> (cutUpperBounds, cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew);

              //if the new coordinate is same as previous.---floating point operations requires this check.
              //choose the coordinate.
              if (fabs(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE){
                isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
                allDone -=1;
                cutCoordinatesWork [i] = cutCoordinates[i];
              } else {
                cutCoordinatesWork [i] = newPivot;
              }
            }
            done = true;
            break;
          }
        }
        //if not found, upper bound will be max coordinate.
        //lower bound will be the last cut line.
        if(!done){
          cutUpperBounds[i] = maxCoordinate;
          cutUpperWeight[i] = totalWeight;

          cutLowerBounds[i] = cutCoordinates [ noCuts -1] + rightClosestDistance[noCuts - 1];
          //cutLowerWeight[i] = totalPartWeights [2 * noCuts - 1];
          cutLowerWeight[i] = totalPartWeights [2 * noCuts - 2];

          scalar_t newPivot = pivotPos<scalar_t> (cutUpperBounds, cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew);
          if (fabs(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE){
            isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
            allDone -=1;
            cutCoordinatesWork [ i] = cutCoordinates[i];
          } else {
            cutCoordinatesWork [ i] = newPivot;
          }
        }

      } else {

        //when moving left, and upper and lower are set.
        //set lower bound to current line.
        cutLowerBounds[i] = cutCoordinates[i] + rightClosestDistance[i];
        cutLowerWeight[i] = seenW;

        //compare the upper bound with the current lines.
        for (size_t ii = i + 1; ii < noCuts ; ++ii){
          scalar_t pw = totalPartWeights[ii * 2];

          if(pw >= ew){
            if(pw == ew){
              cutCoordinatesWork[i] = cutCoordinates[ii];
              isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
              allDone -=1;
            } else if (pw < cutUpperWeight[i]){
              //if a cut line is more strict than the current upper bound,
              //update the upper bound.
              cutUpperBounds[i] = cutCoordinates[ii] - leftClosestDistance[ii];
              //cutUpperWeight[i] = totalPartWeights [2 * ii + 1];
              cutUpperWeight[i] = totalPartWeights [2 * ii ];
            }
            break;
          }
          //if a stricter lower bound is found,
          //update the lower bound.
          if (pw <= ew && pw >= cutLowerWeight[i]){
            cutLowerBounds[i] = cutCoordinates[ii] + rightClosestDistance[ii] ;
            cutLowerWeight[i] = pw;
          }
        }


        scalar_t newPivot = pivotPos<scalar_t> (cutUpperBounds, cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew);
        if (fabs(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE){
          isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          allDone -=1;
          cutCoordinatesWork [ i] = cutCoordinates[i];
        } else {
          cutCoordinatesWork [ i] = newPivot;
        }
      }

    } else {

      //moving to left.
      scalar_t ew = totalWeight * expected;
      //if first iteration.
      if (cutLowerBounds[i] == -1){

        bool done = false;

        //set the upper and lower bounds according to the part weights of the cut lines on the left side.
        for (int ii = i - 1; ii >= 0; --ii){
          scalar_t pw = totalPartWeights[ii * 2];
          if(pw <= ew){
            if(pw == ew){
              cutCoordinatesWork[i] = cutCoordinates[ii];
              isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
              allDone -=1;
            } else {
              cutUpperBounds[i] = cutCoordinates[ii + 1] - leftClosestDistance[ii + 1];
              //cutUpperWeight[i] = totalPartWeights [2 * ii + 3];
              cutUpperWeight[i] = totalPartWeights [2 * ii + 2];


              cutLowerBounds[i] = cutCoordinates[ii] + rightClosestDistance[ii];
              //cutLowerWeight[i] = totalPartWeights [2 * ii + 1];
              cutLowerWeight[i] = totalPartWeights [2 * ii];

              scalar_t newPivot = pivotPos<scalar_t> (cutUpperBounds, cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew);
              /*
if(i == 2) {
cout << "pw:" << pw << " ew:" << ew << endl;
cout << "new:" << newPivot << " old:" << cutCoordinates[i] << endl;
cout <<"((" << cutUpperBounds[i] << "-" <<  cutLowerBounds[i] << " ) /" <<
"(" << cutUpperWeight[i] << " -" <<  cutLowerWeight[i] << " ))  * ( " << ew << "-" << cutLowerWeight[i] << ") + " << cutLowerBounds[i] << endl;

}
               */
              if (fabs(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE){
                isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
                allDone -=1;
                cutCoordinatesWork [ i] = cutCoordinates[i];
              } else {
                cutCoordinatesWork [ i] = newPivot;
              }

            }
            done = true;
            break;
          }

        }

        //if not found a part weight lower than the expected one,
        //set the lower bound to min coordinate.
        if(!done){

          cutUpperBounds[i] = cutCoordinates[0] - leftClosestDistance[0];
          cutLowerBounds[i] = minCoordinate;

          cutUpperWeight[i] = totalPartWeights [1];
          cutLowerWeight[i] = 0;

          scalar_t newPivot = pivotPos<scalar_t> (cutUpperBounds, cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew);
          if (fabs(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE){
            isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
            allDone -=1;
            cutCoordinatesWork [ i] = cutCoordinates[i];
          } else {
            cutCoordinatesWork [ i] = newPivot;
          }

        }
      } else {
        //moving left, both upper and lower is set.
        //set upper to current line.
        cutUpperBounds[i] = cutCoordinates[i] - leftClosestDistance[i];
        cutUpperWeight[i] = seenW;

        // compare the current cut line weights with previous upper and lower bounds.
        for (int ii = i - 1; ii >= 0; --ii){
          scalar_t pw = totalPartWeights[ii * 2];

          if(pw <= ew){
            if(pw == ew){
              cutCoordinatesWork[i] = cutCoordinates[ii];
              isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
              allDone -=1;


            } else if (pw > cutLowerWeight[i]){
              cutLowerBounds[i] = cutCoordinates[ii] + rightClosestDistance[ii];
              //cutLowerWeight[i] = totalPartWeights [2 * ii + 1];
              cutLowerWeight[i] = totalPartWeights [2 * ii ];
            }
            break;
          }
          if (pw >= ew && pw < cutUpperWeight[i]){
            cutUpperBounds[i] = cutCoordinates[ii] - leftClosestDistance[ii] ;
            cutUpperWeight[i] = pw;
          }

        }

        scalar_t newPivot = pivotPos<scalar_t> (cutUpperBounds, cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew);
        if (fabs(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE){
          isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          allDone -=1;
          cutCoordinatesWork [ i] = cutCoordinates[i];
        } else {
          cutCoordinatesWork [ i] = newPivot;
        }
      }

    }
  }

  //swap the work with cut coordinates.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#pragma omp single
#endif
  {
    scalar_t *t = cutCoordinatesWork;
    cutCoordinatesWork = cutCoordinates;
    cutCoordinates = t;
  }
}

/*! \brief Function that is responsible from 1D partitioning of the given range of coordinates.
 * \param pqJagged_coordinates is 1 dimensional array holding coordinate values.
 * \param pqJagged_weights is 1 dimensional array holding the weights of points.
 * \param pqJagged_uniformWeights is a boolean value if the points have uniform weights.
 * \param numLocalCoords is the number of local coordinates.
 * \param numGlobalCoords is the number of global coordinates.
 * \param minCoordinate is the minimum coordinate in the given dimension and in the given range of points.
 * \param maxCoordinate is the maximum coordinate in the given dimension and in the given range of points.
 * \param pqJagged_partSizes is the array holding the cumulative weight ratios of each part. Size P.
 * \param partNo is the total number of parts.
 * \param noThreads is the number of available threads by each processor.
 * \param imbalanceTolerance is the maximum allowed imbalance ratio.
 * \param cutCoordinates is the array holding the coordinates of the cut.
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param coordinateBegin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinateEnd is the end index of the given partition on partitionedPointPermutations.
 * \param cutCoordinatesWork is a work array sized P - 1.
 * \param cutPartRatios is the cumulative desired weight ratios for each part. Sized P.
 * \param cutUpperBounds is the array holding the upper bound coordinate for each cut line. Sized P - 1.
 * \param cutLowerBounds is the array holding the lower bound coordinate for each cut line. Sized P - 1.
 * \param cutLowerWeight is the array holding the weight of the parts at the left of lower bound coordinates.
 * \param cutUpperWeight is the array holding the weight of the parts at the left of upper bound coordinates.
 * \param isDone is the boolean array to determine if the correct position for a cut line is found.
 * \param partWeights is the two dimensional array holding the part weights. One dimension for each thread.
 * \param leftClosestDistance is the two dimensional array holding the distances to the closest points to the cut lines from left. One dimension for each thread.
 * \param rightClosestDistance is the two dimensional array holding the distances to the closest points to the cut lines from right. One dimension for each thread.
 * \param totalPartWeights is one dimensional array holding the part weights. Sum of all threads.
 */
template <typename scalar_t, typename lno_t>
void pqJagged_1DPart_simple(const RCP<const Environment> &env,RCP<Comm<int> > &comm,
    scalar_t *pqJagged_coordinates,	scalar_t *pqJagged_weights,	bool pqJagged_uniformWeights,
    const size_t &numLocalCoords, const global_size_t &numGlobalCoords,	scalar_t &minCoordinate, scalar_t &maxCoordinate,
    scalar_t *pqJagged_partSizes, size_t partNo, int noThreads,
    scalar_t imbalanceTolerance, scalar_t *cutCoordinates, lno_t *partitionedPointPermutations,
    lno_t coordinateBegin, lno_t coordinateEnd,

    scalar_t *cutCoordinatesWork, 	// work array to manipulate coordinate of cutlines in different iterations.
    scalar_t *cutPartRatios,		// the weight ratios at left side of the cuts. last is 1.

    scalar_t *cutUpperBounds,  //to determine the next cut line with binary search
    scalar_t *cutLowerBounds,  //to determine the next cut line with binary search

    scalar_t *cutLowerWeight,  //to determine the next cut line with binary search
    scalar_t *cutUpperWeight,   //to determine the next cut line with binary search
    bool *isDone,
    scalar_t **partWeights,
    scalar_t **leftClosestDistance,
    scalar_t **rightClosestDistance,
    scalar_t *totalPartWeights,
    scalar_t &globaltotalWeight,
    float *nonRectelinearPart,
    bool allowNonRectelinearPart,
    scalar_t maxScalar,
    scalar_t minScalar
){


  scalar_t *cutCoordinates_tmp = cutCoordinates;
  size_t noCuts = partNo - 1;



  scalar_t totalWeight = 0;
  size_t total_part_count = partNo + noCuts;
  size_t allDone = noCuts;
  globaltotalWeight = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(allDone, globaltotalWeight)
#endif
  {
    //int iterationCount = 0;
    //calculate total weight
    if (pqJagged_uniformWeights) {
      totalWeight = coordinateEnd - coordinateBegin;
    } else {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for reduction(+:totalWeight)
#endif
      for (size_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
        int i = partitionedPointPermutations[ii];
        totalWeight += pqJagged_weights[i];
      }
    }

    globaltotalWeight = totalWeight;
#ifdef mpi

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#pragma omp single
#endif
    {
      scalar_t tw = 0;
      //TODO: data type.

      try{
        //MPI_Allreduce ( &totalWeight, &tw, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        reduceAll<int,scalar_t>(
            *comm, Teuchos::REDUCE_SUM,
            1, &totalWeight, &tw
        );

      }
      Z2_THROW_OUTSIDE_ERROR(*env)


      globaltotalWeight = tw;
    }

#endif

    if(noCuts == 0){
      totalPartWeights[0] = globaltotalWeight;
    }
    else {


      int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
      me = omp_get_thread_num();
#endif
      scalar_t *myPartWeights = partWeights[me];
      scalar_t *myLeftClosest = leftClosestDistance[me];
      scalar_t *myRightClosest = rightClosestDistance[me];

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(size_t i = 0; i < noCuts; ++i){
        isDone[i] = false;
        cutLowerBounds[i] = -1;
        cutUpperBounds[i] = -1;
        if(allowNonRectelinearPart){
          nonRectelinearPart[i] = 0;
        }
      }


      while (allDone != 0){

        //iterationCount++;

/*
{
cout << endl << endl << "allDone:" << allDone << endl;
for (size_t i = 0; i < noCuts; ++i){
if(isDone[i] == false)
cout << "i:" << i <<  " c:" << cutCoordinates_tmp[i] << " u:" << cutUpperBounds[i] << " l:" << cutLowerBounds[i] << " not done" << endl;
else
cout << "i:" << i <<  " c:" << cutCoordinates_tmp[i] <<  " done" << endl;
}
}
*/

        for (size_t i = 0; i < total_part_count; ++i){
          if(i/2 < noCuts && isDone[i/2]) continue;
          myPartWeights[i] = 0;
        }
        for(size_t i = 0; i < noCuts; ++i){
          if(isDone[i]) continue;
          myLeftClosest[i] = -1;
          myRightClosest[i] = -1;
        }


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for (size_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
          int i = partitionedPointPermutations[ii];

          //get a coordinate and compare it with cut lines from left to right.
          for(size_t j = 0; j < noCuts; ++j){

            if(isDone[j]) continue;
            scalar_t distance = pqJagged_coordinates[i] - cutCoordinates_tmp[j];

            //if it is on the left
            if (distance < 0) {
              distance = -distance;
              if (myLeftClosest[j] < 0 || myLeftClosest[j] > distance){
                myLeftClosest[j] = distance;
              }
              break;
            }


            else if (distance == 0){
              scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
              myPartWeights[j * 2] -=	w;
              myLeftClosest[j] = 0;
              myRightClosest[j] = 0;
              //TODO: check if the next cut coordinates are in same.
              //break;
            } else {
              //if on the right, continue with the next line.
              scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
              myPartWeights[j * 2] -=	w;
              myPartWeights[j * 2 + 1] -=	w;
              if (myRightClosest[j] < 0 || myRightClosest[j] > distance){
                myRightClosest[j] = distance;
              }
            }
          }
        }

        //TODO: if the number of cutlines are small, run it sequential.
        //accumulate left and right closest distances for different threads.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for(size_t i = 0; i < noCuts; ++i){
          if(isDone[i]) continue;
          scalar_t minl = leftClosestDistance[0][i], minr = rightClosestDistance[0][i];

          for (int j = 1; j < noThreads; ++j){
            if ((rightClosestDistance[j][i] < minr || minr < 0) && rightClosestDistance[j][i] >= 0 ){
              minr = rightClosestDistance[j][i];
            }
            if ((leftClosestDistance[j][i] < minl || minl < 0) && leftClosestDistance[j][i] >= 0){
              minl = leftClosestDistance[j][i];
            }
          }
          leftClosestDistance[0][i] = minl;
          rightClosestDistance[0][i] = minr;

          //cout << "left:i:"<< i << ":" << leftClosestDistance[0][i] << " right:" << rightClosestDistance[0][i] <<endl;
        }

        //TODO: if the number of cutlines are small, run it sequential.
        //accumulate part weights for different threads.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for(size_t j = 0; j < total_part_count; ++j){
          if(noCuts != 0 && j/2 < noCuts && isDone[j/2]) continue;
          scalar_t pwj = 0;
          for (int i = 0; i < noThreads; ++i){
            pwj += partWeights[i][j];
          }
          totalPartWeights[j] = pwj + totalWeight;

        }


        //all to all partweight send;
        //get totalpartweights from different nodes sized total_part_count
        //reduce part sizes here within all mpi processes.
#ifdef mpi

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#pragma omp single
#endif
        {

          //TODO: get these as parameter work arrays
          scalar_t *tw = new scalar_t [total_part_count];
          scalar_t *lc = new scalar_t [noCuts];
          scalar_t *rc = new scalar_t [noCuts];

          //if the no point on the left or right of the cut line in the current processor.
          for(size_t ii = 0; ii < noCuts; ++ii){
            if(leftClosestDistance[0][ii] == -1){
              leftClosestDistance[0][ii] = maxScalar;
            }

            if(rightClosestDistance[0][ii] == -1){
              rightClosestDistance[0][ii] = maxScalar;
            }
          }

          //TODO data types. and data operations.
          try{
            //MPI_Allreduce ( totalPartWeights, tw, total_part_count, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            reduceAll<int,scalar_t>(
                *comm, Teuchos::REDUCE_SUM,
                total_part_count, totalPartWeights, tw
            );
            //MPI_Allreduce ( leftClosestDistance[0], lc, noCuts, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            reduceAll<int,scalar_t>(
                *comm, Teuchos::REDUCE_MIN,
                noCuts, leftClosestDistance[0], lc
            );

            //MPI_Allreduce ( rightClosestDistance[0], rc, noCuts ,MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

            reduceAll<int,scalar_t>(
                *comm, Teuchos::REDUCE_MIN,
                noCuts, rightClosestDistance[0], rc
            );
          }
          Z2_THROW_OUTSIDE_ERROR(*env)




          for(size_t ii = 0; ii < total_part_count; ++ii){
            totalPartWeights[ii] = tw[ii];

          }
          for(size_t ii = 0; ii < noCuts; ++ii){
            //reverse operation.
            //these cases should not occur though.
            if(lc[ii] == maxScalar){
              lc[ii] = -1;
            }
            if(lc[ii] == minScalar){
              lc[ii] = -1;
            }
            //totalPartWeights[ii] = tw[ii];
            leftClosestDistance[0][ii] = lc[ii];
            rightClosestDistance[0][ii] = rc[ii];
          }
          delete []tw;
          delete []lc;
          delete [] rc;
        }
#endif

        {
          /*
for (size_t i = 0; i < noCuts+1; ++i){
cout << "pw:" << i << " is:" << totalPartWeights[i*2] << endl;
}
           */
        }


        //get the new cut coordinates.
        getNewCoordinates_simple<scalar_t>(total_part_count, totalPartWeights, isDone, cutPartRatios,
            globaltotalWeight, imbalanceTolerance, allDone, cutUpperBounds, cutLowerBounds,
            cutCoordinates_tmp, noCuts,maxCoordinate, minCoordinate, leftClosestDistance[0],
            rightClosestDistance[0],cutLowerWeight, cutUpperWeight,cutCoordinatesWork,

            nonRectelinearPart,
            allowNonRectelinearPart);
      }


      //we cannot swap the arrays
      if (cutCoordinates != cutCoordinates_tmp){

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for(size_t i = 0; i < noCuts; ++i){
          cutCoordinates[i] = cutCoordinates_tmp[i];
        }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
        {
          cutCoordinatesWork = cutCoordinates_tmp;
        }

      }
      /*
#pragma omp single
{
cout << "Iteration:" << iterationCount << endl;
}
       */
    }
  }
}



/*! \brief Function that determines the permutation indices of the coordinates.
 * \param partNo is the number of parts.
 * \param numCoordsInPart is the number of points in the current part.
 * \param noThreads is the number of threads avaiable for each processor.
 * \param pqJagged_coordinates is 1 dimensional array holding the coordinate values.
 * \param cutCoordinates is 1 dimensional array holding the cut coordinates.
 * \param totalCounts are the number points in each output part.
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param newpartitionedPointPermutations is the indices of coordinates calculated for the partition on next dimension.
 * \param coordinateBegin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinateEnd is the end index of the given partition on partitionedPointPermutations.
 * \param coordinate_linked_list is the array with size as the number of coordinates.
 * \param coordinate_starts is the two dimensional array with size as the number of parts to be divided in current coordinate dimension. 1 dimension for each thread.
 * \param coordinate_ends is the two dimensional array with size as the number of parts to be divided in current coordinate dimension. 1 dimension for each thread.
 * \param numLocalCoord is the number of local coordinates.
 */
template <typename lno_t, typename scalar_t>
void getChunksFromCoordinates(size_t partNo, size_t numCoordsInPart, int noThreads,
    scalar_t *pqJagged_coordinates, scalar_t *cutCoordinates, lno_t *totalCounts,
    lno_t *partitionedPointPermutations,
    lno_t *newpartitionedPointPermutations, lno_t coordinateBegin, lno_t coordinateEnd,
    lno_t *coordinate_linked_list, lno_t **coordinate_starts, lno_t **coordinate_ends, size_t numLocalCoord, float *actual_ratios, bool allowNonRectelinearPart,
    scalar_t *totalPartWeights, scalar_t *coordWeights, bool pqJagged_uniformWeights, int myRank, int worldSize, scalar_t **partWeights, float **nonRectelinearRatios){

  ++myRank;
  //TODO: move this allocation to outside. Do not free at each dimension.
  lno_t **partPointCounts = new lno_t *[noThreads];
  size_t noCuts = partNo - 1;

  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();
  for(int i = 0; i < noThreads; ++i){
    if(i == 0){
      partPointCounts[i] = totalCounts;
    } else {
      partPointCounts[i] = new lno_t[partNo];
    }
  }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
  {
    int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
    me = omp_get_thread_num();
    //cout << "me:" << me << endl;
#endif

    lno_t *myStarts = coordinate_starts[me];
    lno_t *myEnds = coordinate_ends[me];
    lno_t *myPartPointCounts = partPointCounts[me];
    scalar_t *myPartWeights = partWeights[me];
    float *myRatios = NULL;
    if (allowNonRectelinearPart){
      myRatios = nonRectelinearRatios[me];
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for (size_t i = 0; i < noCuts; ++i){
        float r = actual_ratios[i];
        //if(i > 0 && fabs(cutCoordinates[i] - cutCoordinates[i -1]) < _EPSILON){
        //  r -= actual_ratios[i - 1];
        //}
        scalar_t leftWeight = r * (totalPartWeights[i * 2 + 1] - totalPartWeights[i * 2]);
        for(int ii = 0; ii < noThreads; ++ii){
          if(leftWeight > _EPSILON){
            scalar_t ithWeight = partWeights[ii][i * 2 + 1] - partWeights[ii][i * 2 ];
            if(ithWeight < leftWeight){
              //nonRectelinearRatios[ii][i] = 1;
              nonRectelinearRatios[ii][i] = ithWeight;
            }
            else {

              //nonRectelinearRatios[ii][i] = leftWeight / ithWeight;
              nonRectelinearRatios[ii][i] = leftWeight ;
            }
            cout << "l:" << leftWeight << " ith:" << ithWeight << endl;
            leftWeight -= ithWeight;
          }
          else {
            nonRectelinearRatios[ii][i] = 0;
          }
        }
      }


      if(noCuts > 0){
        for (size_t i = noCuts - 1; i > 0 ; --i){

          cout << "i:"<< i<<" r:" << myRatios[i] << " actual:" << actual_ratios[i]<<endl;
          if(fabs(cutCoordinates[i] - cutCoordinates[i -1]) < _EPSILON){
            myRatios[i] -= myRatios[i - 1] ;
            //myRatios[i] /= partWeights[me][i * 2 + 1] - partWeights[me][i * 2 ];
          }

          cout << "ni:" << myRatios[i] << endl;
        }
      }

    }


    pqJagged_PartVertices <lno_t, size_t> pqPV;
    pqPV.set(coordinate_linked_list, myStarts, myEnds);
    memset(myPartPointCounts, 0, sizeof(lno_t) * partNo);


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (size_t i = 0; i < numLocalCoord; ++i){
      coordinate_linked_list[i] = -1;
    }

    for (size_t i = 0; i < partNo; ++i){
      myEnds[i] = -1;
      myStarts[i] = -1;
    }




    /*
for(size_t j = 0; j < noCuts; ++j){
cout << "nonRectelinearPart[" << j << "]:" << nonRectelinearPart[j] << endl;
}
     */

    //determine part of each point

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

      lno_t i = partitionedPointPermutations[ii];

      //cout << "i:" << i << " me:" << me << " pqJagged_coordinates[i]:" << pqJagged_coordinates[i] << " cutCoordinates[j]:" << cutCoordinates[0]<<  endl;
      bool inserted = false;
      for(size_t j = 0; j < noCuts; ++j){
        if (pqJagged_coordinates[i] < cutCoordinates[j]){
          pqPV.inserToPart(j, i);
          //cout << "i:" << i << " c:" << pqJagged_coordinates[i] << " is placed to:" << j << endl;
          inserted = true;
          break;
        }
        else if (allowNonRectelinearPart && myRatios[j] > _EPSILON * EPS_SCALE && fabs(pqJagged_coordinates[i] - cutCoordinates[j]) < _EPSILON ){

          scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];
          //cout << "was:" << nonRectelinearPart[j] << " ";
          //TODO RACE CONDITION HERE!!
          //scalar_t decrease = w / (totalPartWeights[j * 2 + 1] - totalPartWeights[j * 2]) ;
          //scalar_t decrease = w / (myPartWeights[j * 2 + 1] - myPartWeights[j * 2]) ;
          scalar_t decrease = w;

          //cout << "nonRectelinearPart[j]:" << myRatios[j] << " actual_ratios[j]:" << actual_ratios[j] << " me:" << me << " world:" << worldSize << endl;
          //cout << " myPartWeights[j * 2 + 1]:" << myPartWeights[j * 2 + 1] << " myPartWeights[j * 2 ]:" << myPartWeights[j * 2 ] << endl;
          //cout << "decrease:" << decrease << endl;
//#pragma omp atomic
          myRatios[j] -= decrease;
          if(myRatios[j] < 0){
//#pragma omp atomic
            myRatios[j + 1] += myRatios[j];
          }
          inserted = true;
          pqPV.inserToPart(j, i);
          //cout << "i:" << i << " c:" << pqJagged_coordinates[i] << " is placed to:" << j << endl;
          break;
        }
        else {
          myPartPointCounts[j] -=	 1;
        }

      }
      if(!inserted){
        pqPV.inserToPart(noCuts, i);
        //cout << "i:" << i << " c:" << pqJagged_coordinates[i] << " is placed to:" << noCuts << endl;
      }
    }

    //accumulate the starts of each thread.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(size_t i = 0; i < partNo; ++i){
      int head = 0, tail = 0;
      while(head < noThreads && coordinate_ends[head][i] == -1){
        ++head;
      }
      int firstHead = head;
      for(int j = head; j < noThreads; ){
        tail = j+1;
        bool foundTail = false;
        while(tail < noThreads){
          if(coordinate_starts[tail][i] == -1){
            ++tail;
          }
          else {
            foundTail = true;
            break;
          }
        }
        if(foundTail){
          coordinate_linked_list[coordinate_ends[head][i]] = coordinate_starts[tail][i];
          head = tail;
        }
        j = tail;

      }
      coordinate_starts[0][i] = firstHead >= noThreads ? -1:coordinate_starts[firstHead][i];
    }

    //accumulate the count.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(size_t j = 0; j < partNo; ++j){
      scalar_t pwj = 0;
      for (int i = 0; i < noThreads; ++i){
        pwj += partPointCounts[i][j];
      }
      totalCounts[j] = pwj + numCoordsInPart;
    }


    //write new indices

    //cout << "\n\npart no:" << partNo << endl;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(size_t i = 0; i < partNo; ++i){
      //cout << "starting part:" << i << endl;
      lno_t nextPoint = coordinate_starts[0][i];
      lno_t pcnt = 0;

      lno_t prevCount = coordinateBegin;
      if (i > 0) prevCount = totalCounts[i -1] + coordinateBegin;

      while(nextPoint != -1){
        //cout << "next:" << nextPoint << " is in Part:" << i << endl;
        newpartitionedPointPermutations[prevCount + pcnt++] = nextPoint;
        nextPoint = coordinate_linked_list[nextPoint];
      }
    }
  }

  for(int i = 0; i < noThreads; ++i){
    if(i != 0){
      delete [] partPointCounts[i];
    }
  }

  delete []partPointCounts;
}


/*! \brief PQJagged coordinate partitioning algorithm.
 *
 *  \param env   library configuration and problem parameters
 *  \param comm the communicator for the problem
 *  \param coords    a CoordinateModel with user data
 *  \param solution  a PartitioningSolution, on input it
 *      contains part information, on return it also contains
 *      the solution and quality metrics.
 */

template <typename Adapter>
void AlgPQJagged(
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,
    const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &coords,
    RCP<PartitioningSolution<Adapter> > &solution
)
{
  cout << "start" << endl;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;

  typedef typename Adapter::lno_t lno_t;
  const Teuchos::ParameterList &pl = env->getParameters();

  std::bitset<NUM_RCB_PARAMS> params;
  int numTestCuts = 5;
  scalar_t imbalanceTolerance;

  multiCriteriaNorm mcnorm;
  bool ignoreWeights;
  pqJagged_getParameters<scalar_t>(pl, imbalanceTolerance, mcnorm, params, numTestCuts, ignoreWeights);

  int coordDim, weightDim; size_t numLocalCoords; global_size_t numGlobalCoords; int criteriaDim;
  pqJagged_getCoordinateValues<Adapter>( coords, coordDim, weightDim, numLocalCoords, numGlobalCoords, criteriaDim, ignoreWeights);



  scalar_t **pqJagged_coordinates = new scalar_t *[coordDim];
  scalar_t **pqJagged_weights = new scalar_t *[criteriaDim];
  bool *pqJagged_uniformParts = new bool[criteriaDim];
  scalar_t **pqJagged_partSizes = new scalar_t*[criteriaDim];
  bool *pqJagged_uniformWeights = new bool[criteriaDim];

  ArrayView<const gno_t> pqJagged_gnos;
  size_t numGlobalParts;
  int pqJagged_multiVectorDim;


  pqJagged_getInputValues<Adapter, scalar_t, gno_t>(
      env, coords, solution,params,coordDim,weightDim,numLocalCoords,
      numGlobalParts, pqJagged_multiVectorDim,
      pqJagged_coordinates,criteriaDim, pqJagged_weights,pqJagged_gnos, ignoreWeights,
      pqJagged_uniformWeights, pqJagged_uniformParts, pqJagged_partSizes
  );



  /*
pqJagged_printInput<scalar_t, gno_t>(coordDim, weightDim, numLocalCoords, numGlobalCoords,
criteriaDim, pqJagged_coordinates, pqJagged_weights,
pqJagged_uniformParts, pqJagged_uniformWeights, pqJagged_gnos,
ignoreWeights,numGlobalParts, pqJagged_partSizes);
   */

  double start = 0;
#ifdef HAVE_ZOLTAN2_OMP
  start = omp_get_wtime( );
#endif
  int numThreads = pqJagged_getNumThreads();

  if(comm->getRank() == 0){
    cout << "numGlobalParts:" << numGlobalParts << endl;
    cout << "numGlobalCoords:" << numGlobalCoords << endl;
    cout << "numThreads=" << numThreads << endl;
    cout << "numLocalCoord:" << numLocalCoords << endl;
  }
  scalar_t minCoordinate, maxCoordinate;


  size_t totalDimensionCut = 0;
  size_t totalPartCount = 1;
  size_t maxPartNo = 0;

  //TODO: part numbers are given as integers,
  //although we assume they are size_t
  const int *partNo = pl.getPtr<Array <int> >("pqParts")->getRawPtr();
  for (int i = 0; i < coordDim; ++i){
    totalPartCount *= partNo[i];
    if(partNo[i] > maxPartNo) maxPartNo = partNo[i];
  }
  totalDimensionCut = totalPartCount - 1;

  scalar_t *allCutCoordinates = new scalar_t[totalDimensionCut]; // coordinates of the cut lines. First one is the min, last one is max coordinate.
  lno_t *partitionedPointCoordinates = new lno_t [numLocalCoords];
  lno_t *newpartitionedPointCoordinates = new lno_t [numLocalCoords];
  scalar_t *max_min_array = new scalar_t[numThreads * 2];


  //initial configuration
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
  for(lno_t i = 0; i < numLocalCoords; ++i){
    //set each pointer-i to i.
    partitionedPointCoordinates[i] = i;
  }

  //initially there is a single partition
  lno_t currentPartitionCount = 1;

  //single partition starts at index-0, and ends at numLocalCoords
  lno_t coordinateBegin = 0;
  lno_t coordinateEnd = numLocalCoords;

  //inTotalCounts array holds the end points in partitionedPointCoordinates array
  //for each partition. Initially sized 1, and single element is set to numLocalCoords.
  lno_t *inTotalCounts = new lno_t [1];
  inTotalCounts[0] = numLocalCoords;

  //the ends points of the output.
  lno_t *outTotalCounts = NULL;

  //the array holding the points in each part as linked list
  lno_t *coordinate_linked_list = new lno_t[numLocalCoords];
  //the start and end coordinate of  each part.
  lno_t **coordinate_starts = new lno_t *[numThreads];
  lno_t **coordinate_ends = new lno_t *[numThreads];

  //assign the max size to starts, as it will be reused.
  for(int i = 0; i < numThreads; ++i){
    coordinate_starts[i] = new lno_t[maxPartNo];
    coordinate_ends[i] = new lno_t[maxPartNo];
  }

  size_t maxCutNo = maxPartNo - 1;

  bool allowNonRectelinearPart = true;
  float *nonRectelinearPart = NULL;
  float **nonRectRatios = NULL;
  if(allowNonRectelinearPart){
    nonRectelinearPart = new float[maxCutNo];
    nonRectRatios = new float *[numThreads];
    for(int i = 0; i < numThreads; ++i){
      nonRectRatios[i] = new float[maxCutNo];
    }
  }

  scalar_t *cutCoordinatesWork = new scalar_t [maxCutNo]; // work array to manipulate coordinate of cutlines in different iterations.
  scalar_t *cutPartRatios = new scalar_t[maxPartNo]; // the weight ratios at left side of the cuts. First is 0, last is 1.

  scalar_t *cutUpperBounds = new scalar_t [maxCutNo];  //to determine the next cut line with binary search
  scalar_t *cutLowerBounds = new scalar_t [maxCutNo];  //to determine the next cut line with binary search

  scalar_t *cutLowerWeight = new scalar_t [maxCutNo];  //to determine the next cut line with binary search
  scalar_t *cutUpperWeight = new scalar_t [maxCutNo];  //to determine the next cut line with binary search

  //scalar_t *imbalance_tolerances = new scalar_t [maxPartNo];
  //scalar_t *used_imbalance_tolerances = new scalar_t [maxPartNo];
  //imbalance_tolerances[0] = imbalanceTolerance;

  bool *isDone = new bool [maxCutNo];
  scalar_t **partWeights = new scalar_t *[numThreads];

  scalar_t **leftClosestDistance = new scalar_t* [numThreads];
  scalar_t **rightClosestDistance = new scalar_t* [numThreads];

  size_t maxTotalPartCount = maxPartNo + maxCutNo;
  for(int i = 0; i < numThreads; ++i){
    partWeights[i] = new scalar_t[maxTotalPartCount];
    rightClosestDistance[i] = new scalar_t [maxCutNo];
    leftClosestDistance[i] = new scalar_t [maxCutNo];
  }

  scalar_t *totalPartWeights = new scalar_t[maxTotalPartCount];
  scalar_t *cutCoordinates =  allCutCoordinates;


  size_t leftPartitions = totalPartCount;

  scalar_t maxScalar_t = numeric_limits<float>::max();
  scalar_t minScalar_t = -numeric_limits<float>::max();

  for (int i = 0; i < coordDim; ++i){

    lno_t partitionCoordinateBegin = 0;
    outTotalCounts = new lno_t[currentPartitionCount * partNo[i]];

    size_t currentOut = 0;
    size_t currentIn = 0;
    size_t previousEnd = 0;

    //cout << "i:" << i << endl;

    leftPartitions /= partNo[i];
    /*
int readImbalanceIndex = 0;
int writeImbalanceIndex = 0;
     */
    for (int j = 0; j < currentPartitionCount; ++j, ++currentIn){
      //depending on wheter it is a leaf or intermediate node, used imbalance is adapted.
      /*
if(i == coordDim - 1){
used_imbalance_tolerances[readImbalanceIndex] = imbalance_tolerances[readImbalanceIndex];
} else {
used_imbalance_tolerances[readImbalanceIndex] = imbalance_tolerances[readImbalanceIndex] * (LEAF_IMBALANCE_FACTOR  * (1.0 / leftPartitions) + LEAF_IMBALANCE_FACTOR );
}
       */

      //			cout << "imbalance_tolerances[readImbalanceIndex]:" << imbalance_tolerances[readImbalanceIndex] << endl;
      //used_imbalance_tolerances[readImbalanceIndex] = imbalance_tolerances[readImbalanceIndex] * (LEAF_IMBALANCE_FACTOR   / leftPartitions)/* + LEAF_IMBALANCE_FACTOR*/ ;

      scalar_t used_imbalance = imbalanceTolerance * (LEAF_IMBALANCE_FACTOR + (1 - LEAF_IMBALANCE_FACTOR)   / leftPartitions) * 0.7;

      //cout << "used:" << used_imbalance << endl;

      coordinateEnd= inTotalCounts[currentIn];
      coordinateBegin = currentIn==0 ? 0: inTotalCounts[currentIn -1];


      //cout << "myRank:" << comm->getRank()<< " i:" << i << " beg:" << coordinateBegin << " end:" << coordinateEnd << endl;
      pqJagged_getMinMaxCoord<scalar_t, lno_t>(comm, pqJagged_coordinates[i], minCoordinate,maxCoordinate,
          env, numThreads, partitionedPointCoordinates, coordinateBegin, coordinateEnd, max_min_array, maxScalar_t, minScalar_t);

      if(minCoordinate <= maxCoordinate){
        //cout << "myRank:" << comm->getRank()<< " i:" << i << " beg:" << coordinateBegin << " end:" << coordinateEnd  << " min:" << minCoordinate << " max:" << maxCoordinate << endl;
        pqJagged_getCutCoord_Weights<scalar_t>(
            minCoordinate, maxCoordinate,
            pqJagged_uniformParts[0], pqJagged_partSizes[0], partNo[i] - 1,
            cutCoordinates, cutPartRatios
        );

        /*
for (size_t ii = 0;ii < partNo[i] - 1; ++ii){
cout << "i:"<< i<<" initial: " << cutCoordinates[ii]<<  " max:" << maxCoordinate << " min:"<< minCoordinate<< endl;
}
         */


        scalar_t globalTotalWeight = 0;
        pqJagged_1DPart_simple<scalar_t, lno_t>(env,comm,pqJagged_coordinates[i], pqJagged_weights[0], pqJagged_uniformWeights[0],
            numLocalCoords, numGlobalCoords,	minCoordinate,
            maxCoordinate, pqJagged_partSizes[0], partNo[i], numThreads,
            used_imbalance, cutCoordinates,
            partitionedPointCoordinates, coordinateBegin, coordinateEnd,
            cutCoordinatesWork,
            cutPartRatios,
            cutUpperBounds,
            cutLowerBounds,
            cutLowerWeight,
            cutUpperWeight,
            isDone,
            partWeights,
            leftClosestDistance,
            rightClosestDistance,
            totalPartWeights,
            globalTotalWeight,
            nonRectelinearPart,
            allowNonRectelinearPart,
            maxScalar_t,
            minScalar_t);
      }
      /*
cout << "total:" << totalPartWeights[0] << endl;

for (size_t ii = 0;ii < partNo[i] - 1; ++ii){
cout << "nN: " << nonRectelinearPart[ii]<< endl;
}
       */
      //cout << cutCoordinates[0] << endl;
      getChunksFromCoordinates<lno_t,scalar_t>(partNo[i], coordinateEnd - coordinateBegin, numThreads,
          pqJagged_coordinates[i], cutCoordinates, outTotalCounts + currentOut,
          partitionedPointCoordinates, newpartitionedPointCoordinates, coordinateBegin, coordinateEnd,
          coordinate_linked_list, coordinate_starts, coordinate_ends, numLocalCoords,
          nonRectelinearPart, allowNonRectelinearPart, totalPartWeights, pqJagged_weights[0], pqJagged_uniformWeights, comm->getRank(), comm->getSize(), partWeights,nonRectRatios);



      cutCoordinates += partNo[i] - 1;
      partitionCoordinateBegin += coordinateBegin - coordinateEnd;

      for (size_t ii = 0;ii < partNo[i]; ++ii){
        outTotalCounts[currentOut+ii] += previousEnd;
      }

      //cout << "global:" << globalTotalWeight << endl;
      //Imbalance calculation
      /*
if(i != coordDim - 1){
for(int iii = 0; iii < partNo[i]; ++iii){
//cout << "i:"<<iii << " w:" << totalPartWeights[iii * 2] <<  " e:"<< cutPartRatios[iii]<< endl;
float expectedR = 0;
scalar_t achievedW = 0;
if(iii == 0) {
achievedW =  totalPartWeights[0];
expectedR = cutPartRatios[0];
} else {
expectedR = cutPartRatios[iii] - cutPartRatios[iii - 1];
achievedW =  totalPartWeights[iii * 2] -  totalPartWeights[iii * 2  - 2];
}
//scalar_t expectedW = globalTotalWeight * expectedR;
scalar_t futureImbalance = 0;

futureImbalance = used_imbalance * (expectedR * globalTotalWeight) / achievedW;
imbalance_tolerances[writeImbalanceIndex++] = futureImbalance;
//cout << "future imb:" << futureImbalance << endl;
}
}
       */
      previousEnd = outTotalCounts[currentOut + partNo[i] - 1];
      currentOut += partNo[i];

    }
    /*
scalar_t *t_imb = imbalance_tolerances;
imbalance_tolerances = used_imbalance_tolerances;
used_imbalance_tolerances = t_imb;
     */
    lno_t * tmp = newpartitionedPointCoordinates;
    newpartitionedPointCoordinates = partitionedPointCoordinates;
    partitionedPointCoordinates = tmp;

    currentPartitionCount *= partNo[i];
    delete [] inTotalCounts;
    inTotalCounts = outTotalCounts;
    //outTotalCounts = NULL;

  }


  lno_t * tot = new lno_t[totalPartCount];
  partId_t *partIds = new partId_t[numLocalCoords];

  ArrayRCP<partId_t> partId = arcp(partIds, 0, numLocalCoords, true);

  for(size_t i = 0; i < totalPartCount;++i){
    lno_t begin = 0;
    lno_t end = inTotalCounts[i];
    if(i > 0) begin = inTotalCounts[i -1];
    //if(i < totalPartCount - 1)
    //cout << "part:" << i << " cut:" << allCutCoordinates[i] << endl;

    //cout << "part:i - " << end - begin << endl;
    for (lno_t ii = begin; ii < end; ++ii){

      lno_t k = partitionedPointCoordinates[ii];
      partIds[k] = i;
      cout << "part of coordinate:";
      for(int iii = 0; iii < coordDim; ++iii){
        cout <<  pqJagged_coordinates[iii][k] << " ";
      }
      cout << i;
      cout << endl;
    }

    //cout << endl;
  }

  try{

    //MPI_Allreduce(inTotalCounts, tot, totalPartCount, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    reduceAll<int,lno_t>(
        *comm, Teuchos::REDUCE_SUM,
        totalPartCount, inTotalCounts, tot
    );

  }
  Z2_THROW_OUTSIDE_ERROR(*env)
  inTotalCounts = tot;
  double maxImbalance = 0;
  for(size_t i = 0; i < totalPartCount;++i){
    float imbalance = 0;
    if(i > 0)
      imbalance = imbalanceOf(float(inTotalCounts[i] - inTotalCounts[i -1]), float(numGlobalCoords), 1 / float(totalPartCount));
    else
      imbalance = imbalanceOf(float(inTotalCounts[i] - 0), float(numGlobalCoords), 1 / float(totalPartCount) );
    if(imbalance > maxImbalance) maxImbalance = imbalance;
  }
  if(comm->getRank() == 0){
    cout << "Max imbalance:" << maxImbalance << endl;
    cout << "totalPartCount:" << totalPartCount << endl;
    cout << numeric_limits<scalar_t>::min() << " " <<
        numeric_limits<scalar_t>::max() << endl;

    for(size_t i = 0; i < totalPartCount - 1;++i){
      cout << "cut coordinate:" << allCutCoordinates[i] << endl;
    }
    double end = start;
#ifdef HAVE_ZOLTAN2_OMP
    end = omp_get_wtime( );
#endif
    cout << "start = "<< start << " end = " << end << " diff = " << end - start << endl;
  }


  ArrayRCP<const gno_t> gnoList = arcpFromArrayView(pqJagged_gnos);



  solution->setParts(gnoList, partId);





  ////////////////////////////////////////////////////////
  // Done: Compute quality metrics and update the solution
  // TODO: The algorithm will not compute the metrics.
  // It provides the solution to the PartitioningSolution.
  // Metrics can be computed with a method that takes
  // a model and a solution.

  for(int i = 0; i < numThreads; ++i){
    delete []coordinate_starts[i];
    delete []coordinate_ends[i] ;
  }

  //delete [] imbalance_tolerances;
  delete []coordinate_linked_list;
  //the start and end coordinate of  each part.
  delete []coordinate_starts;
  delete []coordinate_ends;

  //assign the max size to starts, as it will be reused.
  if(allowNonRectelinearPart){
    delete []nonRectelinearPart;

    for(int i = 0; i < numThreads; ++i){
      delete [] nonRectRatios[i];
    }
    delete [] nonRectRatios;
  }


  //delete [] partNo;
  delete []max_min_array;
  delete [] outTotalCounts;
  delete []partitionedPointCoordinates ;
  delete []newpartitionedPointCoordinates ;
  delete []allCutCoordinates;
  delete []pqJagged_coordinates;
  delete []pqJagged_weights;
  delete []pqJagged_uniformParts;
  delete []pqJagged_partSizes;
  delete []pqJagged_uniformWeights;

  delete []cutCoordinatesWork; // work array to manipulate coordinate of cutlines in different iterations.
  delete []cutPartRatios; // the weight ratios at left side of the cuts. First is 0, last is 1.

  delete []cutUpperBounds;  //to determine the next cut line with binary search
  delete []cutLowerBounds;  //to determine the next cut line with binary search

  delete []cutLowerWeight;  //to determine the next cut line with binary search
  delete []cutUpperWeight;  //to determine the next cut line with binary search

  delete []isDone;
  delete []totalPartWeights;

  for(int i = 0; i < numThreads; ++i){
    delete [] partWeights[i] ;
    delete [] rightClosestDistance[i];
    delete [] leftClosestDistance[i];
  }

  delete [] partWeights;

  delete [] leftClosestDistance ;
  delete [] rightClosestDistance;



}

} // namespace Zoltan2







#ifdef oldCode
template <typename scalar_t>
void getNewCoordinates(const size_t &total_part_count, const scalar_t * totalPartWeights, bool *isDone, const scalar_t *cutPartRatios,
    const scalar_t &totalWeight, const scalar_t &imbalanceTolerance, size_t &allDone, scalar_t *cutUpperBounds, scalar_t *cutLowerBounds,
    scalar_t *cutCoordinates, const size_t &noCuts, const scalar_t &maxCoordinate, const scalar_t &minCoordinate){
  scalar_t seenW = 0;
  float expected = 0;
  scalar_t leftImbalance = 0, rightImbalance = 0;

#pragma omp for
  for (size_t i = 0; i < total_part_count; i += 2){
    seenW = totalPartWeights[i];
    int currentLine = i/2;
    if(isDone[currentLine]) continue;
    expected = cutPartRatios[i / 2];

    leftImbalance = imbalanceOf(seenW, totalWeight, expected);
    rightImbalance = imbalanceOf(totalWeight - seenW, totalWeight, 1 - expected);

    //cout << "cut:" << currentLine<< " left and right:" << leftImbalance << " " << rightImbalance << endl;
    bool isLeftValid = leftImbalance <= imbalanceTolerance && leftImbalance >= -imbalanceTolerance;
    bool isRightValid = rightImbalance <= imbalanceTolerance && rightImbalance >= -imbalanceTolerance;

    if(isLeftValid && isRightValid){
      isDone[currentLine] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
      allDone -=1;
      //cout << "\tboth valid" << endl;
    } else if(cutUpperBounds[currentLine] != -1 && cutUpperBounds[currentLine] == cutLowerBounds[currentLine]){
      isDone[currentLine] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
      allDone -=1;
      //cout << "\tconverged upper:" <<  cutUpperBounds[currentLine] << " loweR:" << cutLowerBounds[currentLine] << endl;

    } else if(leftImbalance < 0){
      //TODO if boundary is enough binary search on processors
      //if boundary is not enough, binary search
      //cout << "moving to right" << endl;
      if (cutUpperBounds[currentLine] == -1){
        bool done = false;
        for (size_t ii = currentLine + 1; ii < noCuts ; ++ii){
          scalar_t pw = totalPartWeights[ii * 2];
          scalar_t ew = totalWeight * expected;
          if(pw >= ew){
            if(pw == ew){
              cutCoordinates[currentLine] = cutCoordinates[ii];
              isDone[currentLine] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
allDone -=1;
            } else {
              cutUpperBounds[currentLine] = cutCoordinates[ii];
              cutLowerBounds[currentLine] = cutCoordinates [ ii -1];
              cutCoordinates[currentLine] = (cutLowerBounds[currentLine] + cutUpperBounds[currentLine]) / 2;
            }
            done = true;
            break;
          }
        }
        if(!done){
          cutUpperBounds[currentLine] = maxCoordinate;
          cutLowerBounds[currentLine] = cutCoordinates [ noCuts -1];
          cutCoordinates[currentLine] = (cutLowerBounds[currentLine] + cutUpperBounds[currentLine]) / 2;
        }

      } else {
        cutLowerBounds[currentLine] = cutCoordinates[currentLine];
        cutCoordinates[currentLine] = (cutLowerBounds[currentLine] + cutUpperBounds[currentLine]) / 2;
      }

    } else {
      //move left
      //new coordinate binary search
      if (cutLowerBounds[currentLine] == -1){

        bool done = false;
        for (int ii = currentLine - 1; ii >= 0; --ii){
          scalar_t pw = totalPartWeights[ii * 2];
          scalar_t ew = totalWeight * expected;
          if(pw <= ew){
            if(pw == ew){
              cutCoordinates[currentLine] = cutCoordinates[ii];
              isDone[currentLine] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
allDone -=1;
            } else {
              cutLowerBounds[currentLine] = cutCoordinates[ii];
              cutUpperBounds[currentLine] = cutCoordinates[ii + 1];
              cutCoordinates[currentLine] = (cutLowerBounds[currentLine] + cutUpperBounds[currentLine]) / 2;
            }
            done = true;
            break;
          }

        }
        if(!done){
          cutUpperBounds[currentLine] = cutCoordinates[0];
          cutLowerBounds[currentLine] = minCoordinate;
          cutCoordinates[currentLine] = (cutLowerBounds[currentLine] + cutUpperBounds[currentLine]) / 2;

        }
      } else {
        cutUpperBounds[currentLine] = cutCoordinates[currentLine];
        cutCoordinates[currentLine] = (cutLowerBounds[currentLine] + cutUpperBounds[currentLine]) / 2;
      }

    }

  }
}


template <typename scalar_t, typename lno_t>
void pqJagged_1DPart(scalar_t *pqJagged_coordinates,	scalar_t *pqJagged_weights,	bool pqJagged_uniformWeights,
    const size_t &numLocalCoords,const global_size_t &numGlobalCoords,	scalar_t &minCoordinate,
    scalar_t &maxCoordinate, bool pqJagged_uniformParts, scalar_t *pqJagged_partSizes, size_t partNo, int noThreads,
    scalar_t imbalanceTolerance){


  imbalanceTolerance = 0.001;
  size_t noCuts = partNo - 1;
  scalar_t *cutCoordinates = new scalar_t[noCuts]; // coordinates of the cut lines. First one is the min, last one is max coordinate.
  scalar_t *cutPartRatios = new scalar_t[noCuts];  // the weight ratios at left side of the cuts. First is 0, last is 1.

  scalar_t *cutUpperBounds = new scalar_t [noCuts];  //to determine the next cut line with binary search
  scalar_t *cutLowerBounds = new scalar_t [noCuts];  //to determine the next cut line with binary search

  pqJagged_getCutCoord_Weights<scalar_t>(minCoordinate, maxCoordinate,
      pqJagged_uniformParts, pqJagged_partSizes, noCuts,
      cutCoordinates, cutPartRatios);

  //make a function for this.
  //assuming that multiplication via 0 or 1 is faster than branches.
  scalar_t *wghts = NULL;
  scalar_t unit_weight = 1;
  lno_t v_scaler;
  if (pqJagged_uniformWeights){
    wghts = &unit_weight;
    v_scaler = 0;
  }
  else {
    wghts = pqJagged_weights;
    v_scaler = 1;
  }


  //calculate total weight
  scalar_t totalWeight = 0;

  lno_t *coordinate_linked_list = new lno_t[numLocalCoords];
  lno_t **coordinate_starts = new lno_t *[noThreads];
  lno_t **coordinate_ends = new lno_t *[noThreads];
  scalar_t **partWeights = new scalar_t *[noThreads];

  bool *isDone = new bool [noCuts];

  size_t total_part_count = partNo + noCuts;
  for(int i = 0; i < noThreads; ++i){
    coordinate_starts[i] = new lno_t[total_part_count];
    coordinate_ends[i] = new lno_t[total_part_count];
    partWeights[i] = new scalar_t[total_part_count];
  }

  scalar_t *totalPartWeights = new scalar_t[total_part_count];

  size_t allDone = noCuts;
#pragma omp parallel shared(allDone)
  {
#pragma omp for reduction(+:totalWeight)
    for (size_t i = 0; i < numLocalCoords; ++i){
      totalWeight += wghts[i * v_scaler];
    }


    int me = omp_get_thread_num();
    lno_t *myStarts = coordinate_starts[me];
    lno_t *myEnds = coordinate_ends[me];
    scalar_t *myPartWeights = partWeights[me];


    pqJagged_PartVertices <lno_t, size_t> pqPV;
    pqPV.set(coordinate_linked_list, myStarts, myEnds);

#pragma omp for
    for(size_t i = 0; i < noCuts; ++i){
      isDone[i] = false;
      cutLowerBounds[i] = -1;
      cutUpperBounds[i] = -1;
    }

#pragma omp for
    for (size_t i = 0; i < numLocalCoords; ++i){
      coordinate_linked_list[i] = -1;
    }
    //implicit barrier


    while (allDone != 0){
      cout << "allDone:" << allDone << endl;
      /*
for (size_t i = 0; i < noCuts; ++i){
cout << "i:" << i << " coordinate:" << cutCoordinates[i] << endl;
}
       */

      for (size_t i = 0; i < total_part_count; ++i){
        myEnds[i] = -1;
        myStarts[i] = -1;
        myPartWeights[i] = 0;
      }

#pragma omp for
      for (size_t i = 0; i < numLocalCoords; ++i){
        for(size_t j = 0; j < noCuts; ++j){
          if (pqJagged_coordinates[i] <= cutCoordinates[j]){
            if (pqJagged_coordinates[i] == cutCoordinates[j]){
              myPartWeights[j * 2] -=	wghts [i * v_scaler];
              pqPV.inserToPart(j * 2 + 1, i);
            }
            else {
              pqPV.inserToPart(j * 2, i);
            }
            break;
          } else {
            scalar_t w = wghts [i * v_scaler];
            myPartWeights[j * 2] -=	 w ;
            myPartWeights[j * 2 + 1] -=	w;
          }
        }
      }


#pragma omp for
      for(size_t j = 0; j < total_part_count; ++j){
        scalar_t pwj = 0;
        for (int i = 0; i < noThreads; ++i){
          pwj += partWeights[i][j];
        }
        totalPartWeights[j] = pwj + totalWeight;
      }

      //all to all partweight send;
      //get totalpartweights from different nodes sized total_part_count

      getNewCoordinates<scalar_t>(total_part_count, totalPartWeights, isDone, cutPartRatios,
          totalWeight, imbalanceTolerance, allDone, cutUpperBounds, cutLowerBounds,
          cutCoordinates, noCuts,maxCoordinate, minCoordinate);
    }
  }

  for(size_t j = 0; j < total_part_count; j+=2){
    cout << "j:" << j/2 << " pw:" << totalPartWeights[j] << endl;
  }
}

#endif


#endif





