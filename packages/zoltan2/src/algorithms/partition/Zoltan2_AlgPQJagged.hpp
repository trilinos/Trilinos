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
#include "stdio.h"

#ifdef HAVE_ZOLTAN2_OMP
//#define USE_LESS_THREADS
#endif
#define CACHE_LINE_SIZE 64
#define EPS_SCALE 1
#define LEAST_SIGNIFICANCE 0.0001
#define SIGNIFICANCE_MUL 1000

//#undef HAVE_ZOLTAN2_OMP
#ifdef HAVE_ZOLTAN2_OMP
#include <omp.h>
#endif
#define FIRST_TOUCH

//#define RCBCODE
#define LEAF_IMBALANCE_FACTOR 0.1

//imbalance calculation. Wreal / Wexpected - 1
#define imbalanceOf(Wachieved, totalW, expectedRatio) \
    (Wachieved) / ((totalW) * (expectedRatio)) - 1


namespace Teuchos{
template <typename Ordinal, typename T>
class PQJaggedCombinedReductionOp  : public ValueTypeReductionOp<Ordinal,T>
{
private:
  Ordinal numSum_, numMin_1, numMin_2;

public:
  /*! \brief Default Constructor
   */
  PQJaggedCombinedReductionOp ():numSum_(0), numMin_1(0), numMin_2(0){}

  /*! \brief Constructor
   *   \param nsum  the count of how many sums will be computed at the
   *             start of the list.
   *   \param nmin  following the sums, this many minimums will be computed.
   *   \param nmax  following the minimums, this many maximums will be computed.
   */
  PQJaggedCombinedReductionOp (Ordinal nsum, Ordinal nmin1, Ordinal nmin2):
    numSum_(nsum), numMin_1(nmin1), numMin_2(nmin2){}

  /*! \brief Implement Teuchos::ValueTypeReductionOp interface
   */
  void reduce( const Ordinal count, const T inBuffer[], T inoutBuffer[]) const
  {
    Ordinal next=0;
    for (Ordinal i=0; i < numSum_; i++, next++)
      inoutBuffer[next] += inBuffer[next];

    for (Ordinal i=0; i < numMin_1; i++, next++)
      if (inoutBuffer[next] > inBuffer[next])
        inoutBuffer[next] = inBuffer[next];

    for (Ordinal i=0; i < numMin_2; i++, next++)
      if (inoutBuffer[next] > inBuffer[next])
        inoutBuffer[next] = inBuffer[next];
  }


};



template <typename Ordinal, typename T>
class PQJaggedCombinedMinMaxReductionOp  : public ValueTypeReductionOp<Ordinal,T>
{
private:
  Ordinal numMin, numMax;

public:
  /*! \brief Default Constructor
   */
  PQJaggedCombinedMinMaxReductionOp ():numMin(0), numMax(0){}

  /*! \brief Constructor
   *   \param nsum  the count of how many sums will be computed at the
   *             start of the list.
   *   \param nmin  following the sums, this many minimums will be computed.
   *   \param nmax  following the minimums, this many maximums will be computed.
   */
  PQJaggedCombinedMinMaxReductionOp (Ordinal nmin, Ordinal nmax):
    numMin(nmin), numMax(nmax){}

  /*! \brief Implement Teuchos::ValueTypeReductionOp interface
   */
  void reduce( const Ordinal count, const T inBuffer[], T inoutBuffer[]) const
  {
    Ordinal next=0;

    for (Ordinal i=0; i < numMin; i++, next++)
      if (inoutBuffer[next] > inBuffer[next])
        inoutBuffer[next] = inBuffer[next];

    for (Ordinal i=0; i < numMax; i++, next++)
      if (inoutBuffer[next] < inBuffer[next])
        inoutBuffer[next] = inBuffer[next];
  }
};
} // namespace Teuchos

namespace Zoltan2{

//diffclock for temporary timing experiments.



/*! \brief A helper class containing array representation of
 *  coordinate linked lists.
 */

template <typename lno_t, typename partId_t>
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
  void inserToPart (partId_t partNo, lno_t coordinateIndex){

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

template<typename T>
inline void firstTouch(T *arrayName, size_t arraySize){
#pragma omp parallel for
  for(size_t jj = 0; jj < arraySize; ++jj){
    arrayName[jj] = 0;
  }
}

template <class lno_t,class scalar_t>
class PQJagged_MaxPriorityQueue {
private:
    scalar_t *weights; //gains according to which the heap order is determined.
    lno_t *real_indices;   //the vertex indices corresponding to the gains.
    lno_t *inverse_vertex_positions; //points the position of each vertex
    lno_t heap_size;    //the end of heap inclusive
    lno_t max_size;    //max size of heap. --for addition

    void push_down(lno_t i){
      lno_t child1 = 2 * i + 1;
      lno_t child2 = 2 * i + 2;
      if((child1 >= heap_size || this->weights[i] >= this->weights[child1]) && (child2 >= heap_size || this->weights[i] >= this->weights[child2])){
          return;
      } else {
        lno_t biggerchild = 0;
          if(this->weights[child1] >= this->weights[child2]){
              biggerchild = child1;
          } else {
              biggerchild = child2;
          }
          this->swap(i, biggerchild);
          push_down(biggerchild);
      }
    }

    void push_up(lno_t i){
      if(i == 0 ) return;
      lno_t parent = (i - 1) / 2;
      //std::cout << "parent:" << parent << " child" << i << std::endl;
      if(this->weights[parent] >= this->weights[i]){
        return;
      }
      else {
        this->swap(i, parent);
        push_up(parent);
      }
    }

    void swap(lno_t index1, lno_t index2){
        scalar_t wtmp = this->weights[index1];
        this->weights[index1] = this->weights[index2];
        this->weights[index2] = wtmp;

        lno_t itmp = this->inverse_vertex_positions[this->real_indices[index1]];
        this->inverse_vertex_positions[this->real_indices[index1]] = this->inverse_vertex_positions[this->real_indices[index2]];
        this->inverse_vertex_positions[this->real_indices[index2]] = itmp;

        itmp = this->real_indices[index1];
        this->real_indices[index1] = this->real_indices[index2];
        this->real_indices[index2] = itmp;
    }

public:
    //assumes indices starts from 1.
    PQJagged_MaxPriorityQueue(lno_t max_size, scalar_t *gains, lno_t *real_indices, lno_t *positions){
      this->set_Heap(max_size, gains, real_indices, positions);
    }
    PQJagged_MaxPriorityQueue(){}
    ~PQJagged_MaxPriorityQueue(){}
    void set_Heap(lno_t max_size_, scalar_t *weights_, lno_t *real_indices_, lno_t *positions_){
        this->weights = weights_;
        this->max_size = max_size_;
        this->real_indices = real_indices_;
        this->inverse_vertex_positions = positions_;
        this->heap_size = max_size_;
    }

    void update(lno_t real_index, scalar_t weight_increase){
      lno_t heap_index = this->inverse_vertex_positions[real_index];
      if(heap_index == -1){
        return;
      }
      this->weights[heap_index] += weight_increase;
      if(weight_increase > 0){
        this->push_up(heap_index);
      }
      else {
        this->push_down(heap_index);
      }
    }

    void top(lno_t *index, scalar_t *weight){
      if(heap_size >= 0){
          *index = this->real_indices[1];
          *weight = this->weights[1];
      } else {
          *index = -1;
          *weight = -1;
      }
    }

    void pop(){
      if(heap_size > 0){
        this->inverse_vertex_positions[this->real_indices[0]] = -1;
        this->weights[0] = this->weights[heap_size - 1];
        this->real_indices[0] = this->real_indices[heap_size - 1];
        this->inverse_vertex_positions[this->real_indices[heap_size - 1]] = 0;
        --heap_size;
        this->push_down(0);
      }
    }
};

template <typename lno_t, typename gno_t, typename scalar_t, typename partId_t>
void getMigrationGroups(lno_t *totalPartPointCounts, gno_t totalCount, partId_t partNo, int worldSize, scalar_t *proc_work_array,
    partId_t *part_real_indices_work, partId_t *part_positions_work
    ){
  scalar_t idealWeight = totalCount / double(worldSize);
  partId_t maxGroupNo = partNo;

  for (int i = 0; i < worldSize; ++i){
    proc_work_array[i] = idealWeight;
  }
  //PQJagged_MaxPriorityQueue<lno_t, scalar_t> parts(partNo, totalPartPointCounts_work, part_real_indices_work, lno_t *part_positions_work);
  //PQJagged_MaxPriorityQueue<lno_t, scalar_t> procs(worldSize, proc_work_array, part_real_indices_work, lno_t *part_positions_work);

}

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

  if(cutUpperWeight[currentCutIndex] == cutLowerWeight[currentCutIndex]){
    return cutLowerBounds[currentCutIndex];
  }
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
    multiCriteriaNorm &mcnorm, std::bitset<NUM_RCB_PARAMS> &params,  int &numTestCuts, bool &ignoreWeights, bool &allowNonrectelinear){

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
    imbalanceTolerance = 10e-4;

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

  if (isSet && intChoice==1){
    params.set(rcb_rectilinearBlocks);
    allowNonrectelinear = false;
  }

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

  //scalar_t minm = 0;scalar_t maxm = 0;

  {
    //cout << "t:" << omp_get_num_threads();
    scalar_t localMinMax[2], globalMinMax[2];
    localMinMax[0] = minCoordinate;
    localMinMax[1] = maxCoordinate;

    Teuchos::PQJaggedCombinedMinMaxReductionOp<int, scalar_t> reductionOp(
       1,   // number of left closest mins
       1);  // number of right closest max

    try{

      reduceAll<int, scalar_t>(*comm, reductionOp,
          2, localMinMax, globalMinMax
        );
      /*
      reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_MIN,
          1, &minCoordinate, &(globalMinMax[0])
        );
      reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_MAX ,
          1, &maxCoordinate, &(globalMinMax[1])
        );
        */
    }
    Z2_THROW_OUTSIDE_ERROR(*env)

    minCoordinate = globalMinMax[0];
    maxCoordinate = globalMinMax[1];
  }
}


template <typename scalar_t, typename lno_t>
void pqJagged_getMinMaxCoordofAllParts(RCP<Comm<int> > &comm, scalar_t *pqJagged_coordinates, scalar_t *global_min_max_coordinates,
    const RCP<const Environment> &env, int numThreads,
    lno_t *partitionedPointPermutations, lno_t *coordinateBegins,
    scalar_t *max_min_array /*sized nothreads * 2*/,
    scalar_t maxScalar, scalar_t minScalar, partId_t currentPartCount, scalar_t *local_min_max_coordinates /*sized currentPartCount * 2 */){

  for (partId_t kk = 0; kk < currentPartCount; ++kk){
    lno_t coordinateBegin = kk == 0 ? 0 : coordinateBegins[kk - 1];
    lno_t coordinateEnd = coordinateBegins[kk];
    scalar_t minCoordinate = maxScalar;
    scalar_t maxCoordinate = minScalar;
    if(coordinateBegin < coordinateEnd){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
      {
        int myId = 0;
#ifdef HAVE_ZOLTAN2_OMP
        myId = omp_get_thread_num();
#endif
        scalar_t myMin, myMax;
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
    local_min_max_coordinates[kk] = minCoordinate;
    local_min_max_coordinates[kk + currentPartCount] = maxCoordinate;
  }
  //scalar_t minm = 0;scalar_t maxm = 0;

  {
    Teuchos::PQJaggedCombinedMinMaxReductionOp<int, scalar_t> reductionOp(
       currentPartCount,   // number of left closest mins
       currentPartCount);  // number of right closest max

    try{

      reduceAll<int, scalar_t>(*comm, reductionOp,
          currentPartCount * 2, local_min_max_coordinates, global_min_max_coordinates
        );
    }
    Z2_THROW_OUTSIDE_ERROR(*env)
  }
}


template <typename scalar_t>
void pqJagged_getCutCoord_Weights(
    scalar_t minCoordinate, scalar_t maxCoordinate,
    bool pqJagged_uniformParts, scalar_t *pqJagged_partSizes /*p sized, weight ratios of each part*/,
    partId_t noCuts/*p-1*/ ,
    scalar_t *cutCoordinates /*p - 1 sized, coordinate of each cut line*/,
    scalar_t *cutPartRatios /*cumulative weight ratios, at left side of each cut line. p-1 sized*/,
    int numThreads){

  scalar_t coordinateRange = maxCoordinate - minCoordinate;
  if(pqJagged_uniformParts){
    scalar_t uniform = 1. / (noCuts + 1);
    scalar_t slice = uniform * coordinateRange;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
    {


      partId_t myStart = 0;
      partId_t myEnd = noCuts;

#ifdef USE_LESS_THREADS
      int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
      me = omp_get_thread_num();
#endif

      int scalar_t_bytes = sizeof(scalar_t);
      int requiredPull = ceil(noCuts / float(numThreads));
      int minPull = CACHE_LINE_SIZE / scalar_t_bytes;
      if(requiredPull > minPull){
        minPull = requiredPull;
      }
      myStart = me * minPull;
      myEnd = min(myStart + minPull, noCuts);
#endif


#ifndef USE_LESS_THREADS
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
#endif
      for(partId_t i = myStart; i < myEnd; ++i){
        cutPartRatios[i] =  uniform * (i + 1);
        cutCoordinates[i] = minCoordinate + slice * (i + 1);
      }
      cutPartRatios[noCuts] = 1;
    }
  }
  else {
    cutPartRatios[0] = pqJagged_partSizes[0];
    cutCoordinates[0] = coordinateRange * cutPartRatios[0];
    for(partId_t i = 1; i < noCuts; ++i){
      cutPartRatios[i] = pqJagged_partSizes[i] + cutPartRatios[i - 1];
      cutCoordinates[i] = coordinateRange * cutPartRatios[i];
    }
  }
  /*
  for(partId_t i = 0; i < noCuts; ++i){
    cout << "c[" << i << "]:" << cutCoordinates[i] << endl;

  }
  */

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
    partId_t &allDone, scalar_t *cutUpperBounds, scalar_t *cutLowerBounds,
    scalar_t *&cutCoordinates, const partId_t &noCuts,
    const scalar_t &maxCoordinate, const scalar_t &minCoordinate,
    scalar_t *leftClosestDistance, scalar_t *rightClosestDistance,
    scalar_t * cutLowerWeight,scalar_t * cutUpperWeight, scalar_t *&cutCoordinatesWork,
    float *nonRectelinearPart,
    bool allowNonRectelinearPart, scalar_t maxScalar, const scalar_t * localtw,
    const RCP<const Environment> &env,RCP<Comm<int> > &comm,
    partId_t *rectelinearCount, scalar_t *cutWeights, scalar_t *globalCutWeights,
    partId_t myCutStart, partId_t myCutEnd){


  scalar_t seenW = 0;
  float expected = 0;
  scalar_t leftImbalance = 0, rightImbalance = 0;

  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();

#ifndef USE_LESS_THREADS
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
#endif
  for (partId_t i = myCutStart; i < myCutEnd; i++){

    //if a left and right closes point is not found, set the distance to 0.
    if(leftClosestDistance[i] == maxScalar)
      leftClosestDistance[i] = 0;
    if(rightClosestDistance[i] == maxScalar)
      rightClosestDistance[i] = 0;

  }
#ifdef HAVE_ZOLTAN2_OMP
#ifdef USE_LESS_THREADS
#pragma omp barrier
#endif
#endif

#ifndef USE_LESS_THREADS
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
#endif
  for (partId_t i = myCutStart; i < myCutEnd; i++){
    //if already determined at previous iterations, do nothing.

    if(isDone[i]) {
      cutCoordinatesWork[i] = cutCoordinates[i];
      continue;
    }
    //current weight of the part at the left of the cut line.
    seenW = totalPartWeights[i * 2];
    globalCutWeights[i] = 0;
    cutWeights[i] = 0;

    expected = cutPartRatios[i];
    leftImbalance = imbalanceOf(seenW, totalWeight, expected);
    rightImbalance = imbalanceOf(totalWeight - seenW, totalWeight, 1 - expected);

    //bool isLeftValid = leftImbalance <= imbalanceTolerance && leftImbalance >= -imbalanceTolerance;
    //bool isRightValid = rightImbalance <= imbalanceTolerance && rightImbalance >= -imbalanceTolerance;

    bool isLeftValid = fabs(leftImbalance) - imbalanceTolerance < _EPSILON ;//* 10;
    bool isRightValid = fabs(rightImbalance) - imbalanceTolerance < _EPSILON; // * 10;


    //if the cut line reaches to desired imbalance.
    if(isLeftValid && isRightValid){

      isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
      allDone -=1;
      cutCoordinatesWork [i] = cutCoordinates[i];
      continue;
    }
    //if upper bound and lower bound reaches to the same point.
    //the imbalance cannot be satisfied. Accept as it is.
    //if left imbalance is lower than 0, then cut line should be moved to right.
    else if(leftImbalance < 0){
      scalar_t ew = totalWeight * expected;



      if(allowNonRectelinearPart){

        if (totalPartWeights[i * 2 + 1] == ew){
          isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          allDone -=1;
          cutCoordinatesWork [i] = cutCoordinates[i];
          nonRectelinearPart[i] = 1;
          continue;
        }
        else if (totalPartWeights[i * 2 + 1] > ew){
          isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          allDone -=1;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          *rectelinearCount += 1;
          cutCoordinatesWork [ i] = cutCoordinates[i];
          scalar_t myWeightOnLine = localtw[i * 2 + 1] - localtw[i * 2];
          cutWeights[i] = myWeightOnLine;
          continue;
        }
      }

      //when moving left, and upper and lower are set.
      //set lower bound to current line.

      cutLowerBounds[i] = cutCoordinates[i] + rightClosestDistance[i];
      cutLowerWeight[i] = seenW;

      //compare the upper bound with the current lines.
      for (partId_t ii = i + 1; ii < noCuts ; ++ii){
        scalar_t pw = totalPartWeights[ii * 2];
        scalar_t lw = totalPartWeights[ii * 2 + 1];
        if(pw >= ew){
          if(pw == ew){
            cutUpperBounds[i] = cutCoordinates[ii];
            //cutUpperWeight[i] = totalPartWeights [2 * ii + 1];
            cutUpperWeight[i] = pw;
            cutLowerBounds[i] = cutCoordinates[ii];
            cutLowerWeight[i] = pw;
          } else if (pw < cutUpperWeight[i]){
            //if a cut line is more strict than the current upper bound,
            //update the upper bound.
            cutUpperBounds[i] = cutCoordinates[ii] - leftClosestDistance[ii];
            //cutUpperWeight[i] = totalPartWeights [2 * ii + 1];
            cutUpperWeight[i] = pw;
          }
          break;
        }
        //if comes here then pw < ew
        if(lw >= ew){
          cutUpperBounds[i] = cutCoordinates[ii];
          //cutUpperWeight[i] = totalPartWeights [2 * ii + 1];
          cutUpperWeight[i] = lw;
          cutLowerBounds[i] = cutCoordinates[ii];
          cutLowerWeight[i] = pw;
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
    } else {

      //moving to left.
      scalar_t ew = totalWeight * expected;

      //moving left, both upper and lower is set.
      //set upper to current line.
      cutUpperBounds[i] = cutCoordinates[i] - leftClosestDistance[i];
      cutUpperWeight[i] = seenW;

      // compare the current cut line weights with previous upper and lower bounds.
      for (int ii = i - 1; ii >= 0; --ii){
        scalar_t pw = totalPartWeights[ii * 2];
        scalar_t lw = totalPartWeights[ii * 2 + 1];


        if(pw <= ew){
          if(pw == ew){
            cutUpperBounds[i] = cutCoordinates[ii];
            //cutUpperWeight[i] = totalPartWeights [2 * ii + 1];
            cutUpperWeight[i] = pw;
            cutLowerBounds[i] = cutCoordinates[ii];
            cutLowerWeight[i] = pw;
          }
          else if (pw > cutLowerWeight[i]){
            cutLowerBounds[i] = cutCoordinates[ii] + rightClosestDistance[ii];
            //cutLowerWeight[i] = totalPartWeights [2 * ii + 1];
            cutLowerWeight[i] = pw;
            if(lw > ew){
              cutUpperBounds[i] = cutCoordinates[ii] + rightClosestDistance[ii];
              //cutLowerWeight[i] = totalPartWeights [2 * ii + 1];
              cutUpperWeight[i] = lw;
            }
          }
          break;
        }
        if (pw >= ew && (pw < cutUpperWeight[i] || (pw == cutUpperWeight[i] && cutUpperBounds[i] > cutCoordinates[ii] - leftClosestDistance[ii]))){
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


#ifdef HAVE_ZOLTAN2_OMP
#ifdef USE_LESS_THREADS
#pragma omp barrier
#endif
#endif
#pragma omp barrier
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
  {
    if(*rectelinearCount > 0){
      try{
        Teuchos::scan<int,scalar_t>(
            *comm, Teuchos::REDUCE_SUM,
            noCuts,cutWeights, globalCutWeights
        );
      }
      Z2_THROW_OUTSIDE_ERROR(*env)


      for (partId_t i = 0; i < noCuts; ++i){
        //cout << "gw:" << globalCutWeights[i] << endl;
        if(globalCutWeights[i] > 0) {
          scalar_t expected = cutPartRatios[i];
          scalar_t ew = totalWeight * expected;
          scalar_t expectedWeightOnLine = ew - totalPartWeights[i * 2];
          scalar_t myWeightOnLine = cutWeights[i];
          scalar_t weightOnLineBefore = globalCutWeights[i];
          scalar_t incMe = expectedWeightOnLine - weightOnLineBefore;
          scalar_t mine = incMe + myWeightOnLine;
          if(mine < 0){
            nonRectelinearPart[i] = 0;
          }
          else if(mine >= myWeightOnLine){
            nonRectelinearPart[i] = 1;
          }
          else {
            nonRectelinearPart[i] = mine / myWeightOnLine;
          }
          //cout << "nonrec:" << i << " : " << nonRectelinearPart[i] << endl;
        }
      }
      *rectelinearCount = 0;
    }

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
    scalar_t &minCoordinate, scalar_t &maxCoordinate, partId_t partNo, int numThreads,
    scalar_t imbalanceTolerance, scalar_t *cutCoordinates, lno_t *partitionedPointPermutations,
    lno_t coordinateBegin, lno_t coordinateEnd,

    scalar_t *cutCoordinatesWork, 	// work array to manipulate coordinate of cutlines in different iterations.
    scalar_t *cutPartRatios,		// the weight ratios at left side of the cuts. last is 1.

    scalar_t *cutUpperBounds,  //to determine the next cut line with binary search
    scalar_t *cutLowerBounds,  //to determine the next cut line with binary search

    scalar_t *cutLowerWeight,  //to determine the next cut line with binary search
    scalar_t *cutUpperWeight,   //to determine the next cut line with binary search
    bool *isDone,
    double **partWeights,
    scalar_t **leftClosestDistance,
    scalar_t **rightClosestDistance,
    scalar_t *totalPartWeights_leftClosest_rightCloset,
    scalar_t *global_totalPartWeights_leftClosest_rightCloset,
    scalar_t &globaltotalWeight,
    float *nonRectelinearPart,
    bool allowNonRectelinearPart,
    scalar_t maxScalar,
    scalar_t minScalar,
    scalar_t *cutWeights, scalar_t *globalCutWeights
){

  partId_t recteLinearCount = 0;
  scalar_t *cutCoordinates_tmp = cutCoordinates;
  partId_t noCuts = partNo - 1;
  scalar_t totalWeight = 0;
  size_t total_part_count = partNo + size_t (noCuts) ;
  partId_t allDone = noCuts;

  scalar_t *tw = totalPartWeights_leftClosest_rightCloset;
  scalar_t *lc = totalPartWeights_leftClosest_rightCloset + total_part_count;
  scalar_t *rc = totalPartWeights_leftClosest_rightCloset + total_part_count + noCuts;

  scalar_t *glc = global_totalPartWeights_leftClosest_rightCloset + total_part_count;
  scalar_t *grc = global_totalPartWeights_leftClosest_rightCloset + total_part_count + noCuts;
  scalar_t *gtw = global_totalPartWeights_leftClosest_rightCloset;

  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();

  globaltotalWeight = 0;

  Teuchos::PQJaggedCombinedReductionOp<int, scalar_t> reductionOp(
     total_part_count,      // number of sums
     noCuts,   // number of left closest mins
     noCuts);  // number of right closest mins



#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(allDone, globaltotalWeight, recteLinearCount)
#endif
  {
    //int iterationCount = 0;
    //calculate total weight
    if (pqJagged_uniformWeights) {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
      {
        totalWeight = coordinateEnd - coordinateBegin;
      }
    } else {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for reduction(+:totalWeight)
#endif
      for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
        int i = partitionedPointPermutations[ii];
        totalWeight += pqJagged_weights[i];
      }
    }


//#pragma omp barrier
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
    {
      globaltotalWeight = totalWeight;
      scalar_t tw = 0;
      try{
        reduceAll<int,scalar_t>(
            *comm, Teuchos::REDUCE_SUM,
            1, &totalWeight, &tw
        );

      }
      Z2_THROW_OUTSIDE_ERROR(*env)
      globaltotalWeight = tw;
    }

    if(noCuts == 0){
      tw[0] = globaltotalWeight;
    }
    else {


      int me = 0;

      partId_t myCutStart = 0;
      partId_t myCutEnd = noCuts;
#ifdef HAVE_ZOLTAN2_OMP
      me = omp_get_thread_num();
#ifdef USE_LESS_THREADS
      int scalar_t_bytes = sizeof(scalar_t);
      partId_t requiredPull = ceil(noCuts / double(numThreads));
      partId_t minPull = partId_t (ceil(CACHE_LINE_SIZE / double(scalar_t_bytes)));
      if(requiredPull > minPull){
        minPull = requiredPull;
      }
      myCutStart = me * minPull;
      myCutEnd = min(myCutStart + minPull, noCuts);
#endif
#endif
      //scalar_t *myPartWeights = partWeights[me];
      double *myPartWeights = partWeights[me];
      scalar_t *myLeftClosest = leftClosestDistance[me];
      scalar_t *myRightClosest = rightClosestDistance[me];


#ifndef USE_LESS_THREADS
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
#endif
      for(partId_t i = myCutStart; i < myCutEnd; ++i){
        isDone[i] = false;
        //cutLowerBounds[i] = -1;
        //cutUpperBounds[i] = -1;
        //instead of commented.
        cutLowerBounds[i] = minCoordinate;
        cutUpperBounds[i] = maxCoordinate;
        cutUpperWeight[i] = globaltotalWeight;
        cutLowerWeight[i] = 0;
        //instead of commented.

        if(allowNonRectelinearPart){
          nonRectelinearPart[i] = 0;
        }
      }


#ifdef USE_LESS_THREADS
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#endif
#endif
      int iteration = 0;
      while (allDone != 0){
        iteration += 1;
        /*
        if(comm->getRank() == 0)
        {
#pragma omp single
          {
            cout << endl << endl << "allDone:" << allDone << endl;
            for (size_t i = 0; i < noCuts; ++i){

              if(isDone[i] == false)
                cout << "i:" << i <<  " c:" << cutCoordinates_tmp[i] << " u:" << cutUpperBounds[i] << " l:" << cutLowerBounds[i] << " not done" << endl;
              else
                cout << "i:" << i <<  " c:" << cutCoordinates_tmp[i] <<  " done" << endl;
            }
          }
        }
        */

        for (size_t i = 0; i < total_part_count; ++i){
          if(i/2 < size_t(noCuts) && isDone[i/2]) continue;
          //myPartWeights[i] = 0;
          myPartWeights[i] = 0;
        }
        for(partId_t i = 0; i < noCuts; ++i){
          if(isDone[i]) continue;
          myLeftClosest[i] = maxScalar;
          myRightClosest[i] = maxScalar;
        }


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
          int i = partitionedPointPermutations[ii];

          //get a coordinate and compare it with cut lines from left to right.
          for(partId_t j = 0; j < noCuts; ++j){

            if(isDone[j]) continue;
            scalar_t distance = pqJagged_coordinates[i] - cutCoordinates_tmp[j];
            scalar_t absdistance = fabs(distance);
            if(absdistance < _EPSILON){
              scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
              myPartWeights[j * 2] -= w;
              //myPartWeights[j * 2] += w;
              myLeftClosest[j] = 0;
              myRightClosest[j] = 0;
              //break;
            }
            else
              if (distance < 0) {
                distance = -distance;
                if (/*myLeftClosest[j] < 0 ||*/ myLeftClosest[j] > distance){
                  myLeftClosest[j] = distance;
                }
                break;
              }
            //if it is on the left
              else {
                //if on the right, continue with the next line.
                scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
                myPartWeights[j * 2] -=	w;
                myPartWeights[j * 2 + 1] -=	w;
                //myPartWeights[j * 2] += w;
                //myPartWeights[j * 2 + 1] += w;

                if (/*myRightClosest[j] < 0 ||*/ myRightClosest[j] > distance){
                  myRightClosest[j] = distance;
                }
              }
          }
        }

/*
#pragma omp critical
        {
          for(size_t j = 0; j < noCuts; ++j){
            cout << "me:" << me << " j:" << j << " w:" << myPartWeights[j * 2] << endl;
          }

          for(size_t j = 0; j < noCuts; ++j){
            cout << "cut w:" << me << " j:" << j << " w:" << myPartWeights[j * 2 + 1] << endl;
          }
        }
*/
#ifndef USE_LESS_THREADS
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
#endif
        for(partId_t i = myCutStart; i < myCutEnd; ++i){
          if(isDone[i]) continue;
          scalar_t minl = leftClosestDistance[0][i], minr = rightClosestDistance[0][i];

          for (int j = 1; j < numThreads; ++j){
            if (rightClosestDistance[j][i] < minr ){
              minr = rightClosestDistance[j][i];
            }
            if (leftClosestDistance[j][i] < minl ){
              minl = leftClosestDistance[j][i];
            }

          }
          lc[i] = minl;
          rc[i] = minr;
        }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for(size_t j = 0; j < total_part_count; ++j){
          if(noCuts != 0 && partId_t(j/2) < noCuts && isDone[j/2]) continue;
          double pwj = 0;
          for (int i = 0; i < numThreads; ++i){
            pwj += partWeights[i][j];
          }
          tw[j] = totalWeight + pwj;
        }


        //all to all partweight send;
        //get totalpartweights from different nodes sized total_part_count
        //reduce part sizes here within all mpi processes.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
        {
          //if the no point on the left or right of the cut line in the current processor.
          try{
            reduceAll<int, scalar_t>(*comm, reductionOp,
                total_part_count + 2 * noCuts, tw, global_totalPartWeights_leftClosest_rightCloset
              );
          }
          Z2_THROW_OUTSIDE_ERROR(*env)
        }

        getNewCoordinates_simple<scalar_t>(total_part_count, gtw, isDone, cutPartRatios,
            globaltotalWeight, imbalanceTolerance, allDone, cutUpperBounds, cutLowerBounds,
            cutCoordinates_tmp, noCuts,maxCoordinate, minCoordinate, glc,
            grc,cutLowerWeight, cutUpperWeight,cutCoordinatesWork,
            nonRectelinearPart,
            allowNonRectelinearPart, maxScalar, tw,
            env,comm,
            &recteLinearCount, cutWeights, globalCutWeights, myCutStart, myCutEnd);
      }

      if(comm->getRank() == 0)
      if(me == 0) cout << "it:" << iteration << endl;
      //we cannot swap the arrays
      if (cutCoordinates != cutCoordinates_tmp){


#ifndef USE_LESS_THREADS
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
#endif
        for(partId_t i = myCutStart; i < myCutEnd; ++i){
          cutCoordinates[i] = cutCoordinates_tmp[i];
        }



#ifdef HAVE_ZOLTAN2_OMP
#ifdef USE_LESS_THREADS
#pragma omp barrier
#endif
#pragma omp single
#endif
        {
          cutCoordinatesWork = cutCoordinates_tmp;
        }
      }
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
void getChunksFromCoordinates(partId_t partNo, int noThreads,
    scalar_t *pqJagged_coordinates, scalar_t *cutCoordinates, lno_t *totalCounts,
    lno_t *partitionedPointPermutations,
    lno_t *newpartitionedPointPermutations, lno_t coordinateBegin, lno_t coordinateEnd,
    lno_t *coordinate_linked_list, lno_t **coordinate_starts, lno_t **coordinate_ends, lno_t numLocalCoord, float *actual_ratios, bool allowNonRectelinearPart,
    scalar_t *totalPartWeights, scalar_t *coordWeights, bool pqJagged_uniformWeights, int myRank, int worldSize, double **partWeights, float **nonRectelinearRatios,
    lno_t ** partPointCounts){

  lno_t numCoordsInPart =  coordinateEnd - coordinateBegin;
  //++myRank;

  partId_t noCuts = partNo - 1;


  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
  {
    int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
    me = omp_get_thread_num();
#endif

    lno_t *myStarts = coordinate_starts[me];
    lno_t *myEnds = coordinate_ends[me];
    lno_t *myPartPointCounts = partPointCounts[me];
    float *myRatios = NULL;
    if (allowNonRectelinearPart){


      myRatios = nonRectelinearRatios[me];
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for (partId_t i = 0; i < noCuts; ++i){
        float r = actual_ratios[i];
        scalar_t leftWeight = r * (totalPartWeights[i * 2 + 1] - totalPartWeights[i * 2]);
        for(int ii = 0; ii < noThreads; ++ii){
          if(leftWeight > _EPSILON){
            scalar_t ithWeight = partWeights[ii][i * 2 + 1] - partWeights[ii][i * 2 ];
            if(ithWeight < leftWeight){
              nonRectelinearRatios[ii][i] = ithWeight;
            }
            else {
              nonRectelinearRatios[ii][i] = leftWeight ;
            }
            leftWeight -= ithWeight;
          }
          else {
            nonRectelinearRatios[ii][i] = 0;
          }
        }

        //cout << "m:" << i << ":" <<  actual_ratios[i] << endl;
        //cout << "m:" << i << ":" <<  myRatios[i]  << endl;

      }


      if(noCuts > 0){
        for (partId_t i = noCuts - 1; i > 0 ; --i){
          if(fabs(cutCoordinates[i] - cutCoordinates[i -1]) < _EPSILON){
            myRatios[i] -= myRatios[i - 1] ;
          }
          myRatios[i] = int ((myRatios[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL) / scalar_t(SIGNIFICANCE_MUL);

          //cout << "m:" << i << ":" <<  myRatios[i] << endl;
        }
      }

    }


    pqJagged_PartVertices <lno_t, partId_t> pqPV;
    pqPV.set(coordinate_linked_list, myStarts, myEnds);
    memset(myPartPointCounts, 0, sizeof(lno_t) * partNo);


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t i = 0; i < numLocalCoord; ++i){
      coordinate_linked_list[i] = -1;
    }

    for (partId_t i = 0; i < partNo; ++i){
      myEnds[i] = -1;
      myStarts[i] = -1;
    }


    //determine part of each point

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

      lno_t i = partitionedPointPermutations[ii];

      bool inserted = false;
      for(partId_t j = 0; j < noCuts; ++j){
        scalar_t distance = pqJagged_coordinates[i] - cutCoordinates[j];
        scalar_t absdistance = fabs(distance);

        if (allowNonRectelinearPart && myRatios[j] > _EPSILON * EPS_SCALE * 10 && absdistance < _EPSILON ){
          /*
          if(j == 5 || j == 6 || j ==7 ){
            cout << "myRatios["<<j<<"]:" <<  myRatios[j] << endl;
          }
          */
          scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];
          scalar_t decrease = w;
          myRatios[j] -= decrease;

          if(myRatios[j] < 0 && j < noCuts - 1 && fabs(cutCoordinates[j+1] - cutCoordinates[j]) < _EPSILON){
            myRatios[j + 1] += myRatios[j];
          }
          inserted = true;
          pqPV.inserToPart(j, i);
          break;
        }
        else if (distance < 0 && absdistance > _EPSILON){
          pqPV.inserToPart(j, i);
          inserted = true;
          break;
        }
        else {
          myPartPointCounts[j] -=	 1;
        }

      }
      if(!inserted){
        pqPV.inserToPart(noCuts, i);
      }
    }

    //accumulate the starts of each thread.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(partId_t i = 0; i < partNo; ++i){
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
    for(partId_t j = 0; j < partNo; ++j){
      lno_t pwj = 0;
      for (int i = 0; i < noThreads; ++i){
        pwj += partPointCounts[i][j];
      }
      totalCounts[j] = pwj + numCoordsInPart;
    }



#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(partId_t i = 0; i < partNo; ++i){
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

}

template <typename T>
T *allocMemory(size_t size){
  if (size > 0){
    T * a = new T[size];
    if (a == NULL) {
      throw  "cannot allocate memory";
    }
    return a;
  }
  else {
    return NULL;
  }
}

template <typename T>
void freeArray(T *&array){
  if(array != NULL){
    delete [] array;
    array = NULL;
  }
}

template <typename tt>
std::string toString(tt obj){
  std::stringstream ss (std::stringstream::in |std::stringstream::out);
  ss << obj;
  std::string tmp = "";
  ss >> tmp;
  return tmp;
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

  env->timerStart(MACRO_TIMERS, "PQJagged Total");


  env->timerStart(MACRO_TIMERS, "PQJagged Problem_Init");
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;

  typedef typename Adapter::lno_t lno_t;
  const Teuchos::ParameterList &pl = env->getParameters();

  std::bitset<NUM_RCB_PARAMS> params;
  int numTestCuts = 5;
  scalar_t imbalanceTolerance;

  multiCriteriaNorm mcnorm;
  bool ignoreWeights;

  bool allowNonRectelinearPart = true;
  pqJagged_getParameters<scalar_t>(pl, imbalanceTolerance, mcnorm, params, numTestCuts, ignoreWeights,allowNonRectelinearPart);

  int coordDim, weightDim; size_t nlc; global_size_t gnc; int criteriaDim;
  pqJagged_getCoordinateValues<Adapter>( coords, coordDim, weightDim, nlc, gnc, criteriaDim, ignoreWeights);
  lno_t numLocalCoords = nlc;
  gno_t numGlobalCoords = gnc;


  scalar_t **pqJagged_coordinates = allocMemory<scalar_t *>(coordDim);

  scalar_t **pqJagged_weights = allocMemory<scalar_t *>(criteriaDim);

  bool *pqJagged_uniformParts = allocMemory< bool >(criteriaDim);

  scalar_t **pqJagged_partSizes =  allocMemory<scalar_t *>(criteriaDim);

  bool *pqJagged_uniformWeights = allocMemory< bool >(criteriaDim);

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

  int numThreads = pqJagged_getNumThreads();

  if(comm->getRank() == 0){
    cout << "numGlobalParts:" << numGlobalParts << endl;
    cout << "numGlobalCoords:" << numGlobalCoords << endl;
    cout << "numThreads=" << numThreads << endl;
    cout << "numLocalCoord:" << numLocalCoords << endl;
  }
  scalar_t minCoordinate, maxCoordinate;


  partId_t totalDimensionCut = 0;
  partId_t totalPartCount = 1;
  partId_t maxPartNo = 0;

  const partId_t *partNo = pl.getPtr<Array <partId_t> >("pqParts")->getRawPtr();
  for (int i = 0; i < coordDim; ++i){
    totalPartCount *= partNo[i];
    if(partNo[i] > maxPartNo) maxPartNo = partNo[i];
  }
  totalDimensionCut = totalPartCount - 1;

  partId_t maxTotalCumulativePartCount = totalPartCount / partNo[coordDim];


  scalar_t *global_min_max_coordinates = allocMemory< scalar_t>(2 * maxTotalCumulativePartCount);
  scalar_t *local_min_max_coordinates = allocMemory< scalar_t>(2 * maxTotalCumulativePartCount);


  // coordinates of the cut lines. First one is the min, last one is max coordinate.
  scalar_t *allCutCoordinates = allocMemory< scalar_t>(totalDimensionCut);

  for (partId_t i = 0; i < totalDimensionCut; ++i){
    allCutCoordinates[i] = 0;
  }

  lno_t *partitionedPointCoordinates =  allocMemory< lno_t>(numLocalCoords);

  lno_t *newpartitionedPointCoordinates = allocMemory< lno_t>(numLocalCoords);

  scalar_t *max_min_array =  allocMemory< scalar_t>(numThreads * 2);


  //initial configuration
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
  for(lno_t i = 0; i < numLocalCoords; ++i){
    //set each pointer-i to i.
    partitionedPointCoordinates[i] = i;
  }

#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
  firstTouch<lno_t>(newpartitionedPointCoordinates, numLocalCoords);
#endif
#endif

  //initially there is a single partition
  lno_t currentPartitionCount = 1;

  //single partition starts at index-0, and ends at numLocalCoords
  lno_t coordinateBegin = 0;
  lno_t coordinateEnd = numLocalCoords;

  //inTotalCounts array holds the end points in partitionedPointCoordinates array
  //for each partition. Initially sized 1, and single element is set to numLocalCoords.
  lno_t *inTotalCounts = allocMemory<lno_t>(1);
  inTotalCounts[0] = numLocalCoords;

  //the ends points of the output.
  lno_t *outTotalCounts = NULL;

  //the array holding the points in each part as linked list
  lno_t *coordinate_linked_list = allocMemory<lno_t>(numLocalCoords);

  //the start and end coordinate of  each part.
  lno_t **coordinate_starts =  allocMemory<lno_t *>(numThreads);

  lno_t **coordinate_ends = allocMemory<lno_t *>(numThreads);

  //assign the max size to starts, as it will be reused.
  for(int i = 0; i < numThreads; ++i){
    coordinate_starts[i] = allocMemory<lno_t>(maxPartNo);
    coordinate_ends[i] = allocMemory<lno_t>(maxPartNo);
  }

#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
#pragma omp parallel
  {
    int me = omp_get_thread_num();
    for(partId_t ii = 0; ii < maxPartNo; ++ii){
      coordinate_starts[me][ii] = 0;
      coordinate_ends[me][ii] = 0;
    }
  }
#endif
#endif

  partId_t maxCutNo = maxPartNo - 1;


  float *nonRectelinearPart = NULL;
  float **nonRectRatios = NULL;

  if(allowNonRectelinearPart){
    nonRectelinearPart = allocMemory<float>(maxCutNo);
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<float>(nonRectelinearPart, maxCutNo);
#endif
#endif

    nonRectRatios = allocMemory<float *>(numThreads);

    for(int i = 0; i < numThreads; ++i){
      nonRectRatios[i] = allocMemory<float>(maxCutNo);
    }
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
#pragma omp parallel
  {
    int me = omp_get_thread_num();
    for(partId_t ii = 0; ii < maxCutNo; ++ii){
      nonRectRatios[me][ii] = 0;
    }
  }
#endif
#endif

  }

  scalar_t *cutCoordinatesWork = allocMemory<scalar_t>(maxCutNo); // work array to manipulate coordinate of cutlines in different iterations.
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<scalar_t>(cutCoordinatesWork, maxCutNo);
#endif
#endif

  scalar_t *cutPartRatios = allocMemory<scalar_t>(maxPartNo); // the weight ratios at left side of the cuts. First is 0, last is 1.
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<scalar_t>(cutPartRatios, maxCutNo);
#endif
#endif
  scalar_t *cutUpperBounds = allocMemory<scalar_t>(maxCutNo);  //to determine the next cut line with binary search

  scalar_t *cutLowerBounds = allocMemory<scalar_t>(maxCutNo);  //to determine the next cut line with binary search

  scalar_t *cutLowerWeight = allocMemory<scalar_t>(maxCutNo);  //to determine the next cut line with binary search

  scalar_t *cutUpperWeight = allocMemory<scalar_t>(maxCutNo);  //to determine the next cut line with binary search

  bool *isDone = allocMemory<bool>(maxCutNo);

  //scalar_t **partWeights = allocMemory<scalar_t *>(numThreads);
  double **partWeights = allocMemory<double *>(numThreads);


  scalar_t **leftClosestDistance = allocMemory<scalar_t *>(numThreads);

  scalar_t **rightClosestDistance = allocMemory<scalar_t *>(numThreads);

  size_t maxTotalPartCount = maxPartNo + size_t(maxCutNo);


  lno_t **partPointCounts = allocMemory<lno_t *>(numThreads);

  for(int i = 0; i < numThreads; ++i){
    //partWeights[i] = allocMemory<scalar_t>(maxTotalPartCount);
    partWeights[i] = allocMemory < double >(maxTotalPartCount);
    rightClosestDistance[i] = allocMemory<scalar_t>(maxCutNo);
    leftClosestDistance[i] = allocMemory<scalar_t>(maxCutNo);
    partPointCounts[i] =  allocMemory<lno_t>(maxPartNo);
  }
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
#pragma omp parallel
  {
    int me = omp_get_thread_num();
    for(partId_t ii = 0; ii < maxPartNo; ++ii){
      partPointCounts[me][ii] = 0;
    }
    for(size_t ii = 0; ii < maxTotalPartCount; ++ii){
      partWeights[me][ii] = 0;
    }
    for(partId_t ii = 0; ii < maxCutNo; ++ii){
      rightClosestDistance[me][ii] = 0;
      leftClosestDistance[me][ii] = 0;

    }
  }
#endif
#endif

  scalar_t *cutWeights = allocMemory<scalar_t>(maxCutNo);
  scalar_t *globalCutWeights = allocMemory<scalar_t>(maxCutNo);

  //for faster communication concatanation.
  scalar_t *totalPartWeights_leftClosests_rightClosests = allocMemory<scalar_t>(maxTotalPartCount + maxCutNo * 2);
  scalar_t *global_totalPartWeights_leftClosests_rightClosests = allocMemory<scalar_t>(maxTotalPartCount + maxCutNo * 2);

  scalar_t *cutCoordinates =  allCutCoordinates;

  partId_t leftPartitions = totalPartCount;
  scalar_t maxScalar_t = numeric_limits<float>::max();
  scalar_t minScalar_t = -numeric_limits<float>::max();

  env->timerStop(MACRO_TIMERS, "PQJagged Problem_Init");

  env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning");
  for (int i = 0; i < coordDim; ++i){
    if(partNo[i] == 1) continue;
    env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i));
    lno_t partitionCoordinateBegin = 0;
    outTotalCounts = allocMemory<lno_t>(currentPartitionCount * partNo[i]);

    size_t currentOut = 0;
    size_t currentIn = 0;
    size_t previousEnd = 0;

    leftPartitions /= partNo[i];

#ifdef FIRST_TOUCH_2
    lno_t *pc =  allocMemory< lno_t>(numLocalCoords);
    lno_t *npc = allocMemory< lno_t>(numLocalCoords);
    scalar_t * pqCoord = allocMemory< scalar_t>(numLocalCoords);
#ifdef HAVE_ZOLTAN2_OMP
    for (int j = 0; j < currentPartitionCount; ++j){

      coordinateEnd= inTotalCounts[currentIn + j];
      coordinateBegin = currentIn==0 ? 0: inTotalCounts[currentIn + j -1];
#pragma omp parallel for
      for(lno_t jj = coordinateBegin; jj < coordinateEnd; ++jj){
        lno_t ind = pc[jj] = partitionedPointCoordinates[jj];
        pqCoord[ind] = pqJagged_coordinates[i][ind];
        npc[jj] = 0;
      }
    }
    delete [] partitionedPointCoordinates;
    delete [] newpartitionedPointCoordinates;
#endif
#endif

#ifndef FIRST_TOUCH_2
    scalar_t * pqCoord = pqJagged_coordinates[i];
    lno_t *pc = partitionedPointCoordinates;
    lno_t *npc = newpartitionedPointCoordinates;
#endif


    env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) +"_min_max");
    pqJagged_getMinMaxCoordofAllParts(comm, pqCoord, global_min_max_coordinates,
        env, numThreads,
        pc, inTotalCounts,
        max_min_array ,
        maxScalar_t, minScalar_t, currentPartitionCount, local_min_max_coordinates );
    env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) +"_min_max");


    for (int j = 0; j < currentPartitionCount; ++j, ++currentIn){

      if(comm->getRank() == 0)
        cout << "i: " << i << " j:" << j << " ";
      scalar_t used_imbalance = imbalanceTolerance * (LEAF_IMBALANCE_FACTOR + (1 - LEAF_IMBALANCE_FACTOR)   / leftPartitions) * 0.7;
      //scalar_t used_imbalance = 0;

      coordinateEnd= inTotalCounts[currentIn];
      coordinateBegin = currentIn==0 ? 0: inTotalCounts[currentIn -1];
      minCoordinate = global_min_max_coordinates[j];
      maxCoordinate = global_min_max_coordinates[j + currentPartitionCount];

      //cout << "beg:" << coordinateBegin << " end:" << coordinateEnd << endl;
/*
      env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) +"_min_max");


      pqJagged_getMinMaxCoord<scalar_t, lno_t>(comm, pqCoord, minCoordinate,maxCoordinate,
          env, numThreads, pc, coordinateBegin, coordinateEnd, max_min_array, maxScalar_t, minScalar_t);

      env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) +"_min_max");

*/
      //cout << "min:" << minCoordinate << " max:" << maxCoordinate << endl;
      //cout << "min:" << minCoordinate << " max:" << maxCoordinate << endl;

      for (partId_t ii = 0;ii < partNo[i]; ++ii){
        outTotalCounts[currentOut + ii] = 0;
      }

      if(minCoordinate <= maxCoordinate /*&& partNo[i] > 1*/){

        env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) + "_cut_coord");
        pqJagged_getCutCoord_Weights<scalar_t>(
            minCoordinate, maxCoordinate,
            pqJagged_uniformParts[0], pqJagged_partSizes[0], partNo[i] - 1,
            cutCoordinates, cutPartRatios, numThreads
        );


        env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) + "_cut_coord");

        scalar_t globalTotalWeight = 0;

        env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) + "_1d");
        pqJagged_1DPart_simple<scalar_t, lno_t>(env,comm,pqCoord, pqJagged_weights[0], pqJagged_uniformWeights[0],
            minCoordinate,
            maxCoordinate, partNo[i], numThreads,
            used_imbalance, cutCoordinates,
            pc, coordinateBegin, coordinateEnd,
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
            totalPartWeights_leftClosests_rightClosests,
            global_totalPartWeights_leftClosests_rightClosests,
            globalTotalWeight,
            nonRectelinearPart,
            allowNonRectelinearPart,
            maxScalar_t,
            minScalar_t,
            cutWeights, globalCutWeights
        );

        env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) + "_1d");
        env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) + "_chunks");
        getChunksFromCoordinates<lno_t,scalar_t>(partNo[i], numThreads,
            pqCoord, cutCoordinates, outTotalCounts + currentOut,
            pc, npc, coordinateBegin, coordinateEnd,
            coordinate_linked_list, coordinate_starts, coordinate_ends, numLocalCoords,
            nonRectelinearPart, allowNonRectelinearPart, totalPartWeights_leftClosests_rightClosests,
            pqJagged_weights[0], pqJagged_uniformWeights, comm->getRank(), comm->getSize(), partWeights,nonRectRatios,
            partPointCounts);
        env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i) + "_chunks");
      }
      /*
      else {
#pragma omp parallel for
        for(lno_t iii = coordinateBegin; iii < coordinateEnd; ++iii){
          newpartitionedPointCoordinates[iii] = partitionedPointCoordinates[iii];
        }
        outTotalCounts[currentOut] = coordinateEnd - coordinateBegin;
      }
       */



      cutCoordinates += partNo[i] - 1;
      partitionCoordinateBegin += coordinateBegin - coordinateEnd;

      //cout << "current out:" << currentOut << endl;
      for (partId_t ii = 0;ii < partNo[i]; ++ii){

        //cout << "total was:" << outTotalCounts[currentOut + ii] << endl;
        outTotalCounts[currentOut+ii] += previousEnd;
        //cout << "total now:" << outTotalCounts[currentOut + ii] << endl;
      }
      previousEnd = outTotalCounts[currentOut + partNo[i] - 1];
      currentOut += partNo[i];

    }

    //lno_t * tmp = npc;
    //npc = pc;
    //pc = tmp;
    partitionedPointCoordinates = npc;
    newpartitionedPointCoordinates = pc;
    npc = NULL;
    pc = NULL;

    currentPartitionCount *= partNo[i];
    //delete [] inTotalCounts;
    freeArray<lno_t>(inTotalCounts);
    inTotalCounts = outTotalCounts;
#ifdef FIRST_TOUCH_2
    delete [] pqCoord;
    //delete []pc;
    //delete []npc;
#endif

    env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + toString<int>(i));
  }

  partId_t *partIds = NULL;
  ArrayRCP<partId_t> partId;
  if(numLocalCoords > 0){
  /*partId_t **/partIds = allocMemory<partId_t>(numLocalCoords);
  /*ArrayRCP<partId_t> */partId = arcp(partIds, 0, numLocalCoords, true);
  }

  for(partId_t i = 0; i < totalPartCount;++i){
    lno_t begin = 0;
    lno_t end = inTotalCounts[i];
    if(i > 0) begin = inTotalCounts[i -1];
#pragma omp parallel for
    for (lno_t ii = begin; ii < end; ++ii){

      lno_t k = partitionedPointCoordinates[ii];
      partIds[k] = i;
      /*
      cout << "part of coordinate:";
      for(int iii = 0; iii < coordDim; ++iii){
        cout <<  pqJagged_coordinates[iii][k] << " ";
      }
      cout << i;
      cout << endl;
      */
    }
    //cout << "begin:" << begin << " end:" << end << endl;
  }
  ArrayRCP<const gno_t> gnoList;
  if(numLocalCoords > 0){
    gnoList = arcpFromArrayView(pqJagged_gnos);
  }

  solution->setParts(gnoList, partId);

  env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning");

  env->timerStart(MACRO_TIMERS, "PQJagged Problem_Free");

/*
  if(comm->getRank() == 0){
    for(size_t i = 0; i < totalPartCount - 1;++i){
      cout << "cut coordinate:" << allCutCoordinates[i] << endl;
    }
  }
*/



  freeArray<scalar_t> (global_min_max_coordinates);
  freeArray<scalar_t> (local_min_max_coordinates);

  for(int i = 0; i < numThreads; ++i){
    //delete []coordinate_starts[i];
    freeArray<lno_t>(coordinate_starts[i]);
    //delete []coordinate_ends[i] ;
    freeArray<lno_t>(coordinate_ends[i]);
    //delete [] partPointCounts[i];
    freeArray<lno_t>(partPointCounts[i]);
  }

  //delete []partPointCounts;
  freeArray<lno_t *>(partPointCounts);

  //delete [] imbalance_tolerances;
  //delete []coordinate_linked_list;
  freeArray<lno_t>(coordinate_linked_list);
  //the start and end coordinate of  each part.
  //delete []coordinate_starts;
  freeArray<lno_t *>(coordinate_starts);

  //delete []coordinate_ends;
  freeArray<lno_t *>(coordinate_ends);

  //assign the max size to starts, as it will be reused.
  if(allowNonRectelinearPart){
    //delete []nonRectelinearPart;
    freeArray<float>(nonRectelinearPart);

    for(int i = 0; i < numThreads; ++i){
      //delete [] nonRectRatios[i];
      freeArray<float>(nonRectRatios[i]);
    }
    //delete [] nonRectRatios;
    freeArray<float *>(nonRectRatios);
  }

  freeArray<scalar_t>(cutWeights);
  freeArray<scalar_t>(globalCutWeights);


  //delete [] partNo;
  //delete []max_min_array;
  freeArray<scalar_t>(max_min_array);

  //delete [] outTotalCounts;
  freeArray<lno_t>(outTotalCounts);

  //delete []partitionedPointCoordinates ;
  freeArray<lno_t>(partitionedPointCoordinates);

  //delete []newpartitionedPointCoordinates ;
  freeArray<lno_t>(newpartitionedPointCoordinates);

  //delete []allCutCoordinates;
  freeArray<scalar_t>(allCutCoordinates);

  //delete []pqJagged_coordinates;
  freeArray<scalar_t *>(pqJagged_coordinates);

  //delete []pqJagged_weights;
  freeArray<scalar_t *>(pqJagged_weights);

  //delete []pqJagged_uniformParts;
  freeArray<bool>(pqJagged_uniformParts);


  //delete []pqJagged_partSizes;
  freeArray<scalar_t *>(pqJagged_partSizes);

  //delete []pqJagged_uniformWeights;
  freeArray<bool>(pqJagged_uniformWeights);


  //delete []cutCoordinatesWork; // work array to manipulate coordinate of cutlines in different iterations.
  freeArray<scalar_t>(cutCoordinatesWork);

  //delete []cutPartRatios; // the weight ratios at left side of the cuts. First is 0, last is 1.
  freeArray<scalar_t>(cutPartRatios);

  //delete []cutUpperBounds;  //to determine the next cut line with binary search
  freeArray<scalar_t>(cutUpperBounds);

  //delete []cutLowerBounds;  //to determine the next cut line with binary search
  freeArray<scalar_t>(cutLowerBounds);

  //delete []cutLowerWeight;  //to determine the next cut line with binary search
  freeArray<scalar_t>(cutLowerWeight);

  //delete []cutUpperWeight;  //to determine the next cut line with binary search
  freeArray<scalar_t>(cutUpperWeight);

  //delete []isDone;
  freeArray<bool>(isDone);

  //delete []totalPartWeights;
  freeArray<scalar_t>(totalPartWeights_leftClosests_rightClosests);
  freeArray<scalar_t>(global_totalPartWeights_leftClosests_rightClosests);

  for(int i = 0; i < numThreads; ++i){
    //delete [] partWeights[i] ;
    freeArray<double>(partWeights[i]);//freeArray<scalar_t>(partWeights[i]);
    //delete [] rightClosestDistance[i];
    freeArray<scalar_t>(rightClosestDistance[i]);
    //delete [] leftClosestDistance[i];
    freeArray<scalar_t>(leftClosestDistance[i]);
  }

  //delete [] partWeights;
  //freeArray<scalar_t *>(partWeights);
  freeArray<double *>(partWeights);

  //delete [] leftClosestDistance ;
  freeArray<scalar_t *>(leftClosestDistance);

  //delete [] rightClosestDistance;
  freeArray<scalar_t *>(rightClosestDistance);


  env->timerStop(MACRO_TIMERS, "PQJagged Problem_Free");
  env->timerStop(MACRO_TIMERS, "PQJagged Total");
}
} // namespace Zoltan2





#endif





