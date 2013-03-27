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

/*! \file Zoltan2_AlgPQJagged.hpp
\brief Contains the PQ-jagged algorthm.
 */

#ifndef _ZOLTAN2_ALGPQJagged_HPP_
#define _ZOLTAN2_ALGPQJagged_HPP_

#include <Zoltan2_AlgRCB_methods.hpp>
#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_Metric.hpp>             // won't need thiss

#include <Teuchos_ParameterList.hpp>

#include <bitset>

#define EPS_SCALE 1
#define LEAST_SIGNIFICANCE 0.0001
#define SIGNIFICANCE_MUL 1000
//#define INCLUDE_ZOLTAN2_EXPERIMENTAL
#ifdef HAVE_ZOLTAN2_OMP
#include <omp.h>
#endif
//#define FIRST_TOUCH

//#define BINARYCUTSEARCH

#define ABS(x) ((x) >= 0 ? (x) : -(x))

#define LEAF_IMBALANCE_FACTOR 0.1
#define BINARYCUTOFF 32
//imbalance calculation. Wreal / Wexpected - 1
#define imbalanceOf(Wachieved, totalW, expectedRatio) \
    (Wachieved) / ((totalW) * (expectedRatio)) - 1
//#define mpi_communication

#define KCUTOFF 0.80
#define defaultK 16

namespace Teuchos{
template <typename Ordinal, typename T>
class PQJaggedCombinedReductionOp  : public ValueTypeReductionOp<Ordinal,T>
{
private:
  Ordinal numSum_, numMin_1, numMin_2;
  Ordinal k;

public:
  /*! \brief Default Constructor
   */
  PQJaggedCombinedReductionOp ():numSum_(0), numMin_1(0), numMin_2(0), k(0){}

  /*! \brief Constructor
   *   \param nsum  the count of how many sums will be computed at the
   *             start of the list.
   *   \param nmin  following the sums, this many minimums will be computed.
   *   \param nmax  following the minimums, this many maximums will be computed.
   */
  PQJaggedCombinedReductionOp (Ordinal nsum, Ordinal nmin1, Ordinal nmin2, Ordinal k_):
    numSum_(nsum), numMin_1(nmin1), numMin_2(nmin2), k(k_){}

  /*! \brief Implement Teuchos::ValueTypeReductionOp interface
   */
  void reduce( const Ordinal count, const T inBuffer[], T inoutBuffer[]) const
  {
    Ordinal next=0;
    for(Ordinal ii = 0; ii < k ; ++ii){
      for (Ordinal i=0; i < numSum_; i++, next++)
        inoutBuffer[next] += inBuffer[next];

      for (Ordinal i=0; i < numMin_1; i++, next++)
        if (inoutBuffer[next] > inBuffer[next])
          inoutBuffer[next] = inBuffer[next];

      for (Ordinal i=0; i < numMin_2; i++, next++)
        if (inoutBuffer[next] > inBuffer[next])
          inoutBuffer[next] = inBuffer[next];
    }
  }
};








template <typename Ordinal, typename T>
class PQJaggedCombinedMinMaxTotalReductionOp  : public ValueTypeReductionOp<Ordinal,T>
{
private:
  Ordinal numMin, numMax, numTotal;

public:
  /*! \brief Default Constructor
   */
  PQJaggedCombinedMinMaxTotalReductionOp ():numMin(0), numMax(0), numTotal(0){}

  /*! \brief Constructor
   *   \param nsum  the count of how many sums will be computed at the
   *             start of the list.
   *   \param nmin  following the sums, this many minimums will be computed.
   *   \param nmax  following the minimums, this many maximums will be computed.
   */
  PQJaggedCombinedMinMaxTotalReductionOp (Ordinal nmin, Ordinal nmax, Ordinal nTotal):
    numMin(nmin), numMax(nmax), numTotal(nTotal){}

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


    for (Ordinal i=0; i < numTotal; i++, next++)
        inoutBuffer[next] += inBuffer[next];
  }
};
} // namespace Teuchos

namespace Zoltan2{

//diffclock for temporary timing experiments.

#ifdef mpi_communication
  partId_t concurrent = 0;
  void sumMinMin(void *in, void *inout, int *count, MPI_Datatype *type) {

    int k = *count;
    int numCut = (k  - 1) / 4;
    int total_part_count = numCut * 2 + 1;
    int next = 0;

    float *inoutBuffer = (float *) inout;
    float *inBuffer = (float *) in;


    for(partId_t ii = 0; ii < concurrent ; ++ii){
    for (long i=0; i < total_part_count; i++, next++)
      inoutBuffer[next] += inBuffer[next];

    for (long i=0; i < numCut; i++, next++)
      if (inoutBuffer[next] > inBuffer[next])
        inoutBuffer[next] = inBuffer[next];

    for (long i=0; i < numCut; i++, next++)
      if (inoutBuffer[next] > inBuffer[next])
        inoutBuffer[next] = inBuffer[next];
    }
  }


  void minMaxSum(void *in, void *inout, int *count, MPI_Datatype *type) {
//    long k = (*((int *) count));
    int k = *count;
    int num = k / 3;
    int next = 0;
    float *inoutBuffer = (float *) inout;
    float *inBuffer = (float *) in;


    for (long i=0; i < num; i++, next++)
      if (inoutBuffer[next] > inBuffer[next])
        inoutBuffer[next] = inBuffer[next];
    for (long i=0; i < num; i++, next++)
      if (inoutBuffer[next] < inBuffer[next])
        inoutBuffer[next] = inBuffer[next];
    for (long i=0; i < num; i++, next++)
      inoutBuffer[next] += inBuffer[next];

  }
#endif

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
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
  for(size_t jj = 0; jj < arraySize; ++jj){
    arrayName[jj] = 0;
  }
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
 * \param concurrentPartCount is the number of parts whose cut lines will be calculated concurrently.
 */
template <typename T>
void pqJagged_getParameters(const Teuchos::ParameterList &pl, T &imbalanceTolerance,
    multiCriteriaNorm &mcnorm, std::bitset<NUM_RCB_PARAMS> &params,  int &numTestCuts, bool &ignoreWeights, bool &allowNonrectilinear, partId_t &concurrentPartCount,
    bool &force_binary, bool &force_linear){

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

  imbalanceTolerance = .1;
  pe = pl.getEntryPtr("imbalance_tolerance");
  if (pe){
    double tol;
    tol = pe->getValue(&tol);
    imbalanceTolerance = tol - 1.0;
  }

  if (imbalanceTolerance <= 0)
    imbalanceTolerance = 10e-4;

  force_binary = false;
  pe = pl.getEntryPtr("force_binary_search");
  if (pe){
    int val = 0;
    val = pe->getValue(&val);
    if (val == 1)
      force_binary = true;
  }

  force_linear = false;
  pe = pl.getEntryPtr("force_linear_search");
  if (pe){
    int val;
    val = pe->getValue(&val);
    if (val == 1)
      force_linear = true;
  }

  //TODO: FIX ME.
  //double aa = 1;
  pe = pl.getEntryPtr("parallel_part_calculation_count");
  if (pe){
    //aa = pe->getValue(&aa);
    concurrentPartCount = pe->getValue(&concurrentPartCount);
  }else {
    concurrentPartCount = 1;
  //concurrentPartCount = partId_t(aa);
  }
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

  if (val == 1){
    params.set(rcb_rectilinearBlocks);
    allowNonrectilinear = false;
  } else {
    allowNonrectilinear = true;
  }

  numTestCuts = 1;
  pe = pl.getEntryPtr("bisection_num_test_cuts");
  if (pe)
    numTestCuts = pe->getValue(&numTestCuts);

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
    pqJagged_weights[0] = NULL;
  }
  else{
    for (int wdim = 0; wdim < weightDim; wdim++){
      if (wgts[wdim].size() == 0){
        pqJagged_uniformWeights[wdim] = true;
        pqJagged_weights[wdim] = NULL;
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
      pqJagged_partSizes[wdim] = NULL;
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
  for(size_t i = 0; i < numLocalCoords; ++i){
    for (int ii = 0; ii < coordDim; ++ii){
      std::cout <<  pqJagged_values[ii][i] << " ";
    }
    std::cout << std::endl;
  }


  std::cout << "criteriaDim:" << criteriaDim << std::endl;
  std::cout << "weightDim:" << weightDim << std::endl;
  if(weightDim){
    for(size_t i = 0; i < numLocalCoords; ++i){
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
  for(size_t i = 0; i < numLocalCoords; ++i){
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
      for(size_t ii = 0; ii < numGlobalParts; ++ii){
        std::cout << pqJagged_partSizes[i][ii] << " ";
      }
    std::cout << std::endl;
  }
}


/*! \brief Function to determine the local minimum and maximum coordinate, and local total weight
 * in the given set of local points.
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param pqJagged_coordinates float-like array representing the coordinates in a single dimension. Sized as numLocalCoords.
 * \param pqJagged_uniformWeights boolean value whether each coordinate value has a uniform weight or not. If true, pqJagged_weights is not used.
 * \param pqJagged_weights float-like array representing the weights of the coordinates in a single criteria dimension. Sized as numLocalCoords.
 *
 * \param numThreads is the integer value to represent the number of threads available for each processor.
 * \param coordinateBegin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinateEnd is the end index of the given partition on partitionedPointPermutations.
 *
 * \param max_min_array provided work array sized numThreads * 2.
 * \param maxScalar to initialize minimum coordinate if there are no points in the given range.
 * \param minScalar to initialize maximum coordinate if there are no points in the given range.
 *
 * \param minCoordinate is the output to represent the local minimumCoordinate in  given range of coordinates.
 * \param maxCoordinate is the output to represent the local maximum coordinate in the given range of coordinates.
 * \param totalWeight is the output to represent the local total weight in the coordinate in the given range of coordinates.
 *
 */
template <typename scalar_t, typename lno_t>
void pqJagged_getLocalMinMaxTotalCoord(
    lno_t *partitionedPointPermutations,
    scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    scalar_t *pqJagged_weights,


    int numThreads,
    lno_t coordinateBegin,
    lno_t coordinateEnd,

    scalar_t *max_min_array /*sized nothreads * 2*/,
    scalar_t maxScalar,
    scalar_t minScalar,
    scalar_t &minCoordinate,
    scalar_t &maxCoordinate,
    scalar_t &totalWeight

    ){


  if(coordinateBegin >= coordinateEnd)
  {
    minCoordinate = maxScalar;
    maxCoordinate = minScalar;
    totalWeight = 0;
  }
  else {
    scalar_t mytotalWeight = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
    {
      if (pqJagged_uniformWeights) {
  #ifdef HAVE_ZOLTAN2_OMP
  #pragma omp single
  #endif
        {
          mytotalWeight = coordinateEnd - coordinateBegin;
        }
      } else {
  #ifdef HAVE_ZOLTAN2_OMP
  #pragma omp for reduction(+:mytotalWeight)
  #endif
        for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
          int i = partitionedPointPermutations[ii];
          mytotalWeight += pqJagged_weights[i];
        }
      }

      int myId = 0;
#ifdef HAVE_ZOLTAN2_OMP
      myId = omp_get_thread_num();
#endif
      scalar_t myMin, myMax;

      myMin=myMax
           =pqJagged_coordinates[partitionedPointPermutations[coordinateBegin]];

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(lno_t j = coordinateBegin + 1; j < coordinateEnd; ++j){
        int i = partitionedPointPermutations[j];
        if(pqJagged_coordinates[i] > myMax) myMax = pqJagged_coordinates[i];
        if(pqJagged_coordinates[i] < myMin) myMin = pqJagged_coordinates[i];
      }
      max_min_array[myId] = myMin;
      max_min_array[myId + numThreads] = myMax;

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
    totalWeight = mytotalWeight;
  }
}



/*! \brief Function that reduces global minimum and maximum coordinates with global total weight from given local arrays.
 * \param comm the communicator for the problem
 * \param env   library configuration and problem parameters
 * \param concurrentPartCount is the number of parts whose cut lines will be calculated concurrently.
 * \param localMinMaxTotal is the array holding local min and max coordinate values with local total weight.
 * First concurrentPartCount entries are minimums of the parts, next concurrentPartCount entries are max, and then the total weights.
 * \param localMinMaxTotal is the output array holding global min and global coordinate values with global total weight.
 * The structure is same as localMinMaxTotal.
 */
template <typename scalar_t>
void pqJagged_getGlobalMinMaxTotalCoord(
    RCP<Comm<int> > &comm,
    const RCP<const Environment> &env,
    partId_t concurrentPartCount,
    scalar_t *localMinMaxTotal,
    scalar_t *globalMinMaxTotal){


  //reduce min for first concurrentPartCount elements, reduce max for next concurrentPartCount elements,
  //reduce sum for the last concurrentPartCount elements.
  if(comm->getSize()  > 1){

#ifndef mpi_communication
  Teuchos::PQJaggedCombinedMinMaxTotalReductionOp<int, scalar_t> reductionOp(
     concurrentPartCount,
     concurrentPartCount,
     concurrentPartCount);
#endif

#ifdef mpi_communication
    MPI_Op myop;
    MPI_Op_create(minMaxSum, 0, &myop);   /* step 3 */
#endif

  try{


#ifdef mpi_communication

    MPI_Allreduce(localMinMaxTotal, globalMinMaxTotal, 3 * concurrentPartCount, MPI_FLOAT, myop,MPI_COMM_WORLD);
#endif
#ifndef mpi_communication
    reduceAll<int, scalar_t>(*comm, reductionOp,
        3 * concurrentPartCount, localMinMaxTotal, globalMinMaxTotal
      );
#endif

  }
  Z2_THROW_OUTSIDE_ERROR(*env)
  }
  else {
    partId_t s = 3 * concurrentPartCount;
    for (partId_t i = 0; i < s; ++i){
      globalMinMaxTotal[i] = localMinMaxTotal[i];
    }
  }
}


/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 * \param minCoordinate minimum coordinate in the range.
 * \param maxCoordinate maximum coordinate in the range.
 * \param pqJagged_uniformParts is a boolean value holding whether the desired partitioning is uniform.
 * \param pqJagged_partSizes holds the desired parts sizes if desired partitioning is not uniform.
 * \param cutCoordinates is the output array for the initial cut lines.
 * \param cutPartRatios is the output array holding the cumulative ratios of parts in current partitioning.
 * For partitioning to 4 uniformly, cutPartRatios will be (0.25, 0.5 , 0.75, 1).
 * \param numThreads hold the number of threads available per mpi.
 */
template <typename scalar_t>
void pqJagged_getCutCoord_Weights(
    scalar_t minCoordinate,
    scalar_t maxCoordinate,
    bool pqJagged_uniformParts,
    scalar_t *pqJagged_partSizes /*p sized, weight ratios of each part*/,
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
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(partId_t i = myStart; i < myEnd; ++i){
        cutPartRatios[i] =  uniform * (i + 1);
        cutCoordinates[i] = minCoordinate + slice * (i + 1);
      }
      cutPartRatios[noCuts] = 1;
    }


  }
  else {
    //TODO fix here!!
    cutPartRatios[0] = pqJagged_partSizes[0];
    cutCoordinates[0] = coordinateRange * cutPartRatios[0];
    for(partId_t i = 1; i < noCuts; ++i){
      cutPartRatios[i] = pqJagged_partSizes[i] + cutPartRatios[i - 1];
      cutCoordinates[i] = coordinateRange * cutPartRatios[i];
    }
  }
}


/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 *
 * \param env library configuration and problem parameters
 * \param comm the communicator for the problem
 * \param total_part_count is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
 * \param noCuts is the number of cut lines. P - 1.
 * \param maxCoordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
 * \param minCoordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
 * \param globalTotalWeight is the global total weight in the current range of coordinates.
 * \param imbalanceTolerance is the maximum allowed imbalance ratio.
 * \param maxScalar is the maximum value that scalar_t can represent.
 *
 *
 * \param globalPartWeights is the array holding the weight of parts. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param localPartWeights is the local totalweight of the processor.
 * \param targetPartWeightRatios are the desired cumulative part ratios, sized P.
 * \param isDone is the boolean array to determine if the correct position for a cut line is found.
 *
 * \param cutCoordinates is the array holding the coordinates of each cut line. Sized P - 1.
 * \param cutUpperBounds is the array holding the upper bound coordinate for each cut line. Sized P - 1.
 * \param cutLowerBounds is the array holding the lower bound coordinate for each cut line. Sized P - 1.
 * \param leftClosestDistance is the array holding the distances to the closest points to the cut lines from left.
 * \param rightClosestDistance is the array holding the distances to the closest points to the cut lines from right.
 * \param cutLowerWeight is the array holding the weight of the parts at the left of lower bound coordinates.
 * \param cutUpperWeight is the array holding the weight of the parts at the left of upper bound coordinates.
 * \param newCutCoordinates is the work array, sized P - 1.
 *
 * \param allowNonRectelinearPart is the boolean value whether partitioning should allow distributing the points on same coordinate to different parts.
 * \param nonRectelinearPartRatios holds how much percentage of the coordinates on the cutline should be put on left side.
 * \param rectilinearCutCount is the count of cut lines whose balance can be achived via distributing the points in same coordinate to different parts.
 * \param localCutWeights is the local weight of coordinates on cut lines.
 * \param globalCutWeights is the global weight of coordinates on cut lines.
 *
 * \param myNoneDoneCount is the number of cutlines whose position has not been determined yet. For K > 1 it is the count in a single part (whose cut lines are determined).
 */
template <typename scalar_t>
void getNewCoordinates(
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,
    const size_t &total_part_count,
    const partId_t &noCuts,
    const scalar_t &maxCoordinate,
    const scalar_t &minCoordinate,
    const scalar_t &globalTotalWeight,
    const scalar_t &imbalanceTolerance,
    scalar_t maxScalar,

    const scalar_t * globalPartWeights,
    const scalar_t * localPartWeights,
    const scalar_t *targetPartWeightRatios,
    bool *isDone,

    scalar_t *cutCoordinates,
    scalar_t *cutUpperBounds,
    scalar_t *cutLowerBounds,
    scalar_t *leftClosestDistance,
    scalar_t *rightClosestDistance,
    scalar_t * cutLowerWeight,
    scalar_t * cutUpperWeight,
    scalar_t *newCutCoordinates,

    bool allowNonRectelinearPart,
    float *nonRectelinearPartRatios,
    partId_t *rectilinearCutCount,
    scalar_t *localCutWeights,
    scalar_t *globalCutWeights,

    partId_t &myNoneDoneCount
)
    {


  scalar_t seenW = 0;
  float expected = 0;
  scalar_t leftImbalance = 0, rightImbalance = 0;
  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();
  //scalar_t _EPSILON = numeric_limits<float>::epsilon();

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
  for (partId_t i = 0; i < noCuts; i++){

    //if a left and right closes point is not found, set the distance to 0.
    if(leftClosestDistance[i] == maxScalar)
      leftClosestDistance[i] = 0;
    if(rightClosestDistance[i] == maxScalar)
      rightClosestDistance[i] = 0;

  }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
  for (partId_t i = 0; i < noCuts; i++){
    globalCutWeights[i] = 0;
    localCutWeights[i] = 0;
    //if already determined at previous iterations, do nothing.
    if(isDone[i]) {
      newCutCoordinates[i] = cutCoordinates[i];
      continue;
    }
    //current weight of the part at the left of the cut line.
    seenW = globalPartWeights[i * 2];


    //expected ratio
    expected = targetPartWeightRatios[i];
    leftImbalance = imbalanceOf(seenW, globalTotalWeight, expected);
    rightImbalance = imbalanceOf(globalTotalWeight - seenW, globalTotalWeight, 1 - expected);

    bool isLeftValid = ABS(leftImbalance) - imbalanceTolerance < _EPSILON ;
    bool isRightValid = ABS(rightImbalance) - imbalanceTolerance < _EPSILON;

    //if the cut line reaches to desired imbalance.
    if(isLeftValid && isRightValid){

      isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
      myNoneDoneCount -= 1;
      newCutCoordinates [i] = cutCoordinates[i];
      continue;
    }
    else if(leftImbalance < 0){

      scalar_t ew = globalTotalWeight * expected;
      if(allowNonRectelinearPart){

        if (globalPartWeights[i * 2 + 1] == ew){
          isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          myNoneDoneCount -= 1;
          newCutCoordinates [i] = cutCoordinates[i];
          nonRectelinearPartRatios[i] = 1;
          continue;
        }
        else if (globalPartWeights[i * 2 + 1] > ew){
          isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          *rectilinearCutCount += 1;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
          myNoneDoneCount -= 1;
          newCutCoordinates [i] = cutCoordinates[i];
          scalar_t myWeightOnLine = localPartWeights[i * 2 + 1] - localPartWeights[i * 2];
          localCutWeights[i] = myWeightOnLine;
          continue;
        }
      }
      //when moving right, set lower bound to current line.
      cutLowerBounds[i] = cutCoordinates[i] + rightClosestDistance[i];
      cutLowerWeight[i] = seenW;

      //compare the upper bound with the current lines.
      for (partId_t ii = i + 1; ii < noCuts ; ++ii){
        scalar_t pw = globalPartWeights[ii * 2];
        scalar_t lw = globalPartWeights[ii * 2 + 1];
        if(pw >= ew){
          if(pw == ew){
            cutUpperBounds[i] = cutCoordinates[ii];
            cutUpperWeight[i] = pw;
            cutLowerBounds[i] = cutCoordinates[ii];
            cutLowerWeight[i] = pw;
          } else if (pw < cutUpperWeight[i]){
            //if a cut line is more strict than the current upper bound,
            //update the upper bound.
            cutUpperBounds[i] = cutCoordinates[ii] - leftClosestDistance[ii];
            cutUpperWeight[i] = pw;
          }
          break;
        }
        //if comes here then pw < ew
        if(lw >= ew){
          cutUpperBounds[i] = cutCoordinates[ii];
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
      //if cut line does not move significantly.
      if (ABS(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE || cutUpperBounds[i] < cutLowerBounds[i]){
        isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
        myNoneDoneCount -= 1;
        newCutCoordinates [i] = cutCoordinates[i];
      } else {
        newCutCoordinates [i] = newPivot;
      }
    } else {
      //moving to left.
      scalar_t ew = globalTotalWeight * expected;
      //moving left, set upper to current line.
      cutUpperBounds[i] = cutCoordinates[i] - leftClosestDistance[i];
      cutUpperWeight[i] = seenW;

      // compare the current cut line weights with previous upper and lower bounds.
      for (int ii = i - 1; ii >= 0; --ii){
        scalar_t pw = globalPartWeights[ii * 2];
        scalar_t lw = globalPartWeights[ii * 2 + 1];
        if(pw <= ew){
          if(pw == ew){
            cutUpperBounds[i] = cutCoordinates[ii];
            cutUpperWeight[i] = pw;
            cutLowerBounds[i] = cutCoordinates[ii];
            cutLowerWeight[i] = pw;
          }
          else if (pw > cutLowerWeight[i]){
            cutLowerBounds[i] = cutCoordinates[ii] + rightClosestDistance[ii];
            cutLowerWeight[i] = pw;
            if(lw > ew){
              cutUpperBounds[i] = cutCoordinates[ii] + rightClosestDistance[ii];

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
      //if cut line does not move significantly.
      if (ABS(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE  || cutUpperBounds[i] < cutLowerBounds[i]){
        isDone[i] = true;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
        myNoneDoneCount -= 1;
        newCutCoordinates [ i] = cutCoordinates[i];
      } else {
        newCutCoordinates [ i] = newPivot;
      }
    }
  }



  //communication to determine the ratios of processors for the distribution of coordinates on the cut lines.
#ifdef HAVE_ZOLTAN2_OMP
//#pragma omp barrier
#pragma omp single
#endif
  {
    if(*rectilinearCutCount > 0){
      try{
        Teuchos::scan<int,scalar_t>(
            *comm, Teuchos::REDUCE_SUM,
            noCuts,localCutWeights, globalCutWeights
        );
      }
      Z2_THROW_OUTSIDE_ERROR(*env)


      for (partId_t i = 0; i < noCuts; ++i){
        //cout << "gw:" << globalCutWeights[i] << endl;
        if(globalCutWeights[i] > 0) {
          scalar_t ew = globalTotalWeight * targetPartWeightRatios[i];
          scalar_t expectedWeightOnLine = ew - globalPartWeights[i * 2];
          scalar_t myWeightOnLine = localCutWeights[i];
          scalar_t weightOnLineBefore = globalCutWeights[i];
          scalar_t incMe = expectedWeightOnLine - weightOnLineBefore;
          scalar_t mine = incMe + myWeightOnLine;
          if(mine < 0){
            nonRectelinearPartRatios[i] = 0;
          }
          else if(mine >= myWeightOnLine){
            nonRectelinearPartRatios[i] = 1;
          }
          else {
            nonRectelinearPartRatios[i] = mine / myWeightOnLine;
          }
        }
      }
      *rectilinearCutCount = 0;
    }
  }
}


/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 *
 * \param total_part_count is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
 * \param noCuts is the number of cut lines. P - 1.
 * \param maxScalar is the maximum value that scalar_t can represent.
 * \param _EPSILON is the smallest error value for scalar_t.
 * \param numThreads hold the number of threads available per processor.
 *
 * \param coordinateBegin is the index of the first coordinate in current part.
 * \param coordinateEnd is the index of the last coordinate in current part.
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param pqJagged_coordinates is 1 dimensional array holding coordinate values.
 * \param pqJagged_uniformWeights is a boolean value if the points have uniform weights.
 * \param pqJagged_weights is 1 dimensional array holding the weights of points.
 *
 * \param globalPartWeights is the array holding the weight of parts. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param localPartWeights is the local totalweight of the processor.
 * \param targetPartWeightRatios are the desired cumulative part ratios, sized P.
 *
 * \param cutCoordinates_tmp is the array holding the coordinates of each cut line. Sized P - 1.
 * \param isDone is the boolean array to determine if the correct position for a cut line is found.
 * \param myPartWeights is the array holding the part weights for the calling thread.
 * \param myLeftClosest is the array holding the distances to the closest points to the cut lines from left for the calling thread..
 * \param myRightClosest is the array holding the distances to the closest points to the cut lines from right for the calling thread.
 * \param useBinarySearch is boolean parameter whether to search for cut lines with binary search of linear search.
 * \param partIds is the array that holds the part ids of the coordinates
 * kddnote:  The output of this function myPartWeights differs depending on whether
 * kddnote:  binary or linear search is used.  Values in myPartWeights should be the
 * kddnote:  same only after accumulateThreadResults is done.
 */
template <typename scalar_t, typename lno_t>
void pqJagged_1DPart_getPartWeights(
    size_t total_part_count,
    partId_t noCuts,
    scalar_t maxScalar,
    scalar_t _EPSILON,
    int numThreads,
    lno_t coordinateBegin,
    lno_t coordinateEnd,
    lno_t *partitionedPointPermutations,
    scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    scalar_t *pqJagged_weights,

    scalar_t *cutCoordinates_tmp, //TODO change name
    bool *isDone,
    double *myPartWeights,
    scalar_t *myLeftClosest,
    scalar_t *myRightClosest,
    bool useBinarySearch,
    partId_t *partIds
){

  // initializations for part weights, left/right closest
    for (size_t i = 0; i < total_part_count; ++i){
      myPartWeights[i] = 0;
    }


  for(partId_t i = 0; i < noCuts; ++i){
    //if(isDone[i]) continue;
    myLeftClosest[i] = maxScalar;
    myRightClosest[i] = maxScalar;
  }
  //cout << "iter:" << endl;
  if(useBinarySearch){

	  //lno_t comparison_count = 0;
    scalar_t minus_EPSILON = -_EPSILON;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
      int i = partitionedPointPermutations[ii];
      //partId_t j = (uc + lc) / 2;
      partId_t j = partIds[i] / 2;

      if(j >= noCuts){
    	  j = noCuts - 1;
      }

      partId_t lc = 0;
      partId_t uc = noCuts - 1;

      scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
      bool isInserted = false;
      bool onLeft = false;
      bool onRight = false;
      partId_t lastPart = -1;
      
      scalar_t coord = pqJagged_coordinates[i];

      while(uc >= lc)
      {
    	//comparison_count++;
        lastPart = -1;
        onLeft = false;
        onRight = false;
        scalar_t cut = cutCoordinates_tmp[j];
        scalar_t distance = coord - cut;
        scalar_t absdistance = ABS(distance);

        if(absdistance < _EPSILON){
          myPartWeights[j * 2 + 1] += w;
          partIds[i] = j * 2 + 1;

          myLeftClosest[j] = 0;
          myRightClosest[j] = 0;
          partId_t kk = j + 1;
          while(kk < noCuts){  // Needed when cuts shared the same position
                               // kddnote Can this loop be disabled for RECTILINEAR BLOCKS?
                               // kddnote Mehmet says it is probably needed anyway.
            distance =ABS(cutCoordinates_tmp[kk] - cut);
            if(distance < _EPSILON){
              myPartWeights[2 * kk + 1] += w;

              myLeftClosest[kk] = 0;
              myRightClosest[kk] = 0;
              kk++;
            }
            else{
              if(myLeftClosest[kk] > distance){
                myLeftClosest[kk] = distance;
              }
              break;
            }
          }

          kk = j - 1;
          while(kk >= 0){
            distance =ABS(cutCoordinates_tmp[kk] - cut);
            if(distance < _EPSILON){
              myPartWeights[2 * kk + 1] += w;

              myLeftClosest[kk] = 0;
              myRightClosest[kk] = 0;
              kk--;
            }
            else{
              if(myRightClosest[kk] > distance){
                myRightClosest[kk] = distance;
              }
              break;
            }
          }
          isInserted = true;
          break;
        }
        else {
          if (distance < 0) {
            //TODO fix abs
            distance = absdistance;
            /*
            if (myLeftClosest[j] > distance){
              myLeftClosest[j] = distance;
            }
             */
            bool _break = false;
            if(j > 0){
              distance = coord - cutCoordinates_tmp[j - 1];
              if(distance > _EPSILON){
                /*
                if (myRightClosest[j - 1] > distance){
                  myRightClosest[j - 1] = distance;
                }
                 */
                _break = true;
              }
              /*
              else if(distance < minus_EPSILON){
                distance = -distance;
                if (myLeftClosest[j - 1] > distance){
                  myLeftClosest[j - 1] = distance;
                }
              }
              else {
                myLeftClosest[j - 1] = 0;
                myRightClosest[j - 1 ] = 0;
              }*/
            }
            uc = j - 1;
            onLeft = true;
            lastPart = j;
            if(_break) break;
          }
          else {
            /*
            if (myRightClosest[j] > distance){
              myRightClosest[j] = distance;
            }
             */
            bool _break = false;
            if(j < noCuts - 1){
              scalar_t distance_ = coord - cutCoordinates_tmp[j + 1];
              /*
              if(distance > _EPSILON){
                if (myRightClosest[j + 1] > distance){
                  myRightClosest[j + 1] = distance;
                }
              } else */if(distance_ < minus_EPSILON){
                /*
                distance = -distance;
                if (myLeftClosest[j + 1] > distance){
                  myLeftClosest[j + 1] = distance;
                }*/
                _break = true;
              }
              /*
              else {
                myLeftClosest[j + 1] = 0;
                myRightClosest[j + 1 ] = 0;
              }
               */
            }
            lc = j + 1;
            onRight = true;
            lastPart = j;
            if(_break) break;
          }
        }

        j = (uc + lc) / 2;
      }
      if(!isInserted){
        if(onRight){
          
          
          myPartWeights[2 * lastPart + 2] += w;
          partIds[i] = 2 * lastPart + 2;
          scalar_t distance = coord - cutCoordinates_tmp[lastPart];
          if(myRightClosest[lastPart] > distance){
            myRightClosest[lastPart] = distance;
          }
          if(lastPart+1 < noCuts){
            scalar_t distance_ = cutCoordinates_tmp[lastPart + 1] - coord;
            if(myLeftClosest[lastPart + 1] > distance_){
              myLeftClosest[lastPart + 1] = distance_;
            }
          }
          
        }
        else if(onLeft){
          myPartWeights[2 * lastPart] += w;
          partIds[i] = 2 * lastPart;
          scalar_t distance = cutCoordinates_tmp[lastPart ] - coord;
          if(myLeftClosest[lastPart] > distance){
            myLeftClosest[lastPart] = distance;
          }
          
          if(lastPart-1 >= 0){
            scalar_t distance_ = coord - cutCoordinates_tmp[lastPart - 1];
            if(myRightClosest[lastPart -1] > distance_){
              myRightClosest[lastPart -1] = distance_;
            }
          }
        }
      }
    }
/*
    //cout << "comp:" << comparison_count << " size:" << coordinateEnd- coordinateBegin << endl;
    // prefix sum computation.
    for (size_t i = 1; i < total_part_count; ++i){
      // if check for cuts sharing the same position; all cuts sharing a position
      // have the same weight == total weight for all cuts sharing the position.
      // don't want to accumulate that total weight more than once.
      if(i % 2 == 0 && i > 1 && i < total_part_count - 1 &&
         ABS(cutCoordinates_tmp[i / 2] - cutCoordinates_tmp[i /2 - 1]) < _EPSILON){
        myPartWeights[i] = myPartWeights[i-2];
        continue;
      }
      myPartWeights[i] += myPartWeights[i-1];
    }

*/
  }
  else {

    //lno_t comp = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
      int i = partitionedPointPermutations[ii];

      scalar_t w = pqJagged_uniformWeights ? 1 : pqJagged_weights[i];

      scalar_t coord = pqJagged_coordinates[i];

      partId_t j = partIds[i] / 2;

      if(j >= noCuts){
    	  j = noCuts - 1;
      }
      scalar_t cut = cutCoordinates_tmp[j];
      scalar_t distance = coord - cut;
      scalar_t absdistance = ABS(distance);

      //comp++;
      if(absdistance < _EPSILON){
    	  myPartWeights[j * 2 + 1] += w;
        //cout << "1to part:" << 2*j+1 << " coord:" << coord << endl;
    	  myLeftClosest[j] = 0;
    	  myRightClosest[j] = 0;
        partIds[i] = j * 2 + 1;

        //bas
        partId_t kk = j + 1;
        while(kk < noCuts){  // Needed when cuts shared the same position
          // kddnote Can this loop be disabled for RECTILINEAR BLOCKS?
          // kddnote Mehmet says it is probably needed anyway.
          distance =ABS(cutCoordinates_tmp[kk] - cut);
          if(distance < _EPSILON){
            //cout << "yo" << endl;
            myPartWeights[2 * kk + 1] += w;

            //cout << "2to part:" << 2*kk+1 << " coord:" << coord << endl;
            myLeftClosest[kk] = 0;
            myRightClosest[kk] = 0;
            kk++;
          }
          else{
            if(myLeftClosest[kk] > distance){
              myLeftClosest[kk] = distance;
            }
            break;
          }
        }

        kk = j - 1;
        while(kk >= 0){
          distance =ABS(cutCoordinates_tmp[kk] - cut);
          if(distance < _EPSILON){
            myPartWeights[2 * kk + 1] += w;
            //cout << "3to part:" << 2*kk+1 << " coord:" << coord << endl;

            myLeftClosest[kk] = 0;
            myRightClosest[kk] = 0;
            kk--;
          }
          else{
            if(myRightClosest[kk] > distance){
              myRightClosest[kk] = distance;
            }
            break;
          }
        }
      }
      else if (distance < 0) {
    	  while((absdistance > _EPSILON) && distance < 0){
          //comp++;
    		  if (myLeftClosest[j] > absdistance){
    			  myLeftClosest[j] = absdistance;
    		  }
    		  --j;
    		  if(j  < 0) break;
    		  distance = coord - cutCoordinates_tmp[j];

    		  absdistance = ABS(distance);
      	  }
        if(absdistance < _EPSILON)
        {
          myPartWeights[j * 2 + 1] += w;
          myLeftClosest[j] = 0;
          myRightClosest[j] = 0;
          cut = cutCoordinates_tmp[j];
          //cout << "4to part:" << 2*j+1 <<" j:" << j <<  " coord:" << coord << endl;
          //cout << "cut:" << cutCoordinates_tmp[j] << " dis:" << distance << endl;
          partIds[i] = j * 2 + 1;

          partId_t kk = j + 1;
          while(kk < noCuts){  // Needed when cuts shared the same position
            // kddnote Can this loop be disabled for RECTILINEAR BLOCKS?
            // kddnote Mehmet says it is probably needed anyway.
            distance =ABS(cutCoordinates_tmp[kk] - cut);
            //cout << "distance:" << distance << endl;
            if(distance < _EPSILON){
              myPartWeights[2 * kk + 1] += w;
              //cout << "5to part:" << 2*kk+1 << " kk:" << kk << " coord:" << coord << endl;
              //cout << "cut:" << cutCoordinates_tmp[kk] << " dis:" << distance << endl;
              myLeftClosest[kk] = 0;
              myRightClosest[kk] = 0;
              kk++;
            }
            else{
              if(myLeftClosest[kk] > distance){
                myLeftClosest[kk] = distance;
              }
              break;
            }
          }

          kk = j - 1;
          while(kk >= 0){
            distance =ABS(cutCoordinates_tmp[kk] - cut);
            if(distance < _EPSILON){
              myPartWeights[2 * kk + 1] += w;
              //cout << "6to part:" << 2*kk+1 << " coord:" << coord << endl;
              //cout << "cut:" << cutCoordinates_tmp[kk] << " dis:" << distance << endl;

              myLeftClosest[kk] = 0;
              myRightClosest[kk] = 0;
              kk--;
            }
            else{
              if(myRightClosest[kk] > distance){
                myRightClosest[kk] = distance;
              }
              break;
            }
          }
        }
        else {
          myPartWeights[j * 2 + 2] += w;
          if (j >= 0 && myRightClosest[j] > absdistance){
    			  myRightClosest[j] = absdistance;
    		  }
          partIds[i] = j * 2 + 2;
        }

      }
      //if it is on the left
      else {

    	  while((absdistance > _EPSILON) && distance > 0){
          //comp++;
        	  if ( myRightClosest[j] > absdistance){
        		  myRightClosest[j] = absdistance;
        	  }
    		  ++j;
    		  if(j  >= noCuts) break;
    		  distance = coord - cutCoordinates_tmp[j];
    		  absdistance = ABS(distance);
      	  }



        if(absdistance < _EPSILON)
        {
          myPartWeights[j * 2 + 1] += w;
          myLeftClosest[j] = 0;
          myRightClosest[j] = 0;
          partIds[i] = j * 2 + 1;
          cut = cutCoordinates_tmp[j];
          partId_t kk = j + 1;
          while(kk < noCuts){  // Needed when cuts shared the same position
            // kddnote Can this loop be disabled for RECTILINEAR BLOCKS?
            // kddnote Mehmet says it is probably needed anyway.
            distance =ABS(cutCoordinates_tmp[kk] - cut);
            if(distance < _EPSILON){
              myPartWeights[2 * kk + 1] += w;
              myLeftClosest[kk] = 0;
              myRightClosest[kk] = 0;
              kk++;
            }
            else{
              if(myLeftClosest[kk] > distance){
                myLeftClosest[kk] = distance;
              }
              break;
            }
          }

          kk = j - 1;
          while(kk >= 0){
            distance =ABS(cutCoordinates_tmp[kk] - cut);
            if(distance < _EPSILON){
              myPartWeights[2 * kk + 1] += w;

              myLeftClosest[kk] = 0;
              myRightClosest[kk] = 0;
              kk--;
            }
            else{
              if(myRightClosest[kk] > distance){
                myRightClosest[kk] = distance;
              }
              break;
            }
          }
        }
        else {
    	  myPartWeights[j * 2] += w;
        if (j < noCuts && myLeftClosest[j] > absdistance){
          myLeftClosest[j] = absdistance;
        }
        partIds[i] = j * 2 ;
        }
      }
    }

    /*
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
      int i = partitionedPointPermutations[ii];

      scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
      //get a coordinate and compare it with cut lines from left to right.
      for(partId_t j = 0; j < noCuts; ++j){

        if(isDone[j]) continue;
        scalar_t distance = pqJagged_coordinates[i] - cutCoordinates_tmp[j];
        scalar_t absdistance = ABS(distance);
        if(absdistance < _EPSILON){
          myPartWeights[j * 2] -= w;
          //myPartWeights[j * 2] += w;
          myLeftClosest[j] = 0;
          myRightClosest[j] = 0;
          //break;
        }
        else
          if (distance < 0) {
            distance = -distance;
            if (myLeftClosest[j] > distance){
              myLeftClosest[j] = distance;
            }
            break;
          }
        //if it is on the left
          else {
            //if on the right, continue with the next line.
            myPartWeights[j * 2] -= w;
            myPartWeights[j * 2 + 1] -= w;
            //myPartWeights[j * 2] += w;
            //myPartWeights[j * 2 + 1] += w;

            if ( myRightClosest[j] > distance){
              myRightClosest[j] = distance;
            }
          }
      }
    }
    */

  }

    //cout << "comp:" << comp <<endl;
    // prefix sum computation.
    for (size_t i = 1; i < total_part_count; ++i){
      // if check for cuts sharing the same position; all cuts sharing a position
      // have the same weight == total weight for all cuts sharing the position.
      // don't want to accumulate that total weight more than once.
      if(i % 2 == 0 && i > 1 && i < total_part_count - 1 &&
         ABS(cutCoordinates_tmp[i / 2] - cutCoordinates_tmp[i /2 - 1]) < _EPSILON){
        myPartWeights[i] = myPartWeights[i-2];
        continue;
      }
      myPartWeights[i] += myPartWeights[i-1];
      //cout << "p:" << "i:" << i<< " :" <<myPartWeights[i]  << endl;
    }
  /*
  for (size_t i = 0; i < total_part_count; ++i){
    cout << "p:" << i << ":" << myPartWeights[i] << endl;
  }
   */


}


/*! \brief Function that reduces the result of multiple threads for left and right closest points and part weights in a single mpi.
 *
 * \param noCuts is the number of cut lines. P - 1.
 * \param total_part_count is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
 * \param concurrentPartCount is the number of parts whose cut lines will be calculated concurrently.
 * \param numThreads hold the number of threads available per processor.
 * \param isDone is the boolean array to determine if the correct position for a cut line is found.
 * \param leftClosestDistance is the 2 dimensional array holding the distances to the closest points to the cut lines from left for each thread.
 * \param rightClosestDistance is the 2 dimensional array holding the distances to the closest points to the cut lines from right for each thread.
 * \param partWeights is the array holding the weight of parts for each thread. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param localMinMaxTotal is the array holding the local minimum and maximum coordinate and local total weight of each part.
 * \param totalPartWeights_leftClosest_rightCloset is the output array of accumulation, where total part weights (2P - 1),
 * then left closest distances (P -1), then right closest distance (P -1) are stored.
 */
template <typename scalar_t>
void accumulateThreadResults(
    partId_t noCuts,
    size_t total_part_count,
    partId_t concurrentPartCount,
    int numThreads,
    bool *isDone,
    scalar_t **leftClosestDistance,
    scalar_t **rightClosestDistance,
    double **partWeights,
    scalar_t *localMinMaxTotal,
    scalar_t *totalPartWeights_leftClosest_rightCloset//,
    //bool useBinarySearch
    ){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
  for(partId_t i = 0; i < noCuts * concurrentPartCount; ++i){
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
    size_t ishift = i % noCuts + (i / noCuts) * (total_part_count + 2 * noCuts);
    totalPartWeights_leftClosest_rightCloset[total_part_count + ishift] = minl;
    totalPartWeights_leftClosest_rightCloset[total_part_count + noCuts + ishift] = minr;
  }

 // if(1 || useBinarySearch){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(size_t j = 0; j < total_part_count * concurrentPartCount; ++j){
      size_t actualCutInd = (j % total_part_count);
      partId_t cutInd = actualCutInd / 2 + (j / total_part_count) * noCuts;

      if(actualCutInd !=  total_part_count - 1 && isDone[cutInd]) continue;
      double pwj = 0;
      for (int i = 0; i < numThreads; ++i){
        pwj += partWeights[i][j];
      }
      size_t jshift = j % total_part_count + (j / total_part_count) * (total_part_count + 2 * noCuts);
      totalPartWeights_leftClosest_rightCloset[jshift] = pwj;
    }
/*
  } else {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(size_t j = 0; j < total_part_count * concurrentPartCount; ++j){
      size_t actualCutInd = (j % total_part_count);
      partId_t cutInd = actualCutInd / 2 + (j / total_part_count) * noCuts;

      if(actualCutInd !=  total_part_count - 1 && isDone[cutInd]) continue;
      double pwj = 0;
      for (int i = 0; i < numThreads; ++i){
        pwj += partWeights[i][j];
      }
      scalar_t totalWeight = localMinMaxTotal[j / total_part_count + concurrentPartCount * 2];
      size_t jshift = j % total_part_count + (j / total_part_count) * (total_part_count + 2 * noCuts);
      totalPartWeights_leftClosest_rightCloset[jshift] = totalWeight + pwj;

    }
  }
*/
}




/*! \brief Function that is responsible from 1D partitioning of the given range of coordinates.
 * \param env library configuration and problem parameters
 * \param comm the communicator for the problem
 *
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param pqJagged_coordinates is 1 dimensional array holding coordinate values.
 * \param pqJagged_uniformWeights is a boolean value if the points have uniform weights.
 * \param pqJagged_weights is 1 dimensional array holding the weights of points.
 *
 * \param targetPartWeightRatios is the cumulative desired weight ratios for each part. Sized P.
 * \param globalMinMaxTotal is the array holding the global min,max and total weight.
 * minimum of part k is on k, maximum is on k + concurrentPartcount, totalweight is on k + 2*concurrentPartCount
 * \param localMinMaxTotal is the array holding the local min,max and total weight. Structure is same with globalMinMaxTotal array.
 *
 * \param partNo is the total number of parts.
 * \param numThreads is the number of available threads by each processor.
 * \param maxScalar is the maximum value that scalar_t can represent.
 * \param minScalar is the minimum value that scalar_t can represent.
 * \param imbalanceTolerance is the maximum allowed imbalance ratio.
 * \param currentPartBeginIndex is the beginning index of concurrentPartCount parts on inTotalCounts array.
 * \param concurrentPartCount is the number of parts whose cut lines will be calculated concurrently.
 * \param inTotalCounts is the array holding the beginning indices of the parts from previous dimension partitioning.
 *
 * \param cutCoordinates is the array holding the coordinates of the cut.
 * \param cutCoordinatesWork is a work array sized P - 1.
 * \param leftClosestDistance is the two dimensional array holding the distances to the closest points to the cut lines from left. One dimension for each thread.
 * \param rightClosestDistance is the two dimensional array holding the distances to the closest points to the cut lines from right. One dimension for each thread.
 * \param cutUpperBounds is the array holding the upper bound coordinate for each cut line. Sized P - 1.
 * \param cutLowerBounds is the array holding the lower bound coordinate for each cut line. Sized P - 1.
 * \param cutUpperWeight is the array holding the weight of the parts at the left of upper bound coordinates.
 * \param cutLowerWeight is the array holding the weight of the parts at the left of lower bound coordinates.
 * \param isDone is the boolean array to determine if the correct position for a cut line is found.
 * \param partWeights is the two dimensional array holding the part weights. One dimension for each thread.
 * \param local_totalPartWeights_leftClosest_rightCloset is one dimensional array holding the local part weights, left and right closest points. Convenient for mpi-reduction.
 * \param global_totalPartWeights_leftClosest_rightCloset is one dimensional array holding the local part weights, left and right closest points. Convenient for mpi-reduction.
 * \param allowNonRectelinearPart is whether to allow distribution of points on same coordinate to different parts or not.
 * \param nonRectelinearPartRatios is one dimensional array holding the ratio of the points on cut lines representing how many of them will be put to left of the line.
 * \param localCutWeights is one dimensional array holding the local part weights only on cut lines.
 * \param globalCutWeights is one dimensional array holding the global part weights only on cut lines.
 * \param allDone is the number of cut lines whose positions should be calculated.
 * \param myNonDoneCounts is one dimensional array holding the number of cut lines in each part whose positions should be calculated.
 * \param useBinarySearch is boolean parameter whether to search for cut lines with binary search of linear search.
 * \param partIds is the array that holds the part ids of the coordinates
 *
 */
template <typename scalar_t, typename lno_t>
void pqJagged_1D_Partition(
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,

    lno_t *partitionedPointPermutations,
    scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    scalar_t *pqJagged_weights,

    scalar_t *targetPartWeightRatios,   // the weight ratios at left side of the cuts. last is 1.
    scalar_t *globalMinMaxTotal,
    scalar_t *localMinMaxTotal,

    partId_t partNo,
    int numThreads,
    scalar_t maxScalar,
    scalar_t minScalar,
    scalar_t imbalanceTolerance,
    partId_t currentPartBeginIndex,
    partId_t concurrentPartCount,
    lno_t *inTotalCounts,

    scalar_t *cutCoordinates,
    scalar_t *cutCoordinatesWork, 	// work array to manipulate coordinate of cutlines in different iterations.
    scalar_t **leftClosestDistance,
    scalar_t **rightClosestDistance,
    scalar_t *cutUpperBounds,  //to determine the next cut line with binary search
    scalar_t *cutLowerBounds,  //to determine the next cut line with binary search
    scalar_t *cutUpperWeight,   //to determine the next cut line with binary search
    scalar_t *cutLowerWeight,  //to determine the next cut line with binary search
    bool *isDone,
    double **partWeights,
    scalar_t *local_totalPartWeights_leftClosest_rightCloset,
    scalar_t *global_totalPartWeights_leftClosest_rightCloset,
    bool allowNonRectelinearPart,
    float *nonRectelinearPartRatios,
    scalar_t *localCutWeights,
    scalar_t *globalCutWeights,

    partId_t allDone,
    partId_t *myNonDoneCounts,
    bool useBinarySearch,
//    string dimension,
    partId_t * partIds
){

  partId_t recteLinearCutCount = 0;
  scalar_t *cutCoordinates_tmp = cutCoordinates;
  partId_t noCuts = partNo - 1;
  size_t total_part_count = partNo + size_t (noCuts) ;


#ifdef mpi_communication
  MPI_Op myop;
  MPI_Op_create(sumMinMin, 0, &myop);   /* step 3 */
#endif


  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();


//  env->timerStart(MACRO_TIMERS, "PQJagged creation_" + dimension);
  Teuchos::PQJaggedCombinedReductionOp<int, scalar_t> reductionOp(
      total_part_count , noCuts  , noCuts , concurrentPartCount);

 // env->timerStop(MACRO_TIMERS, "PQJagged creation_" + dimension);

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(allDone,  recteLinearCutCount)
#endif
  {
    //env->timerStart(MACRO_TIMERS, "PQJagged thread_init_" + dimension);
    int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
    me = omp_get_thread_num();
#endif
    double *myPartWeights = partWeights[me];
    scalar_t *myLeftClosest = leftClosestDistance[me];
    scalar_t *myRightClosest = rightClosestDistance[me];

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(partId_t i = 0; i < noCuts * concurrentPartCount; ++i){
      isDone[i] = false;
      partId_t ind = i / noCuts;
      cutLowerBounds[i] = globalMinMaxTotal[ind];
      cutUpperBounds[i] = globalMinMaxTotal[ind + concurrentPartCount];
      cutUpperWeight[i] = globalMinMaxTotal[ind + 2 * concurrentPartCount];
      cutLowerWeight[i] = 0;
      if(allowNonRectelinearPart){
        nonRectelinearPartRatios[i] = 0;
      }
    }

    //env->timerStop(MACRO_TIMERS, "PQJagged thread_init_" + dimension);


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#endif
    int iteration = 0;
    while (allDone != 0){
      iteration += 1;
/*
      if(comm->getRank() == 0)
        {
#pragma omp single
          { cout << endl << endl << "allDone:" << allDone << endl;

            for (size_t i = 0; i < noCuts * concurrentPartCount; ++i){

              if(isDone[i] == false)
            	  cout << "i:" << i <<  " c:" << cutCoordinates_tmp[i] << " u:" << cutUpperBounds[i] << " l:" << cutLowerBounds[i] <<
            	  " uw:" << cutUpperWeight[i] << " lw:" << cutLowerWeight[i] << " not done" << endl;
              else cout << "i:" << i <<  " c:" << cutCoordinates_tmp[i] <<  " done" << endl; }

          }
        }
*/

/*
#pragma omp single
      {
      env->timerStart(MACRO_TIMERS, "PQJagged 1D_" + dimension);
      }
      */
      for (partId_t kk = 0; kk < concurrentPartCount; ++kk){
        if (/*globalMinMax[kk] > globalMinMax[kk + k] ||*/ myNonDoneCounts[kk] == 0) continue;
        partId_t cutShift = noCuts * kk;
        size_t totalPartShift = total_part_count * kk;

        bool *currentDone = isDone + cutShift;
        double *myCurrentPartWeights = myPartWeights + totalPartShift;
        scalar_t *myCurrentLeftClosest = myLeftClosest + cutShift;
        scalar_t *myCurrentRightClosest = myRightClosest + cutShift;

        partId_t current = currentPartBeginIndex + kk;
        lno_t coordinateBegin = current ==0 ? 0: inTotalCounts[current -1];
        lno_t coordinateEnd = inTotalCounts[current];
        scalar_t *cutCoordinates_ = cutCoordinates_tmp + kk * noCuts;

        // compute part weights using existing cuts
        pqJagged_1DPart_getPartWeights<scalar_t, lno_t>(
            total_part_count,
            noCuts,
            maxScalar,
            _EPSILON,
            numThreads,
            coordinateBegin,
            coordinateEnd,
            partitionedPointPermutations,
            pqJagged_coordinates,
            pqJagged_uniformWeights,
            pqJagged_weights,
            cutCoordinates_,
            currentDone,
            myCurrentPartWeights,
            myCurrentLeftClosest,
            myCurrentRightClosest,
            useBinarySearch,
            partIds);
      }
      /*
#pragma omp single
      {
      env->timerStop(MACRO_TIMERS, "PQJagged 1D_" + dimension);
      env->timerStart(MACRO_TIMERS, "PQJagged 1D_accumulation_" + dimension);
      }
      */
      accumulateThreadResults<scalar_t>(
          noCuts, total_part_count, concurrentPartCount,
          numThreads, isDone,
          leftClosestDistance, rightClosestDistance, partWeights,
          localMinMaxTotal,
          local_totalPartWeights_leftClosest_rightCloset//,
//          useBinarySearch
      );
      /*
#pragma omp single
      {
      env->timerStop(MACRO_TIMERS, "PQJagged 1D_accumulation_" + dimension);

      env->timerStart(MACRO_TIMERS, "PQJagged Communication_" + dimension);
      }
      */
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
      {
        if(comm->getSize() > 1){
        try{



#ifdef mpi_communication

            MPI_Allreduce(local_totalPartWeights_leftClosest_rightCloset, global_totalPartWeights_leftClosest_rightCloset,

                          (total_part_count + 2 * noCuts) * concurrentPartCount, MPI_FLOAT, myop,MPI_COMM_WORLD);
#endif
#ifndef mpi_communication

          reduceAll<int, scalar_t>(*comm, reductionOp,
              (total_part_count + 2 * noCuts) * concurrentPartCount,
              local_totalPartWeights_leftClosest_rightCloset,
              global_totalPartWeights_leftClosest_rightCloset
          );
#endif

          //cout << "reducing" << endl;
          /*
          MPI_Allreduce (
              local_totalPartWeights_leftClosest_rightCloset,
              global_totalPartWeights_leftClosest_rightCloset,
              total_part_count ,
              MPI_DOUBLE,
              MPI_SUM,
              MPI_COMM_WORLD);

          MPI_Allreduce (
              local_totalPartWeights_leftClosest_rightCloset + total_part_count ,
              global_totalPartWeights_leftClosest_rightCloset + total_part_count ,
              noCuts,
              MPI_DOUBLE,
              MPI_MIN,
              MPI_COMM_WORLD);

          MPI_Allreduce (
              local_totalPartWeights_leftClosest_rightCloset + total_part_count + noCuts,
              global_totalPartWeights_leftClosest_rightCloset + total_part_count + noCuts,
              noCuts,
              MPI_DOUBLE,
              MPI_MIN,
              MPI_COMM_WORLD);
        */
        }

        Z2_THROW_OUTSIDE_ERROR(*env)
        }
        else {
          size_t s = (total_part_count + 2 * noCuts) * concurrentPartCount;
          for(size_t i = 0; i < s; ++i){
            global_totalPartWeights_leftClosest_rightCloset[i] = local_totalPartWeights_leftClosest_rightCloset[i];
          }
        }
      }
      /*
#pragma omp single
      {
      env->timerStop(MACRO_TIMERS, "PQJagged Communication_" + dimension);

      env->timerStart(MACRO_TIMERS, "PQJagged new_cut_calculation_" + dimension);
      }
      */
      for (partId_t kk = 0; kk < concurrentPartCount; ++kk){


        if (/*globalMinMax[kk] > globalMinMax[kk + k] || */myNonDoneCounts[kk] == 0) continue;
        partId_t cutShift = noCuts * kk;
        size_t tlrShift = (total_part_count + 2 * noCuts) * kk;

        scalar_t *localPartWeights = local_totalPartWeights_leftClosest_rightCloset  + tlrShift ;
        scalar_t *gtlr = global_totalPartWeights_leftClosest_rightCloset + tlrShift;
        scalar_t *glc = gtlr + total_part_count; //left closest points
        scalar_t *grc = gtlr + total_part_count + noCuts; //right closest points
        scalar_t *globalPartWeights = gtlr;
        bool *currentDone = isDone + cutShift;
        scalar_t *currentTargetPartWeightRatios = targetPartWeightRatios + (noCuts + 1) * kk;
        float *currentnonRectelinearPartRatios = nonRectelinearPartRatios + cutShift;

        scalar_t minCoordinate = globalMinMaxTotal[kk];
        scalar_t maxCoordinate = globalMinMaxTotal[kk + concurrentPartCount];
        scalar_t globalTotalWeight = globalMinMaxTotal[kk + concurrentPartCount * 2];
        scalar_t *currentcutLowerWeight = cutLowerWeight + cutShift;
        scalar_t *currentcutUpperWeight = cutUpperWeight + cutShift;
        scalar_t *currentcutUpperBounds = cutUpperBounds + cutShift;
        scalar_t *currentcutLowerBounds = cutLowerBounds + cutShift;

        partId_t prevDoneCount = myNonDoneCounts[kk];
        // Now compute the new cut coordinates.
        getNewCoordinates<scalar_t>(
            env, comm,
            total_part_count, noCuts, maxCoordinate, minCoordinate, globalTotalWeight, imbalanceTolerance, maxScalar,
            globalPartWeights, localPartWeights, currentTargetPartWeightRatios, currentDone,

            cutCoordinates_tmp + cutShift, currentcutUpperBounds, currentcutLowerBounds,
            glc, grc, currentcutLowerWeight, currentcutUpperWeight,
            cutCoordinatesWork +cutShift, //new cut coordinates

            allowNonRectelinearPart,
            currentnonRectelinearPartRatios,
            &recteLinearCutCount,
            localCutWeights,
            globalCutWeights,
            myNonDoneCounts[kk]
        );

        partId_t reduction = prevDoneCount - myNonDoneCounts[kk];
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
        {
          allDone -= reduction;
        }
      }
      /*
#pragma omp single
      {
      env->timerStop(MACRO_TIMERS, "PQJagged new_cut_calculation_" + dimension);
      }
      */
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
      {
      scalar_t *t = cutCoordinates_tmp;
      cutCoordinates_tmp = cutCoordinatesWork;
      cutCoordinatesWork = t;
      }
    }
/*
    if(comm->getRank() == 0)
      if(me == 0) cout << "it:" << iteration << endl;
*/
    // Needed only if keep_cuts; otherwise can simply swap array pointers
    // cutCoordinates and cutCoordinatesWork.
    // (at first iteration, cutCoordinates == cutCoorindates_tmp).
    // computed cuts must be in cutCoordinates.
    if (cutCoordinates != cutCoordinates_tmp){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(partId_t i = 0; i < noCuts *concurrentPartCount; ++i){
        cutCoordinates[i] = cutCoordinates_tmp[i];
      }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
      {
        cutCoordinatesWork = cutCoordinates_tmp;
      }
    }

  }
}





/*! \brief Function that determines the permutation indices of the coordinates.
 * \param partNo is the number of parts.
 * \param noThreads is the number of threads avaiable for each processor.
 *
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param pqJagged_coordinates is 1 dimensional array holding the coordinate values.
 * \param pqJagged_uniformWeights is a boolean value if the points have uniform weights.
 * \param pqJagged_weights is 1 dimensional array holding the weights of points.
 *
 * \param cutCoordinates is 1 dimensional array holding the cut coordinates.
 * \param coordinateBegin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinateEnd is the end index of the given partition on partitionedPointPermutations.
 * \param numLocalCoord is the number of local coordinates.
 *
 * \param allowNonRectelinearPart is the boolean value whether partitioning should allow distributing the points on same coordinate to different parts.
 * \param actual_ratios holds how much percentage of the coordinates on the cutline should be put on left side.
 * \param localPartWeights is the local totalweight of the processor.
 * \param partWeights is the two dimensional array holding the weight of parts for each thread. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param nonRectelinearRatios is the two dimensional work array holding ratios of weights to be put left and right of the cut line.
 * \param partPointCounts is the two dimensional array holding the number of points in each part for each thread.
 *
 * \param newpartitionedPointPermutations is the indices of coordinates calculated for the partition on next dimension.
 * \param totalCounts are the number points in each output part.
 * \param partIds is the array that holds the part ids of the coordinates
 */
template <typename lno_t, typename scalar_t>
void getChunksFromCoordinates(
    partId_t partNo,
    int noThreads,

    lno_t *partitionedPointPermutations,
    scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    scalar_t *coordWeights,
    scalar_t *cutCoordinates,
    lno_t coordinateBegin,
    lno_t coordinateEnd,

    bool allowNonRectelinearPart,
    float *actual_ratios,
    scalar_t *localPartWeights,
    double **partWeights,
    float **nonRectelinearRatios,

    //lno_t *coordinate_linked_list,
    //lno_t **coordinate_starts,
    //lno_t **coordinate_ends,
    lno_t ** partPointCounts,

    lno_t *newpartitionedPointPermutations,
    lno_t *totalCounts,
    partId_t *partIds //,
    //bool useBinarySearch
    ){

  //lno_t numCoordsInPart =  coordinateEnd - coordinateBegin;
  partId_t noCuts = partNo - 1;
  //size_t total_part_count = noCuts + partNo;
  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
  {
    int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
    me = omp_get_thread_num();
#endif

    //lno_t *myStarts = coordinate_starts[me];
    //lno_t *myEnds = coordinate_ends[me];
    lno_t *myPartPointCounts = partPointCounts[me];
    float *myRatios = NULL;
    if (allowNonRectelinearPart){


      myRatios = nonRectelinearRatios[me];
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for (partId_t i = 0; i < noCuts; ++i){
        float r = actual_ratios[i];

        //cout << "real i:" << i << " :" << r << " " << endl;
        scalar_t leftWeight = r * (localPartWeights[i * 2 + 1] - localPartWeights[i * 2]);
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
      }


      if(noCuts > 0){
        for (partId_t i = noCuts - 1; i > 0 ; --i){          if(ABS(cutCoordinates[i] - cutCoordinates[i -1]) < _EPSILON){
            myRatios[i] -= myRatios[i - 1] ;
          }
          //cout << "i:" << i << " :" << myRatios[i] << " ";
          myRatios[i] = int ((myRatios[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL) / scalar_t(SIGNIFICANCE_MUL);
        }
      }
      /*

      for (partId_t i = 0; i < noCuts; ++i){
      cout << "r i:" << i << " :" <<  myRatios[i] << " " << endl;
      }
       */

    }

    for(partId_t ii = 0; ii < partNo; ++ii){
      myPartPointCounts[ii] = 0;
    }


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#pragma omp for
#endif
      for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

        lno_t i = partitionedPointPermutations[ii];
        partId_t pp = partIds[i];
        partId_t p = pp / 2;
        if(pp % 2 == 1){
          //cout << "on: " << pp << ":" << p << " cut:" << p << endl;

          if(allowNonRectelinearPart && myRatios[p] > _EPSILON * EPS_SCALE){
            //cout << "p:" << p << endl;
            scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];
            myRatios[p] -= w;
            if(myRatios[p] < 0 && p < noCuts - 1 && ABS(cutCoordinates[p+1] - cutCoordinates[p]) < _EPSILON){
              myRatios[p + 1] += myRatios[p];
            }
            ++myPartPointCounts[p];
            partIds[i] = p;
          }
          else{
            //scalar_t currentCut = cutCoordinates[p];
            //TODO:currently cannot divide 1 line more than 2 parts.
            //bug cannot be divided, therefore this part should change.
            //cout << "p:" << p+1 << endl;
            ++myPartPointCounts[p + 1];
            partIds[i] = p + 1;
          }
      }
      else {
        ++myPartPointCounts[p];

        partIds[i] = p;
      }
    }



#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(partId_t j = 0; j < partNo; ++j){
      lno_t pwj = 0;//, prev = 0;
      //if(j == 0) prev = 0;
      //else prev = totalCounts[j - 1];
      //lno_t prev2 = prev;
      for (int i = 0; i < noThreads; ++i){
        lno_t threadPartPointCount = partPointCounts[i][j];
        partPointCounts[i][j] = pwj;
        pwj += threadPartPointCount;

        }
      totalCounts[j] = pwj;// + prev2; //+ coordinateBegin;
    }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
    {
      for(partId_t j = 1; j < partNo; ++j){
        totalCounts[j] += totalCounts[j - 1];
      }
    }

    for(partId_t j = 1; j < partNo; ++j){
      myPartPointCounts[j] += totalCounts[j - 1] ;
    }


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
      lno_t i = partitionedPointPermutations[ii];
      partId_t p =  partIds[i];
      //cout << "p:" << p << " my:" << myPartPointCounts[p] << " begin:" << coordinateBegin << endl;
      newpartitionedPointPermutations[coordinateBegin + myPartPointCounts[p]++] = i;
    }

    /*
    pqJagged_PartVertices <lno_t, partId_t> pqPV;
    pqPV.set(coordinate_linked_list, myStarts, myEnds);
    memset(myPartPointCounts, 0, sizeof(lno_t) * partNo);


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (lno_t i = coordinateBegin; i < coordinateEnd; ++i){
      int ii = partitionedPointPermutations[i];
      coordinate_linked_list[ii] = -1;
    }

    for (partId_t i = 0; i < partNo; ++i){
      myEnds[i] = -1;
      myStarts[i] = -1;
    }


    //determine part of each point
    if(1 ||useBinarySearch){

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

        int i = partitionedPointPermutations[ii];
        partId_t lc = 0;
        partId_t uc = noCuts - 1;

        //bool isInserted = false;
        bool onLeft = false;
        bool onRight = false;
        partId_t lastPart = -1;
        while(uc >= lc)
        {
          lastPart = -1;
          onLeft = false;
          onRight = false;
          partId_t j = (uc + lc) / 2;
          scalar_t distance = pqJagged_coordinates[i] - cutCoordinates[j];
          scalar_t absdistance = ABS(distance);

          if(allowNonRectelinearPart  && absdistance < _EPSILON ){
            scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];
            partId_t jj = j -1;
            while(jj >= 0 && ABS(cutCoordinates[jj] - cutCoordinates[j]) < _EPSILON){
              --jj;
            }
            ++jj;
            for (;jj < noCuts && ABS(cutCoordinates[jj] - cutCoordinates[j]) < _EPSILON; ++jj){
              if(myRatios[jj] > _EPSILON * EPS_SCALE ){
                myRatios[jj] -= w;
                if(myRatios[jj] < 0 && jj < noCuts - 1 && ABS(cutCoordinates[jj+1] - cutCoordinates[jj]) < _EPSILON){
                  myRatios[jj + 1] += myRatios[jj];
                }
                onLeft = true;
                lastPart = jj;
                break;
              }
            }
            if(!onLeft){
              onRight= true;
              lastPart = jj -1;
            }

            break;
          }
          else {
            if (distance < 0) {
              uc = j - 1;
              onLeft = true;
              lastPart = j;
            }
            else {
              lc = j + 1;
              onRight = true;
              lastPart = j;
            }
          }
        }
        if(onRight){
          pqPV.inserToPart(lastPart + 1, i);
          myPartPointCounts[lastPart + 1] +=  1;
        }
        else if(onLeft){
          pqPV.inserToPart(lastPart , i);
          myPartPointCounts[lastPart] +=  1;
        }
      }
      for(partId_t i = 1; i < partNo; ++i){
        myPartPointCounts[i] += myPartPointCounts[i-1];
      }
    } else {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for (lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

        lno_t i = partitionedPointPermutations[ii];

        bool inserted = false;
        for(partId_t j = 0; j < noCuts; ++j){
          scalar_t distance = pqJagged_coordinates[i] - cutCoordinates[j];
          scalar_t absdistance = ABS(distance);

          if (allowNonRectelinearPart && myRatios[j] > _EPSILON * EPS_SCALE && absdistance < _EPSILON ){
            scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];
            scalar_t decrease = w;
            myRatios[j] -= decrease;

            if(myRatios[j] < 0 && j < noCuts - 1 && ABS(cutCoordinates[j+1] - cutCoordinates[j]) < _EPSILON){
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


    if(1 || useBinarySearch){
      //accumulate the count.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(partId_t j = 0; j < partNo; ++j){
        lno_t pwj = 0;
        for (int i = 0; i < noThreads; ++i){
          pwj += partPointCounts[i][j];
        }
        totalCounts[j] = pwj;
      }
    } else {
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
    }
*/

/*
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for(partId_t i = 0; i < partNo; ++i){
      lno_t nextPoint = coordinate_starts[0][i];
      lno_t pcnt = 0;

      lno_t prevCount = coordinateBegin;
      if (i > 0) prevCount = totalCounts[i -1] + coordinateBegin;

      while(nextPoint != -1){
        newpartitionedPointPermutations[prevCount + pcnt++] = nextPoint;
        nextPoint = coordinate_linked_list[nextPoint];
      }
    }
    */
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
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL

  Z2_THROW_EXPERIMENTAL("Zoltan2 PQJagged is strictly experimental software "
                        "while it is being developed and tested.")

#else
/*
 *   typedef typename Adapter::scalar_t scalar_t;
 *     typedef typename Adapter::gno_t gno_t;
 *       typedef typename Adapter::lno_t lno_t;
  if(comm->getRank() == 0){
    cout << "size of gno:" << sizeof(gno_t) << endl;
    cout << "size of lno:" << sizeof(lno_t) << endl;
    cout << "size of scalar_t:" << sizeof(scalar_t) << endl;
  }
 */
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
  bool ignoreWeights=false;

  bool allowNonRectelinearPart = false;
  int concurrentPartCount = 1;
  bool force_binary = false, force_linear = false;
  pqJagged_getParameters<scalar_t>(pl, imbalanceTolerance, mcnorm, params, numTestCuts, ignoreWeights,allowNonRectelinearPart,  concurrentPartCount,
      force_binary, force_linear);

  int coordDim, weightDim; size_t nlc; global_size_t gnc; int criteriaDim;
  pqJagged_getCoordinateValues<Adapter>( coords, coordDim, weightDim, nlc, gnc, criteriaDim, ignoreWeights);
  size_t numLocalCoords = nlc;

  //allocate only two dimensional pointer.
  //raw pointer addresess will be obtained from multivector.
  scalar_t **pqJagged_coordinates = allocMemory<scalar_t *>(coordDim);
  scalar_t **pqJagged_weights = allocMemory<scalar_t *>(criteriaDim);
  bool *pqJagged_uniformParts = allocMemory< bool >(criteriaDim); //if the partitioning results wanted to be uniform.
  scalar_t **pqJagged_partSizes =  allocMemory<scalar_t *>(criteriaDim); //if in a criteria dimension, uniform part is false this shows ratios of the target part weights.
  bool *pqJagged_uniformWeights = allocMemory< bool >(criteriaDim); //if the weights of coordinates are uniform in a criteria dimension.

  ArrayView<const gno_t> pqJagged_gnos;
  size_t numGlobalParts;
  int pqJagged_multiVectorDim;

  pqJagged_getInputValues<Adapter, scalar_t, gno_t>(
      env, coords, solution,params,coordDim,weightDim,numLocalCoords,
      numGlobalParts, pqJagged_multiVectorDim,
      pqJagged_coordinates,criteriaDim, pqJagged_weights,pqJagged_gnos, ignoreWeights,
      pqJagged_uniformWeights, pqJagged_uniformParts, pqJagged_partSizes
  );

  int numThreads = 1;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(numThreads)
  {
    numThreads = omp_get_num_threads();
  }

#endif


  partId_t totalDimensionCut = 0;
  partId_t totalPartCount = 1;
  partId_t maxPartNo = 0;

  const partId_t *partNo = pl.getPtr<Array <partId_t> >("pqParts")->getRawPtr();
  int partArraySize = pl.getPtr<Array <partId_t> >("pqParts")->size() - 1;

  for (int i = 0; i < partArraySize; ++i){
    totalPartCount *= partNo[i];
    if(partNo[i] > maxPartNo) maxPartNo = partNo[i];
  }
  totalDimensionCut = totalPartCount - 1;
  partId_t maxCutNo = maxPartNo - 1;
  partId_t maxTotalCumulativePartCount = totalPartCount / partNo[partArraySize - 1];
  size_t maxTotalPartCount = maxPartNo + size_t(maxCutNo);
  //maxPartNo is P, maxCutNo = P-1, matTotalPartcount = 2P-1

  if(concurrentPartCount > maxTotalCumulativePartCount){
    if(comm->getRank() == 0){
      cerr << "Warning: Concurrent part calculation count ("<< concurrentPartCount << ") has been set bigger than maximum amount that can be used." << " Setting to:" << maxTotalCumulativePartCount << "." << endl;
    }
    concurrentPartCount = maxTotalCumulativePartCount;
  }

  // coordinates of the cut lines. First one is the min, last one is max coordinate.
  // kddnote if (keep_cuts)
  // coordinates of the cut lines.
  //only store this much if cuts are needed to be stored.
  scalar_t *allCutCoordinates = allocMemory< scalar_t>(totalDimensionCut);
  // kddnote else
  //scalar_t *allCutCoordinates = allocMemory< scalar_t>(maxCutNo * concurrentPartCount);

  //as input indices.
  lno_t *partitionedPointCoordinates =  allocMemory< lno_t>(numLocalCoords);
  //as output indices
  lno_t *newpartitionedPointCoordinates = allocMemory< lno_t>(numLocalCoords);
  scalar_t *max_min_array =  allocMemory< scalar_t>(numThreads * 2);

  //initial configuration
  //set each pointer-i to i.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < numLocalCoords; ++i){
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

  //inTotalCounts array holds the end points in partitionedPointCoordinates array
  //for each partition. Initially sized 1, and single element is set to numLocalCoords.
  lno_t *inTotalCounts = allocMemory<lno_t>(1);
  inTotalCounts[0] = static_cast<lno_t>(numLocalCoords);//the end of the initial partition is the end of coordinates.

  //the ends points of the output.
  lno_t *outTotalCounts = NULL;

  //the array holding the points in each part as linked list
  //lno_t *coordinate_linked_list = allocMemory<lno_t>(numLocalCoords);
  //the start and end coordinate of  each part.
  //lno_t **coordinate_starts =  allocMemory<lno_t *>(numThreads);
  //lno_t **coordinate_ends = allocMemory<lno_t *>(numThreads);

  //assign the max size to starts, as it will be reused.
  /*
  for(int i = 0; i < numThreads; ++i){
    coordinate_starts[i] = allocMemory<lno_t>(maxPartNo);
    coordinate_ends[i] = allocMemory<lno_t>(maxPartNo);
  }
   */


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




  float *nonRectelinearPart = NULL; //how much weight percentage should a MPI put left side of the each cutline
  float **nonRectRatios = NULL; //how much weight percentage should each thread in MPI put left side of the each cutline

  if(allowNonRectelinearPart){
    //cout << "allowing" << endl;
    nonRectelinearPart = allocMemory<float>(maxCutNo * concurrentPartCount);
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

  // work array to manipulate coordinate of cutlines in different iterations.
  //necessary because previous cut line information is used for determining the next cutline information.
  //therefore, cannot update the cut work array until all cutlines are determined.
  scalar_t *cutCoordinatesWork = allocMemory<scalar_t>(maxCutNo * concurrentPartCount);
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
  firstTouch<scalar_t>(cutCoordinatesWork, maxCutNo);
#endif
#endif

  //cumulative part weight ratio array.
  scalar_t *targetPartWeightRatios = allocMemory<scalar_t>(maxPartNo * concurrentPartCount); // the weight ratios at left side of the cuts. First is 0, last is 1.
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
  firstTouch<scalar_t>(cutPartRatios, maxCutNo);
#endif
#endif


  scalar_t *cutUpperBounds = allocMemory<scalar_t>(maxCutNo * concurrentPartCount);  //upper bound coordinate of a cut line
  scalar_t *cutLowerBounds = allocMemory<scalar_t>(maxCutNo* concurrentPartCount);  //lower bound coordinate of a cut line
  scalar_t *cutLowerWeight = allocMemory<scalar_t>(maxCutNo* concurrentPartCount);  //lower bound weight of a cut line
  scalar_t *cutUpperWeight = allocMemory<scalar_t>(maxCutNo* concurrentPartCount);  //upper bound weight of a cut line

  scalar_t *localMinMaxTotal = allocMemory<scalar_t>(3 * concurrentPartCount); //combined array to exchange the min and max coordinate, and total weight of part.
  scalar_t *globalMinMaxTotal = allocMemory<scalar_t>(3 * concurrentPartCount);//global combined array with the results for min, max and total weight.

  //isDone is used to determine if a cutline is determined already.
  //If a cut line is already determined, the next iterations will skip this cut line.
  bool *isDone = allocMemory<bool>(maxCutNo * concurrentPartCount);
  //myNonDoneCount count holds the number of cutlines that have not been finalized for each part
  //when concurrentPartCount>1, using this information, if myNonDoneCount[x]==0, then no work is done for this part.
  partId_t *myNonDoneCount =  allocMemory<partId_t>(concurrentPartCount);
  //local part weights of each thread.
  double **partWeights = allocMemory<double *>(numThreads);
  //the work manupulation array for partweights.
  double **pws = allocMemory<double *>(numThreads);

  //leftClosesDistance to hold the min distance of a coordinate to a cutline from left (for each thread).
  scalar_t **leftClosestDistance = allocMemory<scalar_t *>(numThreads);
  //leftClosesDistance to hold the min distance of a coordinate to a cutline from right (for each thread)
  scalar_t **rightClosestDistance = allocMemory<scalar_t *>(numThreads);

  //to store how many points in each part a thread has.
  lno_t **partPointCounts = allocMemory<lno_t *>(numThreads);

  for(int i = 0; i < numThreads; ++i){
    //partWeights[i] = allocMemory<scalar_t>(maxTotalPartCount);
    partWeights[i] = allocMemory < double >(maxTotalPartCount * concurrentPartCount);
    rightClosestDistance[i] = allocMemory<scalar_t>(maxCutNo * concurrentPartCount);
    leftClosestDistance[i] = allocMemory<scalar_t>(maxCutNo * concurrentPartCount);
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
  // kddnote  Needed only when non-rectilinear parts.
  scalar_t *cutWeights = allocMemory<scalar_t>(maxCutNo);
  scalar_t *globalCutWeights = allocMemory<scalar_t>(maxCutNo);

  //for faster communication, concatanation of
  //totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
  //leftClosest distances sized P-1, since P-1 cut lines
  //rightClosest distances size P-1, since P-1 cut lines.
  scalar_t *totalPartWeights_leftClosests_rightClosests = allocMemory<scalar_t>((maxTotalPartCount + maxCutNo * 2) * concurrentPartCount);
  scalar_t *global_totalPartWeights_leftClosests_rightClosests = allocMemory<scalar_t>((maxTotalPartCount + maxCutNo * 2) * concurrentPartCount);

  partId_t *partIds = NULL;
  ArrayRCP<partId_t> partId;
  if(numLocalCoords > 0){
    partIds = allocMemory<partId_t>(numLocalCoords);
    partId = arcp(partIds, 0, numLocalCoords, true);
  }

  scalar_t *cutCoordinates =  allCutCoordinates;

  //partId_t leftPartitions = totalPartCount;
  scalar_t maxScalar_t = numeric_limits<scalar_t>::max();
  scalar_t minScalar_t = -numeric_limits<scalar_t>::max();

  env->timerStop(MACRO_TIMERS, "PQJagged Problem_Init");

  env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning");


  //int myRank = comm->getRank();
  //int worldSize = comm->getSize();

  scalar_t _EPSILON = numeric_limits<scalar_t>::epsilon();
  for (int i = 0; i < partArraySize; ++i){
    if(partNo[i] == 1) continue;
    int coordInd = i % coordDim;
    string istring = toString<int>(i);

    env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring);
    outTotalCounts = allocMemory<lno_t>(currentPartitionCount * partNo[i]);

    //the index where in the outtotalCounts will be written.
    partId_t currentOut = 0;
    //whatever is written to outTotalCounts will be added with previousEnd so that the points will be shifted.
    partId_t previousEnd = 0;
    scalar_t * pqCoord = pqJagged_coordinates[coordInd];

    partId_t currentWorkPart = 0;
    partId_t concurrentPart = min(currentPartitionCount - currentWorkPart, concurrentPartCount);
    bool useBinarySearch = false;
    if(partNo[i] > BINARYCUTOFF){
      useBinarySearch = true;
    }
    if(force_binary){
      useBinarySearch = true;
    }
    if (force_linear){
      useBinarySearch = false;
    }

//    size_t total_part_count = partNo[i] * 2 + 1;

    //run for all available parts.
    for (; currentWorkPart < currentPartitionCount; currentWorkPart += concurrentPart){

      concurrentPart = min(currentPartitionCount - currentWorkPart, concurrentPartCount);
#ifdef mpi_communication
      concurrent = concurrentPart;
#endif
      /*
      if(myRank == 0)
        cout << "i: " << i << " j:" << currentWorkPart << " ";
       */
      //scalar_t used_imbalance = imbalanceTolerance * (LEAF_IMBALANCE_FACTOR + (1 - LEAF_IMBALANCE_FACTOR)   / leftPartitions) * 0.7;
      scalar_t used_imbalance = 0;
      //env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring +"_min_max");

      for(int kk = 0; kk < concurrentPart; ++kk){

        partId_t currentPart = currentWorkPart + kk;
        lno_t coordinateEnd= inTotalCounts[currentPart];
        lno_t coordinateBegin = currentPart==0 ? 0: inTotalCounts[currentPart -1];

        pqJagged_getLocalMinMaxTotalCoord<scalar_t, lno_t>(
            partitionedPointCoordinates,
            pqCoord,
            pqJagged_uniformWeights[0],
            pqJagged_weights[0],
            numThreads,
            coordinateBegin,
            coordinateEnd,
            max_min_array,
            maxScalar_t,
            minScalar_t,
            localMinMaxTotal[kk], //min coordinate
            localMinMaxTotal[kk + concurrentPart], //max coordinate
            localMinMaxTotal[kk + 2*concurrentPart] //total weight);
                          );
      }
      pqJagged_getGlobalMinMaxTotalCoord<scalar_t>(comm,env, concurrentPart, localMinMaxTotal, globalMinMaxTotal);

      //env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring +"_min_max");




      partId_t allDone = 0;
      // Compute weight ratios for parts & cuts:
      //  e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
      //       part0  cut0  part1 cut1 part2 cut2 part3
      for(int kk = 0; kk < concurrentPart; ++kk){
        scalar_t minCoordinate = globalMinMaxTotal[kk];
        scalar_t maxCoordinate = globalMinMaxTotal[kk + concurrentPart];
        scalar_t *usedCutCoordinate = cutCoordinates + (partNo[i] - 1) * kk;
        scalar_t *usedCutPartRatios = targetPartWeightRatios + (partNo[i]) * kk;

        if(minCoordinate <= maxCoordinate){
          allDone += partNo[i] - 1;
          myNonDoneCount[kk] = partNo[i] - 1;
          //env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring + "_cut_coord");
          pqJagged_getCutCoord_Weights<scalar_t>(
              minCoordinate, maxCoordinate,
              pqJagged_uniformParts[0], pqJagged_partSizes[0], partNo[i] - 1,
              usedCutCoordinate, usedCutPartRatios, numThreads
          );
          //env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring + "_cut_coord");
          scalar_t coordinate_range = maxCoordinate - minCoordinate;

          partId_t currentPart = currentWorkPart + kk;
          lno_t coordinateEnd= inTotalCounts[currentPart];
          lno_t coordinateBegin = currentPart==0 ? 0: inTotalCounts[currentPart -1];

          if(ABS(coordinate_range) < _EPSILON ){
            for(lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
              partIds[partitionedPointCoordinates[ii]] = 0;
            }
          }
          else{

            scalar_t slice = coordinate_range / partNo[i];

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
            for(lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

              lno_t iii = partitionedPointCoordinates[ii];
              partId_t pp = partId_t((pqCoord[iii] - minCoordinate) / slice);
              //if( pp >= partNo[iii])
              //{
              //  partIds[iii] = 0;
              //} else {
              //  partIds[iii] = total_part_count - 2 * pp;
              //}
              partIds[iii] = 2 * pp;
            }
          }
        }
        else {
          // e.g., if have fewer coordinates than parts, don't need to do next dim.
          myNonDoneCount[kk] = 0;
        }
      }

      //cout << "partId:" << partIds[0] << endl;
      //env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring + "_1d");
      // Determine cut lines for k parts here.
      pqJagged_1D_Partition<scalar_t, lno_t>(
          env, comm,
          partitionedPointCoordinates,
          pqCoord, pqJagged_uniformWeights[0], pqJagged_weights[0],
          targetPartWeightRatios,
          globalMinMaxTotal,localMinMaxTotal,
          partNo[i], numThreads,
          maxScalar_t,
          minScalar_t,
          used_imbalance,
          currentWorkPart,
          concurrentPart,
          inTotalCounts,
          cutCoordinates,
          cutCoordinatesWork,
          leftClosestDistance,
          rightClosestDistance,
          cutUpperBounds,
          cutLowerBounds,
          cutUpperWeight,
          cutLowerWeight,
          isDone,
          partWeights,
          totalPartWeights_leftClosests_rightClosests,
          global_totalPartWeights_leftClosests_rightClosests,
          allowNonRectelinearPart,
          nonRectelinearPart,
          cutWeights,
          globalCutWeights,
          allDone,
          myNonDoneCount,
          useBinarySearch, // istring,
          partIds
          );


      //env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring + "_1d");

      //partId_t pEnd = currentOut + partNo[i] * concurrentPart;
      //for (partId_t ii = currentOut; ii < pEnd; ++ii){
      //  outTotalCounts[ii] = 0;
      //}

      //env->timerStart(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring + "_chunks");
      for(int kk = 0; kk < concurrentPart; ++kk){

        if(globalMinMaxTotal[kk] > globalMinMaxTotal[kk + concurrentPart]) {
        	for(partId_t jj = 0; jj < partNo[i]; ++jj){
        		outTotalCounts[currentOut + kk * (partNo[i]) + jj] = 0;
        	}

        	continue;
        }
        partId_t curr = currentWorkPart + kk;
        lno_t coordinateEnd= inTotalCounts[curr];
        lno_t coordinateBegin = curr==0 ? 0: inTotalCounts[curr -1];
        partId_t cutShift = (partNo[i] - 1) * kk;
        scalar_t *usedCutCoordinate = cutCoordinates + cutShift;
        float *usednonRectelinearPart = nonRectelinearPart + cutShift;

        scalar_t *tlr =  totalPartWeights_leftClosests_rightClosests + (4 *(partNo[i] - 1) + 1) * kk;

        for(int ii = 0; ii < numThreads; ++ii){
          pws[ii] = partWeights[ii] +  (2 * (partNo[i] - 1) + 1) * kk;
        }

        // Rewrite the indices based on the computed cuts.
        getChunksFromCoordinates<lno_t,scalar_t>(
            partNo[i],
            numThreads,

            partitionedPointCoordinates,
            pqCoord,
            pqJagged_uniformWeights[0],
            pqJagged_weights[0],

            usedCutCoordinate,
            coordinateBegin,
            coordinateEnd,

            allowNonRectelinearPart,
            usednonRectelinearPart,
            tlr,
            pws,
            nonRectRatios,

            //coordinate_linked_list,
            //coordinate_starts,
            //coordinate_ends,
            partPointCounts,

            newpartitionedPointCoordinates,
            outTotalCounts + currentOut + kk * (partNo[i]),partIds//,
            //useBinarySearch
        );
      }

      //env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring + "_chunks");
      cutCoordinates += (partNo[i] - 1) * concurrentPart;

      for(partId_t kk = 0; kk < concurrentPart; ++kk){


        for (partId_t ii = 0;ii < partNo[i] ; ++ii){
          outTotalCounts[currentOut+ii] += previousEnd;
          //cout << "out:" << outTotalCounts[currentOut+ii] << " prev:" << previousEnd<< endl;
        }
        previousEnd = outTotalCounts[currentOut + partNo[i] - 1];

        currentOut += partNo[i] ;
      }
/*
      if(myRank == 0)
        cout << endl;
*/
    } // end of this partitioning dimension

    // swap the indices' memory
    lno_t * tmp = partitionedPointCoordinates;
    partitionedPointCoordinates = newpartitionedPointCoordinates;
    newpartitionedPointCoordinates = tmp;

    currentPartitionCount *= partNo[i];
    freeArray<lno_t>(inTotalCounts);
    inTotalCounts = outTotalCounts;
    env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning_" + istring);
  } // Partitioning is done

  env->timerStop(MACRO_TIMERS, "PQJagged Problem_Partitioning");

  //cout << "biter" << endl;
  env->timerStart(MACRO_TIMERS, "PQJagged Part_Assignment");
/*
  partId_t *partIds = NULL;
  ArrayRCP<partId_t> partId;
  if(numLocalCoords > 0){
    partIds = allocMemory<partId_t>(numLocalCoords);
    partId = arcp(partIds, 0, numLocalCoords, true);
  }
*/

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
  for(partId_t i = 0; i < totalPartCount;++i){
    lno_t begin = 0;
    lno_t end = inTotalCounts[i];
    if(i > 0) begin = inTotalCounts[i -1];
    //cout << "part:" << i << " begin:" << begin << " end:" << end << " count:" << end - begin << endl;
    /*
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
     */
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
  }
  env->timerStop(MACRO_TIMERS, "PQJagged Part_Assignment");
  ArrayRCP<const gno_t> gnoList;
  if(numLocalCoords > 0){
    gnoList = arcpFromArrayView(pqJagged_gnos);
  }

  env->timerStart(MACRO_TIMERS, "PQJagged Solution_Part_Assignment");
  solution->setParts(gnoList, partId, true);
  env->timerStop(MACRO_TIMERS, "PQJagged Solution_Part_Assignment");


  env->timerStart(MACRO_TIMERS, "PQJagged Problem_Free");

/*
  if(comm->getRank() == 0){
    for(partId_t i = 0; i < totalPartCount - 1;++i){
      cout << "i:"<< i<<" cut coordinate:" << allCutCoordinates[i] << endl;
    }
  }
*/


  for(int i = 0; i < numThreads; ++i){
//    freeArray<lno_t>(coordinate_starts[i]);
//    freeArray<lno_t>(coordinate_ends[i]);
    freeArray<lno_t>(partPointCounts[i]);
  }

  freeArray<lno_t *>(partPointCounts);
  //freeArray<lno_t>(coordinate_linked_list);
  //freeArray<lno_t *>(coordinate_starts);
  //freeArray<lno_t *>(coordinate_ends);
  freeArray<double *> (pws);

  if(allowNonRectelinearPart){
    freeArray<float>(nonRectelinearPart);
    for(int i = 0; i < numThreads; ++i){
      freeArray<float>(nonRectRatios[i]);
    }
    freeArray<float *>(nonRectRatios);
  }

  freeArray<partId_t>(myNonDoneCount);
  freeArray<scalar_t>(cutWeights);
  freeArray<scalar_t>(globalCutWeights);
  freeArray<scalar_t>(max_min_array);
  freeArray<lno_t>(outTotalCounts);
  freeArray<lno_t>(partitionedPointCoordinates);
  freeArray<lno_t>(newpartitionedPointCoordinates);
  freeArray<scalar_t>(allCutCoordinates);
  freeArray<scalar_t *>(pqJagged_coordinates);
  freeArray<scalar_t *>(pqJagged_weights);
  freeArray<bool>(pqJagged_uniformParts);
  freeArray<scalar_t> (localMinMaxTotal);
  freeArray<scalar_t> (globalMinMaxTotal);
  freeArray<scalar_t *>(pqJagged_partSizes);
  freeArray<bool>(pqJagged_uniformWeights);
  freeArray<scalar_t>(cutCoordinatesWork);
  freeArray<scalar_t>(targetPartWeightRatios);
  freeArray<scalar_t>(cutUpperBounds);
  freeArray<scalar_t>(cutLowerBounds);
  freeArray<scalar_t>(cutLowerWeight);
  freeArray<scalar_t>(cutUpperWeight);
  freeArray<bool>(isDone);
  freeArray<scalar_t>(totalPartWeights_leftClosests_rightClosests);
  freeArray<scalar_t>(global_totalPartWeights_leftClosests_rightClosests);

  for(int i = 0; i < numThreads; ++i){
    freeArray<double>(partWeights[i]);//freeArray<scalar_t>(partWeights[i]);
    freeArray<scalar_t>(rightClosestDistance[i]);
    freeArray<scalar_t>(leftClosestDistance[i]);
  }

  //freeArray<scalar_t *>(partWeights);
  freeArray<double *>(partWeights);
  freeArray<scalar_t *>(leftClosestDistance);
  freeArray<scalar_t *>(rightClosestDistance);


  env->timerStop(MACRO_TIMERS, "PQJagged Problem_Free");
  env->timerStop(MACRO_TIMERS, "PQJagged Total");

#endif // INCLUDE_ZOLTAN2_EXPERIMENTAL
}
} // namespace Zoltan2





#endif
