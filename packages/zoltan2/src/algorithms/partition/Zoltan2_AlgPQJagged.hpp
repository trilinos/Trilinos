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

#include <Zoltan2_PQJagged_ReductionOps.hpp>
#include <Zoltan2_AlgRCB_methods.hpp> // TODO: Needed for RCB params, not sure
                                      // why they are needed here.
#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_Metric.hpp>             // won't need thiss
#include <Zoltan2_Parameters.hpp>
#include <Tpetra_Distributor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Zoltan2_CoordinatePartitioningGraph.hpp>
#include <new>          // ::operator new[]
#include <algorithm>    // std::sort
#include <Zoltan2_Util.hpp>


#ifdef HAVE_ZOLTAN2_ZOLTAN
#ifdef HAVE_ZOLTAN2_MPI
#define enable_migration2
#endif
#endif

#ifdef enable_migration2
#include "zoltan_comm_cpp.h"
#include "zoltan_types.h" // for error codes
#endif
//#define FIRST_TOUCH
//#define BINARYCUTSEARCH
//#define Zoltan_Comm
//
// TODO list: NUM_RCB_PARAMS cannot be used here. Should be NUM_MJ_PARAMS if 
// needed.
//  

#include <bitset>

#define EPS_SCALE 1
#define LEAST_SIGNIFICANCE 0.0001
#define SIGNIFICANCE_MUL 1000
#define FUTURE_REDUCEALL_CUTOFF 1500000
#define MIN_WORK_LAST_DIM 1000
//#define INCLUDE_ZOLTAN2_EXPERIMENTAL
#ifdef HAVE_ZOLTAN2_OMP
#include <omp.h>
#endif
#define ABS(x) ((x) >= 0 ? (x) : -(x))

#define LEAF_IMBALANCE_FACTOR 0.1
#define BINARYCUTOFF 32

//imbalance calculation. Wreal / Wexpected - 1
#define imbalanceOf(Wachieved, totalW, expectedRatio) \
        (Wachieved) / ((totalW) * (expectedRatio)) - 1
#define imbalanceOf2(Wachieved, wExpected) \
        (Wachieved) / (wExpected) - 1
#define MAXOF (a,b) (((a)>(b))?(a):(b))

#define KCUTOFF 0.80
#define forceMigration 1500000
#define Z2_DEFAULT_CON_PART_COUNT 16

using std::vector;

namespace Zoltan2{

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

#if 0
/*! \brief A helper class containing array representation of
 *  coordinate linked lists.
 */

template <typename pq_lno_t, typename partId_t>
class pqJagged_PartVertices{
private:
    pq_lno_t *linkedList; //initially filled with -1's.
    pq_lno_t *partBegins; //initially filled with -1's.
    pq_lno_t *partEnds; //initially filled with -1's.
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

    void set(pq_lno_t *linkedList_, pq_lno_t *partBegins_, pq_lno_t *partEnds_){
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
    void inserToPart (partId_t partNo, pq_lno_t coordinateIndex){

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
    pq_lno_t *getLinkedList(){ return linkedList;}

    /*! \brief
     * partBegins getter function.
     */
    pq_lno_t *getPartBegins(){ return partBegins;}
    /*! \brief
     * partEnds getter function.
     */
    pq_lno_t *getPartEnds(){ return partEnds;}

};
#endif

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
template <typename pq_scalar_t>
inline pq_scalar_t pivotPos (
    pq_scalar_t * cutUpperBounds,
    pq_scalar_t *cutLowerBounds,
    size_t currentCutIndex,
    pq_scalar_t *cutUpperWeight,
    pq_scalar_t *cutLowerWeight,
    pq_scalar_t ew,
    pq_scalar_t _EPSILON){

    if(ABS(cutUpperWeight[currentCutIndex] - cutLowerWeight[currentCutIndex])
         < _EPSILON){
        return cutLowerBounds[currentCutIndex];
    }

    pq_scalar_t coordinate_range = (cutUpperBounds[currentCutIndex] - cutLowerBounds[currentCutIndex]);
    //pq_scalar_t weight_range = (cutUpperWeight[currentCutIndex] - cutLowerWeight[currentCutIndex]);

    //pq_scalar_t myWeightDif = (ew - cutLowerWeight[currentCutIndex]);
    pq_scalar_t newCut =(coordinate_range ) / 2 /** (myWeightDif / weight_range)*/ + cutLowerBounds[currentCutIndex];
    /*
    cout << "cutIndex:" << currentCutIndex <<
            " upper:" << cutUpperBounds[currentCutIndex] << " uw:" << cutUpperWeight[currentCutIndex] <<
            " lower:" << cutLowerBounds[currentCutIndex] << " lw:" << cutLowerWeight[currentCutIndex] <<
            " found:" << newCut << endl;
            */
    return newCut;
}

template <typename T>
void get_partitioning_params(
   const Teuchos::ParameterList &pl,
   T &imbalanceTolerance,
   multiCriteriaNorm &mcnorm,
   std::bitset<NUM_RCB_PARAMS> &params,
   int &numTestCuts,
   bool &ignoreWeights){

    string obj;
    const Teuchos::ParameterEntry *pe = pl.getEntryPtr
                                                ("partitioning_objective");
    if (pe)
        obj = pe->getValue(&obj);

    if (!pe){
        params.set(rcb_balanceWeight);
        mcnorm = normBalanceTotalMaximum;
    }
    else if (obj == string("balance_object_count")){
        params.set(rcb_balanceCount);
    }
#if 0
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
#endif
    else{
        params.set(rcb_balanceWeight);
        mcnorm = normBalanceTotalMaximum;
    }

#if 0
    int val = 0;
    pe = pl.getEntryPtr("average_cuts");
    if (pe)
        val = pe->getValue(&val);

    if (val == 1)
        params.set(rcb_averageCuts);
#endif

    // TODO: Imbalance tolerance and ignoreWeights are the two parameters
    // that are used. ignoreWeights derives from rcb_balanceCount which 
    // never seem to be set, making me wonder whether it is used at all.
    // MD: ignoreWeights is not used. I copied this part 
    // from Lee Ann's code. Left it in this way, in case needed in future.
    imbalanceTolerance = .1;
    pe = pl.getEntryPtr("imbalance_tolerance");
    if (pe){
        double tol;
        tol = pe->getValue(&tol);
        imbalanceTolerance = tol - 1.0;
    }

    if (imbalanceTolerance <= 0)
        imbalanceTolerance = 10e-4;

#if 0
    numTestCuts = 1;
    pe = pl.getEntryPtr("bisection_num_test_cuts");
    if (pe)
        numTestCuts = pe->getValue(&numTestCuts);
#endif

    ignoreWeights = params.test(rcb_balanceCount);
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
void pqJagged_getParameters(
        const Teuchos::ParameterList &pl,
        bool &allowNonrectilinear, 
        partId_t &concurrentPartCount,
        int &migration_actualMigration_option,
        int &migration_check_option,
        int &migration_all2all_option,
        T &migration_imbalance_cut_off, 
        int &migration_assignment_type,
        int &keep_part_boxes,
        int &enable_rcb,
        int &recursion_depth){


    // TODO: Append all the parameters with mj_
    const Teuchos::ParameterEntry *pe = pl.getEntryPtr("partitioning_objective");
    migration_imbalance_cut_off = 0.35;
    pe = pl.getEntryPtr("migration_imbalance_cut_off");
    if (pe){
        double tol;
        tol = pe->getValue(&tol);
        migration_imbalance_cut_off = tol - 1.0;
    }

    pe = pl.getEntryPtr("migration_all_to_all_type");
    if (pe){
        migration_all2all_option = pe->getValue(&concurrentPartCount);
    }else {
        migration_all2all_option = 1;
    }

    pe = pl.getEntryPtr("migration_check_option");
    if (pe){
        migration_check_option = pe->getValue(&concurrentPartCount);
    }else {
        migration_check_option = 0;
    }
    if (migration_check_option > 1) migration_check_option = -1;

    pe = pl.getEntryPtr("migration_processor_assignment_type");
    if (pe){
        migration_assignment_type = pe->getValue(&concurrentPartCount);
    }else {
        migration_assignment_type = 0;
    }
    pe = pl.getEntryPtr("migration_doMigration_type");
    if (pe){
        migration_actualMigration_option = pe->getValue(&migration_actualMigration_option);
    }else {
        migration_actualMigration_option = 1;
    }

    pe = pl.getEntryPtr("parallel_part_calculation_count");
    if (pe){
        concurrentPartCount = pe->getValue(&concurrentPartCount);
    }else {
        concurrentPartCount = 0; // Set to invalid value
    }

    pe = pl.getEntryPtr("keep_part_boxes");
    if (pe){
        keep_part_boxes = pe->getValue(&keep_part_boxes);
    }else {
        keep_part_boxes = 0; // Set to invalid value
    }

    if (keep_part_boxes == 0){
        pe = pl.getEntryPtr("mapping_type");
        if (pe){
            int mapping_type = -1;
            mapping_type = pe->getValue(&mapping_type);
            if (mapping_type == 0){
                keep_part_boxes  = 1;
            }
        }
    }

    pe = pl.getEntryPtr("mj_enable_rcb");
    if (pe){
        enable_rcb = pe->getValue(&enable_rcb);
    }else {
        enable_rcb = 0; // Set to invalid value
    }

    pe = pl.getEntryPtr("recursion_depth");
    if (pe){
        recursion_depth = pe->getValue(&recursion_depth);
    }else {
        recursion_depth = -1; // Set to invalid value
    }

    int val = 0;
    pe = pl.getEntryPtr("rectilinear_blocks");
    if (pe)
        val = pe->getValue(&val);

    if (val == 1){
        allowNonrectilinear = false;
    } else {
        allowNonrectilinear = true;
    }


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
        int &weightDim, size_t &numLocalCoords, global_size_t &numGlobalCoords,
        int &criteriaDim, const bool &ignoreWeights){
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
template <typename Adapter, typename pq_scalar_t, typename pq_gno_t>
void pqJagged_getInputValues(
    const RCP<const Environment> &env, const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &coords,
    RCP<PartitioningSolution<Adapter> > &solution,
    std::bitset<NUM_RCB_PARAMS> &params,
    const int &coordDim,
    const int &weightDim,
    const size_t &numLocalCoords,
    size_t &numGlobalParts,
    int &pqJagged_multiVectorDim,
    pq_scalar_t **pqJagged_values,
    const int &criteriaDim,
    pq_scalar_t **pqJagged_weights,
    ArrayView<const pq_gno_t> &pqJagged_gnos,
    bool &ignoreWeights,
    bool *pqJagged_uniformWeights,
    bool *pqJagged_uniformParts,
    pq_scalar_t **pqJagged_partSizes
){
    //typedef typename Adapter::node_t pq_node_t;
    typedef typename Adapter::lno_t pq_lno_t;
    typedef StridedData<pq_lno_t, pq_scalar_t> input_t;

    ArrayView<const pq_gno_t> gnos;
    ArrayView<input_t>     xyz;
    ArrayView<input_t>     wgts;

    coords->getCoordinates(gnos, xyz, wgts);
    pqJagged_gnos = gnos;

    for (int dim=0; dim < coordDim; dim++){
        ArrayRCP<const pq_scalar_t> ar;
        xyz[dim].getInputArray(ar);
        //pqJagged coordinate values assignment
        pqJagged_values[dim] =  (pq_scalar_t *)ar.getRawPtr();
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
                ArrayRCP<const pq_scalar_t> ar;
                wgts[wdim].getInputArray(ar);
                pqJagged_uniformWeights[wdim] = false;
                pqJagged_weights[wdim] = (pq_scalar_t *) ar.getRawPtr();
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
            //TODO : Need a test for this, 
            pq_scalar_t *tmp = new pq_scalar_t [numGlobalParts];
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

    // TODO: Do we handle multidimensional weights.
    // MD: no. We can only handle single dimensional weights.
    if (criteriaDim > 1){
        for (int wdim1 = 0; wdim1 < criteriaDim; wdim1++)
            for (int wdim2 = wdim1+1; wdim2 < criteriaDim; wdim2++)
                if (!solution->criteriaHaveSamePartSizes(wdim1, wdim2)){
                    multiplePartSizeSpecs = true;
                    break;
                }
    }

    // TODO: We do not use this at all, should be removed. 
    //MD: Yes. Again, I left this in case used in future.
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
template <typename pq_scalar_t, typename pq_gno_t>
void pqJagged_printInput(
    int coordDim,
    int weightDim,
    size_t numLocalCoords,
    global_size_t numGlobalCoords,
    int criteriaDim,
    pq_scalar_t **pqJagged_values,
    pq_scalar_t **pqJagged_weights,
    bool *pqJagged_uniformParts,
    bool *pqJagged_uniformWeights,
    pq_gno_t *pqJagged_gnos,
    bool ignoreWeights,
    size_t numGlobalParts,
    pq_scalar_t **pqJagged_partSizes
){

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

    std::cout << "pqJagged_uniformWeights:" << pqJagged_uniformWeights[0] <<
        std::endl;
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

    std::cout << "pqJagged_uniformParts:" << pqJagged_uniformParts[0] <<
        std::endl;
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
template <typename pq_scalar_t, typename pq_lno_t>
void pqJagged_getLocalMinMaxTotalCoord(
    pq_lno_t *partitionedPointPermutations,
    pq_scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    pq_scalar_t *pqJagged_weights,
    int numThreads,
    pq_lno_t coordinateBegin,
    pq_lno_t coordinateEnd,
    pq_scalar_t *max_min_array /*sized nothreads * 2*/,
    pq_scalar_t maxScalar,
    pq_scalar_t minScalar,
    pq_scalar_t &minCoordinate,
    pq_scalar_t &maxCoordinate,
    pq_scalar_t &totalWeight
){

    //if the part is empty.
    //set the min and max coordinates as reverse.
    if(coordinateBegin >= coordinateEnd)
    {
        minCoordinate = maxScalar;
        maxCoordinate = minScalar;
        totalWeight = 0;
    }
    else {
        pq_scalar_t mytotalWeight = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
        {
            //if uniform weights are used, then weight is equal to count.
            if (pqJagged_uniformWeights) {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
                {
                    mytotalWeight = coordinateEnd - coordinateBegin;
                }

            }
            else {
                //if not uniform, then weights are reducted from threads.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for reduction(+:mytotalWeight)
#endif
                for (pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
                    int i = partitionedPointPermutations[ii];
                    mytotalWeight += pqJagged_weights[i];
                }
            }

            int myId = 0;
#ifdef HAVE_ZOLTAN2_OMP
            myId = omp_get_thread_num();
#endif
            pq_scalar_t myMin, myMax;
/*
            problemComm->barrier();
            cout << "initial me:" << problemComm->getRank()
                    << " coordinateBegin:" << coordinateBegin
                    << " ind:" << partitionedPointPermutations[coordinateBegin] << endl;
*/
            myMin=myMax
                =pqJagged_coordinates[partitionedPointPermutations[coordinateBegin]];
//            problemComm->barrier();
 //           cout << "initial me:" << problemComm->getRank() << " myMin:" << myMin << " myMax:" << myMax << endl;


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
            for(pq_lno_t j = coordinateBegin + 1; j < coordinateEnd; ++j){
                int i = partitionedPointPermutations[j];
                if(pqJagged_coordinates[i] > myMax)
                    myMax = pqJagged_coordinates[i];
                if(pqJagged_coordinates[i] < myMin)
                    myMin = pqJagged_coordinates[i];
            }
            max_min_array[myId] = myMin;
            max_min_array[myId + numThreads] = myMax;
/*
            problemComm->barrier();
            cout << "after me:" << problemComm->getRank() << " myMin:" << myMin << " myMax:" << myMax << endl;
*/

#ifdef HAVE_ZOLTAN2_OMP
//we need a barrier here, because max_min_array might not be filled by some of the threads.
#pragma omp barrier
#pragma omp single nowait
#endif
            {
                minCoordinate = max_min_array[0];
                for(int i = 1; i < numThreads; ++i){
                    if(max_min_array[i] < minCoordinate)
                        minCoordinate = max_min_array[i];
                }
            }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single nowait
#endif
            {
                maxCoordinate = max_min_array[numThreads];
                for(int i = numThreads + 1; i < numThreads * 2; ++i){
                    if(max_min_array[i] > maxCoordinate)
                        maxCoordinate = max_min_array[i];
                }
            }
        }
        totalWeight = mytotalWeight;
    }

}

// fEpsilon is different for double and float. It is found once
// and passed multiple times.
template <typename partId_t>
inline partId_t getPartCount(partId_t numFuture, double root, double fEpsilon){
    double fp = pow(numFuture, root);
    partId_t ip = partId_t (fp);
    if (fp - ip < fEpsilon * 100){
        return ip;
    }
    else {
        return ip  + 1;
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
template <typename pq_scalar_t>
void pqJagged_getGlobalMinMaxTotalCoord(
    RCP<Comm<int> > &comm,
    const RCP<const Environment> &env,
    partId_t concurrentPartCount,
    pq_scalar_t *localMinMaxTotal,
    pq_scalar_t *globalMinMaxTotal){

    //reduce min for first concurrentPartCount elements, reduce max for next
    //concurrentPartCount elements,
    //reduce sum for the last concurrentPartCount elements.
    if(comm->getSize()  > 1){

#ifndef mpi_communication
        Teuchos::PQJaggedCombinedMinMaxTotalReductionOp<int, pq_scalar_t>
         reductionOp( concurrentPartCount, concurrentPartCount,
                concurrentPartCount);
#endif

#ifdef mpi_communication
        MPI_Op myop;
        MPI_Op_create(minMaxSum, 0, &myop);   /* step 3 */
#endif

        try{
#ifdef mpi_communication
	    //TODO: The parts insides the mpi_communication are not run, and should
        //be removed.
            MPI_Allreduce(localMinMaxTotal, globalMinMaxTotal,
            3 * concurrentPartCount, MPI_FLOAT, myop,MPI_COMM_WORLD);
#endif
#ifndef mpi_communication
            reduceAll<int, pq_scalar_t>(*comm, reductionOp,
                    3 * concurrentPartCount, localMinMaxTotal,
                    globalMinMaxTotal);
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
 *
 * \param pqJagged_uniformParts is a boolean value holding whether the desired partitioning is uniform.
 * \param pqJagged_uniformWeight is a boolean value holding whether the coordinates have unit weights.
 *
 * \param pqJagged_partSizes holds the desired parts sizes if desired partitioning is not uniform.
 * \param noCuts holds the number of cuts in the current partitioning dimension.
 * \param globalTotalWeight holds the global total weight in the current part.
 *
 * \param cutCoordinates is the output array for the initial cut lines.
 * \param cutPartRatios is the output array holding the cumulative ratios of parts in current partitioning.
 * For partitioning to 4 uniformly, cutPartRatios will be (0.25 * globalTotalWeight, 0.5 *globalTotalWeight , 0.75 * globalTotalWeight, globalTotalWeight).
 * \param numThreads hold the number of threads available per mpi.
 *
 * \param currentPartitions is the vector that holds how many more parts each part will be divided into more
 * for the parts at the beginning of this coordinate partitioning
 * \param futurePartitions is the vector that holds how many more parts each part will be divided into more
 * for the parts that will be obtained at the end of this coordinate partitioning.
 * \param partIndex is the index of the part in the currentPartitions vector.
 * \param futureArrayIndex holds the amount of shift in the futurePartitions for the output parts.
 */
template <typename pq_scalar_t>
void pqJagged_getCutCoord_Weights(
    pq_scalar_t minCoordinate,
    pq_scalar_t maxCoordinate,
    bool pqJagged_uniformParts,
    bool pqJagged_uniformWeights,
    pq_scalar_t *pqJagged_partSizes /*p sized, weight ratios of each part*/,
    partId_t noCuts/*p-1*/ ,
    pq_scalar_t globalTotalWeight,
    pq_scalar_t *cutCoordinates /*p - 1 sized, coordinate of each cut line*/,
    pq_scalar_t *cutPartRatios /*cumulative weight ratios, at left side of each cut line. p-1 sized*/,
    int numThreads,
    vector <partId_t> *currentPartitions, //the vecto
    vector <partId_t> *futurePartitions,
    partId_t partIndex,
    partId_t futureArrayIndex
){

    pq_scalar_t coordinateRange = maxCoordinate - minCoordinate;
    if(pqJagged_uniformParts){
        {
            partId_t cumulative = 0;
            pq_scalar_t totalInnerPartCount = pq_scalar_t((*currentPartitions)[partIndex]);
            pq_scalar_t unitWeight = globalTotalWeight / totalInnerPartCount;
            for(partId_t i = 0; i < noCuts; ++i){
                cumulative += (*futurePartitions)[i + futureArrayIndex];
                //cutPartRatios[i] = (cumulative /*+  (*futurePartitions)[i + futureArrayIndex]*/) / (totalInnerPartCount);
                cutPartRatios[i] = cumulative * unitWeight;
                cutCoordinates[i] = minCoordinate + (coordinateRange *
                                         cumulative) / totalInnerPartCount;
            }
            cutPartRatios[noCuts] = 1;
        }
        if (pqJagged_uniformWeights){
            for(partId_t i = 0; i < noCuts + 1; ++i){
                cutPartRatios[i] = long(cutPartRatios[i] + 0.5);
            }
        }
    }
    else {
        /*
        //TODO fix here!!
        cutPartRatios[0] = pqJagged_partSizes[0];
        cutCoordinates[0] = coordinateRange * cutPartRatios[0];
        for(partId_t i = 1; i < noCuts; ++i){
            cutPartRatios[i] = pqJagged_partSizes[i] + cutPartRatios[i - 1];
            cutCoordinates[i] = coordinateRange * cutPartRatios[i];
        }
        */
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
 * \param maxCoordinate is the maximum coordinate in the current part.
 * \param minCoordinate is the minumum coordinate in the current part.
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
template <typename pq_scalar_t>
void getNewCoordinates(
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,
    const size_t &total_part_count,
    const partId_t &noCuts,
    const pq_scalar_t &maxCoordinate,
    const pq_scalar_t &minCoordinate,
    const pq_scalar_t &globalTotalWeight,
    const pq_scalar_t &imbalanceTolerance,
    pq_scalar_t * globalPartWeights,
    const pq_scalar_t * localPartWeights,
    const pq_scalar_t *targetPartWeightRatios,
    bool *isDone,
    pq_scalar_t *cutCoordinates,
    pq_scalar_t *cutUpperBounds,
    pq_scalar_t *cutLowerBounds,
    pq_scalar_t *leftClosestDistance,
    pq_scalar_t *rightClosestDistance,
    pq_scalar_t * cutLowerWeight,
    pq_scalar_t * cutUpperWeight,
    pq_scalar_t *newCutCoordinates,
    bool allowNonRectelinearPart,
    float *nonRectelinearPartRatios,
    partId_t *rectilinearCutCount,
    pq_scalar_t *localCutWeights,
    pq_scalar_t *globalCutWeights,
    partId_t &myNoneDoneCount
)
{
    pq_scalar_t seenW = 0;
    float expected = 0;
    pq_scalar_t leftImbalance = 0, rightImbalance = 0;
    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon();
    //pq_scalar_t _EPSILON = numeric_limits<float>::epsilon();

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
    for (partId_t i = 0; i < noCuts; i++){
        //if a left and right closes point is not found, set the distance to 0.
        if(minCoordinate - leftClosestDistance[i] > _EPSILON)
            leftClosestDistance[i] = cutCoordinates[i];
        if(rightClosestDistance[i] - maxCoordinate > _EPSILON)
            rightClosestDistance[i] = cutCoordinates[i];

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
        //leftImbalance = imbalanceOf(seenW, globalTotalWeight, expected);
        leftImbalance = imbalanceOf2(seenW, expected);
        //rightImbalance = imbalanceOf(globalTotalWeight - seenW, globalTotalWeight, 1 - expected);
        rightImbalance = imbalanceOf2(globalTotalWeight - seenW, globalTotalWeight - expected);

        bool isLeftValid = ABS(leftImbalance) - imbalanceTolerance < _EPSILON ;
        bool isRightValid = ABS(rightImbalance) - imbalanceTolerance < _EPSILON;

/*
        cout << "\t\tc:" << i << "leftImbalance:" << leftImbalance <<
                " seenW:" << seenW <<
                " lineW:" <<  globalPartWeights[i * 2 + 1] <<
                " globalTotalWeight:" << globalTotalWeight <<
                " expected:" << expected <<
                " r:" << rightImbalance <<
                " upper:" <<cutUpperBounds[i] <<
                " lower:" << cutLowerBounds[i] <<
                " leftValid:" << isLeftValid <<
                " rightValid:" << isRightValid <<
                " imbalanceTolerance:" << imbalanceTolerance <<
                endl;
*/
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

            pq_scalar_t ew = expected; //globalTotalWeight * expected;
            if(allowNonRectelinearPart){

                if (globalPartWeights[i * 2 + 1] == ew){

                    isDone[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
                    myNoneDoneCount -= 1;
                    newCutCoordinates [i] = cutCoordinates[i];

                    //nonRectelinearPartRatios[i] = 1;
                    nonRectelinearPartRatios[i] = localPartWeights[i * 2 + 1] - localPartWeights[i * 2];
                    /*
                    cout << "setting i" << i <<
                            " to:" << localPartWeights[i * 2 + 1] <<
                            " - " << localPartWeights[i * 2] <<
                            " = " <<  nonRectelinearPartRatios[i] <<
                            endl;
                            */
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
                    pq_scalar_t myWeightOnLine = localPartWeights[i * 2 + 1] -
                                                     localPartWeights[i * 2];
                    localCutWeights[i] = myWeightOnLine;
                    continue;
                }
            }
            //when moving right, set lower bound to current line.
            cutLowerBounds[i] = /*cutCoordinates[i] + */rightClosestDistance[i];
            cutLowerWeight[i] = seenW;

            //compare the upper bound with the current lines.
            for (partId_t ii = i + 1; ii < noCuts ; ++ii){
                pq_scalar_t pw = globalPartWeights[ii * 2];
                pq_scalar_t lw = globalPartWeights[ii * 2 + 1];
                if(pw >= ew){
                    if(pw == ew){
                        cutUpperBounds[i] = cutCoordinates[ii];
                        cutUpperWeight[i] = pw;
                        cutLowerBounds[i] = cutCoordinates[ii];
                        cutLowerWeight[i] = pw;
                    } else if (pw < cutUpperWeight[i]){
                        //if a cut line is more strict than the current upper
                        //bound update the upper bound.
                        cutUpperBounds[i] = /*cutCoordinates[ii] - */leftClosestDistance[ii];
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
                    cutLowerBounds[i] = /*cutCoordinates[ii] +*/ rightClosestDistance[ii] ;
                    cutLowerWeight[i] = pw;
                }
            }


            pq_scalar_t newPivot = pivotPos<pq_scalar_t> (cutUpperBounds,
             cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew, _EPSILON);

            //if cut line does not move significantly.
            if (ABS(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE ||
                 cutLowerBounds[i] - cutUpperBounds[i] > _EPSILON
                 /*cutUpperBounds[i] < cutLowerBounds[i]*/){
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
            pq_scalar_t ew = expected; //globalTotalWeight * expected;
            //moving left, set upper to current line.
            //cout << "setting upper bound to:" << leftClosestDistance[i] << endl;
            cutUpperBounds[i] = /*cutCoordinates[i] -*/ leftClosestDistance[i];
            cutUpperWeight[i] = seenW;

            // compare the current cut line weights with previous upper and lower bounds.
            for (int ii = i - 1; ii >= 0; --ii){
                pq_scalar_t pw = globalPartWeights[ii * 2];
                pq_scalar_t lw = globalPartWeights[ii * 2 + 1];
                if(pw <= ew){
                    if(pw == ew){
                        cutUpperBounds[i] = cutCoordinates[ii];
                        cutUpperWeight[i] = pw;
                        cutLowerBounds[i] = cutCoordinates[ii];
                        cutLowerWeight[i] = pw;
                    }
                    else if (pw > cutLowerWeight[i]){
                        cutLowerBounds[i] = /*cutCoordinates[ii] +*/ rightClosestDistance[ii];
                        cutLowerWeight[i] = pw;
                        if(lw > ew){
                            cutUpperBounds[i] = /*cutCoordinates[ii] +*/ rightClosestDistance[ii];

                            cutUpperWeight[i] = lw;
                        }
                    }
                    break;
                }
                if (pw >= ew && (pw < cutUpperWeight[i] || (pw == cutUpperWeight[i] && cutUpperBounds[i] > /*cutCoordinates[ii] - */leftClosestDistance[ii]))){
                    cutUpperBounds[i] = /*cutCoordinates[ii] - */leftClosestDistance[ii] ;

                    cutUpperWeight[i] = pw;
                }
            }

            pq_scalar_t newPivot = pivotPos<pq_scalar_t> (cutUpperBounds, cutLowerBounds,i, cutUpperWeight, cutLowerWeight, ew, _EPSILON);
            //if cut line does not move significantly.
            if (ABS(cutCoordinates[i] - newPivot) < _EPSILON * EPS_SCALE  || cutLowerBounds[i] - cutUpperBounds[i] > _EPSILON ){
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

    //communication to determine the ratios of processors for the distribution
    //of coordinates on the cut lines.
#ifdef HAVE_ZOLTAN2_OMP
    //#pragma omp barrier
#pragma omp single
#endif
    {
        if(*rectilinearCutCount > 0){
            try{
                Teuchos::scan<int,pq_scalar_t>(
                        *comm, Teuchos::REDUCE_SUM,
                        noCuts,localCutWeights, globalCutWeights
                );
            }
            Z2_THROW_OUTSIDE_ERROR(*env)

            for (partId_t i = 0; i < noCuts; ++i){
                //cout << "gw:" << globalCutWeights[i] << endl;
                if(globalCutWeights[i] > 0) {
                    //pq_scalar_t ew = globalTotalWeight * targetPartWeightRatios[i];
                    pq_scalar_t ew = targetPartWeightRatios[i];
                    //ew = long (ew * long(1000)) / pq_scalar_t(1000.0);
                    //globalPartWeights[i * 2] = long(globalPartWeights[i * 2] * long(1000)) / pq_scalar_t (1000);
                    pq_scalar_t expectedWeightOnLine = ew - globalPartWeights[i * 2];
                    pq_scalar_t myWeightOnLine = localCutWeights[i];
                    pq_scalar_t weightOnLineBefore = globalCutWeights[i];
                    pq_scalar_t incMe = expectedWeightOnLine - weightOnLineBefore;
                    pq_scalar_t mine = incMe + myWeightOnLine;
                    if(mine < 0){
                        nonRectelinearPartRatios[i] = 0;
/*
                        cout << "setting i" << i <<
                                " to:" <<  nonRectelinearPartRatios[i] <<
                                endl;
*/
                    }
                    else if(mine >= myWeightOnLine){
                        //nonRectelinearPartRatios[i] = 1;
                        nonRectelinearPartRatios[i] = myWeightOnLine;
/*
                        cout << "setting i" << i <<
                                " to:" << myWeightOnLine <<
                                " myWeight= " <<  nonRectelinearPartRatios[i] <<
                                endl;
*/
                    }
                    else {
                        //nonRectelinearPartRatios[i] = mine / myWeightOnLine;
                        nonRectelinearPartRatios[i] = mine ;
/*
                        cout << "setting i" << i <<
                                " to:" << mine <<
                                " mine= " <<  nonRectelinearPartRatios[i] <<
                                " expectedWeightOnLine:" << expectedWeightOnLine <<
                                " ratio:" << targetPartWeightRatios[i] <<
                                " ew:" << ew << " globalPartWeights[i * 2]:" << 
                                globalPartWeights[i * 2] <<
                                " weightOnLineBefore:" << weightOnLineBefore <<
                                " incMe:" << incMe <<
                                " mine:" << mine <<
                                endl;
*/
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
 * \param maxCoordinate is the maximum coordinate in the part.
 * \param minCoordinate is the min coordinate in the part.
 * \param _EPSILON is the smallest error value for pq_scalar_t.
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
 * \param myLeftClosest is the array holding the coordinate of the closest points to the cut lines from left for the calling thread..
 * \param myRightClosest is the array holding the coordinate of the closest points to the cut lines from right for the calling thread.
 * \param useBinarySearch is boolean parameter whether to search for cut lines with binary search of linear search.
 * \param partIds is the array that holds the part ids of the coordinates
 */
template <typename pq_scalar_t, typename pq_lno_t>
void pqJagged_1DPart_getPartWeights(
    size_t total_part_count,
    partId_t noCuts,
    pq_scalar_t maxCoordinate,
    pq_scalar_t minCoordinate,
    pq_scalar_t _EPSILON,
    int numThreads,
    pq_lno_t coordinateBegin,
    pq_lno_t coordinateEnd,
    pq_lno_t *partitionedPointPermutations,
    pq_scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    pq_scalar_t *pqJagged_weights,
    pq_scalar_t *cutCoordinates_tmp, //TODO change name
    bool *isDone,
    double *myPartWeights,
    pq_scalar_t *myLeftClosest,
    pq_scalar_t *myRightClosest,
    bool useBinarySearch,
    partId_t *partIds
){

    // initializations for part weights, left/right closest
    for (size_t i = 0; i < total_part_count; ++i){
        myPartWeights[i] = 0;
    }

    for(partId_t i = 0; i < noCuts; ++i){
        //if(isDone[i]) continue;
        //myLeftClosest[i] = maxScalar;
        myLeftClosest[i] = minCoordinate - 1;
        myRightClosest[i] = maxCoordinate + 1;
    }
    if(useBinarySearch){
        //pq_lno_t comparison_count = 0;
        pq_scalar_t minus_EPSILON = -_EPSILON;
#ifdef HAVE_ZOLTAN2_OMP
//no need for the barrier as all threads uses their local memories.
#pragma omp for
#endif
        for (pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
            int i = partitionedPointPermutations[ii];
            partId_t j = partIds[i] / 2;

            if(j >= noCuts){
                j = noCuts - 1;
            }

            partId_t lc = 0;
            partId_t uc = noCuts - 1;

            pq_scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
            bool isInserted = false;
            bool onLeft = false;
            bool onRight = false;
            partId_t lastPart = -1;

            pq_scalar_t coord = pqJagged_coordinates[i];

            while(uc >= lc)
            {
                //comparison_count++;
                lastPart = -1;
                onLeft = false;
                onRight = false;
                pq_scalar_t cut = cutCoordinates_tmp[j];
                pq_scalar_t distance = coord - cut;
                pq_scalar_t absdistance = ABS(distance);

                if(absdistance < _EPSILON){

                    myPartWeights[j * 2 + 1] += w;
                    partIds[i] = j * 2 + 1;
                    //cout << "\ti:" << i << " assigned:" << j * 2 + 1 << endl;

                    myLeftClosest[j] = coord;
                    myRightClosest[j] = coord;
                    partId_t kk = j + 1;
                    while(kk < noCuts){
                        // Needed when cuts shared the same position
                        // kddnote Can this loop be disabled for RECTILINEAR BLOCKS?
                        // kddnote Mehmet says it is probably needed anyway.
                        distance =ABS(cutCoordinates_tmp[kk] - cut);
                        if(distance < _EPSILON){
                            myPartWeights[2 * kk + 1] += w;

                            myLeftClosest[kk] = coord;
                            myRightClosest[kk] = coord;
                            kk++;
                        }
                        else{
                            if(coord - myLeftClosest[kk] > _EPSILON){
                                myLeftClosest[kk] = coord;
                            }
                            break;
                        }
                    }

                    kk = j - 1;
                    while(kk >= 0){
                        distance =ABS(cutCoordinates_tmp[kk] - cut);
                        if(distance < _EPSILON){
                            myPartWeights[2 * kk + 1] += w;

                            //try to write the partId as the leftmost cut.
                            partIds[i] = kk * 2 + 1;
                            myLeftClosest[kk] = coord;
                            myRightClosest[kk] = coord;
                            kk--;
                        }
                        else{
                            if(myRightClosest[kk] - coord > _EPSILON){
                                myRightClosest[kk] = coord;
                            }
                            break;
                        }
                    }

                    isInserted = true;
                    break;
                }
                else {
                    if (distance < 0) {
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
                            pq_scalar_t distance_ = coord -
                                                     cutCoordinates_tmp[j + 1];
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
                    //cout << "\ti:" << i << " assigned:" << 2 * lastPart + 2
                    //  << endl;
                    myPartWeights[2 * lastPart + 2] += w;
                    partIds[i] = 2 * lastPart + 2;
                    //pq_scalar_t distance = coord - cutCoordinates_tmp[lastPart];
                    if(myRightClosest[lastPart] - coord > _EPSILON){
                        myRightClosest[lastPart] = coord;
                    }
                    if(lastPart+1 < noCuts){
                        //pq_scalar_t distance_ = cutCoordinates_tmp[lastPart + 1] - coord;
                        if(coord - myLeftClosest[lastPart + 1] > _EPSILON){
                            myLeftClosest[lastPart + 1] = coord;
                        }
                    }

                }
                else if(onLeft){
                    //cout << "\ti:" << i << " assigned:" <<
                    //  2 * lastPart << endl;
                    myPartWeights[2 * lastPart] += w;
                    partIds[i] = 2 * lastPart;
                    //pq_scalar_t distance = cutCoordinates_tmp[lastPart ] - coord;
                    if(coord - myLeftClosest[lastPart] > _EPSILON){
                        myLeftClosest[lastPart] = coord;
                    }

                    if(lastPart-1 >= 0){
                        //pq_scalar_t distance_ = coord - cutCoordinates_tmp[lastPart - 1];
                        if(myRightClosest[lastPart -1] - coord > _EPSILON){
                            myRightClosest[lastPart -1] = coord;
                        }
                    }
                }
            }
        }
    }
    else {

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for (pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
            int i = partitionedPointPermutations[ii];

            pq_scalar_t w = pqJagged_uniformWeights ? 1 : pqJagged_weights[i];

            pq_scalar_t coord = pqJagged_coordinates[i];

            partId_t j = partIds[i] / 2;

            if(j >= noCuts){
                j = noCuts - 1;
            }
            pq_scalar_t cut = cutCoordinates_tmp[j];
            pq_scalar_t distance = coord - cut;
            pq_scalar_t absdistance = ABS(distance);

            //comp++;
            if(absdistance < _EPSILON){
                myPartWeights[j * 2 + 1] += w;

                myLeftClosest[j] = 0;
                myRightClosest[j] = 0;
                partIds[i] = j * 2 + 1;

                //bas
                partId_t kk = j + 1;
                while(kk < noCuts){ // Needed when cuts shared the same position
                    // kddnote Can this loop be disabled for RECTILINEAR BLOCKS?
                    // kddnote Mehmet says it is probably needed anyway.
                    distance =ABS(cutCoordinates_tmp[kk] - cut);
                    if(distance < _EPSILON){
                        myPartWeights[2 * kk + 1] += w;
                        //cout << "2to part:" << 2*kk+1 << " coord:" << coord
                        //<< endl;
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
                        //cout << "3to part:" << 2*kk+1 << " coord:" << coord
                        //<< endl;

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
                    //cout << "4to part:" << 2*j+1 <<" j:" << j <<  " coord:"
                    //<< coord << endl;
                    //cout << "cut:" << cutCoordinates_tmp[j] << " dis:" <<i
                    //distance << endl;
                    partIds[i] = j * 2 + 1;

                    partId_t kk = j + 1;
                    while(kk < noCuts){//Needed when cuts share the same position
                        // kddnote Can this loop be disabled for RECTILINEAR BLOCKS?
                        // kddnote Mehmet says it is probably needed anyway.
                        distance =ABS(cutCoordinates_tmp[kk] - cut);
                        //cout << "distance:" << distance << endl;
                        if(distance < _EPSILON){
                            myPartWeights[2 * kk + 1] += w;
                            //cout << "5to part:" << 2*kk+1 << " kk:" << kk <<
                            //" coord:" << coord << endl;
                            //cout << "cut:" << cutCoordinates_tmp[kk] <<
                            //" dis:" << distance << endl;
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
                            //cout << "6to part:" << 2*kk+1 << " coord:" <<
                            //coord << endl;
                            //cout << "cut:" << cutCoordinates_tmp[kk] <<i
                            //" dis:" << distance << endl;

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
                    while(kk < noCuts){//Needed when cuts share the same position
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
for (pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
int i = partitionedPointPermutations[ii];

pq_scalar_t w = pqJagged_uniformWeights? 1:pqJagged_weights[i];
//get a coordinate and compare it with cut lines from left to right.
for(partId_t j = 0; j < noCuts; ++j){

if(isDone[j]) continue;
pq_scalar_t distance = pqJagged_coordinates[i] - cutCoordinates_tmp[j];
pq_scalar_t absdistance = ABS(distance);
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

    // prefix sum computation.
    for (size_t i = 1; i < total_part_count; ++i){
        // check for cuts sharing the same position; all cuts sharing a position
        // have the same weight == total weight for all cuts sharing the position.
        // don't want to accumulate that total weight more than once.
        if(i % 2 == 0 && i > 1 && i < total_part_count - 1 &&
                ABS(cutCoordinates_tmp[i / 2] - cutCoordinates_tmp[i /2 - 1])
                < _EPSILON){
            myPartWeights[i] = myPartWeights[i-2];
            continue;
        }
        myPartWeights[i] += myPartWeights[i-1];
    }
}


/*! \brief Function that reduces the result of multiple threads for left and right closest points and part weights in a single mpi.
 *
 * \param pVector is the vector that holds the number of cut lines for each part.
 * \param vBegin holds the index of the first part (important when concurrent parts are used.)
 * \param concurrentPartCount is the number of parts whose cut lines will be calculated concurrently.
 * \param numThreads hold the number of threads available per processor.
 * \param isDone is the boolean array to determine if the correct position for a cut line is found.
 * \param leftClosestDistance is the 2 dimensional array holding the coordinates of the closest points to the cut lines from left for each thread.
 * \param rightClosestDistance is the 2 dimensional array holding the coordinates of the closest points to the cut lines from right for each thread.
 * \param partWeights is the array holding the weight of parts for each thread. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param localMinMaxTotal is the array holding the local minimum and maximum coordinate and local total weight of each part.
 * \param totalPartWeights_leftClosest_rightCloset is the output array of accumulation, where total part weights (2P - 1),
 * then left closest distances (P -1), then right closest distance (P -1) are stored.
 */
template <typename pq_scalar_t>
void accumulateThreadResults(
    const vector <partId_t> &pVector,
    partId_t vBegin,
    partId_t concurrentPartCount,
    int numThreads,
    bool *isDone,
    pq_scalar_t **leftClosestPoints,
    pq_scalar_t **rightClosestPoints,
    double **partWeights,
    pq_scalar_t *totalPartWeights_leftClosest_rightCloset
){

#ifdef HAVE_ZOLTAN2_OMP
//needs barrier here, as it requires all threads to finish pqJagged_1DPart_getPartWeights
#pragma omp barrier
#pragma omp single
#endif
        {
            size_t tlr_shift = 0;
            partId_t cut_shift = 0;
            for(partId_t i = 0; i < concurrentPartCount; ++i){

                partId_t partNo =  pVector[vBegin + i];
                partId_t noCuts = partNo - 1;
                size_t total_part_count = partNo + size_t (noCuts) ;

                for(partId_t ii = 0; ii < noCuts ; ++ii){
                    partId_t next = tlr_shift + ii;
                    partId_t nCut = cut_shift + ii;
                    if(isDone[nCut]) continue;
                    pq_scalar_t minl = leftClosestPoints[0][nCut],
                                minr = rightClosestPoints[0][nCut];

                    for (int j = 1; j < numThreads; ++j){
                        if (rightClosestPoints[j][nCut] < minr ){
                            minr = rightClosestPoints[j][nCut];
                        }
                        if (leftClosestPoints[j][nCut] > minl ){
                            minl = leftClosestPoints[j][nCut];
                        }
                    }
                    totalPartWeights_leftClosest_rightCloset[total_part_count +
                                next] = minl;
                    totalPartWeights_leftClosest_rightCloset[total_part_count +
                                noCuts + next] = minr;
                }
                tlr_shift += (total_part_count + 2 * noCuts);
                cut_shift += noCuts;
            }

            tlr_shift = 0;
            cut_shift = 0;
            size_t totalPartShift = 0;

            for(partId_t i = 0; i < concurrentPartCount; ++i){

                partId_t partNo =  pVector[vBegin + i];
                partId_t noCuts = partNo - 1;
                size_t total_part_count = partNo + size_t (noCuts) ;

                for(size_t j = 0; j < total_part_count; ++j){

                    partId_t cutInd = j / 2 + cut_shift;

                    if(j !=  total_part_count - 1 && isDone[cutInd]) continue;
                    double pwj = 0;
                    for (int k = 0; k < numThreads; ++k){
                        pwj += partWeights[k][totalPartShift + j];
                    }
                    //size_t jshift = j % total_part_count + i * (total_part_count + 2 * noCuts);
                    totalPartWeights_leftClosest_rightCloset[tlr_shift + j] = pwj;
                }
                cut_shift += noCuts;
                tlr_shift += total_part_count + 2 * noCuts;
                totalPartShift += total_part_count;
            }
        }
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
 * \param pVector is the vector that holds how many parts each part will be divided into.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t>
void pqJagged_1D_Partition(
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,

    pq_lno_t *partitionedPointPermutations,
    pq_scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    pq_scalar_t *pqJagged_weights,

    pq_scalar_t *targetPartWeightRatios,   // the weight ratios at left side of the cuts. last is 1.
    pq_scalar_t *globalMinMaxTotal,
    pq_scalar_t *localMinMaxTotal,

    int numThreads,
    //pq_scalar_t maxScalar,
    //pq_scalar_t minScalar,
    pq_scalar_t imbalanceTolerance,
    partId_t currentPartBeginIndex,
    partId_t concurrentPartCount,
    pq_lno_t *inTotalCounts,

    pq_scalar_t *cutCoordinates,
    pq_scalar_t *cutCoordinatesWork, 	// work array to manipulate coordinate of cutlines in different iterations.
    pq_scalar_t **leftClosestDistance,
    pq_scalar_t **rightClosestDistance,
    pq_scalar_t *cutUpperBounds,  //to determine the next cut line with binary search
    pq_scalar_t *cutLowerBounds,  //to determine the next cut line with binary search
    pq_scalar_t *cutUpperWeight,   //to determine the next cut line with binary search
    pq_scalar_t *cutLowerWeight,  //to determine the next cut line with binary search
    bool *isDone,
    double **partWeights,
    pq_scalar_t *local_totalPartWeights_leftClosest_rightCloset,
    pq_scalar_t *global_totalPartWeights_leftClosest_rightCloset,
    bool allowNonRectelinearPart,
    float *nonRectelinearPartRatios,
    pq_scalar_t *localCutWeights,
    pq_scalar_t *globalCutWeights,

    partId_t allDone,
    partId_t *myNonDoneCounts,
    bool useBinarySearch,

    partId_t * partIds,
    vector <partId_t> &pVector
){

    partId_t recteLinearCutCount = 0;
    pq_scalar_t *cutCoordinates_tmp = cutCoordinates;

#ifdef mpi_communication
    MPI_Op myop;
    MPI_Op_create(sumMinMin, 0, &myop);   /* step 3 */
#endif

    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon();

    Teuchos::PQJaggedCombinedReductionOp<partId_t, pq_scalar_t>
                 *reductionOp = NULL;
    reductionOp = new Teuchos::PQJaggedCombinedReductionOp
                     <partId_t, pq_scalar_t>(&pVector , currentPartBeginIndex ,
                      concurrentPartCount);

    size_t totalReductionSize = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(allDone,  recteLinearCutCount)
#endif
    {
        int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
        me = omp_get_thread_num();
#endif
        double *myPartWeights = partWeights[me];
        pq_scalar_t *myLeftClosest = leftClosestDistance[me];
        pq_scalar_t *myRightClosest = rightClosestDistance[me];

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            {
                //initialize the lower and upper bounds of the cuts.
                partId_t next = 0;
                for(partId_t i = 0; i < concurrentPartCount; ++i){

                    partId_t partNo =  pVector[currentPartBeginIndex + i];
                    partId_t noCuts = partNo - 1;
                    totalReductionSize += (4 * noCuts + 1);

                    for(partId_t ii = 0; ii < noCuts; ++ii){
                        isDone[next] = false;
                        cutLowerBounds[next] = globalMinMaxTotal[i]; //min coordinate
                        cutUpperBounds[next] = globalMinMaxTotal[i + concurrentPartCount]; //max coordinate

                        cutUpperWeight[next] = globalMinMaxTotal[i + 2 * concurrentPartCount]; //total weight
                        cutLowerWeight[next] = 0;

                        if(allowNonRectelinearPart){
                            nonRectelinearPartRatios[next] = 0;
                        }
                        ++next;
                    }
                }
            }

        //no need to have barrier here.
        //pragma omp single have implicit barrier.

        int iteration = 0;
        while (allDone != 0){
            iteration += 1;
            partId_t cutShifts = 0;
            size_t totalPartShift = 0;

            for (partId_t kk = 0; kk < concurrentPartCount; ++kk){
                partId_t partNo =  -1;
                partNo =  pVector[currentPartBeginIndex + kk];
                //cout << "p:" << partNo << endl;
                partId_t noCuts = partNo - 1;
                size_t total_part_count = partNo + size_t (noCuts) ;
                if (myNonDoneCounts[kk] > 0){

                    //although isDone shared, currentDone is private and same for all.
                    bool *currentDone = isDone + cutShifts;
                    double *myCurrentPartWeights = myPartWeights + totalPartShift;
                    pq_scalar_t *myCurrentLeftClosest = myLeftClosest + cutShifts;
                    pq_scalar_t *myCurrentRightClosest = myRightClosest + cutShifts;

                    partId_t current = currentPartBeginIndex + kk;
                    pq_lno_t coordinateBegin = current ==0 ? 0: inTotalCounts[current -1];
                    pq_lno_t coordinateEnd = inTotalCounts[current];
                    pq_scalar_t *cutCoordinates_ = cutCoordinates_tmp + cutShifts;

                    pq_scalar_t minCoordinate = globalMinMaxTotal[kk];
                    pq_scalar_t maxCoordinate = globalMinMaxTotal[kk + concurrentPartCount];

                    // compute part weights using existing cuts
                    pqJagged_1DPart_getPartWeights<pq_scalar_t, pq_lno_t>(
                        total_part_count, noCuts,
                        maxCoordinate,//globalMinMaxTotal[kk + concurrentPartCount],//maxScalar,
                        minCoordinate,//globalMinMaxTotal[kk]//minScalar,
                        _EPSILON, numThreads, coordinateBegin, coordinateEnd,
                        partitionedPointPermutations, pqJagged_coordinates,
                        pqJagged_uniformWeights, pqJagged_weights,
                        cutCoordinates_, currentDone, myCurrentPartWeights,
                        myCurrentLeftClosest, myCurrentRightClosest,
                        useBinarySearch, partIds);

                }

                cutShifts += noCuts;
                totalPartShift += total_part_count;
            }

            //sum up the results of threads
            accumulateThreadResults<pq_scalar_t>(
                pVector, currentPartBeginIndex, concurrentPartCount, numThreads,
                isDone, leftClosestDistance, rightClosestDistance, partWeights,
                local_totalPartWeights_leftClosest_rightCloset
            );
            /*
#pragma omp single
            {
                int cutShifts = 0;
                int partNo = pVector[currentPartBeginIndex];
                cout << "iteration:" << iteration << endl;

                pq_scalar_t *cutCoordinates_ = cutCoordinates_tmp + cutShifts;
                pq_scalar_t *myCurrentLeftClosest = myLeftClosest + cutShifts;
                pq_scalar_t *myCurrentRightClosest = myRightClosest + cutShifts;
                bool *currentDone = isDone + cutShifts;

                for (int i = 0; i < partNo - 1; ++i){
                    if (currentDone[i]){
                        cout << "\tcut:" << i << " c:" << cutCoordinates_[i] <<
                                " lc:" << local_totalPartWeights_leftClosest_rightCloset[2 * partNo - 1 + i] <<
                                " rc:" << local_totalPartWeights_leftClosest_rightCloset[partNo - 1 + 2 * partNo - 1 + i]<<
                                " w:" << local_totalPartWeights_leftClosest_rightCloset[2 * i]<<
                                " done" <<endl;
                    }else {

                        cout << "\tcut:" << i << " c:" << cutCoordinates_[i] <<
                                " lc:" << local_totalPartWeights_leftClosest_rightCloset[2 * partNo - 1 + i] <<
                                " rc:" << local_totalPartWeights_leftClosest_rightCloset[partNo - 1 + 2 * partNo - 1 + i]<<
                                " w:" << local_totalPartWeights_leftClosest_rightCloset[2 * i]<<
                                " not done" << endl;
                    }
                }
            }
            */
            //now sum up the results of mpi processors.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            {
                if(comm->getSize() > 1){
                    try{
#ifdef mpi_communication

                        MPI_Allreduce(local_totalPartWeights_leftClosest_rightCloset,
                            global_totalPartWeights_leftClosest_rightCloset,
                            (total_part_count + 2 * noCuts) *
                             concurrentPartCount, MPI_FLOAT, myop,
                             MPI_COMM_WORLD);
#endif
#ifndef mpi_communication

                        reduceAll<int, pq_scalar_t>( *comm, *reductionOp,
                            totalReductionSize,
                            local_totalPartWeights_leftClosest_rightCloset,
                            global_totalPartWeights_leftClosest_rightCloset);
#endif
                    }
                    Z2_THROW_OUTSIDE_ERROR(*env)
                }
                else {
                        memcpy(global_totalPartWeights_leftClosest_rightCloset,
                            local_totalPartWeights_leftClosest_rightCloset,
                            totalReductionSize * sizeof(pq_scalar_t));
                }
            }
            partId_t cutShift = 0;
            size_t tlrShift = 0;
            for (partId_t kk = 0; kk < concurrentPartCount; ++kk){
                partId_t partNo =  -1;
                    partNo = pVector[currentPartBeginIndex + kk];

                partId_t noCuts = partNo - 1;
                size_t total_part_count = partNo + size_t (noCuts) ;

                if (myNonDoneCounts[kk] == 0) {
                    cutShift += noCuts;
                    tlrShift += (total_part_count + 2 * noCuts);
                    continue;
                }

                pq_scalar_t *localPartWeights = local_totalPartWeights_leftClosest_rightCloset  + tlrShift ;
                pq_scalar_t *gtlr = global_totalPartWeights_leftClosest_rightCloset + tlrShift;
                pq_scalar_t *glc = gtlr + total_part_count; //left closest points
                pq_scalar_t *grc = gtlr + total_part_count + noCuts; //right closest points
                pq_scalar_t *globalPartWeights = gtlr;
                bool *currentDone = isDone + cutShift;

                pq_scalar_t *currentTargetPartWeightRatios = targetPartWeightRatios + cutShift + kk;
                float *currentnonRectelinearPartRatios = nonRectelinearPartRatios + cutShift;

                pq_scalar_t minCoordinate = globalMinMaxTotal[kk];
                pq_scalar_t maxCoordinate = globalMinMaxTotal[kk + concurrentPartCount];
                pq_scalar_t globalTotalWeight = globalMinMaxTotal[kk + concurrentPartCount * 2];
                pq_scalar_t *currentcutLowerWeight = cutLowerWeight + cutShift;
                pq_scalar_t *currentcutUpperWeight = cutUpperWeight + cutShift;
                pq_scalar_t *currentcutUpperBounds = cutUpperBounds + cutShift;
                pq_scalar_t *currentcutLowerBounds = cutLowerBounds + cutShift;

                partId_t prevDoneCount = myNonDoneCounts[kk];

                // Now compute the new cut coordinates.
                getNewCoordinates<pq_scalar_t>(
                    env, comm, total_part_count, noCuts, maxCoordinate,
                    minCoordinate, globalTotalWeight, imbalanceTolerance,
                    //maxCoordinate,//globalMinMaxTotal[kk + concurrentPartCount],//maxScalar,
                    //minCoordinate,//globalMinMaxTotal[kk],//minScalar,
                    globalPartWeights, localPartWeights,
                    currentTargetPartWeightRatios, currentDone, 
                    cutCoordinates_tmp + cutShift,
                    currentcutUpperBounds, currentcutLowerBounds, glc, grc,
                    currentcutLowerWeight, currentcutUpperWeight,
                    cutCoordinatesWork +cutShift, //new cut coordinates 
                    allowNonRectelinearPart, currentnonRectelinearPartRatios,
                    &recteLinearCutCount, localCutWeights, globalCutWeights,
                    myNonDoneCounts[kk]);

                cutShift += noCuts;
                tlrShift += (total_part_count + 2 * noCuts);
                partId_t reduction = prevDoneCount - myNonDoneCounts[kk];
#ifdef HAVE_ZOLTAN2_OMP
//#pragma omp barrier
//ne need for the barrier as getNewCoordinates function have a implicit barrier at the end.
#pragma omp single
#endif
                {
                    allDone -= reduction;
                }

            }
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#pragma omp single
#endif
            {
                pq_scalar_t *t = cutCoordinates_tmp;
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
#pragma omp single
#endif
            {
                partId_t next = 0;
                for(partId_t i = 0; i < concurrentPartCount; ++i){
                    partId_t partNo = -1;
                        partNo = pVector[currentPartBeginIndex + i];
                    //cout << "kk:" << " partition - 5:" << partNo << endl;
                    partId_t noCuts = partNo - 1;
                    for(partId_t ii = 0; ii < noCuts; ++ii){
                        cutCoordinates[next + ii] = cutCoordinates_tmp[next + ii];
                    }
                    next += noCuts;
                }
            }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            {
                cutCoordinatesWork = cutCoordinates_tmp;
            }
        }
    }
    delete reductionOp;
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


template <typename IT, typename CT, typename WT>
class uMultiSortItem
{
public:
    IT index;
    CT count;
    //unsigned int val;
    WT *val;
    WT _EPSILON;

    uMultiSortItem(){
        this->index = 0;
        this->count = 0;
        this->val = NULL;
        this->_EPSILON = numeric_limits<WT>::epsilon();
    }


    uMultiSortItem(IT index_ ,CT count_, WT *vals_){
        this->index = index_;
        this->count = count_;
        this->val = vals_;
        this->_EPSILON = numeric_limits<WT>::epsilon();
    }

    uMultiSortItem( const uMultiSortItem<IT,CT,WT>& other ){
        this->index = other.index;
        this->count = other.count;
        this->val = other.val;
        this->_EPSILON = other._EPSILON;
    }

    ~uMultiSortItem(){
        //freeArray<WT>(this->val);
    }

    void set(IT index_ ,CT count_, WT *vals_){
        this->index = index_;
        this->count = count_;
        this->val = vals_;
    }


    uMultiSortItem<IT,CT,WT> operator=(const uMultiSortItem<IT,CT,WT>& other){
        this->index = other.index;
        this->count = other.count;
        this->val = other.val;
        return *(this);
    }

    bool operator<(const uMultiSortItem<IT,CT,WT>& other) const{
        assert (this->count == other.count);
        for(CT i = 0; i < this->count; ++i){
            //if the values are equal go to next one.
            if (ABS(this->val[i] - other.val[i]) < this->_EPSILON){
                continue;
            }
            //if next value is smaller return true;
            if(this->val[i] < other.val[i]){
                return true;
            }
            //if next value is bigger return false;
            else {
                return false;
            }
        }
        //if they are totally equal.
        return this->index < other.index;
    }
    bool operator>(const uMultiSortItem<IT,CT,WT>& other) const{
        assert (this->count == other.count);
        for(CT i = 0; i < this->count; ++i){
            //if the values are equal go to next one.
            if (ABS(this->val[i] - other.val[i]) < this->_EPSILON){
                continue;
            }
            //if next value is bigger return true;
            if(this->val[i] > other.val[i]){
                return true;
            }
            //if next value is smaller return false;
            else //(this->val[i] > other.val[i])
            {
                return false;
            }
        }
        //if they are totally equal.
        return this->index > other.index;
    }
};// uSortItem;


template <class IT, class WT>
struct uSortItem
{
    IT id;
    //unsigned int val;
    WT val;
};// uSortItem;


template <class IT, class WT>
void uqsort(IT n, uSortItem<IT, WT> * arr)
{
#define SWAP(a,b,temp) temp=(a);(a)=(b);(b)=temp;
    int NSTACK = 50;
    int M = 7;
    IT         i, ir=n, j, k, l=1;
    IT         jstack=0, istack[50];
    WT aval;
    uSortItem<IT,WT>    a, temp;

    --arr;
    for (;;)
    {
        if (ir-l < M)
        {
            for (j=l+1;j<=ir;j++)
            {
                a=arr[j];
                aval = a.val;
                for (i=j-1;i>=1;i--)
                {
                    if (arr[i].val <= aval)
                        break;
                    arr[i+1] = arr[i];
                }
                arr[i+1]=a;
            }
            if (jstack == 0)
                break;
            ir=istack[jstack--];
            l=istack[jstack--];
        }
        else
        {
            k=(l+ir) >> 1;
            SWAP(arr[k],arr[l+1], temp)
            if (arr[l+1].val > arr[ir].val)
            {
                SWAP(arr[l+1],arr[ir],temp)
            }
            if (arr[l].val > arr[ir].val)
            {
                SWAP(arr[l],arr[ir],temp)
            }
            if (arr[l+1].val > arr[l].val)
            {
                SWAP(arr[l+1],arr[l],temp)
            }
            i=l+1;
            j=ir;
            a=arr[l];
            aval = a.val;
            for (;;)
            {
                do i++; while (arr[i].val < aval);
                do j--; while (arr[j].val > aval);
                if (j < i) break;
                SWAP(arr[i],arr[j],temp);
            }
            arr[l]=arr[j];
            arr[j]=a;
            jstack += 2;
            if (jstack > NSTACK){
                cout << "uqsort: NSTACK too small in sort." << endl;
                exit(1);
            }
            if (ir-i+1 >= j-l)
            {
                istack[jstack]=ir;
                istack[jstack-1]=i;
                ir=j-1;
            }
            else
            {
                istack[jstack]=j-1;
                istack[jstack-1]=l;
                l=i;
            }
        }
    }
}

#ifdef enable_migration2
template <typename pq_gno_t, typename pq_lno_t,typename pq_scalar_t, typename pq_node_t>
RCP<const Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> > 
create_initial_multi_vector(
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,
    pq_gno_t numGlobalPoints,
    pq_lno_t numLocalPoints,
    int coord_dim,
    pq_scalar_t **coords,
    int weight_dim,
    pq_scalar_t **weight,
    int pqJagged_multiVectorDim
){

    typedef Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> tMVector_t;
#ifdef memory_debug
    sleep(1); env->memory("initial multivector before map");
#endif
    RCP<Tpetra::Map<pq_lno_t, pq_gno_t, pq_node_t> > mp = rcp(
            new Tpetra::Map<pq_lno_t, pq_gno_t, pq_node_t> (numGlobalPoints, numLocalPoints, 0, comm));


#ifdef memory_debug
    if(comm->getRank() == 0) cout << "nn:"<< numLocalPoints<<"created map:" << mp.getRawPtr() << endl;
    sleep(1); env->memory("initial multivector after map");
#endif
    Teuchos::Array<Teuchos::ArrayView<const pq_scalar_t> > coordView(pqJagged_multiVectorDim);
    for (int i=0; i < coord_dim; i++){
        //cout << "setting i:" << i << endl;
        if(numLocalPoints > 0){
            Teuchos::ArrayView<const pq_scalar_t> a(coords[i], numLocalPoints);
            coordView[i] = a;
        } else{
            Teuchos::ArrayView<const pq_scalar_t> a;
            coordView[i] = a;
        }
#ifdef memory_debug
        sleep(1); env->memory("initial multivector after coordinates -" + toString<int>(i));
#endif
    }
#ifdef memory_debug
    sleep(1); env->memory("initial multivector after coordinates");
#endif
    for (int i=0;  i < weight_dim; ++i){
        int j = i + coord_dim;
        //if (j >= pqJagged_multiVectorDim - 1) break;
        //cout << "setting j:" << j << endl;
        if(numLocalPoints > 0){
            Teuchos::ArrayView<const pq_scalar_t> a(weight[i], numLocalPoints);
            coordView[j] = a;
        } else{
            Teuchos::ArrayView<const pq_scalar_t> a;
            coordView[j] = a;
        }
    }

#ifdef memory_debug
    sleep(1); env->memory("initial multivector after weights");
#endif
#ifdef migrate_gid
    if(numLocalPoints > 0){
        //cout << "writing :" << weight_dim + coord_dim << endl;
        Teuchos::ArrayView<const pq_scalar_t> a(mappedGnos, numLocalPoints);
        coordView[pqJagged_multiVectorDim - 1] = a;
    } else{
        Teuchos::ArrayView<const pq_scalar_t> a;
        coordView[pqJagged_multiVectorDim - 1] = a;
    }
#endif
    RCP< Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> >tmVector =
            RCP< Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> >(
                    new Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t>( mp, coordView.view(0, pqJagged_multiVectorDim), pqJagged_multiVectorDim));

#ifdef memory_debug
    sleep(1); env->memory("initial multivector after create");
#endif
    RCP<const tMVector_t> coordsConst = Teuchos::rcp_const_cast<const tMVector_t>(tmVector);

#ifdef memory_debug
    sleep(1); env->memory("initial multivector after cast");
#endif
    return coordsConst;
}


//fills up the p_gno_np_global_num_coord_each_part_actual array and
//returns the allocation size of arrays.
template <typename pq_gno_t,
          typename pq_lno_t,
          typename partId_t>
void getProcessorCoordinatePartCounts(
    RCP<Comm<int> > &pcomm,
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,
    int migration_proc_assignment_type,
    partId_t nprocs,
    partId_t myRank,
    partId_t num_parts,
    pq_lno_t *partBeginArray,
    pq_gno_t *&p_gno_np_global_num_coord_each_part_actual){

    //if nprocs is less than numparts we still need to allocate more memory.
    if (nprocs <= num_parts) migration_proc_assignment_type = 1;

    //initially allocation_size is num_parts
    size_t allocation_size = num_parts;
    if (migration_proc_assignment_type == 1){
        //increase it by nprocs + 1 so that everyone will
        //know how many coordinate each processors has in each part.
        allocation_size = num_parts * (nprocs + 1);
    }


    //this will be output
    p_gno_np_global_num_coord_each_part_actual  = allocMemory<pq_gno_t>(allocation_size);

    //allocate memory for the local num coordinates in each part.
    pq_gno_t *p_gno_np_local_num_coord_each_part_actual = allocMemory<pq_gno_t>(allocation_size);


    pq_gno_t *p_gno_np_local_num_coord_each_part = p_gno_np_local_num_coord_each_part_actual;
    pq_gno_t *p_gno_np_local_num_coord_each_part_mypart = p_gno_np_local_num_coord_each_part_actual;

    //each processor will write a certain index of the local arrays.
    //then this arrays will be used in reduceAll function.
    if (migration_proc_assignment_type == 1){
        partId_t shift_amount = nprocs * num_parts;
        partId_t my_part_shift = myRank * num_parts;
        p_gno_np_local_num_coord_each_part += shift_amount;
        p_gno_np_local_num_coord_each_part_mypart += my_part_shift;
    }

    memset(p_gno_np_local_num_coord_each_part_actual, 0, sizeof(pq_gno_t)*allocation_size);

    //write the number of coordinates in each part.
    for (partId_t i = 0; i < num_parts; ++i){
        pq_lno_t pBegin = 0;
        if (i > 0){
            pBegin = partBeginArray[i - 1];
        }
        pq_lno_t pEnd = partBeginArray[i];
        p_gno_np_local_num_coord_each_part[i] = pEnd - pBegin;
    }

    //copy the local num parts to the last portion of array,
    //so that this portion will represent the global num points in each part after the reduction.
    if (migration_proc_assignment_type == 1){
        memcpy (p_gno_np_local_num_coord_each_part_mypart,
                p_gno_np_local_num_coord_each_part,
                sizeof(pq_gno_t) * (num_parts) );
    }

    //reduceAll operation.
    //when allocation_size = num_parts * (nprocs + 1),
    //the portion that belongs to a processor with index p
    //will start from myRank * num_parts.
    //the global number of points will be held at the index
    //nprocs * num_parts size
    try{
        reduceAll<int, pq_gno_t>( *comm, Teuchos::REDUCE_SUM, allocation_size,
                p_gno_np_local_num_coord_each_part_actual,
                p_gno_np_global_num_coord_each_part_actual);
    }
    Z2_THROW_OUTSIDE_ERROR(*env)
    freeArray<pq_gno_t>(p_gno_np_local_num_coord_each_part_actual);
    //free local num coordinates array, as it is no use after this point.

}


//checks if should do migration or not.
//It returns true to point that migration should be done when
// -futureReduceAlls are higher than a predetermined value
// -numCoordinates that left for the last dimension partitioning is less than a predetermined value
// -the imbalance of the processors on the parts are higher than given threshold.
//returns false otherwise, and migration is not performed.
template <typename pq_gno_t, typename pq_lno_t, typename partId_t>
bool checkMigration(
    RCP<Comm<int> > &pcomm,
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,
    int migration_check_option,
    size_t futureReduceAll,
    pq_lno_t numCoordinatesForLastDimPartitioning,
    float migration_imbalance_cut_off,
    int migration_proc_assignment_type,
    partId_t nprocs,
    partId_t myRank,
    partId_t num_parts,
    pq_gno_t *p_gno_np_global_num_coord_each_part_actual,
    pq_lno_t *partBeginArray){

    /*
    if (myRank == 0){
        cout << "futureReduceAll:" << futureReduceAll <<
                " FUTURE_REDUCEALL_CUTOFF:" << FUTURE_REDUCEALL_CUTOFF <<
                " numCoordinatesForLastDimPartitioning:" <<
                 numCoordinatesForLastDimPartitioning <<
                " MIN_WORK_LAST_DIM:" << MIN_WORK_LAST_DIM << endl;
    }
    */

    if (futureReduceAll > FUTURE_REDUCEALL_CUTOFF) return true;
    if (numCoordinatesForLastDimPartitioning < MIN_WORK_LAST_DIM) return true;

    if (migration_check_option == 0){
        double diff = 0, global_diff = 0;
        if (migration_proc_assignment_type == 0 )
        {

            for (partId_t i = 0; i < num_parts; ++i){
                double ideal_num = p_gno_np_global_num_coord_each_part_actual[i] /  double(nprocs);

                pq_lno_t pBegin = 0;
                if (i > 0){
                    pBegin = partBeginArray[i - 1];
                }
                pq_lno_t pEnd = partBeginArray[i];
                diff += ABS(ideal_num -
                        (pEnd - pBegin)) /  (ideal_num);
            }
            diff /= num_parts;
            reduceAll<int, double>( *comm, Teuchos::REDUCE_SUM, 1, &diff,
                    &global_diff);
        }
        else {
            size_t global_shift = nprocs * num_parts;
            for (partId_t ii = 0; ii < nprocs; ++ii){
                for (partId_t i = 0; i < num_parts; ++i){
                    double ideal_num = p_gno_np_global_num_coord_each_part_actual[global_shift + i] /  double(nprocs);
                    global_diff += ABS(ideal_num -
                            p_gno_np_global_num_coord_each_part_actual[ii * num_parts + i]) /  (ideal_num);
                    /*
                    if (myRank == 0&& nprocs == 96){
                        cout << "i:" << i << " ii:" << ii << " idealNum:" << ideal_num << " procHas:" << p_gno_np_global_num_coord_each_part_actual[ii * num_parts + i] << endl;
                    }
                    */
                }
            }
            global_diff /= num_parts;

            /*
            if (myRank == 0&& nprocs == 96){
                cout << "Global Diff:" << global_diff << " nprocs:" << nprocs << endl;
            }
            */

        }
        global_diff /= nprocs;
        if (myRank == 0) {
            cout << "imbalance for next iteration:" << global_diff << endl;
        }

        if(global_diff <= migration_imbalance_cut_off){
            return false;
        }
        else {
            return true;
        }
    }
    else {
        return true;
    }
}


//creates the new multivector in which the partId's are also added.
//this function is only called when nprocs < num_parts, in this case
//a processor will have coordinates from multiple parts.
//the partId's are migrated as well to separate the points from different parts.
template <typename pq_gno_t, typename pq_lno_t,typename pq_scalar_t, typename pq_node_t>
RCP< Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> > createNewMultivector(
    RCP<Comm<int> > &comm,
    pq_gno_t numGlobalPoints,
    pq_lno_t numLocalPoints,
    int coord_dim,
    pq_scalar_t **coords,
    int weight_dim,
    pq_scalar_t **weight,
    partId_t num_parts,
    partId_t *partBeginArray,
    pq_lno_t *permutationArray,
    int &multiVectorDim,
    RCP<const Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> > old_mvector){


    typedef Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> tMVector_t;
    RCP<const Tpetra::Map<pq_lno_t, pq_gno_t, pq_node_t> > mp = old_mvector->getMap();

    Teuchos::Array<Teuchos::ArrayView<const pq_scalar_t> > coordView(multiVectorDim + 1);

    for (int i=0; i < coord_dim; i++){
        if(numLocalPoints > 0){
            Teuchos::ArrayView<const pq_scalar_t> a(coords[i], numLocalPoints);
            coordView[i] = a;
        } else{
            Teuchos::ArrayView<const pq_scalar_t> a;
            coordView[i] = a;
        }
    }

    for (int i=0; i < weight_dim; i++){
        int j = i + coord_dim;
        if(numLocalPoints > 0){
            Teuchos::ArrayView<const pq_scalar_t> a(weight[i], numLocalPoints);
            coordView[j] = a;
        } else{
            Teuchos::ArrayView<const pq_scalar_t> a;
            coordView[j] = a;
        }
    }
    pq_scalar_t *assigned_parts = allocMemory<pq_scalar_t>(numLocalPoints);


    for (partId_t i = 0; i < num_parts; ++i){
        pq_lno_t pBegin = 0;
        if (i > 0) pBegin = partBeginArray[i - 1];
        pq_lno_t pEnd = partBeginArray[i];

        //cout << "when setting me:" << comm->getRank() << " p:" << i << " c:" << pEnd - pBegin  <<  " assigned_parts:" << assigned_parts<< endl;
        pq_scalar_t p = pq_scalar_t (i);
        for (pq_lno_t j = pBegin; j < pEnd; ++j){
            pq_lno_t ind = permutationArray[j];
            assigned_parts[ind] = p;
        }
    }
    if(numLocalPoints > 0){
        Teuchos::ArrayView<const pq_scalar_t> a(assigned_parts, numLocalPoints);
        coordView[multiVectorDim] = a;
    } else{
        Teuchos::ArrayView<const pq_scalar_t> a;
        coordView[multiVectorDim] = a;
    }
    multiVectorDim += 1;

    RCP< Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> >tmVector = RCP< Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> >(
            new Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t>( mp, coordView.view(0, multiVectorDim), multiVectorDim));

    RCP<const tMVector_t> coordsConst = Teuchos::rcp_const_cast<const tMVector_t>(tmVector);

    return tmVector;
}


template <typename partId_t, typename pq_lno_t, typename pq_gno_t>
void fillContinousSendBuffer2(
    partId_t num_parts,
    pq_lno_t *partBegins,
    pq_lno_t *permutation,
    const pq_gno_t *gno_list,
    uSortItem<partId_t, partId_t> * part_assignment, //input sorted wrt processors
    pq_gno_t *sendBuf,
    partId_t &partBeginIndex,
    vector<partId_t> *newFuturePartitions,
    partId_t myRank){

    pq_lno_t nextInd = 0;
    partId_t partShift = partBeginIndex;
    partId_t previous_processor = -1;
    for(partId_t i = 0; i < num_parts; ++i){
        partId_t p = part_assignment[i].id;
        //assigned processors are sorted.
        //partId_t assigned_proc = part_assignment[i].val;

        pq_lno_t pBegin = 0;
        if (p > 0) pBegin = partBegins[p - 1];
        pq_lno_t pEnd = partBegins[p];

        partId_t assigned_proc = part_assignment[i].val;
        if (myRank == assigned_proc && previous_processor != assigned_proc){
            partBeginIndex =  partShift;
        }
        previous_processor = assigned_proc;
        partShift += (*newFuturePartitions)[p];

        for (pq_lno_t j=pBegin; j < pEnd; j++){
            pq_lno_t localInd = permutation[j];
            pq_gno_t to_send = gno_list[localInd];
            sendBuf[nextInd++] = to_send;
        }
        //cout << "sent to:" << part_assignment[i].val << " this much:" << pEnd - pBegin << " because of part:" << p << endl;
    }
}


template <typename partId_t, typename pq_lno_t>
void Z1fillSendBuffer2(
    partId_t num_parts,
    pq_lno_t *partBegins,
    pq_lno_t *permutation,
    uSortItem<partId_t, partId_t> * part_assignment, //input sorted wrt processors
    int *coordinate_destionations,
    partId_t &partBeginIndex,
    vector<partId_t> *newFuturePartitions,
    partId_t myRank){

    partId_t partShift = partBeginIndex;
    partId_t previous_processor = -1;
    for(partId_t i = 0; i < num_parts; ++i){
        partId_t p = part_assignment[i].id;
        //assigned processors are sorted.
        //partId_t assigned_proc = part_assignment[i].val;

        pq_lno_t pBegin = 0;
        if (p > 0) pBegin = partBegins[p - 1];
        pq_lno_t pEnd = partBegins[p];

        partId_t assigned_proc = part_assignment[i].val;
        if (myRank == assigned_proc && previous_processor != assigned_proc){
            partBeginIndex =  partShift;
        }
        previous_processor = assigned_proc;
        partShift += (*newFuturePartitions)[p];

        for (pq_lno_t j=pBegin; j < pEnd; j++){
            pq_lno_t localInd = permutation[j];
            coordinate_destionations[localInd] = assigned_proc;
/*
            cout << " i:" << localInd <<
                    " sending to:" << coordinate_destionations[localInd] <<
                    endl;
                    */

        }
    }
}

template <typename partId_t, typename pq_lno_t, typename pq_gno_t>
void fillContinousSendBuffer1(
    RCP<Comm<int> > &pcomm, //original communication.
    partId_t numParts,
    pq_lno_t *partBegins,
    pq_lno_t *permutation,
    const pq_gno_t *gno_list,
    partId_t nprocs,
    partId_t *part_assign_begins,
    partId_t *proc_chains,
    pq_lno_t *sendCount,
    pq_gno_t *sendBuf){
    //function will fill the sendBuf array in a consecutive way.

    //initially allocate array to store prefixSum of sendCounts,
    //so that we know where to write the send gnos for other processors.
    pq_lno_t *_sendCount_psum = allocMemory<pq_lno_t>(nprocs);
    pq_lno_t prefixsum = 0;
    for (int i = 0; i < nprocs; ++i ){
        _sendCount_psum[i] = prefixsum;
        prefixsum += sendCount[i];
    }

    /*
    for (int i = 0; i < nprocs; ++i ){
        cout << "me:" << pcomm->getRank() << " i:" << i << " _sendCount_psum:" << _sendCount_psum[i] << endl;
    }
    pcomm->barrier();
    cout << "me:" << pcomm->getRank() << " reached -1" << endl;
    */
    //now distribute all parts.
    for (partId_t p = 0; p < numParts; ++p){
        pq_lno_t pBegin = 0;
        if (p > 0) pBegin = partBegins[p - 1];
        pq_lno_t pEnd = partBegins[p];

        //get the first part that current processor will send its part-p.
        partId_t proc_to_sent = part_assign_begins[p];
        //initialize how many point I sent to this processor.
        pq_lno_t totalSend = 0;
        for (pq_lno_t j=pBegin; j < pEnd; j++){
            pq_lno_t localInd = permutation[j];
            pq_gno_t to_send = gno_list[localInd];
            //if the sendCount to this processor
            //reached its limit then start sending to the next processor.
            while (totalSend >= sendCount[proc_to_sent]){
                totalSend = 0;
                //assign new processor to part_assign_begin[p]
                part_assign_begins[p] = proc_chains[proc_to_sent];
                //remove the previous processor
                proc_chains[proc_to_sent] = -1;
                //choose the next processor as the next one to send.
                proc_to_sent = part_assign_begins[p];
            }
            //write the gno index to corresponding position in sendBuf.
            sendBuf[_sendCount_psum[proc_to_sent] + totalSend++] = to_send;
        }

    }
    freeArray<pq_lno_t>(_sendCount_psum);
}


template <typename partId_t, typename pq_lno_t>
void Z1fillSendBuffer1(
    RCP<Comm<int> > &pcomm, //original communication.
    partId_t numParts,
    pq_lno_t *partBegins,
    pq_lno_t *permutation,
    partId_t nprocs,
    partId_t *part_assign_begins,
    partId_t *proc_chains,
    pq_lno_t *sendCount,
    int *coordinate_destionations){

    for (partId_t p = 0; p < numParts; ++p){
        pq_lno_t pBegin = 0;
        if (p > 0) pBegin = partBegins[p - 1];
        pq_lno_t pEnd = partBegins[p];

        //get the first part that current processor will send its part-p.
        partId_t proc_to_sent = part_assign_begins[p];
        //initialize how many point I sent to this processor.
        pq_lno_t totalSend = 0;
        for (pq_lno_t j=pBegin; j < pEnd; j++){
            pq_lno_t localInd = permutation[j];
            while (totalSend >= sendCount[proc_to_sent]){
                totalSend = 0;
                //assign new processor to part_assign_begin[p]
                part_assign_begins[p] = proc_chains[proc_to_sent];
                //remove the previous processor
                proc_chains[proc_to_sent] = -1;
                //choose the next processor as the next one to send.
                proc_to_sent = part_assign_begins[p];
            }
            //write the gno index to corresponding position in sendBuf.
            coordinate_destionations[localInd] = proc_to_sent;
            ++totalSend;
        }

    }
}


template <typename partId_t, typename pq_lno_t, typename pq_gno_t>
void procAssignment2(
    int assignment_type, //either assign to minimize migration, or assign to increase locality.
    pq_gno_t * p_gno_np_global_num_coord_each_part_actual,
    pq_gno_t nGlobalObj,
    //pq_lno_t nLocal,
    partId_t num_parts,
    partId_t nprocs,
    partId_t myRank,
    pq_lno_t *partBegins, //holds the beginning of each part.
    pq_lno_t *permutation, //the permutation array ordered wrt partBegins array.
    const pq_gno_t *gnoList, //gno array

    pq_gno_t *sendBuf, //output: sized nLocal, the buffer is filled by the function with gnos.
    pq_lno_t *sendCount, //output: sized nprocs, show the number of send point counts to each proc.

    //TODO futurePartIndex might need to be changed.sendBuf,
    vector<partId_t> *newFuturePartitions,//input how many more partitions the part will be partitioned into.
    partId_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
    vector<partId_t> &out_part_index, //output: the part indices which the processor is assigned to.
    partId_t &partIndexBegin, //output: how much the part number should be shifted when setting the solution
    int doMigrationType,
    int *coordinate_destionations){

    out_num_part = 0;

    pq_gno_t *p_gno_np_global_num_coord_each_part = p_gno_np_global_num_coord_each_part_actual + nprocs * num_parts;
    //pq_gno_t *p_gno_np_global_num_coord_each_part_mypart = p_gno_np_global_num_coord_each_part_actual + myRank * num_parts;

    out_part_index.clear();

    //to sort the parts that is assigned to the processors.
    //id is the part number, sort value is the assigned processor id.
    uSortItem<partId_t, partId_t> * part_assignment  = allocMemory <uSortItem<partId_t, partId_t> >(num_parts);
    uSortItem<partId_t, pq_gno_t> * proc_load_sort = allocMemory <uSortItem<partId_t, pq_gno_t> >(nprocs);


    //calculate the optimal number of coordinates
    //that should be assigned to each processor.
    pq_lno_t work_each = nGlobalObj / (float (nprocs)) + 0.5f;
    //to hold the left space as the number of coordinates to the optimal number in each proc.
    pq_lno_t *space_in_each_processor = allocMemory <pq_lno_t>(nprocs);
    //initialize left space in each.
    for (partId_t i = 0; i < nprocs; ++i){
        space_in_each_processor[i] = work_each;
    }


    //to sort the parts with decreasing order of their coordiantes.
    //id are the part numbers, sort value is the number of points in each.
    uSortItem<partId_t, pq_gno_t> * part_loads  = allocMemory <uSortItem<partId_t, pq_gno_t> >(num_parts);

    //initially we will sort the parts according to the number of coordinates they have.
    //so that we will start assigning with the part that has the most number of coordinates.
    for (partId_t i = 0; i < num_parts; ++i){
        part_loads[i].id = i;
        part_loads[i].val = p_gno_np_global_num_coord_each_part[i];
    }
    //sort parts with increasing order of loads.
    uqsort<partId_t, pq_gno_t>(num_parts, part_loads);


    //assigning parts to the processors
    //traverse the part win decreasing order of load.
    //first assign the heaviest part.
    for (partId_t j = 0; j < num_parts; ++j){
        //sorted with increasing order, traverse inverse.
        partId_t i = part_loads[num_parts - 1 - j].id;
        //load of the part
        pq_gno_t load = p_gno_np_global_num_coord_each_part[i];

        //assigned processors
        partId_t assigned_proc = -1;
        //if not fit best processor.
        partId_t best_proc_to_assign = 0;


        //sort processors with increasing number of points in this part.
        for (partId_t ii = 0; ii < nprocs; ++ii){
            proc_load_sort[ii].id = ii;
            //how many points processor ii has in part i?
            proc_load_sort[ii].val =  p_gno_np_global_num_coord_each_part_actual[ii * num_parts + i];
        }
        uqsort<partId_t, pq_gno_t>(nprocs, proc_load_sort);

        //traverse all processors with increasing load.
        //TODO there should be a mistake here. (inefficiency)
        //should traverse from end to beginning.
        //currently gets the processors with least number of coordinates,
        //and assings the part to this part.
	//MD: uqsort sorts it in increasing order. 
	//We traverse it from end to beginning to get the highest number of coordinates.
        for (partId_t iii = nprocs - 1; iii >= 0; --iii){

            partId_t ii = proc_load_sort[iii].id;
            pq_lno_t left_space = space_in_each_processor[ii] - load;
            //if enought space, assign to this part.
            if(left_space >= 0 ){
                assigned_proc = ii;
                break;
            }
            //if space is not enough, store the best candidate part.
            if (space_in_each_processor[best_proc_to_assign] < space_in_each_processor[ii]){
                best_proc_to_assign = ii;
            }
        }
        if (assigned_proc == -1){
            assigned_proc = best_proc_to_assign;
        }
        space_in_each_processor[assigned_proc] -= load;

        //to sort later, part-i is assigned to the proccessor - assignment.
        part_assignment[j].id = i; //part i
        part_assignment[j].val = assigned_proc; //assigned to processor - assignment.


        //if assigned processor is me, increase the number.
        if (assigned_proc == myRank){
            out_num_part++;//assigned_part_count;
            out_part_index.push_back(i);
        }
        //increase the send to that processor by the number of points in that part.
        sendCount[assigned_proc] += p_gno_np_global_num_coord_each_part_actual[myRank * num_parts + i];
    }
    freeArray< uSortItem<partId_t, pq_gno_t> > (proc_load_sort);
    freeArray<uSortItem<partId_t, pq_gno_t> >(part_loads);
    freeArray<pq_lno_t >(space_in_each_processor);


    //sort assignments with respect to the assigned processors.
    uqsort<partId_t, partId_t>(num_parts, part_assignment);
    //fill sendBuf.

    if (doMigrationType == 0){
        fillContinousSendBuffer2< partId_t,  pq_lno_t,  pq_gno_t>(
            num_parts, partBegins, permutation, gnoList,
            part_assignment, //input sorted wrt processors
            sendBuf, partIndexBegin, newFuturePartitions, myRank);
    }
    else {
        Z1fillSendBuffer2< partId_t,  pq_lno_t>(
            num_parts, partBegins, permutation,
            part_assignment, //input sorted wrt processors
            coordinate_destionations, partIndexBegin, newFuturePartitions,
            myRank);
    }

    freeArray<uSortItem<partId_t, partId_t> >(part_assignment);
}

template <typename partId_t, typename pq_lno_t, typename pq_gno_t>
void procAssignment1(
    RCP<Comm<int> > &pcomm, //original communication.
    int assignment_type, //either assign to minimize migration, or assign to increase locality.
    pq_gno_t * p_gno_np_global_num_coord_each_part_actual,
    pq_gno_t nGlobalObj,
    pq_lno_t nLocal,
    partId_t num_parts,
    partId_t nprocs,
    partId_t myRank,
    pq_lno_t *partBegins, //holds the beginning of each part.
    pq_lno_t *permutation, //the permutation array ordered wrt partBegins array.
    const pq_gno_t *gnoList, //gno array
    pq_gno_t *sendBuf, //output: sized nLocal, the buffer is filled by the function with gnos.
    pq_lno_t *sendCount, //output: sized nprocs, show the number of send point counts to each proc.
    vector<partId_t> &ids, //output: this holds the id of the processors for the next subcommunicatior.
    //TODO futurePartIndex might need to be changed.
vector<partId_t> *newFuturePartitions,//input how many more partitions the part will be partitioned into.
    //partId_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
    partId_t &out_part_index, //output: the part index which the processor is assigned to.
    partId_t &partIndexBegin, //output: how much the part number should be shifted when setting the solution
    int doMigrationType,
    int *coordinate_destionations){

    pq_gno_t *p_gno_np_global_num_coord_each_part = p_gno_np_global_num_coord_each_part_actual + nprocs * num_parts;
    //pq_gno_t *p_gno_np_global_num_coord_each_part_mypart = p_gno_np_global_num_coord_each_part_actual + myRank * num_parts;


    partId_t *p_pid_np_num_procs_each_part = allocMemory<partId_t>(num_parts);

    bool did_i_find_my_group = false;

    partId_t left_proc = nprocs;
    partId_t min_required_for_rest = num_parts - 1;

    float max_difference = 0;
    partId_t max_differ_part = 0;

    //find how many processor each part requires.
    for (partId_t i=0; i < num_parts; i++){

        //scalar portion of the required processors
        float scalar_required_proc = nprocs *
                (float (p_gno_np_global_num_coord_each_part[i])
                        / float(nGlobalObj));

        //round it to the integer.
        partId_t required_proc = static_cast<partId_t> (
                floor (0.5f + scalar_required_proc));

        if (left_proc - required_proc < min_required_for_rest){
            required_proc = left_proc - (min_required_for_rest);
        }
        left_proc -= required_proc;
        --min_required_for_rest;
        p_pid_np_num_procs_each_part[i] = required_proc;

        //because of the roundings some processors might be left as unassigned.
        //we want to assign those processors to the part with most imbalance.
        float imbalance_wrt_ideal = (scalar_required_proc - required_proc) /  required_proc;
        if (imbalance_wrt_ideal > max_difference){
            max_difference = imbalance_wrt_ideal;
            max_differ_part = i;
        }
    }

    /*
    pcomm->barrier();
    cout << "me:" << pcomm->getRank() << " reached 1" << endl;
    */
    //assign extra processors to those parts.
    if (left_proc > 0){
        p_pid_np_num_procs_each_part[max_differ_part] +=  left_proc;
    }

    //now find what are the best processors with least migration for each part.
    //or assign the parts according to processors locality.

    partId_t *part_assign_begins = allocMemory<partId_t>(num_parts);
    partId_t *proc_chains = allocMemory<partId_t>(nprocs);
    partId_t *proc_part_assignments = allocMemory<partId_t>(nprocs);

    //initialize the assignment of each processor.
    //this has a linked list implementation.
    //the beginning of processors assigned
    //to each part is hold at  part_assign_begins[part].
    //then the next processor assigned to that part is located at
    //proc_part_assignments[part_assign_begins[part]], this is a chain
    //until the value of -1 is reached.
    for (int i = 0; i < nprocs; ++i ){
        proc_part_assignments[i] = -1;
        proc_chains[i] = -1;
    }
    for (int i = 0; i < num_parts; ++i ){
        part_assign_begins[i] = -1;
    }

    partId_t next_proc_to_assign = 0;

    /*
    pcomm->barrier();
    cout << "me:" << pcomm->getRank() << " reached 2" << endl;
     */
    //Allocate memory for sorting data structure.
    uSortItem<partId_t, pq_lno_t> * proc_points_in_part = allocMemory <uSortItem<partId_t, partId_t> > (nprocs);
    for(partId_t i = 0; i < num_parts; ++i){

        //if the assignment type is 0, then the algorithm tries to minimize the cost
        //of migration, by assigning the processors with highest number of coordinates on that part.
        if(assignment_type == 0){
            for(partId_t ii = 0; ii < nprocs; ++ii){

                proc_points_in_part[ii].id = ii;

                if (proc_part_assignments[ii] == -1){
                    proc_points_in_part[ii].val = p_gno_np_global_num_coord_each_part_actual[ii * num_parts + i];
                }
                else {
                    //if processor is already assigned, insert -nLocal - 1.
                    proc_points_in_part[ii].val = -p_gno_np_global_num_coord_each_part_actual[ii * num_parts + i] - 1;
                }
            }
        }
        else {

            //in this case we assign only the parts according to their closeness,
            //assuming that consecutive ranks are close to each other.
            partId_t required_proc = p_pid_np_num_procs_each_part[i];
            partId_t last_proc_to_assign = next_proc_to_assign + required_proc;
            for(partId_t ii = 0; ii < nprocs; ++ii){
                proc_points_in_part[ii].id = ii;
                if (ii >= next_proc_to_assign && ii < last_proc_to_assign){
                    proc_points_in_part[ii].val = p_gno_np_global_num_coord_each_part_actual[ii * num_parts + i];
                }
                else {
                    proc_points_in_part[ii].val = -p_gno_np_global_num_coord_each_part_actual[ii * num_parts + i] - 1;
                }
            }
            next_proc_to_assign = last_proc_to_assign;
        }

        //sort the processors in the part.
        uqsort<partId_t, partId_t>(nprocs, proc_points_in_part);

        partId_t required_proc_count =  p_pid_np_num_procs_each_part[i];
        pq_gno_t total_num_points_in_part = p_gno_np_global_num_coord_each_part[i];
        pq_gno_t ideal_num_points_in_procs =
                ceil (total_num_points_in_part / float (required_proc_count));

        //starts sending to heaviest part.
        partId_t next_part_to_send = nprocs - required_proc_count;
        partId_t next_part_to_send_id = proc_points_in_part[next_part_to_send].id;
        pq_lno_t space_left =  ideal_num_points_in_procs - proc_points_in_part[next_part_to_send].val;

        //find the assigned processors.
        for(partId_t ii = nprocs - 1; ii >= nprocs - required_proc_count; --ii){

            partId_t partid = proc_points_in_part[ii].id;
            proc_part_assignments[partid] = i;

            /*
            if(myRank == 0){
                cout << "proc:" << partid <<
                        " assigned to part:" << i <<
                        " as it has:" << proc_points_in_part[ii].val <<
                        " points in this part." << endl;
            }
            */
        }

        bool did_change_anything = false;
        for(partId_t ii = 0; ii < nprocs; ++ii){
            if (proc_points_in_part[ii].val < 0){
                did_change_anything = true;
                proc_points_in_part[ii].val = -proc_points_in_part[ii].val - 1;
            }
            else {
                break;
            }
        }

        if(did_change_anything){
            //resort the processors in the part.
            uqsort<partId_t, partId_t>(nprocs - required_proc_count, proc_points_in_part);
        }

        /*
        for(partId_t ii = nprocs - required_proc_count - 1; ii >= 0; --ii){

            partId_t partid = proc_points_in_part[ii].id;

            if(myRank == 0){
                cout << "proc:" << partid <<
                        " is not assigned to part:" << i <<
                        " as it has:" << proc_points_in_part[ii].val <<
                        " points in this part." << endl;
            }
        }
        */

        //check if this processors is one of the procs assigned to this part.
        //if it is, then get the group.
        if (!did_i_find_my_group){
            for(partId_t ii = nprocs - 1; ii >= nprocs - required_proc_count;
                                                             --ii){
                partId_t partid = proc_points_in_part[ii].id;
                //ids[nprocs - 1 - ii] = partid;
                ids.push_back(partid);
                //proc_part_assignments[partid] = i;

                if(partid == myRank){
                    did_i_find_my_group = true;
                    part_assign_begins[i] = myRank;
                    proc_chains[myRank] = -1;
                    sendCount[myRank] = proc_points_in_part[ii].val;
                    //TODO futurePartIndex might need to be changed.
                    //partIndexBegin += i * futurePartIndex;
                    for (partId_t in = 0; in < i; ++in){
                        partIndexBegin += (*newFuturePartitions)[in];
                    }

                    out_part_index = i;
                }
            }
            /*
            if (did_i_find_my_group){
                groupSize = required_proc_count;
            }
            */
            if (!did_i_find_my_group){
                ids.clear();
            }
        }

        //send points of the nonassigned coordinates to the assigned coordinates.
        //TODO play with this part, to get a better communication imbalance.
        //starts from the heaviest nonassigned processor.
        for(partId_t ii = nprocs - required_proc_count - 1; ii >= 0; --ii){
            partId_t partid = proc_points_in_part[ii].id;
            pq_lno_t to_sent = proc_points_in_part[ii].val;

            //we set number of points to -to_sent - 1 for the assigned processors.
            //we reverse it here.
            if (to_sent < 0) to_sent = -to_sent - 1;

            //now sends the
            while (to_sent > 0){
                //if the processor has enough space.

                if (to_sent <= space_left){
                    /*
                    if(myRank == 0){
                        cout << " proc:" << partid
                                << " needs to send:"
                                << to_sent << " and left space:"
                                << space_left << " on proc:"
                                << next_part_to_send_id << " because of part:"
                                << i << endl;
                    }
                    */

                    space_left -= to_sent;

                    if (myRank == partid){
                        //set my sent count to the sent processor.
                        sendCount[next_part_to_send_id] = to_sent;
                        //save the processor in the list (proc_chains and
                        //part_assign_begins)
                        //that the processor
                        //will send its point in part-i.
                        partId_t prev_begin = part_assign_begins[i];
                        part_assign_begins[i] = next_part_to_send_id;
                        proc_chains[next_part_to_send_id] = prev_begin;
                    }
                    to_sent = 0;
                }
                else {
                    /*
                    if(myRank == 0){
                        cout << " proc:" << partid
                                << " needs to send:"
                                << to_sent << " but left space:"
                                << space_left << " on proc:"
                                << next_part_to_send_id << " because of part:" << i << endl;
                    }
                    */
                    //there might be no space left in the processor.
                    if(space_left > 0){
                        to_sent -= space_left;

                        //send as the space left in the processor.
                        if (myRank == partid){
                            sendCount[next_part_to_send_id] = space_left;
                            partId_t prev_begin = part_assign_begins[i];
                            part_assign_begins[i] = next_part_to_send_id;
                            proc_chains[next_part_to_send_id] = prev_begin;

                        }
                    }
                    //change the sent part
                    ++next_part_to_send;

#ifdef debug_
                    //TODO remove comment
                    if(next_part_to_send <  nprocs - required_proc_count ){
                        cout << "this should not happen next part to send:" << next_part_to_send << endl;

                    }
#endif
                    next_part_to_send_id =  proc_points_in_part[next_part_to_send].id;
                    space_left = ideal_num_points_in_procs - proc_points_in_part[next_part_to_send].val;
                }
            }
        }

    }

    /*
    pcomm->barrier();
    cout << "me:" << pcomm->getRank() << " reached 3:" << nLocal << endl;
    */
    if (doMigrationType == 0){
        fillContinousSendBuffer1 <partId_t, pq_lno_t, pq_gno_t> ( pcomm,
            num_parts, partBegins, permutation, gnoList, nprocs, 
            part_assign_begins, proc_chains, sendCount, sendBuf);
    } else {
        Z1fillSendBuffer1 <partId_t, pq_lno_t> ( pcomm, num_parts, partBegins,
            permutation, nprocs, part_assign_begins, proc_chains, sendCount,
            coordinate_destionations);
    }
    /*
    pcomm->barrier();
    cout << "me:" << pcomm->getRank() << " reached 4" << endl;
    */
    freeArray<partId_t>(part_assign_begins);
    freeArray<partId_t>(proc_chains);
    freeArray<partId_t>(proc_part_assignments);
    freeArray<uSortItem<partId_t, partId_t> > (proc_points_in_part);
    freeArray<partId_t > (p_pid_np_num_procs_each_part);

}

template <typename partId_t, typename pq_lno_t, typename pq_gno_t>
void getProcGroups_SendCounts_SendBuff(
    RCP<Comm<int> > &pcomm, //original communication.
    int migration_proc_assignment_type,
    int assignment_type, //either assign to minimize migration, or assign to increase locality.
    pq_gno_t * p_gno_np_global_num_coord_each_part_actual,
    pq_gno_t nGlobalObj,
    pq_lno_t nLocal,
    partId_t num_parts,
    partId_t nprocs,
    partId_t myRank,
    pq_lno_t *partBegins, //holds the beginning of each part.
    pq_lno_t *permutation, //the permutation array ordered wrt partBegins array.
    const pq_gno_t *gnoList, //gno array

    pq_gno_t *sendBuf, //output: sized nLocal, the buffer is filled by the function with gnos.
    pq_lno_t *sendCount, //output: sized nprocs, show the number of send point counts to each proc.
    vector<partId_t> &ids, //output: this holds the id of the processors for the next subcommunicatior.
    //TODO futurePartIndex might need to be changed.
        //partId_t futurePartIndex,//input how many more partitions the part will be partitioned into.
    vector<partId_t> *newFuturePartitions,
    partId_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
    vector<partId_t> &out_part_indices, //output: the part index which the processor is assigned to.
    partId_t &partIndexBegin, //output: how much the part number should be shifted when setting the solution
    int doMigrationType,
    int *coordinate_destionations){

    /**fuction calculates
    *how much each processor should send the other processors.
    *fills the sendBuf, whose layout is chosen according to sendCount array.
    *returns the ids object, which is the rank id's that will be used for the next subCommunicator
    *returns groupSize, which is  the size of ids array.
    *returns how many parts the processor is assigned in out_num_part
    *returns out_part_indices, which is the vector of the partId's that the processors is assigned
    *returns partIndexBegin, which will be used when setting the solution.
    */

    ids.clear();
    /*
    if (myRank == 0){
        cout << endl << endl;
    }
    */
    if (1 || nLocal > 0){
        if (nprocs > num_parts){

            partId_t out_part_index = 0;
            /*
            pcomm->barrier();
            cout << "me:" << pcomm->getRank() << endl;
            */
            procAssignment1<partId_t, pq_lno_t, pq_gno_t>(
                pcomm,
                assignment_type, //either assign to minimize migration, or assign to increase locality.
                p_gno_np_global_num_coord_each_part_actual, nGlobalObj, nLocal,
                num_parts, nprocs, myRank, 
                partBegins, //holds the beginning of each part.
                permutation, //the permutation array ordered wrt partBegins array.
                gnoList, //gno array

                sendBuf, //output: sized nLocal, the buffer is filled by the function with gnos.
                sendCount, //output: sized nprocs, show the number of send point counts to each proc.
                ids, //output: this holds the id of the processors for the next subcommunicatior.

                //TODO futurePartIndex might need to be changed.
                newFuturePartitions,//input how many more partitions the part will be partitioned into.
                //partId_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
                out_part_index, //output: the part index which the processor is assigned to.
                partIndexBegin, //output: how much the part number should be shifted when setting the solution
                doMigrationType,
                coordinate_destionations
                );

            /*
            pcomm->barrier();
            cout << "me:" << pcomm->getRank() << " out"<< endl;
            */
            out_num_part = 1;
            out_part_indices.clear();
            out_part_indices.push_back(out_part_index);
            /*
            pcomm->barrier();
            cout << "me:" << pcomm->getRank() << " finish out"<< endl;
            */

            /*
            int partsendCount = 0;
            pq_lno_t sc = 0;
            for(int i = 0; i < nprocs; ++i){

                if(sendCount[i] > 0){
                    partsendCount++;
                    sc += sendCount[i];
                }
            }
            cout << "me:" << myRank <<
                    " sending to:" << partsendCount <<
                    " many processors with assignment-1:" <<
                    " I am assigned to:" << out_part_index << endl;

            for(int i = 0; i < num_parts; ++i){
                cout << "me:" << myRank <<
                        " have " << p_gno_np_global_num_coord_each_part_actual[myRank * num_parts + i] <<
                        " in part:" << i << endl;
            }
            */

        }
        else {
            ids.push_back(myRank);

            procAssignment2<partId_t, pq_lno_t, pq_gno_t>(
                assignment_type, //either assign to minimize migration, or assign to increase locality.
                p_gno_np_global_num_coord_each_part_actual,
                nGlobalObj,
                //pq_lno_t nLocal,
                num_parts,
                nprocs,
                myRank,
                partBegins, //holds the beginning of each part.
                permutation, //the permutation array ordered wrt partBegins array.
                gnoList, //gno array

                sendBuf, //output: sized nLocal, the buffer is filled by the function with gnos.
                sendCount, //output: sized nprocs, show the number of send point counts to each proc.

                //TODO futurePartIndex might need to be changed.
                newFuturePartitions,//input how many more partitions the part will be partitioned into.
                out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
                out_part_indices, //output: the part indices which the processor is assigned to.
                partIndexBegin, //output: how much the part number should be shifted when setting the solution
                doMigrationType,
                coordinate_destionations);

            /*
            pq_lno_t sc = 0;
            int partsendCount = 0;
            for(int i = 0; i < nprocs; ++i){
                if(sendCount[i] > 0){
                    partsendCount++;
                    sc += sendCount[i];
                }
            }

            cout << "me:" << myRank <<
                    " sending to:" << partsendCount <<
                    " many processors with assignment-2:" << endl;

            for (int i = 0; i < out_num_part; ++i){
                cout << "me:" << myRank << " I am assigned to:" << out_part_indices[i] << endl;
            }


            for(int i = 0; i < num_parts; ++i){
                cout << "me:" << myRank <<
                        " have " << p_gno_np_global_num_coord_each_part_actual[myRank * num_parts + i] <<
                        " in part:" << i << endl;
            }
            */
        }

    }
}


template <typename partId_t, typename pq_lno_t, typename pq_gno_t, typename pq_scalar_t>
void doAll2All(
    const RCP<const Environment> &env, //environment
    RCP<Comm<int> > &comm, //current communication object.
    int doMigrationType,
    int all2alloption,
    partId_t nprocs,
    pq_lno_t nLocal,
    pq_lno_t *sendCount,
    pq_gno_t *sendBuf,
    ArrayRCP<pq_gno_t> &recvBuf,
    pq_gno_t &numMyNewGnos,
    string iteration,
    int coord_dim, // coordinate dimension
    pq_scalar_t **coords, //coordinates.
    int weight_dim, //weight dimension
    pq_scalar_t **weights, //weights
    pq_gno_t *&coordinate_gnos,
    int *&coordinate_owners,
    int *coordinate_destionations,
    partId_t *&assigned_parts,
    partId_t num_parts){

    //function to obtain recvBuf that holds the new gno's that processor will own.

    if (doMigrationType == 0){
        if (all2alloption == 2) {
            //cout << "ALL2ALL distribotr" << endl;
            //uses distributor object.
            partId_t *partIds = allocMemory< partId_t>(nLocal);
            partId_t *p = partIds;

            //write which processor each point is going.
            for (int i = 0; i < nprocs; ++i){
                pq_lno_t sendC = sendCount[i];
                for (int ii = 0; ii < sendC; ++ii){
                    *(p++) = i;
                }
            }

            env->timerStart(MACRO_TIMERS,
                     "PQJagged - Migration DistPlanCreating-" + iteration);
            Tpetra::Distributor distributor(comm);

            ArrayView<const partId_t> pIds( partIds, nLocal);
            numMyNewGnos = distributor.createFromSends(pIds);
            env->timerStop(MACRO_TIMERS,
                     "PQJagged - Migration DistPlanCreating-" + iteration);

            ArrayRCP<pq_gno_t> recvBuf2(distributor.getTotalReceiveLength());

            env->timerStart(MACRO_TIMERS,
                     "PQJagged - Migration DistPlanCom-" + iteration);
            ArrayView<pq_gno_t> s(sendBuf, nLocal);
            distributor.doPostsAndWaits<pq_gno_t>(s, 1, recvBuf2());
            env->timerStop(MACRO_TIMERS,
                     "PQJagged - Migration DistPlanCom-" + iteration);
            recvBuf = recvBuf2;
            freeArray<partId_t>(partIds);

        } else if (all2alloption == 1){


            partId_t *partIds = allocMemory< partId_t>(nLocal);
            partId_t *p = partIds;

            //write which processor each point is going.

            for (int i = 0; i < nprocs; ++i){
                pq_lno_t sendC = sendCount[i];
                //cout << "me:" << comm->getRank() << " to:" << i <<
                //" sending:" << sendC << endl;
                for (int ii = 0; ii < sendC; ++ii){
                    *(p++) = i;
                }
            }


            ZOLTAN_COMM_OBJ *plan = NULL; /* pointer for communication object */


            MPI_Comm mpi_comm = Teuchos2MPI (comm);
            pq_lno_t incoming = 0;
            int message_tag = 7859;

            env->timerStart(MACRO_TIMERS,
                     "PQJagged - Migration Z1PlanCreating-" + iteration);
            int ierr = Zoltan_Comm_Create(&plan, nLocal, partIds, mpi_comm,
                message_tag, &incoming);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
            env->timerStop(MACRO_TIMERS, "PQJagged - Migration Z1PlanCreating-"
                + iteration);


            ArrayRCP<pq_gno_t> recvBuf2(incoming);
            pq_gno_t *recieves  = recvBuf2.getRawPtr();


            message_tag++;
            env->timerStart(MACRO_TIMERS, "PQJagged - Migration Z1PlanComm-" +
                iteration);
            ierr = Zoltan_Comm_Do(plan, message_tag, (char *) sendBuf,
                    sizeof(pq_gno_t), (char *) recieves);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
            env->timerStop(MACRO_TIMERS, "PQJagged - Migration Z1PlanComm-" +
                iteration);

            ierr = Zoltan_Comm_Destroy(&plan);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
            numMyNewGnos = incoming;
            recvBuf = recvBuf2;
            freeArray<partId_t>(partIds);

        } else {
            ArrayView<pq_gno_t> s (sendBuf, nLocal);
            ArrayView<int> sc (sendCount, nprocs);
            Array<pq_lno_t> recvCount(nprocs, 0);
            try{
                AlltoAllv<pq_gno_t>(*comm, *env,
                        s, sc,
                        recvBuf, recvCount());
            }
            Z2_FORWARD_EXCEPTIONS

            for (int i=0; i < nprocs; i++){
                numMyNewGnos += recvCount[i];
            }
            recvCount.clear();
        }
    }
    else {
/*
        int coord_dim, // coordinate dimension
        pq_scalar_t **coords, //coordinates.
        int weight_dim, //weight dimension
        pq_scalar_t **weights, //weights
        pq_gno_t * coordinate_gnos,
        int * coordinate_owners,
        int *coordinate_destionations,
        partId_t num_parts

*/
        ZOLTAN_COMM_OBJ *plan = NULL;


        MPI_Comm mpi_comm = Teuchos2MPI (comm);
        pq_lno_t incoming = 0;
        int message_tag = 7859;

        /*
        for (int i = 0; i < nLocal; ++i){
            cout << "me:" << comm->getRank() <<
                    " i:" << i <<
                    " before coordinate_destionations:" << coordinate_destionations[i] <<
                    endl;
        }
        */

        env->timerStart(MACRO_TIMERS, "PQJagged - Migration Z1PlanCreating-" +
            iteration);
        int ierr = Zoltan_Comm_Create( &plan, nLocal, coordinate_destionations,
            mpi_comm, message_tag, &incoming);
        Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
        env->timerStop(MACRO_TIMERS, "PQJagged - Migration Z1PlanCreating-" +
            iteration);

        int toMyself = 0;
        for (int i = 0; i < nLocal; ++i){
            if (coordinate_destionations[i] == comm->getRank()){
                ++toMyself;
            }
        }
        /*
        cout << "iteration:" << iteration <<
                " me:" << comm->getRank() <<
                " mySelf:" << toMyself <<
                " nLocal:" << nLocal << endl;
         */
        /*
        cout << "me:" << comm->getRank() << " incoming:" << incoming << endl;
        */

        pq_gno_t *incoming_gnos = allocMemory< pq_gno_t>(incoming);

        message_tag++;
        env->timerStart(MACRO_TIMERS, "PQJagged - Migration Z1PlanComm-" +
            iteration);
        ierr = Zoltan_Comm_Do( plan, message_tag, (char *) coordinate_gnos,
                sizeof(pq_gno_t), (char *) incoming_gnos);
        Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
        env->timerStop(MACRO_TIMERS, "PQJagged - Migration Z1PlanComm-" +
            iteration);
        freeArray<pq_gno_t>(coordinate_gnos);
        coordinate_gnos = incoming_gnos;

        env->timerStart(MACRO_TIMERS, "PQJagged - Migration Z1Migration-" +
            iteration);
        for (int i = 0; i < coord_dim; ++i){
            message_tag++;
            pq_scalar_t *coord = coords[i];

            coords[i] = allocMemory<pq_scalar_t>(incoming);
            ierr = Zoltan_Comm_Do( plan, message_tag, (char *) coord,
                    sizeof(pq_scalar_t), (char *) coords[i]);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
            freeArray<pq_scalar_t>(coord);
        }

        for (int i = 0; i < weight_dim; ++i){
            message_tag++;
            pq_scalar_t *weight = weights[i];

            weights[i] = allocMemory<pq_scalar_t>(incoming);
            ierr = Zoltan_Comm_Do( plan, message_tag, (char *) weight,
                    sizeof(pq_scalar_t), (char *) weights[i]);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
            freeArray<pq_scalar_t>(weight);
        }


        int *coord_own = allocMemory<int>(incoming);
        message_tag++;
        ierr = Zoltan_Comm_Do( plan, message_tag, (char *) coordinate_owners,
                sizeof(int), (char *) coord_own);
        Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
        freeArray<int>(coordinate_owners);
        coordinate_owners = coord_own;

        /*
        for (int i = 0; i < incoming; ++i){
            cout << "me:" << comm->getRank() <<
                    " i:" << i <<
                    " after move:" << coordinate_owners[i] <<
                    endl;
        }
        */

        partId_t *new_parts = allocMemory<int>(incoming);

        if(nprocs < num_parts){
            message_tag++;
            ierr = Zoltan_Comm_Do( plan, message_tag, (char *) assigned_parts,
                    sizeof(partId_t), (char *) new_parts);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
        }

        freeArray<partId_t>(assigned_parts);
        assigned_parts = new_parts;
        ierr = Zoltan_Comm_Destroy(&plan);
        Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
        env->timerStop(MACRO_TIMERS, "PQJagged - Migration Z1Migration-" +
            iteration);
        numMyNewGnos = incoming;
    }
}

template <typename mvector_t, typename pq_gno_t>
void doActualMigration(
    const RCP<const Environment> &env, //environment
    RCP<const mvector_t> &vectors,    // on return is the new data,
    pq_gno_t numMyNewGnos,
    ArrayRCP<pq_gno_t> recvBuf){

    try{
        vectors = XpetraTraits<mvector_t>::doMigration(
                vectors, numMyNewGnos, recvBuf.getRawPtr()/*, env*/);
    }
    Z2_FORWARD_EXCEPTIONS
}

template <typename partId_t>
void createSubCommunicator(
    RCP<Comm<int> > &comm,
    vector<partId_t> &proc_ids){

    partId_t groupSize = proc_ids.size();
    partId_t *ids = allocMemory<partId_t>(groupSize);
    for(partId_t i = 0; i < groupSize; ++i) {
        //cout << "ids:" << i << " is:" << proc_ids[i] << endl;
        ids[i] = proc_ids[i];
    }

    ArrayView<const partId_t> idView(ids, groupSize);
    comm = comm->createSubcommunicator(idView);
    freeArray<partId_t>(ids);
}

template <typename mvector_t, typename pq_lno_t, typename pq_gno_t, typename pq_scalar_t, typename pq_node_t>
void create_new_multi_vector(
    const RCP<const Environment> &env,
    RCP<Comm<int> > &comm,
    RCP <const mvector_t> &mvector,
    int multiVectorDim){

    typedef ArrayView<const pq_scalar_t> coordList_t;
    typedef Tpetra::Map<pq_lno_t, pq_gno_t, pq_node_t> map_t;

    ArrayView<const pq_gno_t> gnoList = mvector->getMap()->getNodeElementList();
    size_t localSize = mvector->getLocalLength();

    // Tpetra will calculate the globalSize.
    size_t globalSize = Teuchos::OrdinalTraits<size_t>::invalid();

    RCP<map_t> subMap;
    try{
        subMap= rcp(new map_t(globalSize, gnoList, 0, comm));
    }
    Z2_THROW_OUTSIDE_ERROR(*env)

    coordList_t *avSubList = allocMemory<coordList_t>(multiVectorDim);

    for (int dim=0; dim < multiVectorDim; dim++)
        avSubList[dim] = mvector->getData(dim).view(0, localSize);

    ArrayRCP<const ArrayView<const pq_scalar_t> > subVectors =
            arcp(avSubList, 0, multiVectorDim);


    try{
        mvector = rcp(new mvector_t(
                subMap, subVectors.view(0, multiVectorDim), multiVectorDim));
    }
    Z2_THROW_OUTSIDE_ERROR(*env)
}

template <typename pq_lno_t, typename partId_t>
void resizeArrays(
    pq_lno_t numLocalPoints,
    pq_lno_t prev_num_local,
    partId_t *& partArray,
    partId_t *&permutation,
    partId_t *&oldpermutation){

    if (prev_num_local != numLocalPoints){
        freeArray<partId_t>(partArray);
        partArray = allocMemory<partId_t>(numLocalPoints);
        freeArray<pq_lno_t>(permutation);
        freeArray<pq_lno_t>(oldpermutation);
        oldpermutation = allocMemory<pq_lno_t>(numLocalPoints);
        permutation = allocMemory<pq_lno_t>(numLocalPoints);
    }

}


template <typename pq_gno_t, typename pq_lno_t,typename pq_scalar_t, typename pq_node_t>
void getNewMultivectorArrays(
    RCP<Comm<int> > &comm,
    pq_lno_t &numLocalPoints,
    int coord_dim,
    pq_scalar_t **coords,
    int weight_dim,
    pq_scalar_t **weight,
    RCP<const Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> > coordsConst,
    partId_t out_num_parts,
    pq_lno_t *&permutations,
    pq_lno_t *part_begins,
    partId_t num_parts,
    int pqJagged_multiVectorDim
){

    for (int i=0; i < coord_dim; i++){
        ArrayRCP< const pq_scalar_t >  coord = coordsConst->getData(i);
        coords[i] =(pq_scalar_t *) coord.getRawPtr();
    }

    for (int i=0; i < weight_dim; i++){
        if (i + coord_dim >= pqJagged_multiVectorDim) break;
        ArrayRCP< const pq_scalar_t >  wgts = coordsConst->getData(i+ coord_dim);
        weight[i] = (pq_scalar_t *) wgts.getRawPtr();
    }
    //cout << "me:" << comm->getRank() << " obtained till weights:" << endl;

    if (out_num_parts == 1){
        //cout << "me:" << comm->getRank() << " filling permutation:" << endl;
        for(pq_lno_t i = 0; i < numLocalPoints; ++i){
            permutations[i] = i;
        }
        part_begins[0] = numLocalPoints;

        //cout << "me:" << comm->getRank() << " filled permutation:" << endl;
    }
    else {
        pq_scalar_t *assigned_parts = (pq_scalar_t *)coordsConst->getData /*getData*/(pqJagged_multiVectorDim).getRawPtr();
        pq_lno_t *counts = allocMemory<pq_lno_t>(num_parts);

        partId_t *part_shifts = allocMemory<partId_t>(num_parts);

        memset(counts, 0, sizeof(pq_lno_t) * num_parts);

        for(pq_lno_t i = 0; i < numLocalPoints; ++i){
            partId_t ii = assigned_parts[i];
            ++counts[ii];
        }
        partId_t p = 0;
        pq_lno_t prev_index = 0;
        for(partId_t i = 0; i < num_parts; ++i){
            //cout << "me:" << comm->getRank() << " p:" << i << " count:" <<  counts[i] << " a:" << assigned_parts << endl;
            if(counts[i] > 0)  {
                part_begins[p] =  prev_index + counts[i];
                prev_index += counts[i];
                part_shifts[i] = p++;
            }
        }
        partId_t assigned_count = p - 1;

        for (;p < num_parts; ++p){
            part_begins[p] =  part_begins[assigned_count];
        }
        for(partId_t i = 0; i < out_num_parts; ++i){
            counts[i] = part_begins[i];
        }
        for(pq_lno_t i = numLocalPoints - 1; i >= 0; --i){
            partId_t part = part_shifts[partId_t(assigned_parts[i])];
            permutations[--counts[part]] = i;
        }

        freeArray<pq_lno_t>(counts);
        freeArray<partId_t>(part_shifts);
    }
}


template <typename pq_lno_t, typename partId_t>
void fillPermutationArrays(
    RCP<Comm<int> > &comm,
    pq_lno_t &numLocalPoints,
    partId_t out_num_parts,
    pq_lno_t *permutations,
    pq_lno_t *part_begins,
    partId_t *assigned_parts,
    partId_t num_parts
){
    if (out_num_parts == 1){
        //cout << "me:" << comm->getRank() << " filling permutation:" << endl;
        for(pq_lno_t i = 0; i < numLocalPoints; ++i){
            permutations[i] = i;
        }
        part_begins[0] = numLocalPoints;
    }
    else {
        pq_lno_t *counts = allocMemory<pq_lno_t>(num_parts);

        partId_t *part_shifts = allocMemory<partId_t>(num_parts);

        memset(counts, 0, sizeof(pq_lno_t) * num_parts);

        for(pq_lno_t i = 0; i < numLocalPoints; ++i){
            partId_t ii = assigned_parts[i];
            ++counts[ii];
        }
/*
        for(partId_t i = 0; i < num_parts; ++i){
            cout << "me:" <<  comm->getRank() << " i:" << i << " c:" << counts[i] << endl;
        }
*/
        partId_t p = 0;
        pq_lno_t prev_index = 0;
        for(partId_t i = 0; i < num_parts; ++i){
            //cout << "me:" << comm->getRank() << " p:" << i << " count:" <<  counts[i] << " a:" << assigned_parts << endl;
            if(counts[i] > 0)  {
                part_begins[p] =  prev_index + counts[i];
                prev_index += counts[i];
                part_shifts[i] = p++;
            }
        }
        partId_t assigned_count = p - 1;

        for (;p < num_parts; ++p){
            part_begins[p] =  part_begins[assigned_count];
        }
        for(partId_t i = 0; i < out_num_parts; ++i){
            counts[i] = part_begins[i];
        }
        for(pq_lno_t i = numLocalPoints - 1; i >= 0; --i){
            partId_t part = part_shifts[partId_t(assigned_parts[i])];
            permutations[--counts[part]] = i;
        }

        freeArray<pq_lno_t>(counts);
        freeArray<partId_t>(part_shifts);
    }
}


template <typename pq_gno_t,
          typename pq_lno_t,
          typename pq_scalar_t,
          typename pq_node_t,
          typename partId_t>
bool migration_refactored(
    RCP<Comm<int> > &pcomm, //original communication.
    const RCP<const Environment> &env, //environment
    RCP<Comm<int> > &comm, //current communication object.

    RCP<const Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> > &vectors, //multivector
    int pqJagged_multiVectorDim, //multivector dimension
    pq_gno_t &numGlobalPoints, //numGlobal points, output
    pq_lno_t &numLocalPoints, //numLocal points, output
    int coord_dim, // coordinate dimension
    pq_scalar_t **coords, //coordinates.
    int weight_dim, //weight dimension
    pq_scalar_t **weight, //weights
    partId_t * &assigned_parts_, //this should not be necessary anymore.
    partId_t num_parts, //current num parts
    partId_t &out_num_part, //output num parts.
    vector<partId_t> *newFuturePartitions,
    pq_lno_t *&permutation,
    pq_lno_t *&old_permutation,
    pq_lno_t *partBeginArray,
    partId_t &partIndexBegin,
    //partId_t futurePartIndex,
    int all2alloption,
    int assignment_type,
    int doMigrationType,
    int migration_check_option,
    pq_scalar_t migration_imbalance_cut_off,
    size_t futureReduceAll,
    pq_lno_t numCoordinatesForLastDimPartitioning,
    pq_gno_t *&coordinate_gnos,
    int *&actual_gno_owner,
    string iteration,
    int keep_part_boxes,
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > &inPartBoxes,
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > &outPartBoxes
){
/*
    if(pcomm->getRank() == 0){
        cout << "doMigrationType:" << doMigrationType << endl;
        cout << "all2alloption:" << all2alloption << endl;
        cout << "assignment_type:" << assignment_type << endl;
        cout << "migration_check_option:" << migration_check_option << endl;
    }
*/
    int migration_proc_assignment_type = 1;
    partId_t nprocs = comm->getSize();
    partId_t myRank = comm->getRank();
    int multiVectorDim = pqJagged_multiVectorDim;


    pq_gno_t *p_gno_np_global_num_coord_each_part_actual = NULL;

    //get the number of coordinates in each part in each processor.
    //p_gno_np_global_num_coord_each_part_actual is allocated in this function.
    //allocation_size is returned, and this is the size of p_gno_np_global_num_coord_each_part_actual array.
    //size_t allocation_size =
    getProcessorCoordinatePartCounts <pq_gno_t, pq_lno_t ,partId_t>(
            pcomm, env, comm, migration_proc_assignment_type, nprocs, myRank,
            num_parts, partBeginArray,
            p_gno_np_global_num_coord_each_part_actual);


    //check if migration will be performed or not.
    if (!checkMigration <pq_gno_t, pq_lno_t ,partId_t>(
            pcomm, env, comm, migration_check_option, futureReduceAll,
            numCoordinatesForLastDimPartitioning, migration_imbalance_cut_off, 
            migration_proc_assignment_type, nprocs, myRank, num_parts, 
            p_gno_np_global_num_coord_each_part_actual, partBeginArray)){
        freeArray<pq_gno_t>(p_gno_np_global_num_coord_each_part_actual);
        return false;
    }

    //TODO only do it when doMigration is used.
    if (nprocs < num_parts) {
        if (doMigrationType == 0){
            vectors = createNewMultivector(
                    comm, numGlobalPoints, numLocalPoints, coord_dim, coords, 
                    weight_dim, weight, num_parts, partBeginArray, permutation, 
                    multiVectorDim, vectors);
        } else {

            for (partId_t i = 0; i < num_parts; ++i){
                pq_lno_t pBegin = 0;
                if (i > 0) pBegin = partBeginArray[i - 1];
                pq_lno_t pEnd = partBeginArray[i];

                for (pq_lno_t j = pBegin; j < pEnd; ++j){
                    pq_lno_t ind = permutation[j];
                    assigned_parts_[ind] = i;
                }
            }
        }
    }

    const pq_gno_t *gnoList = NULL;
    pq_lno_t *sendCount = NULL;
    pq_gno_t *sendBuf = NULL;
    int *coordinate_destionations = NULL;
    sendCount = allocMemory<pq_lno_t>(nprocs);
    for (int i = 0; i < nprocs; ++i) sendCount[i] = 0;

    if (doMigrationType == 0){

        sendBuf = allocMemory<pq_gno_t>(numLocalPoints);
        ArrayView<const pq_gno_t> gno_list = vectors->getMap()->getNodeElementList();
        gnoList = gno_list.getRawPtr();
    }
    else {
        coordinate_destionations = allocMemory<int>(numLocalPoints);
    }

    vector<partId_t> ids;
    vector<partId_t> out_part_indices;

    getProcGroups_SendCounts_SendBuff<partId_t, pq_lno_t, pq_gno_t>(
            pcomm, migration_proc_assignment_type,
            assignment_type, //either assign to minimize migration, or
                             // assign to increase locality.  
            p_gno_np_global_num_coord_each_part_actual, numGlobalPoints,
            numLocalPoints, num_parts, nprocs, myRank, 
            partBeginArray, //holds the beginning of each part.
            permutation, //the permutation array ordered wrt partBegins array.
            gnoList, //gno array 
            sendBuf, //output: sized nLocal, the buffer is filled by the function with gnos.
            sendCount, //output: sized nprocs, show the number of send point counts to each proc.
            ids, //output: this holds the id of the processors for the next subcommunicatior.

            //TODO futurePartIndex might need to be changed.
            //futurePartIndex,//input how many more partitions the part will be partitioned into.
            newFuturePartitions,
            out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
            out_part_indices, //output: the part index which the processor is assigned to.
            partIndexBegin, //output: how much the part number should be shifted when setting the solution
            doMigrationType,
            coordinate_destionations
               );
    vector <partId_t> tmpv;
    std::sort (out_part_indices.begin(), out_part_indices.end());
    partId_t outP = out_part_indices.size();

    pq_gno_t newGlobalCount = 0;
    pq_gno_t *p_gno_np_global_num_coord_each_part = p_gno_np_global_num_coord_each_part_actual + nprocs * num_parts;

    if (keep_part_boxes){
        inPartBoxes->clear();
    }

    for (partId_t i = 0; i < outP; ++i){
        partId_t ind = out_part_indices[i];

        newGlobalCount += p_gno_np_global_num_coord_each_part[ind];

        tmpv.push_back((*newFuturePartitions)[ind]);

        if (keep_part_boxes){
            //cout << "me:" << pcomm->getRank() << " is assigned to:" << ind
            //<< endl;
            //(*outPartBoxes)[ind].print();
            inPartBoxes->push_back((*outPartBoxes)[ind]);
        }
    }

    if (keep_part_boxes){
        RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > tmpPartBoxes = inPartBoxes;
        inPartBoxes = outPartBoxes;
        outPartBoxes = tmpPartBoxes;
    }

    newFuturePartitions->clear();

    for (partId_t i = 0; i < outP; ++i){
        partId_t p = tmpv[i];
        newFuturePartitions->push_back(p);
    }

    ArrayRCP<pq_gno_t> recvBuf;
    pq_gno_t numMyNewGnos = 0;

    /*
    pcomm->barrier();
    cout << "it:" << iteration << " in doAll2All:" << endl;
    */
    env->timerStart(MACRO_TIMERS, "PQJagged - Migration AlltoAll-" + iteration);
    doAll2All<partId_t, pq_lno_t, pq_gno_t, pq_scalar_t>( env, comm,
        doMigrationType, all2alloption, nprocs, numLocalPoints, sendCount,
        sendBuf, recvBuf, numMyNewGnos, iteration, coord_dim, coords,
        weight_dim, weight, coordinate_gnos, actual_gno_owner,
        coordinate_destionations, assigned_parts_, num_parts);
    /*
    pcomm->barrier();
    cout << "it:" << iteration << " out doAll2All:" << endl;
    */

    env->timerStop(MACRO_TIMERS, "PQJagged - Migration AlltoAll-" + iteration);

    if (doMigrationType == 0){
        freeArray<pq_gno_t>(sendBuf);
    }
    else {
        freeArray<int>(coordinate_destionations);

        if(numLocalPoints != numMyNewGnos){
            freeArray<pq_lno_t>(permutation);
            freeArray<pq_lno_t>(old_permutation);

            permutation = allocMemory<pq_lno_t>(numMyNewGnos);
            old_permutation = allocMemory<pq_lno_t>(numMyNewGnos);
        }
        numLocalPoints = numMyNewGnos;
        numGlobalPoints = newGlobalCount;
    }
    freeArray<pq_lno_t>(sendCount);
    freeArray<pq_gno_t>(p_gno_np_global_num_coord_each_part_actual);

    typedef Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> mvector_t;
    if(doMigrationType == 0){
        env->timerStart(MACRO_TIMERS, "PQJagged - Migration doActualMigration-"
            + iteration);
        // vectors on return will have the new data
        doActualMigration<mvector_t, pq_gno_t>( env, vectors, numMyNewGnos,
            recvBuf);

        env->timerStop(MACRO_TIMERS, "PQJagged - Migration doActualMigration-"
            + iteration);
    }
    createSubCommunicator <partId_t>(comm, ids);
    ids.clear();

    if(doMigrationType == 0){
        create_new_multi_vector<mvector_t, pq_lno_t, pq_gno_t,pq_scalar_t, pq_node_t>(
        env, comm, vectors, multiVectorDim);

        pq_lno_t prev_num_local = numLocalPoints;
        numLocalPoints = vectors->getLocalLength();
        numGlobalPoints = vectors->getGlobalLength();

        resizeArrays<pq_lno_t, partId_t>( numLocalPoints, prev_num_local,
            assigned_parts_, permutation, old_permutation);

        getNewMultivectorArrays<pq_gno_t, pq_lno_t,pq_scalar_t, pq_node_t>(
            comm, numLocalPoints, coord_dim, coords, weight_dim, weight,
            vectors, out_num_part, permutation, partBeginArray, num_parts,
            pqJagged_multiVectorDim );
    }
    else {
        fillPermutationArrays<pq_lno_t,partId_t>( comm, numLocalPoints,
            out_num_part, permutation, partBeginArray, assigned_parts_,
            num_parts);
    }

    return true;
}

#endif

template <typename pq_scalar_t, typename partId_t>
partId_t getPartitionArrays(
    const partId_t *partNo,
    vector <partId_t> &pAlongI, //assumes this vector is empty.
    vector<partId_t> *currentPartitions,
    vector<partId_t> *newFuturePartitions, //assumes this vector is empty.
    partId_t &futurePartNumbers,

    partId_t currentPartitionCount,
    int partArraySize,
    int i,
    partId_t maxPartNo,
    int keep_part_boxes,
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > inPartBoxes,
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > outPartBoxes
){
    partId_t outPartCount = 0;
    if(partNo){
        //when the partNo array is provided as input,
        //each current partition will be partition to the same number of parts.
        //we dont need to use the currentPartition vector in this case.

        partId_t p = partNo[i];
        if (p < 1){
            cout << "i:" << i << " p is given as:" << pAlongI[0] << endl;
            exit(1);
        }
        if (p == 1){
            return currentPartitionCount;
        }

        for (partId_t ii = 0; ii < currentPartitionCount; ++ii){
            pAlongI.push_back(p);

        }

        //TODO this should be removed.
        futurePartNumbers /= pAlongI[0];
        outPartCount = currentPartitionCount * pAlongI[0];

        if (keep_part_boxes){
            for (partId_t k = 0; k < currentPartitionCount; ++k){
                for (partId_t j = 0; j < pAlongI[0]; ++j){
                    outPartBoxes->push_back((*inPartBoxes)[k]);
                }
            }
        }

        //set the how many more parts each part will be divided.
        //this is obvious when partNo array is provided as input.
        //however, fill this so that weights will be calculated according to this array.
        for (partId_t ii = 0; ii < outPartCount; ++ii){
            newFuturePartitions->push_back(futurePartNumbers);
        }
    }
    else {
        //if partNo array is not provided as input,
        //currentPartitions  hold how many parts each part should be divided.
        //initially it has single number with total number of global parts.

        //calculate the futurePartNumbers from beginning,
        //since each part might be divided into different number of parts.
        futurePartNumbers = 1; //TODO this should be removed.

        //cout << "i:" << i << endl;
        double fEpsilon = numeric_limits<double>::epsilon();
        for (partId_t ii = 0; ii < currentPartitionCount; ++ii){
            //get how many parts a part should be divided.
            partId_t numFuture = (*currentPartitions)[ii];

            //get the ideal number of parts that is close to the
            //(partArraySize - i) root of the numFuture.
            partId_t numParts = getPartCount<partId_t>( numFuture,
                                 1.0 / (partArraySize - i), fEpsilon);
            //partId_t numParts = ceil( pow(numFuture, 1.0f / (partArraySize - i)));// + 0.5f;

            //cout << "\tii:" << ii << " numParts:" << numParts << endl;
            //cout << "\tii:" << ii << " numFuture:" << numFuture << endl;
            if (numParts > maxPartNo){
                cerr << "ERROR: maxPartNo calculation is wrong." << endl;
                exit(1);
            }
            //add this number to pAlongI vector.
            pAlongI.push_back(numParts);


            //increase the output number of parts.
            outPartCount += numParts;

            //ideal number of future partitions for each part.
            partId_t idealNumFuture = numFuture / numParts;
            for (partId_t iii = 0; iii < numParts; ++iii){
                partId_t fNofCuts = idealNumFuture;

                if (iii < numFuture % numParts){
                    //if not uniform, add 1 for the extra parts.
                    ++fNofCuts;
                }
                newFuturePartitions->push_back(fNofCuts);

                if (keep_part_boxes){
                    //for (partId_t j = 0; j < numParts; ++j){
                        outPartBoxes->push_back((*inPartBoxes)[ii]);
                    //}
                }


                //TODO this should be removed.
                if (fNofCuts > futurePartNumbers) futurePartNumbers = fNofCuts;
            }
        }
    }
    return outPartCount;
}

/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 * \param minCoordinate minimum coordinate in the range.
 * \param maxCoordinate maximum coordinate in the range.
 *
 * \param currentPart is the index of the part in the inTotalCounts vector.
 * \param inTotalCounts holds the beginning and end of the coordinates in each part on partitionedPointCoordinates array.
 * \param partitionedPointCoordinates is the permutation array, holds the real indices of coordinates on pqCoord array.
 * \param pqCoord is the 1D array holding the coordinates.
 * \param partIds is the array holding the partIds of each coordinate.
 *
 * \param _EPSILON is min element can be differentiated by pq_scalar_t
 * \param partition is the number of parts that the current part will be partitioned into.
 */

template <typename pq_scalar_t, typename pq_lno_t, typename partId_t>
void getInitialPartAssignments(
    pq_scalar_t &maxCoordinate,
    pq_scalar_t &minCoordinate,
    partId_t &currentPart,
    pq_lno_t *inTotalCounts,
    pq_lno_t *partitionedPointCoordinates,
    pq_scalar_t *pqCoord,
    partId_t *partIds,
    pq_scalar_t _EPSILON,
    partId_t &partition
){
    pq_scalar_t coordinate_range = maxCoordinate - minCoordinate;
    //partId_t currentPart = currentWorkPart + kk;
    pq_lno_t coordinateEnd= inTotalCounts[currentPart];
    pq_lno_t coordinateBegin = currentPart==0 ? 0: inTotalCounts[currentPart -1];

    //if there is single point, or if all points are along a line.
    //set initial part to 0 for all.
    if(ABS(coordinate_range) < _EPSILON ){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
            partIds[partitionedPointCoordinates[ii]] = 0;
        }
    }
    else{

        //otherwise estimate an initial part for each coordinate.
        //assuming uniform distribution of points.
        pq_scalar_t slice = coordinate_range / partition;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

            pq_lno_t iii = partitionedPointCoordinates[ii];
            partId_t pp = partId_t((pqCoord[iii] - minCoordinate) / slice);
            partIds[iii] = 2 * pp;
        }
    }
}



/* Either the pqParts array (partNo) or numGlobalParts should be provided in
 * the input. partArraySize should be provided. partNo takes
 * precedence if both are provided. The rest of the parameters are all
 * output. Providing pqParts array also requires its size be provided.
 * */
template <typename partId_t>
void getPartSpecifications(
    const partId_t *partNo,     // pqParts array
    int partArraySize,
    size_t &numGlobalParts,
    partId_t &totalDimensionCut, //how many cuts will be totally
    partId_t &totalPartCount ,    //how many parts will be totally
    partId_t &maxPartNo ,         //maximum cut count along a dimension.
    partId_t &reduceAllCount ,    //estimate on #reduceAlls that can be done.
    partId_t &maxTotalCumulativePartCount, //max no of parts that might occur
                                           //during the partition before the
                                           //last partitioning dimension.
    partId_t &maxCutNo,
    size_t &maxTotalPartCount
){
    if (partNo){
        for (int i = 0; i < partArraySize; ++i){
            reduceAllCount += totalPartCount;
            totalPartCount *= partNo[i];
            if(partNo[i] > maxPartNo) maxPartNo = partNo[i];
        }
        maxTotalCumulativePartCount = totalPartCount / partNo[partArraySize-1];
        numGlobalParts = totalPartCount;
    } else {
        float fEpsilon = numeric_limits<float>::epsilon();
        partId_t futureNumParts = numGlobalParts;

        for (int i = 0; i < partArraySize; ++i){
            partId_t maxNoPartAlongI = getPartCount<partId_t>( futureNumParts,
                                        1.0f / (partArraySize - i), fEpsilon);
            //cout << "futureNumParts:" << futureNumParts << "partArraySize:"
            // << partArraySize << " i:" << i << "maxNoPartAlongI:" <<i
            // maxNoPartAlongI<< endl;
            //partId_t maxNoPartAlongI = ceil(pow(futureNumParts,
            //1.0f / (coordDim - i)));// + 0.5f;
            if (maxNoPartAlongI > maxPartNo){
                maxPartNo = maxNoPartAlongI;
            }
            //cout << "i:" << i << "maxPartNo:" << maxPartNo<< endl;

            partId_t nfutureNumParts = futureNumParts / maxNoPartAlongI;
            if (futureNumParts % maxNoPartAlongI){
                ++nfutureNumParts;
            }
            futureNumParts = nfutureNumParts;
            //cout << "i:" << i << " fp:" << futureNumParts << endl;
        }
        totalPartCount = numGlobalParts;
        //estimate reduceAll Count here.
        //we find the upperbound instead.
        partId_t p = 1;
        for (int i = 0; i < partArraySize; ++i){
            reduceAllCount += p;
            p *= maxPartNo;
        }

        //cout << "maxPartNo:" << maxPartNo << endl;
        maxTotalCumulativePartCount  = p / maxPartNo;
    }

    totalDimensionCut = totalPartCount - 1;
    maxCutNo = maxPartNo - 1;
    maxTotalPartCount = maxPartNo + size_t(maxCutNo);
    //maxPartNo is P, maxCutNo = P-1, matTotalPartcount = 2P-1
}



/*! \brief Function that determines the permutation indices of the coordinates.
 * \param partNo is the number of parts.
 * \param noThreads is the number of threads avaiable for each processor.
 *
 * \param partitionedPointPermutations is the indices of coordinates in the given partition.
 * \param pqJagged_coordinates is 1 dimensional array holding the coordinate values.
 * \param pqJagged_uniformWeights is a boolean value if the points have uniform weights.
 * \param coordWeights is 1 dimensional array holding the weights of points.
 *
 * \param cutCoordinates is 1 dimensional array holding the cut coordinates.
 * \param coordinateBegin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinateEnd is the end index of the given partition on partitionedPointPermutations.
 *
 * \param allowNonRectelinearPart is the boolean value whether partitioning should allow distributing the points on same coordinate to different parts.
 * \param actual_ratios holds how much weight of the coordinates on the cutline should be put on left side.
 * \param localPartWeights is the local totalweight of the processor.
 * \param partWeights is the two dimensional array holding the weight of parts for each thread. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param nonRectelinearRatios is the two dimensional work array holding ratios of weights to be put left and right of the cut line.
 * \param partPointCounts is the two dimensional array holding the number of points in each part for each thread.
 *
 * \param newpartitionedPointPermutations is the indices of coordinates calculated for the partition on next dimension.
 * \param totalCounts are the number points in each output part.
 * \param partIds is the array that holds the part ids of the coordinates
 *
 * \param isSequentialAndPricise is the boolean value. This is true when the pricise and deterministic
 * result is required.
 * \param pqJaggedAllCoordinates is the array that holds the coordinates of points.
 * \param pqCoordDim is dimension of the input
 * \param currentCoord is the index according to which the partitioning is done.
 */
template <typename pq_lno_t, typename pq_scalar_t>
void getChunksFromCoordinates(
    partId_t partNo,
    int noThreads,
    pq_lno_t *partitionedPointPermutations,
    pq_scalar_t *pqJagged_coordinates,
    bool pqJagged_uniformWeights,
    pq_scalar_t *coordWeights,
    pq_scalar_t *cutCoordinates,
    pq_lno_t coordinateBegin,
    pq_lno_t coordinateEnd,
    bool allowNonRectelinearPart,
    float *actual_ratios,
    pq_scalar_t *localPartWeights,
    double **partWeights,
    float **nonRectelinearRatios,
    pq_lno_t ** partPointCounts,
    pq_lno_t *newpartitionedPointPermutations,
    pq_lno_t *totalCounts,
    partId_t *partIds,
    bool isSequentialAndPricise,
    pq_scalar_t **pqJaggedAllCoordinates,
    int pqCoordDim,
    int currentCoord
){

    //pq_lno_t numCoordsInPart =  coordinateEnd - coordinateBegin;
    partId_t noCuts = partNo - 1;
    //size_t total_part_count = noCuts + partNo;
    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon();

    //if (migration_check == true) allowNonRectelinearPart = false;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
    {
        int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
        me = omp_get_thread_num();
#endif

        pq_lno_t *myPartPointCounts = partPointCounts[me];
        float *myRatios = NULL;
        if (allowNonRectelinearPart){


            myRatios = nonRectelinearRatios[me];
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
            for (partId_t i = 0; i < noCuts; ++i){
                //float r = actual_ratios[i];
                //pq_scalar_t leftWeight = r * (localPartWeights[i * 2 + 1] - localPartWeights[i * 2]);
                pq_scalar_t leftWeight = actual_ratios[i];
                for(int ii = 0; ii < noThreads; ++ii){
                    if(leftWeight > _EPSILON){

                        pq_scalar_t ithWeight = partWeights[ii][i * 2 + 1] - partWeights[ii][i * 2 ];
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

/*
            for (partId_t i = 0; i < noCuts; ++i){
                cout << "myR:" << i << " is:" <<  myRatios[i] << endl;
            }
*/
            if(noCuts > 0){
                for (partId_t i = noCuts - 1; i > 0 ; --i){
                    if(ABS(cutCoordinates[i] - cutCoordinates[i -1]) < _EPSILON){
                        myRatios[i] -= myRatios[i - 1] ;
                    }
                    myRatios[i] = int ((myRatios[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL) / pq_scalar_t(SIGNIFICANCE_MUL);
                }
            }
/*
            for (partId_t i = 0; i < noCuts; ++i){
                cout << "myR:" << i << " is:" <<  myRatios[i] << endl;
            }
*/


        }

        for(partId_t ii = 0; ii < partNo; ++ii){
            myPartPointCounts[ii] = 0;
        }

        if(!allowNonRectelinearPart || !isSequentialAndPricise) {

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
            for (pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

                pq_lno_t i = partitionedPointPermutations[ii];
                pq_scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];
                partId_t pp = partIds[i];
                partId_t p = pp / 2;
                if(pp % 2 == 1){
                    if(allowNonRectelinearPart && myRatios[p] > _EPSILON * EPS_SCALE){
                        //cout << "pa:" << p << endl;
                        myRatios[p] -= w;
                        if(myRatios[p] < 0 && p < noCuts - 1 && ABS(cutCoordinates[p+1] - cutCoordinates[p]) < _EPSILON){
                            myRatios[p + 1] += myRatios[p];
                        }
                        ++myPartPointCounts[p];
                        partIds[i] = p;
                    }
                    else{
                        ++p;
                        //this while loop is necessary when a line is partitioned into
                        //more than 2 parts.

                        while(allowNonRectelinearPart && p < noCuts){
                        //traverse all the cut lines having the same partitiong
                            if(ABS(cutCoordinates[p] - cutCoordinates[p - 1])
                                 < _EPSILON){
                                //if line has enough space on left, put it there.
                                if(myRatios[p] > _EPSILON * EPS_SCALE &&
                                     myRatios[p] >= ABS(myRatios[p] - w)){
                                    //pq_scalar_t w = pqJagged_uniformWeights?
                                    //1:coordWeights[i];
                                    myRatios[p] -= w;
                                    if(myRatios[p] < 0 && p < noCuts - 1 &&
                                         ABS(cutCoordinates[p+1] -
                                         cutCoordinates[p]) < _EPSILON){
                                        myRatios[p + 1] += myRatios[p];
                                    }
                                    break;
                                }
                            }
                            else {
                                //if cut coordinates are different, put it to
                                //next part.
                                break;
                            }
                            ++p;
                        }

                        //pq_scalar_t currentCut = cutCoordinates[p];
                        //TODO:currently cannot divide 1 line more than 2 parts.
                        //bug cannot be divided, therefore this part should change.
                        //cout << "p:" << p+1 << endl;
                        //++myPartPointCounts[p + 1];
                        //partIds[i] = p + 1;
                        ++myPartPointCounts[p];
                        partIds[i] = p;
                    }
                }
                else {
                    ++myPartPointCounts[p];
                    partIds[i] = p;
                }
            }
/*
            for (partId_t p = 0; p < noCuts; ++p){

                cout << "p:" << p <<
                        " myPartPointCounts[p]:" << myPartPointCounts[p] <<
                        " \t\tcutCoorD:" << cutCoordinates[p] <<
                        " \tmyRatios[i]:" << myRatios[p] <<
                        " \tactual_ratios:" << actual_ratios[p] <<
                        endl;

            }
            cout << "p:" << noCuts << " myPartPointCounts[p]:" << myPartPointCounts[noCuts] << endl;
*/
        } else {
            partId_t *cutMap = allocMemory<partId_t> (noCuts);

            typedef uMultiSortItem<pq_lno_t, int, pq_scalar_t> multiSItem;
            typedef std::vector< multiSItem > multiSVector;
            typedef std::vector<multiSVector> multiS2Vector;

            std::vector<pq_scalar_t *>allocated_memory;
            //pq_scalar_t *currentPartWeight = allocMemory<pq_scalar_t> (partNo);
            //memset(currentPartWeight, 0, sizeof(pq_scalar_t) * partNo);

            multiS2Vector cutPointSortArrays;

            partId_t differentCutCount = 1;
            cutMap[0] = 0;

            multiSVector tmpMultiSVector;
            cutPointSortArrays.push_back(tmpMultiSVector);

            for (partId_t i = 1; i < noCuts ; ++i){
                if(ABS(cutCoordinates[i] - cutCoordinates[i -1]) < _EPSILON){
                    cutMap[i] = cutMap[i-1];
                }
                else {
                    cutMap[i] = differentCutCount++;
                    multiSVector tmp2MultiSVector;
                    cutPointSortArrays.push_back(tmp2MultiSVector);
                }
            }


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
            for (pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){

                pq_lno_t i = partitionedPointPermutations[ii];
                //pq_scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];

                //cout << "index:" << i << " pqJagged_coordinates:" << pqJagged_coordinates[i] << " p:" << partIds[i] << endl;
                partId_t pp = partIds[i];
                partId_t p = pp / 2;
                if(pp % 2 == 1 /*&& ABS(actual_ratios[p] - 1.0f) < _EPSILON*/){
                    pq_scalar_t *vals = allocMemory<pq_scalar_t>(pqCoordDim -1);
                    allocated_memory.push_back(vals);

                    int val_ind = 0;
                    for(int dim = currentCoord + 1; dim < pqCoordDim; ++dim){
                        vals[val_ind++] = pqJaggedAllCoordinates[dim][i];
                    }
                    for(int dim = 0; dim < currentCoord; ++dim){
                        vals[val_ind++] = pqJaggedAllCoordinates[dim][i];
                    }
                    multiSItem tempSortItem(i, pqCoordDim -1, vals);

                    partId_t cmap = cutMap[p];
                    cutPointSortArrays[cmap].push_back(tempSortItem);

                }
                else {
                    ++myPartPointCounts[p];
                    partIds[i] = p;
                    //currentPartWeight[p] += w;
                }
            }

            for (partId_t i = 0; i < differentCutCount; ++i){
                std::sort (cutPointSortArrays[i].begin(), cutPointSortArrays[i].end());
            }


            /*
            cout << "differentCutCount:" << differentCutCount << endl;

            for (partId_t p = 0; p < noCuts; ++p){

                partId_t mappedCut = cutMap[p];
                pq_lno_t cutPointCount = (pq_lno_t)cutPointSortArrays[mappedCut].size();
                cout << "p:" << p <<
                        " mappedCut:" << mappedCut <<
                        " myPartPointCounts[p]:" << myPartPointCounts[p] <<
                        " \t\tcutCoorD:" << cutCoordinates[p] <<
                        " \tcutPointCount[p]:" << cutPointCount <<
                        " \tmyRatios[i]:" << myRatios[p] <<
                        " \tactual_ratios:" << actual_ratios[p] <<
                        endl;

            }
            cout << "p:" << noCuts << " myPartPointCounts[p]:" << myPartPointCounts[noCuts] << endl;
            */
            partId_t prevMap = cutMap[0];
            pq_scalar_t leftOver = 0;
            for (partId_t p = 0; p < noCuts; ++p){
                /*
                {
                    partId_t mappedCut = cutMap[p];
                    pq_lno_t cutPointCount = (pq_lno_t)cutPointSortArrays[mappedCut].size();

                    cout << "p:" << p <<
                            " mappedCut:" << mappedCut <<
                            " myPartPointCounts[p]:" << myPartPointCounts[p] <<
                            " \t\tcutCoorD:" << cutCoordinates[p] <<
                            " \tcutPointCount[p]:" << cutPointCount <<
                            " \tmyRatios[i]:" << myRatios[p] <<
                            " \tactual_ratios:" << actual_ratios[p] <<
                            endl;

                }
                */
                //cout << "p:" << p << " leftOVer:" << leftOver << endl;
                partId_t mappedCut = cutMap[p];
                if (prevMap != mappedCut){
                    //biggestCutForVector[mappedCut] = p;
                    pq_lno_t vEnd = (pq_lno_t)cutPointSortArrays[prevMap].size() - 1;
                    for (; vEnd >= 0; --vEnd){
                        multiSItem t = cutPointSortArrays[prevMap][vEnd];
                        pq_lno_t i = t.index;
                        //pq_scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];
                        //cutPointSortArrays[mappedCut].pop_back();
                        //cout << "putting i : " << i << " on cut:" << prevMap << " to part p:" << p << endl;
                        ++myPartPointCounts[p];
                        partIds[i] = p;
                    }
                    cutPointSortArrays[prevMap].clear();
                }
                //biggestCutForVector[mappedCut] = p;
                pq_lno_t vEnd = (pq_lno_t)cutPointSortArrays[mappedCut].size() - 1;

                //cout << "c:" << p << " r:" << myRatios[p] << " leftOver:" << leftOver << endl;
                //pq_scalar_t leftPartW = currentPartWeight[p];

                for (; vEnd >= 0; --vEnd){
                    multiSItem t = cutPointSortArrays[mappedCut][vEnd];
                    pq_lno_t i = t.index;
                    pq_scalar_t w = pqJagged_uniformWeights? 1:coordWeights[i];

                    if(myRatios[p] + leftOver> _EPSILON * EPS_SCALE &&
                        myRatios[p] + leftOver - ABS(myRatios[p] + leftOver - w)
                         > _EPSILON){

                        myRatios[p] -= w;
                        cutPointSortArrays[mappedCut].pop_back();
                        ++myPartPointCounts[p];
                        //cout << "putting i : " << i << " on cut:" << p <<
                        //" to part p:" << p << endl;
                        partIds[i] = p;
                        if(p < noCuts - 1 && myRatios[p] < _EPSILON){
                            if(mappedCut == cutMap[p + 1] ){
                                if (prevMap != mappedCut){
                                    leftOver = myRatios[p];
                                }
                                else {
                                    leftOver += myRatios[p];
                                }
                            }
                            else{
                                leftOver = -myRatios[p];
                            }
                            break;
                        }
                    } else {
                        if(p < noCuts - 1 && mappedCut == cutMap[p + 1]){
                            if (prevMap != mappedCut){
                                leftOver = myRatios[p];
                            }
                            else {
                                leftOver += myRatios[p];
                            }
                        }
                        else{
                            leftOver = -myRatios[p];
                        }
                        break;
                    }
                }
                prevMap = mappedCut;
            }
            {
                pq_lno_t vEnd = (pq_lno_t)cutPointSortArrays[prevMap].size() - 1;
                for (; vEnd >= 0; --vEnd){
                    multiSItem t = cutPointSortArrays[prevMap][vEnd];
                    pq_lno_t i = t.index;
                    //pq_scalar_t w = pqJagged_uniformWeights? 1:
                    //coordWeights[i];

                    //cutPointSortArrays[mappedCut].pop_back();
                    //cout << "putting i : " << i << " on cut:" << prevMap <<
                    //" to part p:" << noCuts << endl;
                    ++myPartPointCounts[noCuts];
                    partIds[i] = noCuts;
                }
                cutPointSortArrays[prevMap].clear();
            }

            freeArray<partId_t> (cutMap);
            //freeArray<pq_scalar_t> (currentPartWeight);

            pq_lno_t vSize = (pq_lno_t) allocated_memory.size();
            for(pq_lno_t i = 0; i < vSize; ++i){
                freeArray<pq_scalar_t> (allocated_memory[i]);
            }
        }


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
        for(partId_t j = 0; j < partNo; ++j){
            pq_lno_t pwj = 0;
            for (int i = 0; i < noThreads; ++i){
                pq_lno_t threadPartPointCount = partPointCounts[i][j];
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
        for (pq_lno_t ii = coordinateBegin; ii < coordinateEnd; ++ii){
            pq_lno_t i = partitionedPointPermutations[ii];
            partId_t p =  partIds[i];
            newpartitionedPointPermutations[coordinateBegin +
                                     myPartPointCounts[p]++] = i;
        }
    }
}




template <typename pq_scalar_t, typename pq_lno_t, typename partId_t>
void sequentialTaskPartitioning(
    const RCP<const Environment> &env,
    pq_lno_t numLocalCoords,
    pq_lno_t actualNumCoords,
    size_t numGlobalParts,
    int coordDim,
    pq_scalar_t **pqJagged_coordinates,
    pq_lno_t *output_permutation,
    pq_lno_t *output_partIndices,
    int partArraySize,
    const partId_t *partNo
){

    //env->timerStart(MACRO_TIMERS, "PQJagged - " +partitioningName+"-Problem_Partitioning");
#ifdef HAVE_ZOLTAN2_OMP
    int actual_num_threads = omp_get_num_threads();
    omp_set_num_threads(1);
#endif
    const RCP<Comm<int> > commN;
    RCP<Comm<int> >comm =  Teuchos::rcp_const_cast<Comm<int> >
            (Teuchos::DefaultComm<int>::getDefaultSerialComm(commN));

    /*
    const partId_t *partNo = NULL;
    int partArraySize = 0;
    */

    //weights are uniform for task mapping
    bool pqJagged_uniformWeights[1];
    pqJagged_uniformWeights[0] = true;

    //parts are uniform for task mapping
    bool pqJagged_uniformParts[1];
    pqJagged_uniformParts[0] = true;
    pq_scalar_t *pqJagged_partSizes[1];
    pqJagged_partSizes[0] = NULL;
    //no weights
    //int weightDim = 0;
    pq_scalar_t pqJagged_weights[1][1];

    //as input indices.
    pq_lno_t *partitionedPointCoordinates =  allocMemory< pq_lno_t>(numLocalCoords);
    //initial configuration
    //set each pointer-i to i.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < static_cast<size_t>(numLocalCoords); ++i){
        partitionedPointCoordinates[i] = output_permutation[i];
    }

    //as output indices
    pq_lno_t *newpartitionedPointCoordinates = allocMemory< pq_lno_t>(numLocalCoords);
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<pq_lno_t>(newpartitionedPointCoordinates, numLocalCoords);
#endif
#endif

    partId_t *partIds = NULL;
    ArrayRCP<partId_t> partId;
    if(numLocalCoords > 0){
        partIds = allocMemory<partId_t>(numLocalCoords);
    }

    //initially there is a single partition
    partId_t currentPartitionCount = 1;
    //single partition starts at index-0, and ends at numLocalCoords
    //inTotalCounts array holds the end points in partitionedPointCoordinates array
    //for each partition. Initially sized 1, and single element is set to numLocalCoords.
    pq_lno_t *inTotalCounts = allocMemory<pq_lno_t>(1);
    inTotalCounts[0] = static_cast<pq_lno_t>(actualNumCoords);//the end of the initial partition is the end of coordinates.
    //the ends points of the output.
    pq_lno_t *outTotalCounts = NULL;

    //get pqJagged specific parameters.
    bool allowNonRectelinearPart = true;
    int concurrentPartCount = 1; // Set to invalid value

    int numThreads = 1;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(numThreads)
    {
        numThreads = omp_get_num_threads();
    }

#endif

    partId_t totalDimensionCut = 0; //how many cuts will be totally
    partId_t totalPartCount = 1;    //how many parts will be totally
    partId_t maxPartNo = 0;         //maximum cut count along a dimension.
    partId_t reduceAllCount = 0;    //estimate of #reduceAlls that can be done.
    partId_t maxTotalCumulativePartCount = 1; //max no of parts that might occur
                                              //during partitioning before the
                                              //last partitioning dimension.
    partId_t maxCutNo = 0;
    size_t maxTotalPartCount = 0;

    //partArraySize = coordDim * 8;
    getPartSpecifications <partId_t>( partNo, partArraySize, numGlobalParts,
        totalDimensionCut, totalPartCount, maxPartNo, reduceAllCount,
        maxTotalCumulativePartCount, maxCutNo, maxTotalPartCount);

    // coordinates of the cut lines. First one is the min, last one is max coordinate.
    // kddnote if (keep_cuts)
    // coordinates of the cut lines.
    //only store this much if cuts are needed to be stored.
    pq_scalar_t *allCutCoordinates = allocMemory< pq_scalar_t>(totalDimensionCut);
    // kddnote else
    //pq_scalar_t *allCutCoordinates = allocMemory< pq_scalar_t>(maxCutNo * concurrentPartCount);
    pq_scalar_t *max_min_array =  allocMemory< pq_scalar_t>(numThreads * 2);

    float *nonRectelinearPart = NULL; //how much weight percentage should a MPI put left side of the each cutline
    float **nonRectRatios = NULL; //how much weight percentage should each thread in MPI put left side of the each cutline
    if(allowNonRectelinearPart){
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
    pq_scalar_t *cutCoordinatesWork = allocMemory<pq_scalar_t>(maxCutNo * concurrentPartCount);

#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<pq_scalar_t>(cutCoordinatesWork, maxCutNo);
#endif
#endif

    //cumulative part weight ratio array.
    pq_scalar_t *targetPartWeightRatios = allocMemory<pq_scalar_t>(maxPartNo * concurrentPartCount); // the weight ratios at left side of the cuts. First is 0, last is 1.
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<pq_scalar_t>(cutPartRatios, maxCutNo);
#endif
#endif


    pq_scalar_t *cutUpperBounds = allocMemory<pq_scalar_t>(maxCutNo * concurrentPartCount);  //upper bound coordinate of a cut line
    pq_scalar_t *cutLowerBounds = allocMemory<pq_scalar_t>(maxCutNo* concurrentPartCount);  //lower bound coordinate of a cut line
    pq_scalar_t *cutLowerWeight = allocMemory<pq_scalar_t>(maxCutNo* concurrentPartCount);  //lower bound weight of a cut line
    pq_scalar_t *cutUpperWeight = allocMemory<pq_scalar_t>(maxCutNo* concurrentPartCount);  //upper bound weight of a cut line

    pq_scalar_t *localMinMaxTotal = allocMemory<pq_scalar_t>(3 * concurrentPartCount); //combined array to exchange the min and max coordinate, and total weight of part.
    pq_scalar_t *globalMinMaxTotal = allocMemory<pq_scalar_t>(3 * concurrentPartCount);//global combined array with the results for min, max and total weight.

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
    pq_scalar_t **leftClosestDistance = allocMemory<pq_scalar_t *>(numThreads);
    //leftClosesDistance to hold the min distance of a coordinate to a cutline from right (for each thread)
    pq_scalar_t **rightClosestDistance = allocMemory<pq_scalar_t *>(numThreads);

    //to store how many points in each part a thread has.
    pq_lno_t **partPointCounts = allocMemory<pq_lno_t *>(numThreads);

    for(int i = 0; i < numThreads; ++i){
        //partWeights[i] = allocMemory<pq_scalar_t>(maxTotalPartCount);
        partWeights[i] = allocMemory < double >(maxTotalPartCount * concurrentPartCount);
        rightClosestDistance[i] = allocMemory<pq_scalar_t>(maxCutNo * concurrentPartCount);
        leftClosestDistance[i] = allocMemory<pq_scalar_t>(maxCutNo * concurrentPartCount);
        partPointCounts[i] =  allocMemory<pq_lno_t>(maxPartNo);
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
    pq_scalar_t *cutWeights = allocMemory<pq_scalar_t>(maxCutNo);
    pq_scalar_t *globalCutWeights = allocMemory<pq_scalar_t>(maxCutNo);

    //for faster communication, concatanation of
    //totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
    //leftClosest distances sized P-1, since P-1 cut lines
    //rightClosest distances size P-1, since P-1 cut lines.
    pq_scalar_t *totalPartWeights_leftClosests_rightClosests = allocMemory<pq_scalar_t>((maxTotalPartCount + maxCutNo * 2) * concurrentPartCount);
    pq_scalar_t *global_totalPartWeights_leftClosests_rightClosests = allocMemory<pq_scalar_t>((maxTotalPartCount + maxCutNo * 2) * concurrentPartCount);


    pq_scalar_t *cutCoordinates =  allCutCoordinates;

    //partId_t leftPartitions = totalPartCount;
    pq_scalar_t maxScalar_t = numeric_limits<pq_scalar_t>::max();
    pq_scalar_t minScalar_t = -numeric_limits<pq_scalar_t>::max();



    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon();
    //partId_t partIndexBegin = 0;
    partId_t futurePartNumbers = totalPartCount;

    vector<partId_t> *currentPartitions = new vector<partId_t> ();
    vector<partId_t> *newFuturePartitions = new vector<partId_t> ();
    newFuturePartitions->push_back(numGlobalParts);
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > t1;
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > t2;


    for (int i = 0; i < partArraySize; ++i){

        //partitioning array.
        //size will be as the number of current partitions
        //and this hold how many parts each part will be
        //in the current dimension partitioning.
        vector <partId_t> pAlongI;

        //number of parts that will be obtained at the end of this partitioning.
        //currentPartitions is as the size of current number of parts.
        //holds how many more parts each should be divided in the further
        //iterations.
        //this will be used to calculate pAlongI,
        //as the number of parts that the part will be partitioned
        //in the current dimension partitioning.

        //newFuturePartitions will be as the size of outnumParts,
        //and this will hold how many more parts that each output part
        //should be divided.
        //this array will also be used to determine the weight ratios
        //of the parts.
        //swap the arrays.
        vector<partId_t> *tmpPartVect= currentPartitions;
        currentPartitions = newFuturePartitions;
        newFuturePartitions = tmpPartVect;

        //clear newFuturePartitions array as
        //getPartitionArrays expects it to be empty.
        //it also expects pAlongI to be empty as well.
        newFuturePartitions->clear();

        /*
        cout << "i:" << i << " ";
        for (int jj = 0; jj < currentPartitions->size();++jj){
            cout << (*currentPartitions)[jj] << " ";
        }
        cout << endl;
        */
        //returns the total number of output parts for this dimension partitioning.
        partId_t outPartCount = getPartitionArrays<pq_scalar_t, partId_t>(
                partNo, pAlongI, currentPartitions, newFuturePartitions,
                futurePartNumbers, currentPartitionCount, partArraySize, i,
                maxPartNo, 0, t1, t2);

        /*
        cout << "i:" << i << " ";
        for (int jj = 0; jj < pAlongI.size();++jj){
            cout << pAlongI[jj] << " ";
        }
        cout << endl;
        */

        if(outPartCount == currentPartitionCount) {
            tmpPartVect= currentPartitions;
            currentPartitions = newFuturePartitions;
            newFuturePartitions = tmpPartVect;
            continue;
        }

        //get the coordinate axis along which the partitioning will be done.
        int coordInd = i % coordDim;
        pq_scalar_t * pqCoord = pqJagged_coordinates[coordInd];
        //convert i to string to be used for debugging purposes.

        string istring = toString<int>(i);
        //env->timerStart(MACRO_TIMERS, "PQJagged - " +partitioningName+"
        //-Problem_Partitioning_" + istring);

        //alloc Memory to point the indices
        //of the parts in the permutation array.
        outTotalCounts = allocMemory<pq_lno_t>(outPartCount);

        //cout << "outPart:" << outPartCount << " global:" << numGlobalParts
        //<< endl;

        //the index where in the outtotalCounts will be written.
        partId_t currentOut = 0;
        //whatever is written to outTotalCounts will be added with previousEnd
        //so that the points will be shifted.
        partId_t previousEnd = 0;

        partId_t currentWorkPart = 0;
        partId_t concurrentPart = min(currentPartitionCount - currentWorkPart,
                                         concurrentPartCount);

        //always use binary search algorithm.
        bool useBinarySearch = true;
        partId_t obtainedPartCount = 0;

        //run for all available parts.
        for (; currentWorkPart < currentPartitionCount;
                     currentWorkPart += concurrentPart){

            concurrentPart = min(currentPartitionCount - currentWorkPart,
            concurrentPartCount);

            partId_t workPartCount = 0;
            //get the min and max coordinates of each part
            //together with the part weights of each part.
            for(int kk = 0; kk < concurrentPart; ++kk){
                partId_t currentPart = currentWorkPart + kk;

                //if this part wont be partitioned any further
                //dont do any work for this part.
                if (pAlongI[currentPart] == 1){
                    continue;
                }
                ++workPartCount;
                pq_lno_t coordinateEnd= inTotalCounts[currentPart];
                pq_lno_t coordinateBegin = currentPart==0 ? 0: inTotalCounts[currentPart -1];
                //cout << "me:" << problemComm->getRank() << " begin:" <<
                //coordinateBegin  << " end:" << coordinateEnd << endl;
                pqJagged_getLocalMinMaxTotalCoord<pq_scalar_t, pq_lno_t>(
                    partitionedPointCoordinates, pqCoord,
                    pqJagged_uniformWeights[0], pqJagged_weights[0], numThreads,
                    coordinateBegin, coordinateEnd, max_min_array, maxScalar_t,
                    minScalar_t,
                        localMinMaxTotal[kk], //min coordinate
                        localMinMaxTotal[kk + concurrentPart], //max coordinate
                        localMinMaxTotal[kk + 2*concurrentPart] //total weight);
                );


            }


            if (workPartCount > 0){
                //obtain global Min max of the part.
                pqJagged_getGlobalMinMaxTotalCoord<pq_scalar_t>( comm, env,
                concurrentPart, localMinMaxTotal, globalMinMaxTotal);

                //represents the total number of cutlines
                //whose coordinate should be determined.
                partId_t allDone = 0;

                //Compute weight ratios for parts & cuts:
                //e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
                //part0  cut0  part1 cut1 part2 cut2 part3
                partId_t cutShifts = 0;
                partId_t partShift = 0;
                for(int kk = 0; kk < concurrentPart; ++kk){
                    pq_scalar_t minCoordinate = globalMinMaxTotal[kk];
                    pq_scalar_t maxCoordinate = globalMinMaxTotal[kk +
                                                     concurrentPart];
                    pq_scalar_t globalTotalWeight = globalMinMaxTotal[kk +
                                                     2 * concurrentPart];

                    partId_t currentPart = currentWorkPart + kk;

                    partId_t partition = pAlongI[currentPart];

                    pq_scalar_t *usedCutCoordinate = cutCoordinates + cutShifts;
                    pq_scalar_t *usedCutPartRatios = targetPartWeightRatios +
                                                                     partShift;
                    //shift the usedCutCoordinate array as noCuts.
                    cutShifts += partition - 1;
                    //shift the partRatio array as noParts.
                    partShift += partition;

                    //cout << "min:" << minCoordinate << " max:" <<
                    //maxCoordinate << endl;
 
                    //calculate only if part is not empty,
                    //and part will be further partitioend.
                    if(partition > 1 && minCoordinate <= maxCoordinate){

                        //increase allDone by the number of cuts of the current
                        //part's cut line number.
                        allDone += partition - 1;
                        //set the number of cut lines that should be determined
                        //for this part.
                        myNonDoneCount[kk] = partition - 1;

                        //get the target weights of the parts.
                        pqJagged_getCutCoord_Weights<pq_scalar_t>(
                            minCoordinate, maxCoordinate,
                            pqJagged_uniformParts[0],
                            pqJagged_uniformWeights[0], pqJagged_partSizes[0],
                            partition - 1, globalTotalWeight, usedCutCoordinate,
                            usedCutPartRatios, numThreads, currentPartitions,
                            newFuturePartitions, currentPart,
                            obtainedPartCount);

                        //get the initial estimated part assignments of the coordinates.
                        getInitialPartAssignments<pq_scalar_t, pq_lno_t, partId_t>(
                            maxCoordinate, minCoordinate, currentPart,
                            inTotalCounts, partitionedPointCoordinates, pqCoord,
                            partIds, _EPSILON, partition);
                    }
                    else {
                        // e.g., if have fewer coordinates than parts, don't need to do next dim.
                        myNonDoneCount[kk] = 0;
                    }
                    obtainedPartCount += partition;
                }

                //used imbalance, it is always 0, as it is difficult to estimate a range.
                pq_scalar_t used_imbalance = 0;

                // Determine cut lines for k parts here.
                pqJagged_1D_Partition<pq_scalar_t, pq_lno_t>(
                    env, comm, partitionedPointCoordinates, pqCoord,
                    pqJagged_uniformWeights[0], pqJagged_weights[0],
                    targetPartWeightRatios, globalMinMaxTotal, localMinMaxTotal,
                    numThreads, used_imbalance, currentWorkPart, concurrentPart,
                    inTotalCounts, cutCoordinates, cutCoordinatesWork,
                    leftClosestDistance, rightClosestDistance, cutUpperBounds,
                    cutLowerBounds, cutUpperWeight, cutLowerWeight, isDone,
                    partWeights, totalPartWeights_leftClosests_rightClosests,
                    global_totalPartWeights_leftClosests_rightClosests,
                    allowNonRectelinearPart, nonRectelinearPart, cutWeights,
                    globalCutWeights, allDone, myNonDoneCount, useBinarySearch,
                    partIds, pAlongI);
            }

            //create part chunks
            {

                partId_t outShift = 0;
                partId_t cutShift = 0;
                size_t tlrShift = 0;
                size_t pwShift = 0;

                for(int kk = 0; kk < concurrentPart; ++kk){
                    partId_t curr = currentWorkPart + kk;
                    partId_t noParts = pAlongI[curr];

                    //if the part is empty, skip the part.
                    if((noParts != 1  ) && globalMinMaxTotal[kk] >
                             globalMinMaxTotal[kk + concurrentPart]) {

                        for(partId_t jj = 0; jj < noParts; ++jj){
                            outTotalCounts[currentOut + outShift + jj] = 0;
                        }
                        cutShift += noParts - 1;
                        tlrShift += (4 *(noParts - 1) + 1);
                        outShift += noParts;
                        pwShift += (2 * (noParts - 1) + 1);
                        continue;
                    }

                    pq_lno_t coordinateEnd= inTotalCounts[curr];
                    pq_lno_t coordinateBegin = curr==0 ? 0: inTotalCounts[curr
                                                             -1];
                    pq_scalar_t *usedCutCoordinate = cutCoordinates + cutShift;
                    float *usednonRectelinearPart = nonRectelinearPart +
                                                         cutShift;

                    pq_scalar_t *tlr =  totalPartWeights_leftClosests_rightClosests + tlrShift;

                    for(int ii = 0; ii < numThreads; ++ii){
                        pws[ii] = partWeights[ii] +  pwShift;
                    }

                    if(noParts > 1){
                        // Rewrite the indices based on the computed cuts.
                        getChunksFromCoordinates<pq_lno_t,pq_scalar_t>(
                            noParts, numThreads, partitionedPointCoordinates,
                            pqCoord, pqJagged_uniformWeights[0],
                            pqJagged_weights[0], usedCutCoordinate,
                            coordinateBegin, coordinateEnd, 
                            allowNonRectelinearPart, usednonRectelinearPart,
                            tlr, pws, nonRectRatios,
                            partPointCounts, newpartitionedPointCoordinates,
                            outTotalCounts + currentOut + outShift, partIds,
                            true, pqJagged_coordinates, coordDim, coordInd );
                    }
                    else {
                        //if this part is partitioned into 1 then just copy
                        //the old values.
                        pq_lno_t partSize = coordinateEnd - coordinateBegin;
                        *(outTotalCounts + currentOut + outShift) = partSize;
                        memcpy(newpartitionedPointCoordinates + coordinateBegin,
                        partitionedPointCoordinates + coordinateBegin,
                        partSize * sizeof(pq_lno_t));
                    }
                    cutShift += noParts - 1;
                    tlrShift += (4 *(noParts - 1) + 1);
                    outShift += noParts;
                    pwShift += (2 * (noParts - 1) + 1);
                }

                //shift cut coordinates so that all cut coordinates are stored.
                cutCoordinates += cutShift;

                //getChunks from coordinates partitioned the parts and
                //wrote the indices as if there were a single part.
                //now we need to shift the beginning indices.
                for(partId_t kk = 0; kk < concurrentPart; ++kk){
                    partId_t noParts = pAlongI[ currentWorkPart + kk];
                    for (partId_t ii = 0;ii < noParts ; ++ii){
                        //shift it by previousCount
                        outTotalCounts[currentOut+ii] += previousEnd;
                    }
                    //increase the previous count by current end.
                    previousEnd = outTotalCounts[currentOut + noParts - 1];
                    //increase the current out.
                    currentOut += noParts ;
                }
            }
        }
        // end of this partitioning dimension

        currentPartitionCount = outPartCount;
        pq_lno_t * tmp = partitionedPointCoordinates;
        partitionedPointCoordinates = newpartitionedPointCoordinates;
        newpartitionedPointCoordinates = tmp;
        freeArray<pq_lno_t>(inTotalCounts);
        inTotalCounts = outTotalCounts;

        //env->timerStop(MACRO_TIMERS, "PQJagged - " +partitioningName+"-Problem_Partitioning_" + istring);
    }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(pq_lno_t i = 0; i < numLocalCoords; ++i){
        output_permutation[i] = partitionedPointCoordinates[i];
    }

    for(size_t i = 0; i < numGlobalParts ; ++i){
        //cout << "i:" << i << endl;
        output_partIndices[i] = inTotalCounts[i];
    }


    delete currentPartitions;
    delete newFuturePartitions;

    freeArray<partId_t>(partIds);
    freeArray<pq_lno_t>(partitionedPointCoordinates);
    freeArray<pq_lno_t>(newpartitionedPointCoordinates);
    freeArray<pq_lno_t>(inTotalCounts);

    freeArray<pq_scalar_t>(allCutCoordinates);
    freeArray<pq_scalar_t>(cutCoordinatesWork);
    freeArray<pq_scalar_t>(max_min_array);

    if(allowNonRectelinearPart){
        freeArray<float>(nonRectelinearPart);
        for(int i = 0; i < numThreads; ++i){
            freeArray<float>(nonRectRatios[i]);
        }
        freeArray<float *>(nonRectRatios);
    }
    freeArray<pq_scalar_t>(targetPartWeightRatios);

    freeArray<pq_scalar_t>(cutUpperBounds);
    freeArray<pq_scalar_t>(cutLowerBounds);
    freeArray<pq_scalar_t>(cutLowerWeight);
    freeArray<pq_scalar_t>(cutUpperWeight);
    freeArray<pq_scalar_t> (localMinMaxTotal);
    freeArray<pq_scalar_t> (globalMinMaxTotal);
    freeArray<bool>(isDone);
    freeArray<partId_t>(myNonDoneCount);

    for(int i = 0; i < numThreads; ++i){
        freeArray<double>(partWeights[i]);
        freeArray<pq_scalar_t>(rightClosestDistance[i]);
        freeArray<pq_scalar_t>(leftClosestDistance[i]);
    }
    freeArray<double *>(partWeights);
    freeArray<pq_scalar_t *>(leftClosestDistance);
    freeArray<pq_scalar_t *>(rightClosestDistance);

    freeArray<double *> (pws);

    for(int i = 0; i < numThreads; ++i){
        freeArray<pq_lno_t>(partPointCounts[i]);
    }
    freeArray<pq_lno_t *>(partPointCounts);
    freeArray<pq_scalar_t>(cutWeights);

    freeArray<pq_scalar_t>(globalCutWeights);
    freeArray<pq_scalar_t>(totalPartWeights_leftClosests_rightClosests);
    freeArray<pq_scalar_t>(global_totalPartWeights_leftClosests_rightClosests);
    //env->timerStop(MACRO_TIMERS, "PQJagged - " +partitioningName+"-Problem_Partitioning");

#ifdef HAVE_ZOLTAN2_OMP
    omp_set_num_threads(actual_num_threads);
#endif
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
    RCP<Comm<int> > &problemComm,
    const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &coords,
    RCP<PartitioningSolution<Adapter> > &solution
)
{
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL

    Z2_THROW_EXPERIMENTAL("Zoltan2 PQJagged is strictly experimental software "
            "while it is being developed and tested.")

#else

    env->timerStart(MACRO_TIMERS, "PQJagged - Total");
    env->timerStart(MACRO_TIMERS, "PQJagged - Total2");

    typedef typename Adapter::scalar_t pq_scalar_t;
    typedef typename Adapter::gno_t pq_gno_t;
    typedef typename Adapter::lno_t pq_lno_t;
    typedef typename Adapter::node_t pq_node_t;

    env->debug(3, "In PQ Jagged");

    /*
    if(comm->getRank() == 0){
    cout << "size of gno:" << sizeof(pq_gno_t) << endl;
    cout << "size of lno:" << sizeof(pq_lno_t) << endl;
    cout << "size of pq_scalar_t:" << sizeof(pq_scalar_t) << endl;
    }
     */
    const Teuchos::ParameterList &pl = env->getParameters();

    // TODO: This is not used
    std::bitset<NUM_RCB_PARAMS> params;
    int numTestCuts = 5;

    pq_scalar_t imbalanceTolerance;

    multiCriteriaNorm mcnorm;
    bool ignoreWeights=false;

    get_partitioning_params<pq_scalar_t>(pl, imbalanceTolerance, mcnorm, params,
            numTestCuts, ignoreWeights);

    const partId_t *partNo = NULL;
    int partArraySize = 0;

    if (pl.getPtr<Array <partId_t> >("pqParts")){
        partNo = pl.getPtr<Array <partId_t> >("pqParts")->getRawPtr();
        partArraySize = pl.getPtr<Array <partId_t> >("pqParts")->size() - 1;
        env->debug(2, "PQparts provided by user");
    }

    //cout << "partArraySize:" << partArraySize << endl;
    int coordDim, weightDim;
    size_t nlc;
    global_size_t gnc; int criteriaDim;
    pqJagged_getCoordinateValues<Adapter>( coords, coordDim, weightDim, nlc,
            gnc, criteriaDim, ignoreWeights);

    pq_lno_t numLocalCoords = nlc;
#ifdef enable_migration2
    pq_gno_t numGlobalCoords = gnc;
#endif

    //allocate only two dimensional pointer.
    //raw pointer addresess will be obtained from multivector.
    pq_scalar_t **pqJagged_coordinates = allocMemory<pq_scalar_t *>(coordDim);
    pq_scalar_t **pqJagged_weights = allocMemory<pq_scalar_t *>(criteriaDim);
     //if the partitioning results are to be uniform.
    bool *pqJagged_uniformParts = allocMemory< bool >(criteriaDim);

    //if in a criteria dimension, uniform part is false this shows ratios of
    //the target part weights.
    pq_scalar_t **pqJagged_partSizes =  allocMemory<pq_scalar_t *>(criteriaDim);
    //if the weights of coordinates are uniform in a criteria dimension.
    bool *pqJagged_uniformWeights = allocMemory< bool >(criteriaDim);

    ArrayView<const pq_gno_t> pqJagged_gnos;
    size_t numGlobalParts;
    int pqJagged_multiVectorDim;

    pqJagged_getInputValues<Adapter, pq_scalar_t, pq_gno_t>(
            env, coords, solution, params, coordDim, weightDim, numLocalCoords,
            numGlobalParts, pqJagged_multiVectorDim, pqJagged_coordinates,
            criteriaDim, pqJagged_weights, pqJagged_gnos, ignoreWeights,
            pqJagged_uniformWeights, pqJagged_uniformParts, pqJagged_partSizes);

    //////////////////BEGINNING OF THE FUNCTION///////////////////////

#ifdef enable_migration2
    const pq_gno_t *actual_pqgnos = pqJagged_gnos.getRawPtr();
    pq_gno_t *pq_gnos = (pq_gno_t *)actual_pqgnos;
    int *actual_owner_of_coordinate  = NULL;
#endif

    //as input indices.
    pq_lno_t *partitionedPointCoordinates =  allocMemory< pq_lno_t>(numLocalCoords);
    //initial configuration
    //set each pointer-i to i.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < static_cast<size_t>(numLocalCoords); ++i){
        partitionedPointCoordinates[i] = i;
    }

    //as output indices
    pq_lno_t *newpartitionedPointCoordinates = allocMemory< pq_lno_t>(numLocalCoords);
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<pq_lno_t>(newpartitionedPointCoordinates, numLocalCoords);
#endif
#endif

    partId_t *partIds = NULL;
    ArrayRCP<partId_t> partId;
    if(numLocalCoords > 0){
        partIds = allocMemory<partId_t>(numLocalCoords);
    }

    //initially there is a single partition
    partId_t currentPartitionCount = 1;
    //single partition starts at index-0, and ends at numLocalCoords
    //inTotalCounts array holds the end points in partitionedPointCoordinates array
    //for each partition. Initially sized 1, and single element is set to numLocalCoords.
    pq_lno_t *inTotalCounts = allocMemory<pq_lno_t>(1);
    inTotalCounts[0] = static_cast<pq_lno_t>(numLocalCoords);//the end of the initial partition is the end of coordinates.
    //the ends points of the output.
    pq_lno_t *outTotalCounts = NULL;

    //get pqJagged specific parameters.
    bool allowNonRectelinearPart = false;
    int concurrentPartCount = 0; // Set to invalid value
    int migration_actualMigration_option = 1;
    int migration_check_option = 0;
    int migration_all2all_option = 1;
    pq_scalar_t migration_imbalance_cut_off = 0.35;
    int migration_assignment_type = 0;
    int keep_part_boxes = 0;

    int enable_rcb = 0;
    int recursion_depth = -1;

    pqJagged_getParameters<pq_scalar_t>(pl,
        allowNonRectelinearPart, concurrentPartCount,
        migration_actualMigration_option, migration_check_option,
        migration_all2all_option, migration_imbalance_cut_off,
        migration_assignment_type, keep_part_boxes, enable_rcb,
        recursion_depth);

    //cout << "enable_rcb:" << enable_rcb << " recursion_depth:" <<
    //recursion_depth << endl;
    //cout << "partArraySize:" << partArraySize << " recursion_depth:" <<
    //recursion_depth << endl;
    if (enable_rcb){
        recursion_depth = (int)(ceil(log ((numGlobalParts)) / log (2.0)));
    }
    if (partArraySize < 1){
        if (recursion_depth > 0){
            partArraySize = recursion_depth;
        }
        else {
            partArraySize = coordDim;
        }
    }

    //cout << "partArraySize:" << partArraySize << " recursion_depth:" <<
    //recursion_depth << endl;

    int numThreads = 1;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(numThreads)
    {
        numThreads = omp_get_num_threads();
    }

#endif

    partId_t totalDimensionCut = 0; //how many cuts will be totally
    partId_t totalPartCount = 1;    //how many parts will be totally
    partId_t maxPartNo = 0;         //maximum cut count along a dimension.
    partId_t reduceAllCount = 0;    //estimate on #reduceAlls can be done.
    partId_t maxTotalCumulativePartCount = 1; //max no of parts that might occur
                                              //during the partition before the
                                              //last partitioning dimension.
    partId_t maxCutNo = 0;
    size_t maxTotalPartCount = 0;

    getPartSpecifications <partId_t>( partNo, partArraySize, numGlobalParts,
     totalDimensionCut, totalPartCount, maxPartNo, reduceAllCount,
      maxTotalCumulativePartCount, maxCutNo, maxTotalPartCount);


    if (concurrentPartCount == 0)
    {
        if(problemComm->getSize() == 1){
            concurrentPartCount = 1;
        }
        // User did not specify concurrentPartCount parameter, Pick a default.
        // Still a conservative default. We could go as big as
        // maxTotalCumulativePartCount, but trying not to use too much memory.
        // Another assumption the total #parts will be large-very large
        else if (coordDim == partArraySize)
        {
            // partitioning each dimension only once, pick the default as
            // maxPartNo >> default
            concurrentPartCount = min(Z2_DEFAULT_CON_PART_COUNT, maxPartNo);
        }
        else
        {
            // partitioning each dimension more than once, pick the max
            concurrentPartCount = max(Z2_DEFAULT_CON_PART_COUNT, maxPartNo);
        }
    }

    if(concurrentPartCount > maxTotalCumulativePartCount){
        if(problemComm->getRank() == 0){
            cerr << "Warning: Concurrent part count ("<< concurrentPartCount <<
            ") has been set bigger than maximum amount that can be used." <<
            " Setting to:" << maxTotalCumulativePartCount << "." << endl;
        }
        concurrentPartCount = maxTotalCumulativePartCount;
    }

    //We duplicate the comm as we create subcommunicators during migration.
    //We keep the problemComm as it is, while comm changes after each migration.
    RCP<Comm<int> > comm = problemComm->duplicate();

    // coordinates of the cut lines. First one is the min, last one is max
    // coordinate.
    // kddnote if (keep_cuts)
    // coordinates of the cut lines.
    //only store this much if cuts are needed to be stored.
    pq_scalar_t *allCutCoordinates = allocMemory< pq_scalar_t>(totalDimensionCut);
    // kddnote else
    //pq_scalar_t *allCutCoordinates = allocMemory< pq_scalar_t>(maxCutNo * concurrentPartCount);
    pq_scalar_t *max_min_array =  allocMemory< pq_scalar_t>(numThreads * 2);

    float *nonRectelinearPart = NULL; //how much weight percentage should a MPI put left side of the each cutline
    float **nonRectRatios = NULL; //how much weight percentage should each thread in MPI put left side of the each cutline
    if(allowNonRectelinearPart){
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
    //necessary because previous cut line information is used for determining
    //the next cutline information. therefore, cannot update the cut work array
    //until all cutlines are determined.
    pq_scalar_t *cutCoordinatesWork = allocMemory<pq_scalar_t>(maxCutNo *
    concurrentPartCount);

#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<pq_scalar_t>(cutCoordinatesWork, maxCutNo);
#endif
#endif

    //cumulative part weight ratio array.
    pq_scalar_t *targetPartWeightRatios = allocMemory<pq_scalar_t>(maxPartNo *
    concurrentPartCount); // the weight ratios at left side of the cuts. First is 0, last is 1.
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<pq_scalar_t>(cutPartRatios, maxCutNo);
#endif
#endif


    pq_scalar_t *cutUpperBounds = allocMemory<pq_scalar_t>(maxCutNo * concurrentPartCount);  //upper bound coordinate of a cut line
    pq_scalar_t *cutLowerBounds = allocMemory<pq_scalar_t>(maxCutNo* concurrentPartCount);  //lower bound coordinate of a cut line
    pq_scalar_t *cutLowerWeight = allocMemory<pq_scalar_t>(maxCutNo* concurrentPartCount);  //lower bound weight of a cut line
    pq_scalar_t *cutUpperWeight = allocMemory<pq_scalar_t>(maxCutNo* concurrentPartCount);  //upper bound weight of a cut line

    pq_scalar_t *localMinMaxTotal = allocMemory<pq_scalar_t>(3 * concurrentPartCount); //combined array to exchange the min and max coordinate, and total weight of part.
    pq_scalar_t *globalMinMaxTotal = allocMemory<pq_scalar_t>(3 * concurrentPartCount);//global combined array with the results for min, max and total weight.

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
    pq_scalar_t **leftClosestDistance = allocMemory<pq_scalar_t *>(numThreads);
    //leftClosesDistance to hold the min distance of a coordinate to a cutline from right (for each thread)
    pq_scalar_t **rightClosestDistance = allocMemory<pq_scalar_t *>(numThreads);

    //to store how many points in each part a thread has.
    pq_lno_t **partPointCounts = allocMemory<pq_lno_t *>(numThreads);

    for(int i = 0; i < numThreads; ++i){
        //partWeights[i] = allocMemory<pq_scalar_t>(maxTotalPartCount);
        partWeights[i] = allocMemory < double >(maxTotalPartCount * concurrentPartCount);
        rightClosestDistance[i] = allocMemory<pq_scalar_t>(maxCutNo * concurrentPartCount);
        leftClosestDistance[i] = allocMemory<pq_scalar_t>(maxCutNo * concurrentPartCount);
        partPointCounts[i] =  allocMemory<pq_lno_t>(maxPartNo);
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
    pq_scalar_t *cutWeights = allocMemory<pq_scalar_t>(maxCutNo);
    pq_scalar_t *globalCutWeights = allocMemory<pq_scalar_t>(maxCutNo);

    //for faster communication, concatanation of
    //totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
    //leftClosest distances sized P-1, since P-1 cut lines
    //rightClosest distances size P-1, since P-1 cut lines.
    pq_scalar_t *totalPartWeights_leftClosests_rightClosests = allocMemory<pq_scalar_t>((maxTotalPartCount + maxCutNo * 2) * concurrentPartCount);
    pq_scalar_t *global_totalPartWeights_leftClosests_rightClosests = allocMemory<pq_scalar_t>((maxTotalPartCount + maxCutNo * 2) * concurrentPartCount);


    pq_scalar_t *cutCoordinates =  allCutCoordinates;

    //partId_t leftPartitions = totalPartCount;
    pq_scalar_t maxScalar_t = numeric_limits<pq_scalar_t>::max();
    pq_scalar_t minScalar_t = -numeric_limits<pq_scalar_t>::max();


    env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Partitioning");

    //create multivector based on the coordinates, weights, and multivectordim
#ifdef enable_migration2
    typedef Tpetra::MultiVector<pq_scalar_t, pq_lno_t, pq_gno_t, pq_node_t> mvector_t;
    RCP<const mvector_t> mvector;
    if(migration_actualMigration_option == 0){
            mvector =  create_initial_multi_vector
             <pq_gno_t,pq_lno_t,pq_scalar_t,pq_node_t>( env, comm,
             numGlobalCoords, numLocalCoords, coordDim, pqJagged_coordinates,
             weightDim, pqJagged_weights, pqJagged_multiVectorDim);
    }
    else {
        for (int i=0; i < coordDim; i++){
            pq_scalar_t *coord = allocMemory<pq_scalar_t>(numLocalCoords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
            for (pq_lno_t j=0; j < numLocalCoords; j++)
                coord[j] = pqJagged_coordinates[i][j];

            pqJagged_coordinates[i] = coord;
        }
        for (int i=0; i < weightDim; i++){
            pq_scalar_t *w = allocMemory<pq_scalar_t>(numLocalCoords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
            for (pq_lno_t j=0; j < numLocalCoords; j++)
                w[j] = pqJagged_weights[i][j];

            pqJagged_weights[i] = w;
        }
        pq_gnos = allocMemory<pq_gno_t>(numLocalCoords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (pq_lno_t j=0; j < numLocalCoords; j++)
            pq_gnos[j] = actual_pqgnos[j];

        actual_owner_of_coordinate = allocMemory<int>(numLocalCoords);
        int me = comm->getRank();
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (pq_lno_t j=0; j < numLocalCoords; j++)
            actual_owner_of_coordinate[j] = me;
    }
#endif

    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon();
    partId_t partIndexBegin = 0;
    partId_t futurePartNumbers = totalPartCount;
    bool is_data_ever_migrated = false;

    vector<partId_t> *currentPartitions = new vector<partId_t> ();
    vector<partId_t> *newFuturePartitions = new vector<partId_t> ();
    newFuturePartitions->push_back(numGlobalParts);

    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
     inPartBoxes(new vector <coordinateModelPartBox <pq_scalar_t, partId_t> >
      (), true) ;
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
     outPartBoxes(new vector <coordinateModelPartBox <pq_scalar_t, partId_t> >
      (), true);

    if(keep_part_boxes){
        pq_scalar_t *mins = allocMemory<pq_scalar_t>(coordDim);
        pq_scalar_t *gmins = allocMemory<pq_scalar_t>(coordDim);
        pq_scalar_t *maxs = allocMemory<pq_scalar_t>(coordDim);
        pq_scalar_t *gmaxs = allocMemory<pq_scalar_t>(coordDim);
        for (int i = 0; i < coordDim; ++i){
            //cout << " pqJagged_coordinates[i][0]:" <<
            //pqJagged_coordinates[i][0] << endl;
            pq_scalar_t localMin = pqJagged_coordinates[i][0];
            pq_scalar_t localMax = pqJagged_coordinates[i][0];
            for (pq_lno_t j = 1; j < numLocalCoords; ++j){
                //cout << " pqJagged_coordinates[i][i]:" <<
                //pqJagged_coordinates[i][j] << endl;
                if (pqJagged_coordinates[i][j] < localMin){
                    localMin = pqJagged_coordinates[i][j];
                }
                if (pqJagged_coordinates[i][j] > localMax){
                    localMax = pqJagged_coordinates[i][j];
                }
            }
            //cout << " localMin:" << localMin << endl;
            //cout << " localMax:" << localMax << endl;
            mins[i] = localMin;
            maxs[i] = localMax;
        }
        reduceAll<int, pq_scalar_t>(*comm, Teuchos::REDUCE_MIN,
                coordDim, mins, gmins
        );

        reduceAll<int, pq_scalar_t>(*comm, Teuchos::REDUCE_MAX,
                coordDim, maxs, gmaxs
        );
        coordinateModelPartBox <pq_scalar_t, partId_t> tmpBox (0, coordDim,
                                                            gmins, gmaxs);
        //coordinateModelPartBox <pq_scalar_t, partId_t> tmpBox (0, coordDim);
        freeArray<pq_scalar_t>(mins);
        freeArray<pq_scalar_t>(gmins);
        freeArray<pq_scalar_t>(maxs);
        freeArray<pq_scalar_t>(gmaxs);
        outPartBoxes->push_back(tmpBox);
    }

    //cout << "partArraySize:" << partArraySize << endl;
    for (int i = 0; i < partArraySize; ++i){
        //partitioning array.
        //size will be as the number of current partitions
        //and this hold how many parts each part will be
        //in the current dimension partitioning.
        vector <partId_t> pAlongI;

        //number of parts that will be obtained at the end of this partitioning.
        //currentPartitions is as the size of current number of parts.
        //holds how many more parts each should be divided in the further
        //iterations.
        //this will be used to calculate pAlongI,
        //as the number of parts that the part will be partitioned
        //in the current dimension partitioning.

        //newFuturePartitions will be as the size of outnumParts,
        //and this will hold how many more parts that each output part
        //should be divided.
        //this array will also be used to determine the weight ratios
        //of the parts.
        //swap the arrays.
        vector<partId_t> *tmpPartVect= currentPartitions;
        currentPartitions = newFuturePartitions;
        newFuturePartitions = tmpPartVect;

        //clear newFuturePartitions array as
        //getPartitionArrays expects it to be empty.
        //it also expects pAlongI to be empty as well.
        newFuturePartitions->clear();

        if(keep_part_boxes){
            RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
                                                 tmpPartBoxes = inPartBoxes;
            inPartBoxes = outPartBoxes;
            outPartBoxes = tmpPartBoxes;
            outPartBoxes->clear();

/*
            for (partId_t j = 0; j < inPartBoxes->size(); ++j){

                cout << "me:" << j << endl;
                (*inPartBoxes)[j].print();
            }
            */
        }

        //returns the total no. of output parts for this dimension partitioning.
        partId_t outPartCount = getPartitionArrays<pq_scalar_t, partId_t>(
                partNo, pAlongI, currentPartitions, newFuturePartitions,
                futurePartNumbers, currentPartitionCount, partArraySize, i,
                maxPartNo, keep_part_boxes, inPartBoxes, outPartBoxes);

        if(outPartCount == currentPartitionCount) {
            tmpPartVect= currentPartitions;
            currentPartitions = newFuturePartitions;
            newFuturePartitions = tmpPartVect;

            if(keep_part_boxes){
                RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
                                             tmpPartBoxes = inPartBoxes;
                inPartBoxes = outPartBoxes;
                outPartBoxes = tmpPartBoxes;
            }
            continue;
        }

#ifdef enable_migration2
        int worldSize = comm->getSize();
        long migration_reduceAllPop = reduceAllCount * worldSize;
#endif

        //get the coordinate axis along which the partitioning will be done.
        int coordInd = i % coordDim;
        pq_scalar_t * pqCoord = pqJagged_coordinates[coordInd];
        //convert i to string to be used for debugging purposes.
        string istring = toString<int>(i);

        env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Partitioning_" +
                         istring);

        //alloc Memory to point the indices
        //of the parts in the permutation array.
        outTotalCounts = allocMemory<pq_lno_t>(outPartCount);

        //the index where in the outtotalCounts will be written.
        partId_t currentOut = 0;
        //whatever is written to outTotalCounts will be added with previousEndi
        //so that the points will be shifted.
        partId_t previousEnd = 0;

        partId_t currentWorkPart = 0;
        partId_t concurrentPart = min(currentPartitionCount - currentWorkPart,
                            concurrentPartCount);

        //always use binary search algorithm.
        bool useBinarySearch = true;
        partId_t obtainedPartCount = 0;

        //run for all available parts.
        for (; currentWorkPart < currentPartitionCount;
                 currentWorkPart += concurrentPart){

            concurrentPart = min(currentPartitionCount - currentWorkPart,
                                 concurrentPartCount);
#ifdef mpi_communication
            concurrent = concurrentPart;
#endif

            partId_t workPartCount = 0;
            //get the min and max coordinates of each part
            //together with the part weights of each part.
            for(int kk = 0; kk < concurrentPart; ++kk){
                partId_t currentPart = currentWorkPart + kk;

                //if this part wont be partitioned any further
                //dont do any work for this part.
                if (pAlongI[currentPart] == 1){
                    continue;
                }
                ++workPartCount;
                pq_lno_t coordinateEnd= inTotalCounts[currentPart];
                pq_lno_t coordinateBegin = currentPart==0 ? 0: inTotalCounts[currentPart -1];
                //cout << "me:" << problemComm->getRank() << " begin:" <<i
                //coordinateBegin  << " end:" << coordinateEnd << endl;
                pqJagged_getLocalMinMaxTotalCoord<pq_scalar_t, pq_lno_t>(
                    partitionedPointCoordinates, pqCoord,
                    pqJagged_uniformWeights[0], pqJagged_weights[0], numThreads,
                    coordinateBegin, coordinateEnd, max_min_array, maxScalar_t,
                    minScalar_t, localMinMaxTotal[kk],
                    localMinMaxTotal[kk + concurrentPart],
                    localMinMaxTotal[kk + 2*concurrentPart]);
            }

            if (workPartCount > 0){
                //obtain global Min max of the part.
                pqJagged_getGlobalMinMaxTotalCoord<pq_scalar_t>( comm, env,
                 concurrentPart, localMinMaxTotal, globalMinMaxTotal);

                //represents the total number of cutlines
                //whose coordinate should be determined.
                partId_t allDone = 0;

                //Compute weight ratios for parts & cuts:
                //e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
                //part0  cut0  part1 cut1 part2 cut2 part3
                partId_t cutShifts = 0;
                partId_t partShift = 0;
                for(int kk = 0; kk < concurrentPart; ++kk){
                    pq_scalar_t minCoordinate = globalMinMaxTotal[kk];
                    pq_scalar_t maxCoordinate = globalMinMaxTotal[kk +
                                                     concurrentPart];

                    pq_scalar_t globalTotalWeight = globalMinMaxTotal[kk +
                                                        2 * concurrentPart];

                    partId_t currentPart = currentWorkPart + kk;

                    partId_t partition = pAlongI[currentPart];

                    pq_scalar_t *usedCutCoordinate = cutCoordinates + cutShifts;
                    pq_scalar_t *usedCutPartRatios = targetPartWeightRatios +
                                                        partShift;
                    //shift the usedCutCoordinate array as noCuts.
                    cutShifts += partition - 1;
                    //shift the partRatio array as noParts.
                    partShift += partition;

                    //cout << "min:" << minCoordinate << " max:" <<
                    //maxCoordinate << endl;
 
                    //calculate only if part is not empty,
                    //and part will be further partitioend.
                    if(partition > 1 && minCoordinate <= maxCoordinate){

                        //increase allDone by the number of cuts of the current
                        //part's cut line number.
                        allDone += partition - 1;
                        //set the number of cut lines that should be determined
                        //for this part.
                        myNonDoneCount[kk] = partition - 1;

                        //get the target weights of the parts.
                        pqJagged_getCutCoord_Weights<pq_scalar_t>(minCoordinate,
                         maxCoordinate, pqJagged_uniformParts[0],
                         pqJagged_uniformWeights[0], pqJagged_partSizes[0],
                         partition - 1, globalTotalWeight, usedCutCoordinate,
                         usedCutPartRatios, numThreads, currentPartitions,
                         newFuturePartitions, currentPart, obtainedPartCount);

                        //get the initial estimated part assignments of the
                        //coordinates.
                        getInitialPartAssignments<pq_scalar_t, pq_lno_t, partId_t>(
                            maxCoordinate, minCoordinate, currentPart,
                            inTotalCounts, partitionedPointCoordinates,
                            pqCoord, partIds, _EPSILON, partition);
                    }
                    else {
                        // e.g., if have fewer coordinates than parts, don't
                        // need to do next dim.
                        myNonDoneCount[kk] = 0;
                    }
                    obtainedPartCount += partition;
                }



                //used imbalance, it is always 0, as it is difficult to
                //estimate a range.
                pq_scalar_t used_imbalance = 0;


                // Determine cut lines for k parts here.
                pqJagged_1D_Partition<pq_scalar_t, pq_lno_t>(
                    env, comm, partitionedPointCoordinates, pqCoord,
                    pqJagged_uniformWeights[0], pqJagged_weights[0],
                    targetPartWeightRatios, globalMinMaxTotal, localMinMaxTotal,
                    numThreads, used_imbalance, currentWorkPart, concurrentPart,
                    inTotalCounts, cutCoordinates, cutCoordinatesWork,
                    leftClosestDistance, rightClosestDistance, cutUpperBounds,
                    cutLowerBounds, cutUpperWeight, cutLowerWeight, isDone,
                    partWeights, totalPartWeights_leftClosests_rightClosests,
                    global_totalPartWeights_leftClosests_rightClosests,
                    allowNonRectelinearPart, nonRectelinearPart, cutWeights,
                    globalCutWeights, allDone, myNonDoneCount, useBinarySearch,
                    partIds, pAlongI);
            }

            //create part chunks
            {

                partId_t outShift = 0;
                partId_t cutShift = 0;
                size_t tlrShift = 0;
                size_t pwShift = 0;

                for(int kk = 0; kk < concurrentPart; ++kk){
                    partId_t curr = currentWorkPart + kk;
                    partId_t noParts = pAlongI[curr];

                    //if the part is empty, skip the part.
                    if((noParts != 1  )&& globalMinMaxTotal[kk] >
                             globalMinMaxTotal[kk + concurrentPart]) {

                        for(partId_t jj = 0; jj < noParts; ++jj){
                            outTotalCounts[currentOut + outShift + jj] = 0;
                        }
                        cutShift += noParts - 1;
                        tlrShift += (4 *(noParts - 1) + 1);
                        outShift += noParts;
                        pwShift += (2 * (noParts - 1) + 1);
                        continue;
                    }

                    pq_lno_t coordinateEnd= inTotalCounts[curr];
                    pq_lno_t coordinateBegin = curr==0 ? 0: inTotalCounts[
                                                                curr -1];
                    pq_scalar_t *usedCutCoordinate = cutCoordinates + cutShift;
                    float *usednonRectelinearPart = nonRectelinearPart +
                                                            cutShift;

                    pq_scalar_t *tlr =  totalPartWeights_leftClosests_rightClosests + tlrShift;

                    for(int ii = 0; ii < numThreads; ++ii){
                        pws[ii] = partWeights[ii] +  pwShift;
                    }
/*
                    cout << endl;
                    pq_scalar_t prev = 0;
                    for (int i = 0; i < noParts * 2 - 1 ; ++i){
                        pq_scalar_t w = 0;
                        for(int ii = 0; ii < numThreads; ++ii){
                            w += pws[ii][i];
                        }

                        cout << "\ti:" << i  << " w:" << w - prev << endl;
                        prev = w;
                    }

*/
                    if(noParts > 1){
                        if(keep_part_boxes){
                            for (partId_t j = 0; j < noParts - 1; ++j){
/*
                                cout << " outShift + currentOut + j: " <<
                                 outShift + currentOut + j << endl;
                                cout << " coordInd:" << coordInd << " cut:" <<
                                 "usedCutCoordinate["<< j<< "]" <<
                                  usedCutCoordinate[j] << endl;

                                cout << "me:" << problemComm->getRank() << "
                                 before update:" << endl;

                                (*outPartBoxes)[outShift + currentOut +
                                 j].print();
*/
                                (*outPartBoxes)[outShift + currentOut +
                                 j].updateMinMax(usedCutCoordinate[j], 1
                                  /*update max*/, coordInd);
/*
                                cout << "me:" << problemComm->getRank() <<
                                  " after update:"<< endl;
                                (*outPartBoxes)[outShift + currentOut +
                                 j].print();

                                cout <<  "me:" << problemComm->getRank() <<
                                " before update:"<< endl;

                                (*outPartBoxes)[outShift + currentOut + j +
                                1].print();
*/
                                (*outPartBoxes)[outShift + currentOut + j +
                                 1].updateMinMax(usedCutCoordinate[j], 0
                                  /*update min*/, coordInd);
/*
                                cout <<  "me:" << problemComm->getRank() <<
                                 " after update:"<< endl;
                                (*outPartBoxes)[outShift + currentOut + j +
                                 1].print();
*/

                            }
                        }

                        // Rewrite the indices based on the computed cuts.
                        getChunksFromCoordinates<pq_lno_t,pq_scalar_t>(
                            noParts, numThreads, partitionedPointCoordinates,
                            pqCoord, pqJagged_uniformWeights[0],
                            pqJagged_weights[0], usedCutCoordinate,
                            coordinateBegin, coordinateEnd,
                            allowNonRectelinearPart, usednonRectelinearPart,
                            tlr, pws, nonRectRatios, partPointCounts,
                            newpartitionedPointCoordinates,
                            outTotalCounts + currentOut + outShift, partIds,
                            false, pqJagged_coordinates, coordDim, coordInd );
                        /*
                        pq_lno_t *mm = outTotalCounts + currentOut + outShift;
                        for (int i = 0; i < noParts; ++i){
                            pq_lno_t pbeg = 0;
                            pq_lno_t pend = mm[i];
                            if (i > 0) pbeg = mm[i - 1];


                            cout << "\ti:" << i << " w:" << pend - pbeg << endl;
                        }
                        */
                    }
                    else {
                        //if this part is partitioned into 1 then just copy
                        //the old values.
                        pq_lno_t partSize = coordinateEnd - coordinateBegin;
                        *(outTotalCounts + currentOut + outShift) = partSize;
                        memcpy(newpartitionedPointCoordinates + coordinateBegin,
                             partitionedPointCoordinates + coordinateBegin,
                             partSize * sizeof(pq_lno_t));
                    }
                    cutShift += noParts - 1;
                    tlrShift += (4 *(noParts - 1) + 1);
                    outShift += noParts;
                    pwShift += (2 * (noParts - 1) + 1);
                }

                //shift cut coordinates so that all cut coordinates are stored.
                cutCoordinates += cutShift;

                //getChunks from coordinates partitioned the parts and
                //wrote the indices as if there were a single part.
                //now we need to shift the beginning indices.
                for(partId_t kk = 0; kk < concurrentPart; ++kk){
                    partId_t noParts = pAlongI[ currentWorkPart + kk];
                    for (partId_t ii = 0;ii < noParts ; ++ii){
                        //shift it by previousCount
                        outTotalCounts[currentOut+ii] += previousEnd;
                    }
                    //increase the previous count by current end.
                    previousEnd = outTotalCounts[currentOut + noParts - 1];
                    //increase the current out.
                    currentOut += noParts ;
                }
            }
        }
        // end of this partitioning dimension


        bool is_migrated_in_current = false;
        /*
        //now check for the migration.
        cout << "futurePartNumbers :" << futurePartNumbers  <<
                " migration_check_option:" << migration_check_option <<
                " worldSize:" << worldSize  << endl;
                */

#ifdef enable_migration2
        if (
                futurePartNumbers > 1 &&
                migration_check_option >= 0 &&
                worldSize > 1){
            env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Migration-" +
             istring);


            partId_t num_parts = outPartCount;
            if (
                    migration_refactored<pq_gno_t, pq_lno_t, pq_scalar_t, pq_node_t, partId_t>(
                        problemComm, env, comm, mvector, 
                        pqJagged_multiVectorDim, numGlobalCoords,
                        numLocalCoords, coordDim,
                        pqJagged_coordinates, //outout will be modified.
                        weightDim,
                        pqJagged_weights,// output will be modified.
                        partIds, num_parts,
                        currentPartitionCount, //output
                        newFuturePartitions, //output
                        newpartitionedPointCoordinates, //output
                        partitionedPointCoordinates,
                        outTotalCounts //output
                        ,partIndexBegin,
                        migration_all2all_option, migration_assignment_type,
                        migration_actualMigration_option, //migration_proc_assignment_type
                        migration_check_option, migration_imbalance_cut_off,
                        migration_reduceAllPop,
                        numLocalCoords / (futurePartNumbers *
                         currentPartitionCount) ,
                        //used when z1 migration is used.
                        pq_gnos, actual_owner_of_coordinate, istring,
                        keep_part_boxes, inPartBoxes, outPartBoxes)
            )
            {
                is_migrated_in_current = true;
                is_data_ever_migrated = true;
                env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Migration-" +
                                 istring);
                reduceAllCount /= num_parts;

            }
            else {
                is_migrated_in_current = false;
                env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Migration-" +
                             istring);
            }
        }
#endif

        pq_lno_t * tmp = partitionedPointCoordinates;
        partitionedPointCoordinates = newpartitionedPointCoordinates;
        newpartitionedPointCoordinates = tmp;

        if(!is_migrated_in_current){
            reduceAllCount -= currentPartitionCount;
            currentPartitionCount = outPartCount;
        }
        freeArray<pq_lno_t>(inTotalCounts);
        inTotalCounts = outTotalCounts;

        env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Partitioning_" +
                                     istring);
    }

    //cout << "me:" << problemComm->getRank() << " done:" << endl;

    // Partitioning is done
    delete currentPartitions;
    delete newFuturePartitions;

    env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Partitioning");
    /////////////////////////////End of the function////////////////////////

    env->timerStart(MACRO_TIMERS, "PQJagged - Part_Assignment");

#ifdef debug_setparts
    string fname = toString<int>(problemComm->getRank());
    fname = fname   + ".print";
    FILE *f = fopen(fname.c_str(), "w");
#endif
/*
    problemComm->barrier();
    cout    << "me:" << problemComm->getRank()
            << " partIndexBegin:" << partIndexBegin
            << " currentPartitionCount:" << currentPartitionCount
            << endl;
    problemComm->barrier();
    */
#ifdef writeParts
    ofstream *ss = new ofstream[currentPartitionCount];

    for(partId_t i = 0; i < currentPartitionCount;++i){
        string a = toString<partId_t>(i)+".part";
        ss[i].open (a.c_str(), std::ofstream::out);
    }
#endif

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(partId_t i = 0; i < currentPartitionCount;++i){

        pq_lno_t begin = 0;
        pq_lno_t end = inTotalCounts[i];

        if(i > 0) begin = inTotalCounts[i -1];
        partId_t pToSet = i + partIndexBegin;
        if (keep_part_boxes){
            (*outPartBoxes)[i].setpId(pToSet);
        }
        /*
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
         */
        //cout <<"me:" << problemComm->getRank()<< " p:" << i + partIndexBegin
        //<< " count:" << end - begin << endl;
        for (pq_lno_t ii = begin; ii < end; ++ii){
            pq_lno_t k = partitionedPointCoordinates[ii];
            partIds[k] = pToSet;


#ifdef writeParts
            for (int d = 0; d < coordDim; ++d){
                ss[pToSet] << pqJagged_coordinates[d][k] << " ";
            }
            ss[pToSet] << endl;
#endif

#ifdef debug_setparts
            fprintf(f, "setting %d with coords: %lf %lf %lf to part %d\n", k,
             pqJagged_coordinates[0][k],  pqJagged_coordinates[1][k],
             pqJagged_coordinates[2][k], partIds[k]);
#endif

        }
    }

#ifdef writeParts
    for(partId_t i = 0; i < currentPartitionCount;++i){
        ss[i].close();
    }
    if(keep_part_boxes){
        ofstream arrowLines("arrows.gnuplot");

        arrowLines << "set nokey" << endl;
        ofstream corners("corners.gnuplot");
        for (partId_t i = 0; i < partId_t (outPartBoxes->size()); ++i){
            (*outPartBoxes)[i].writeGnuPlot(arrowLines, corners);

        }
        for(partId_t i = 0; i < currentPartitionCount;++i){
            string plotStr = "";
            if (i == 0){
                if (coordDim == 2){
                    plotStr = "plot ";
                }
                else {
                    plotStr = "splot ";
                }
            }
            else {
                plotStr = "replot ";
            }
            string a = toString<partId_t>(i)+".part";
            plotStr += "\""+a + "\"";
            arrowLines << plotStr << endl;
        }

        arrowLines << "set terminal png" << endl;
        arrowLines << "set output 'output.png'\nreplot" << endl;

        //arrowLines << "pause -1" << endl;
        arrowLines.close();
        corners.close();

    }
#endif
    //cout << "me:" << problemComm->getRank() << " currentPartitionCount:"
    //<< currentPartitionCount << endl;

#ifdef debug_setparts
    fclose(f);
#endif
    env->timerStop(MACRO_TIMERS, "PQJagged - Part_Assignment");
    ArrayRCP<const pq_gno_t> gnoList;
    if(!is_data_ever_migrated){
#ifdef enable_migration2
        if(migration_actualMigration_option != 0){
            freeArray<pq_gno_t>(pq_gnos);
        }
#endif
        if(numLocalCoords > 0){
            gnoList = arcpFromArrayView(pqJagged_gnos);
        }
    }
#ifdef enable_migration2
    else {
        if(migration_actualMigration_option == 0){
            gnoList = arcpFromArrayView(mvector->getMap()->getNodeElementList());
        }
        else {

            ZOLTAN_COMM_OBJ *plan = NULL;
            MPI_Comm mpi_comm = Teuchos2MPI (problemComm);

            pq_lno_t incoming = 0;
            int message_tag = 7856;

            env->timerStart(MACRO_TIMERS, "PQJagged - Final Z1PlanCreating");
            int ierr = Zoltan_Comm_Create( &plan, numLocalCoords,
                    actual_owner_of_coordinate, mpi_comm, message_tag,
                    &incoming);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
            env->timerStop(MACRO_TIMERS, "PQJagged - Final Z1PlanCreating" );

            pq_gno_t *incoming_gnos = allocMemory< pq_gno_t>(incoming);

            message_tag++;
            env->timerStart(MACRO_TIMERS, "PQJagged - Final Z1PlanComm");
            ierr = Zoltan_Comm_Do( plan, message_tag, (char *) pq_gnos,
                    sizeof(pq_gno_t), (char *) incoming_gnos);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

            freeArray<pq_gno_t>(pq_gnos);
            pq_gnos = incoming_gnos;

            partId_t *incoming_partIds = allocMemory< partId_t>(incoming);

            message_tag++;
            ierr = Zoltan_Comm_Do( plan, message_tag, (char *) partIds,
                    sizeof(partId_t), (char *) incoming_partIds);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
            freeArray<partId_t>(partIds);
            partIds = incoming_partIds;

            env->timerStop(MACRO_TIMERS, "PQJagged - Final Z1PlanComm");
            ierr = Zoltan_Comm_Destroy(&plan);
            Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

            numLocalCoords = incoming;
            is_data_ever_migrated = false;

            gnoList = arcp(pq_gnos, 0, numLocalCoords, true);
            //cout << " me:" << problemComm->getRank() << " gnoList:" <<
            //gnoList() << endl;
        }

    }
#endif
    env->timerStop(MACRO_TIMERS, "PQJagged - Total2");

    env->timerStart(MACRO_TIMERS, "PQJagged - Solution_Part_Assignment");
    partId = arcp(partIds, 0, numLocalCoords, true);

    // SR TODO: migrate_gid never seem to be defined why is this needed ?
    // Looks like old code when we did not do migration within pqjagged but
    // let setParts take care of it. Now as we do migration above
    // we will never pass true to setParts. The next three lines should be
    // removed in that case.
#ifdef migrate_gid
    solution->setParts(gnoList, partId,true);
#endif

    if (keep_part_boxes){
        solution->setPartBoxes(outPartBoxes);
    }
#ifndef migrate_gid
    solution->setParts(gnoList, partId,!is_data_ever_migrated);
#endif
    env->timerStop(MACRO_TIMERS, "PQJagged - Solution_Part_Assignment");

    env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Free");

    /*
    if(comm->getRank() == 0){
    for(partId_t i = 0; i < totalPartCount - 1;++i){
    cout << "i:"<< i<<" cut coordinate:" << allCutCoordinates[i] << endl;
    }
    }
    */

#ifdef enable_migration2
    if (migration_actualMigration_option != 0){

        for (int i=0; i < coordDim; i++){

            freeArray<pq_scalar_t>(pqJagged_coordinates[i]);
        }
        for (int i=0; i < weightDim; i++){

            freeArray<pq_scalar_t>(pqJagged_weights[i]);
        }

        freeArray<int>(actual_owner_of_coordinate);
    }
#endif

    for(int i = 0; i < numThreads; ++i){
        freeArray<pq_lno_t>(partPointCounts[i]);
    }

    freeArray<pq_lno_t *>(partPointCounts);
    freeArray<double *> (pws);

    if(allowNonRectelinearPart){
        freeArray<float>(nonRectelinearPart);
        for(int i = 0; i < numThreads; ++i){
            freeArray<float>(nonRectRatios[i]);
        }
        freeArray<float *>(nonRectRatios);
    }

    freeArray<partId_t>(myNonDoneCount);

    freeArray<pq_scalar_t>(cutWeights);

    freeArray<pq_scalar_t>(globalCutWeights);

    freeArray<pq_scalar_t>(max_min_array);

    freeArray<pq_lno_t>(outTotalCounts);

    freeArray<pq_lno_t>(partitionedPointCoordinates);

    freeArray<pq_lno_t>(newpartitionedPointCoordinates);

    freeArray<pq_scalar_t>(allCutCoordinates);

    freeArray<pq_scalar_t *>(pqJagged_coordinates);

    freeArray<pq_scalar_t *>(pqJagged_weights);

    freeArray<bool>(pqJagged_uniformParts);

    freeArray<pq_scalar_t> (localMinMaxTotal);

    freeArray<pq_scalar_t> (globalMinMaxTotal);

    freeArray<pq_scalar_t *>(pqJagged_partSizes);

    freeArray<bool>(pqJagged_uniformWeights);

    freeArray<pq_scalar_t>(cutCoordinatesWork);

    freeArray<pq_scalar_t>(targetPartWeightRatios);

    freeArray<pq_scalar_t>(cutUpperBounds);

    freeArray<pq_scalar_t>(cutLowerBounds);


    freeArray<pq_scalar_t>(cutLowerWeight);
    freeArray<pq_scalar_t>(cutUpperWeight);
    freeArray<bool>(isDone);
    freeArray<pq_scalar_t>(totalPartWeights_leftClosests_rightClosests);
    freeArray<pq_scalar_t>(global_totalPartWeights_leftClosests_rightClosests);

    for(int i = 0; i < numThreads; ++i){
        freeArray<double>(partWeights[i]);
        freeArray<pq_scalar_t>(rightClosestDistance[i]);
        freeArray<pq_scalar_t>(leftClosestDistance[i]);
    }

    freeArray<double *>(partWeights);
    freeArray<pq_scalar_t *>(leftClosestDistance);
    freeArray<pq_scalar_t *>(rightClosestDistance);

    env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Free");
    env->timerStop(MACRO_TIMERS, "PQJagged - Total");
    env->debug(3, "Out of PQ Jagged");
#endif // INCLUDE_ZOLTAN2_EXPERIMENTAL
}

} // namespace Zoltan2

#endif
