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
#include <vector>

#define ECLIPSE
#ifdef ECLIPSE
#define INCLUDE_ZOLTAN2_EXPERIMENTAL
#define HAVE_ZOLTAN2_ZOLTAN
#define HAVE_ZOLTAN2_MPI
#define enable_migration2
#define HAVE_ZOLTAN2_OMP
#endif

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
        bool &distribute_points_on_cut_lines,
        partId_t &max_concurrent_part_calculation,
        int &migration_actualMigration_option,
        int &check_migrate_avoid_migration_option,
        int &migration_all2all_option,
        T &minimum_migration_imbalance,
        int &migration_assignment_type,
        int &mj_keep_part_boxes,
        int &mj_run_as_rcb,
        int &mj_user_recursion_depth){


    // TODO: Append all the parameters with mj_
    const Teuchos::ParameterEntry *pe = pl.getEntryPtr("migration_imbalance_cut_off");
    if (pe){
        double tol;
        tol = pe->getValue(&tol);
        minimum_migration_imbalance = tol - 1.0;
    }

    pe = pl.getEntryPtr("migration_all_to_all_type");
    if (pe){
        migration_all2all_option = pe->getValue(&max_concurrent_part_calculation);
    }else {
        migration_all2all_option = 1;
    }

    pe = pl.getEntryPtr("migration_check_option");
    if (pe){
        check_migrate_avoid_migration_option = pe->getValue(&max_concurrent_part_calculation);
    }else {
        check_migrate_avoid_migration_option = 0;
    }
    if (check_migrate_avoid_migration_option > 1) check_migrate_avoid_migration_option = -1;

    pe = pl.getEntryPtr("migration_processor_assignment_type");
    if (pe){
        migration_assignment_type = pe->getValue(&max_concurrent_part_calculation);
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
        max_concurrent_part_calculation = pe->getValue(&max_concurrent_part_calculation);
    }else {
        max_concurrent_part_calculation = 0; // Set to invalid value
    }

    pe = pl.getEntryPtr("keep_part_boxes");
    if (pe){
        mj_keep_part_boxes = pe->getValue(&mj_keep_part_boxes);
    }else {
        mj_keep_part_boxes = 0; // Set to invalid value
    }

    if (mj_keep_part_boxes == 0){
        pe = pl.getEntryPtr("mapping_type");
        if (pe){
            int mapping_type = -1;
            mapping_type = pe->getValue(&mapping_type);
            if (mapping_type == 0){
                mj_keep_part_boxes  = 1;
            }
        }
    }

    pe = pl.getEntryPtr("mj_enable_rcb");
    if (pe){
        mj_run_as_rcb = pe->getValue(&mj_run_as_rcb);
    }else {
        mj_run_as_rcb = 0; // Set to invalid value
    }

    pe = pl.getEntryPtr("recursion_depth");
    if (pe){
        mj_user_recursion_depth = pe->getValue(&mj_user_recursion_depth);
    }else {
        mj_user_recursion_depth = -1; // Set to invalid value
    }

    int val = 0;
    pe = pl.getEntryPtr("rectilinear_blocks");
    if (pe)
        val = pe->getValue(&val);

    if (val == 1){
        distribute_points_on_cut_lines = false;
    } else {
        distribute_points_on_cut_lines = true;
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
        typename Adapter::base_adapter_t> > &coords, int &coord_dim,
        int &weight_dim, size_t &initial_num_loc_coords, global_size_t &initial_num_glob_coords,
        int &criteria_dim, const bool &ignoreWeights){}


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
    const RCP<const Environment> &mj_env, const RCP<const CoordinateModel<
    typename Adapter::base_adapter_t> > &mj_coords,
    RCP<PartitioningSolution<Adapter> > &mj_solution,
    std::bitset<NUM_RCB_PARAMS> &params,
    const int &coord_dim,
    const int &weight_dim,
    const size_t &num_local_coords,
    size_t &num_global_parts,
    int &mj_multivector_dim,
    pq_scalar_t **mj_coordinates,
    const int &criteria_dim,
    pq_scalar_t **mj_weights,
    ArrayView<const pq_gno_t> &mj_gnos,
    bool &ignoreWeights,
    bool *mj_uniform_weights,
    bool *mj_uniform_parts,
    pq_scalar_t **mj_part_sizes
){
    //typedef typename Adapter::node_t pq_node_t;
    typedef typename Adapter::lno_t pq_lno_t;
    typedef StridedData<pq_lno_t, pq_scalar_t> input_t;

    ArrayView<const pq_gno_t> gnos;
    ArrayView<input_t>     xyz;
    ArrayView<input_t>     wgts;

    mj_coords->getCoordinates(gnos, xyz, wgts);
     mj_gnos = gnos;

    for (int dim=0; dim < coord_dim; dim++){
        ArrayRCP<const pq_scalar_t> ar;
        xyz[dim].getInputArray(ar);
        //pqJagged coordinate values assignment
        mj_coordinates[dim] =  (pq_scalar_t *)ar.getRawPtr();
    }

    if (weight_dim == 0 || ignoreWeights){
        mj_uniform_weights[0] = true;
        mj_weights[0] = NULL;
    }
    else{
        for (int wdim = 0; wdim < weight_dim; wdim++){
            if (wgts[wdim].size() == 0){
                mj_uniform_weights[wdim] = true;
                mj_weights[wdim] = NULL;
            }
            else{
                ArrayRCP<const pq_scalar_t> ar;
                wgts[wdim].getInputArray(ar);
                mj_uniform_weights[wdim] = false;
                mj_weights[wdim] = (pq_scalar_t *) ar.getRawPtr();
            }
        }
    }

    ////////////////////////////////////////////////////////
    // From the Solution we get part information.
    // If the part sizes for a given criteria are not uniform,
    // then they are values that sum to 1.0.

    num_global_parts = mj_solution->getTargetGlobalNumberOfParts();

    for (int wdim = 0; wdim < criteria_dim; wdim++){
        if (mj_solution->criteriaHasUniformPartSizes(wdim)){
            mj_uniform_parts[wdim] = true;
            mj_part_sizes[wdim] = NULL;
        }
        else{
            //TODO : Need a test for this,
            pq_scalar_t *tmp = new pq_scalar_t [num_global_parts];
            mj_env->localMemoryAssertion(__FILE__, __LINE__, num_global_parts, tmp) ;
            for (size_t i=0; i < num_global_parts; i++){
                tmp[i] = mj_solution->getCriteriaPartSize(wdim, i);
            }
            mj_part_sizes[wdim] = tmp;
        }
    }

    // It may not be possible to solve the partitioning problem
    // if we have multiple weight dimensions with part size
    // arrays that differ. So let's be aware of this possibility.

    bool multiplePartSizeSpecs = false;

    // TODO: Do we handle multidimensional weights.
    // MD: no. We can only handle single dimensional weights.
    if (criteria_dim > 1){
        for (int wdim1 = 0; wdim1 < criteria_dim; wdim1++)
            for (int wdim2 = wdim1+1; wdim2 < criteria_dim; wdim2++)
                if (!mj_solution->criteriaHaveSamePartSizes(wdim1, wdim2)){
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

    mj_multivector_dim = coord_dim;
    for (int wdim = 0; wdim < criteria_dim; wdim++){
        if (!mj_uniform_weights[wdim]) mj_multivector_dim++;
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



        try{

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
    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon() * 10;
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



    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon() * 10;

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
        this->_EPSILON = numeric_limits<WT>::epsilon() * 10;
    }


    uMultiSortItem(IT index_ ,CT count_, WT *vals_){
        this->index = index_;
        this->count = count_;
        this->val = vals_;
        this->_EPSILON = numeric_limits<WT>::epsilon() * 10;
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
        double fEpsilon = numeric_limits<double>::epsilon() * 10;
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
    const partId_t *part_no_array,     // pqParts array
    int recursion_depth,
    size_t &num_global_parts,
    partId_t &total_num_cut, //how many cuts will be totally
    partId_t &total_num_part ,    //how many parts will be totally
    partId_t &max_num_part_along_dim ,         //maximum part count along a dimension.
    partId_t &last_dim_num_reduce_all ,    //estimate on #reduceAlls that can be done.
    partId_t &last_dim_num_part, //max no of parts that might occur
                                           //during the partition before the
                                           //last partitioning dimension.
    partId_t &max_num_cut_along_dim,
    size_t &max_num_total_part_along_dim
){
    if (part_no_array){
        for (int i = 0; i < recursion_depth; ++i){
            last_dim_num_reduce_all += total_num_part;
            total_num_part *= part_no_array[i];
            if(part_no_array[i] > max_num_part_along_dim) max_num_part_along_dim = part_no_array[i];
        }
        last_dim_num_part = total_num_part / part_no_array[recursion_depth-1];
        num_global_parts = total_num_part;
    } else {
        float fEpsilon = numeric_limits<float>::epsilon();
        partId_t futureNumParts = num_global_parts;

        for (int i = 0; i < recursion_depth; ++i){
            partId_t maxNoPartAlongI = getPartCount<partId_t>( futureNumParts,
                                        1.0f / (recursion_depth - i), fEpsilon);
            //cout << "futureNumParts:" << futureNumParts << "partArraySize:"
            // << partArraySize << " i:" << i << "maxNoPartAlongI:" <<i
            // maxNoPartAlongI<< endl;
            //partId_t maxNoPartAlongI = ceil(pow(futureNumParts,
            //1.0f / (coordDim - i)));// + 0.5f;
            if (maxNoPartAlongI > max_num_part_along_dim){
                max_num_part_along_dim = maxNoPartAlongI;
            }
            //cout << "i:" << i << "maxPartNo:" << maxPartNo<< endl;

            partId_t nfutureNumParts = futureNumParts / maxNoPartAlongI;
            if (futureNumParts % maxNoPartAlongI){
                ++nfutureNumParts;
            }
            futureNumParts = nfutureNumParts;
            //cout << "i:" << i << " fp:" << futureNumParts << endl;
        }
        total_num_part = num_global_parts;
        //estimate reduceAll Count here.
        //we find the upperbound instead.
        partId_t p = 1;
        for (int i = 0; i < recursion_depth; ++i){
            last_dim_num_reduce_all += p;
            p *= max_num_part_along_dim;
        }

        //cout << "maxPartNo:" << maxPartNo << endl;
        last_dim_num_part  = p / max_num_part_along_dim;
    }

    total_num_cut = total_num_part - 1;
    max_num_cut_along_dim = max_num_part_along_dim - 1;
    max_num_total_part_along_dim = max_num_part_along_dim + size_t(max_num_cut_along_dim);
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
    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon() * 10;

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
                    myRatios[i] = int ((myRatios[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL)
                    		/ pq_scalar_t(SIGNIFICANCE_MUL);
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
    pq_lno_t num_total_coords,
    pq_lno_t num_selected_coords,
    size_t num_target_part,
    int coord_dim,
    pq_scalar_t **mj_coordinates,
    pq_lno_t *initial_selected_coords_output_permutation,
    pq_lno_t *output_xadj,
    int recursion_depth,
    const partId_t *part_no_array
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
    pq_lno_t *partitionedPointCoordinates =  allocMemory< pq_lno_t>(num_total_coords);
    //initial configuration
    //set each pointer-i to i.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < static_cast<size_t>(num_total_coords); ++i){
        partitionedPointCoordinates[i] = initial_selected_coords_output_permutation[i];
    }

    //as output indices
    pq_lno_t *newpartitionedPointCoordinates = allocMemory< pq_lno_t>(num_total_coords);
#ifdef HAVE_ZOLTAN2_OMP
#ifdef FIRST_TOUCH
    firstTouch<pq_lno_t>(newpartitionedPointCoordinates, numLocalCoords);
#endif
#endif

    partId_t *partIds = NULL;
    ArrayRCP<partId_t> partId;
    if(num_total_coords > 0){
        partIds = allocMemory<partId_t>(num_total_coords);
    }

    //initially there is a single partition
    partId_t currentPartitionCount = 1;
    //single partition starts at index-0, and ends at numLocalCoords
    //inTotalCounts array holds the end points in partitionedPointCoordinates array
    //for each partition. Initially sized 1, and single element is set to numLocalCoords.
    pq_lno_t *inTotalCounts = allocMemory<pq_lno_t>(1);
    inTotalCounts[0] = static_cast<pq_lno_t>(num_selected_coords);//the end of the initial partition is the end of coordinates.
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
    getPartSpecifications <partId_t>( part_no_array, recursion_depth, num_target_part,
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



    pq_scalar_t _EPSILON = numeric_limits<pq_scalar_t>::epsilon() * 10;
    //partId_t partIndexBegin = 0;
    partId_t futurePartNumbers = totalPartCount;

    vector<partId_t> *currentPartitions = new vector<partId_t> ();
    vector<partId_t> *newFuturePartitions = new vector<partId_t> ();
    newFuturePartitions->push_back(num_target_part);
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > t1;
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > t2;


    for (int i = 0; i < recursion_depth; ++i){

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
                part_no_array, pAlongI, currentPartitions, newFuturePartitions,
                futurePartNumbers, currentPartitionCount, recursion_depth, i,
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
        int coordInd = i % coord_dim;
        pq_scalar_t * pqCoord = mj_coordinates[coordInd];
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
                            true, mj_coordinates, coord_dim, coordInd );
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
    for(pq_lno_t i = 0; i < num_total_coords; ++i){
        initial_selected_coords_output_permutation[i] = partitionedPointCoordinates[i];
    }

    for(size_t i = 0; i < num_target_part ; ++i){
        //cout << "i:" << i << endl;
        output_xadj[i] = inTotalCounts[i];
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
{}



/*! \brief Multi Jagged coordinate partitioning algorithm.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
class AlgMJ{
private:

    RCP<const Environment> mj_env; //the environment object
    RCP<Comm<int> > mj_problemComm; //initial comm object
    //RCP<const CoordinateModel<typename Adapter::base_adapter_t> > mj_coords; //coordinate adapter
    //RCP<PartitioningSolution<Adapter> > mj_solution; //solution object

    double imbalance_tolerance; //input imbalance tolerance.
    partId_t *part_no_array; //input part array specifying num part to divide along each dim.
    int recursion_depth; //the number of steps that partitioning will be solved in.
    int coord_dim, weight_dim; //coordinate and weight dimensions.

    size_t initial_num_loc_coords; //initial num local coords.
    global_size_t initial_num_glob_coords; //initial num global coords.

    pq_lno_t num_local_coords; //number of local coords.
    pq_gno_t num_global_coords; //number of global coords.

    pq_scalar_t **mj_coordinates; //two dimension coordinate array
    pq_scalar_t **mj_weights; //two dimension weight array
    bool *mj_uniform_parts; //if the target parts are uniform
    pq_scalar_t **mj_part_sizes; //target part weight sizes.
    bool *mj_uniform_weights; //if the coordinates have uniform weights.

    ArrayView<const pq_gno_t> mj_gnos; //global ids of the coordinates, comes from the input
    size_t num_global_parts; //the targeted number of parts


    pq_gno_t *initial_mj_gnos; //initial global ids of the coordinates.
    pq_gno_t *current_mj_gnos; //current global ids of the coordinates, might change during migration.
    int *owner_of_coordinate; //the actual processor owner of the coordinate, to track after migrations.

    pq_lno_t *coordinate_permutations; //permutation of coordinates, for partitioning.
    pq_lno_t *new_coordinate_permutations; //permutation work array.
    partId_t *assigned_part_ids; //the part ids assigned to coordinates.

    pq_lno_t *part_xadj; //beginning and end of each part.
    pq_lno_t *new_part_xadj; // work array for beginning and end of each part.

    //get mj specific parameters.
    bool distribute_points_on_cut_lines; //if partitioning can distribute points on same coordiante to different parts.
    int max_concurrent_part_calculation; // how many parts we can calculate concurrently.

    int mj_run_as_rcb; //if this is set, then recursion depth is adjusted to its maximum value.
    int mj_user_recursion_depth; //the recursion depth value provided by user.
    int mj_keep_part_boxes; //if the boxes need to be kept.

    int check_migrate_avoid_migration_option; //whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
    pq_scalar_t minimum_migration_imbalance; //when MJ decides whether to migrate, the minimum imbalance for migration.
    int num_threads; //num threads

    partId_t total_num_cut ; //how many cuts will be totally
    partId_t total_num_part;    //how many parts will be totally

    partId_t max_num_part_along_dim ;         //maximum part count along a dimension.
    partId_t max_num_cut_along_dim; //maximum cut count along a dimension.
    size_t max_num_total_part_along_dim; //maximum part+cut count along a dimension.


    partId_t total_dim_num_reduce_all;    //estimate on #reduceAlls can be done.
    partId_t last_dim_num_part; //max no of parts that might occur
                                              //during the partition before the
                                              //last partitioning dimension.

    RCP<Comm<int> > comm; //comm object than can be altered during execution
    float fEpsilon; //epsilon for float
    pq_scalar_t sEpsilon; //epsilon for pq_scalar_t


    pq_scalar_t maxScalar_t; //max possible scalar
    pq_scalar_t minScalar_t; //min scalar

    pq_scalar_t *all_cut_coordinates;
    pq_scalar_t *max_min_coords;
    pq_scalar_t *process_cut_line_weight_to_put_left; //how much weight should a MPI put left side of the each cutline
    pq_scalar_t **thread_cut_line_weight_to_put_left; //how much weight percentage should each thread in MPI put left side of the each outline

    // work array to manipulate coordinate of cutlines in different iterations.
    //necessary because previous cut line information is used for determining
    //the next cutline information. therefore, cannot update the cut work array
    //until all cutlines are determined.
    pq_scalar_t *cut_coordinates_work_array;

    //cumulative part weight array.
    pq_scalar_t *target_part_weights;

    pq_scalar_t *cut_upper_bound_coordinates ;  //upper bound coordinate of a cut line
    pq_scalar_t *cut_lower_bound_coordinates ;  //lower bound coordinate of a cut line
    pq_scalar_t *cut_lower_bound_weights ;  //lower bound weight of a cut line
    pq_scalar_t *cut_upper_bound_weights ;  //upper bound weight of a cut line

    pq_scalar_t *process_local_min_max_coord_total_weight ; //combined array to exchange the min and max coordinate, and total weight of part.
    pq_scalar_t *global_min_max_coord_total_weight ;//global combined array with the results for min, max and total weight.

    //isDone is used to determine if a cutline is determined already.
    //If a cut line is already determined, the next iterations will skip this cut line.
    bool *is_cut_line_determined;
    //my_incomplete_cut_count count holds the number of cutlines that have not been finalized for each part
    //when concurrentPartCount>1, using this information, if my_incomplete_cut_count[x]==0, then no work is done for this part.
    partId_t *my_incomplete_cut_count;
    //local part weights of each thread.
    double **thread_part_weights;
    //the work manupulation array for partweights.
    double **thread_part_weight_work;

    //thread_cut_left_closest_point to hold the closest coordinate to a cutline from left (for each thread).
    pq_scalar_t **thread_cut_left_closest_point;
    //thread_cut_right_closest_point to hold the closest coordinate to a cutline from right (for each thread)
    pq_scalar_t **thread_cut_right_closest_point;

    //to store how many points in each part a thread has.
    pq_lno_t **thread_point_counts;

    pq_scalar_t *process_rectelinear_cut_weight;
    pq_scalar_t *global_rectelinear_cut_weight;


    //for faster communication, concatanation of
    //totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
    //leftClosest distances sized P-1, since P-1 cut lines
    //rightClosest distances size P-1, since P-1 cut lines.
    pq_scalar_t *total_part_weight_left_right_closests ;
    pq_scalar_t *global_total_part_weight_left_right_closests;

    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > kept_boxes;
    int myRank, myActualRank; //processor rank, and initial rank

    /* \brief Either the mj array (part_no_array) or num_global_parts should be provided in
     * the input. part_no_array takes
     * precedence if both are provided.
     * Depending on these parameters, total cut/part number,
     * maximum part/cut number along a dimension, estimated number of reduceAlls,
     * and the number of parts before the last dimension is calculated.
     * */
    void set_part_specifications();

    /* \brief Tries to determine the part number for current dimension,
     * by trying to make the partitioning as square as possible.
     * \param num_total_future how many more partitionings are required.
     * \param root how many more recursion depth is left.
     */
    inline partId_t get_part_count(
    		partId_t num_total_future,
    		double root);

    /* \brief Allocates the all required memory for the mj partitioning algorithm.
     *
     */
    void allocate_set_work_memory();


    /* \brief for part communication we keep track of the box boundaries.
     * This is performed when either asked specifically, or when geometric mapping is performed afterwards.
     * This function initializes a single box with all global min and max coordinates.
     * \param initial_partitioning_boxes the input and output vector for boxes.
     */
    void init_part_boxes(
    		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > & outPartBoxes);

    /* \brief Function returns how many parts that will be obtained after this dimension partitioning.
     * It sets how many parts each current part will be partitioned into in this dimension to num_partitioning_in_current_dim vector,
     * sets how many total future parts each obtained part will be partitioned into in next_future_num_parts_in_parts vector,
     * If part boxes are kept, then sets initializes the output_part_boxes as its ancestor.
     *
     *  \param num_partitioning_in_current_dim: output. How many parts each current part will be partitioned into.
     *  \param future_num_part_in_parts: input, how many future parts each current part will be partitioned into.
     *  \param next_future_num_parts_in_parts: output, how many future parts each obtained part will be partitioned into.
     *  \param future_num_parts: output, max number of future parts that will be obtained from a single
     *  \param current_num_parts: input, how many parts are there currently.
     *  \param current_iteration: input, current dimension iteration number.
     *  \param input_part_boxes: input, if boxes are kept, current boxes.
     *  \param output_part_boxes: output, if boxes are kept, the initial box boundaries for obtained parts.
     */
    partId_t update_part_num_arrays(
    		vector <partId_t> &num_partitioning_in_current_dim, //assumes this vector is empty.
    		vector<partId_t> *future_num_part_in_parts,
    		vector<partId_t> *next_future_num_parts_in_parts, //assumes this vector is empty.
    		partId_t &future_num_parts,
    		partId_t current_num_parts,
    		int current_iteration,
    		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > input_part_boxes,
    		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > output_part_boxes);

    /*! \brief Function to determine the local minimum and maximum coordinate, and local total weight
     * in the given set of local points.
     * \param coordinate_begin_index is the start index of the given partition on partitionedPointPermutations.
     * \param coordinate_end_index is the end index of the given partition on partitionedPointPermutations.
     * \param mj_current_coordinate_permutations is the permutation array that point to the actual coordinate index. Sized as numLocalCoords.
     * \param mj_current_dim_coords float-like array representing the coordinates in a single dimension. Sized as numLocalCoords.
     * \param min_coordinate is the output to represent the local minimumCoordinate in  given range of coordinates.
     * \param max_coordinate is the output to represent the local maximum coordinate in the given range of coordinates.
     * \param total_weight is the output to represent the local total weight in the coordinate in the given range of coordinates.
     *
     */
    void mj_get_local_min_max_coord_totW(
    		pq_lno_t coordinate_begin_index,
    		pq_lno_t coordinate_end_index,
    		pq_lno_t *mj_current_coordinate_permutations,
    		pq_scalar_t *mj_current_dim_coords,
    		pq_scalar_t &min_coordinate,
    		pq_scalar_t &max_coordinate,
    		pq_scalar_t &total_weight);

    /*! \brief Function that reduces global minimum and maximum coordinates with global total weight from given local arrays.
     * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
     * \param local_min_max_total is the array holding local min and max coordinate values with local total weight.
     * First current_concurrent_num_parts entries are minimums of the parts, next current_concurrent_num_parts entries are max, and then the total weights.
     * \param global_min_max_total is the output array holding global min and global coordinate values with global total weight.
     * The structure is same as local_min_max_total.
     */
    void mj_get_global_min_max_coord_totW(
        partId_t current_concurrent_num_parts,
        pq_scalar_t *local_min_max_total,
        pq_scalar_t *global_min_max_total);

    /*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
     * \param min_coord minimum coordinate in the range.
     * \param max_coord maximum coordinate in the range.
     *
     * \param num_cuts holds the number of cuts in the current partitioning dimension.
     * \param global_weight holds the global total weight in the current part.
     *
     * \param initial_cut_coords is the output array for the initial cut lines.
     * \param target_part_weights is the output array holding the cumulative ratios of parts in current partitioning.
     * For partitioning to 4 uniformly, target_part_weights will be (0.25 * globalTotalWeight, 0.5 *globalTotalWeight , 0.75 * globalTotalWeight, globalTotalWeight).
     *
     * \param future_num_part_in_parts is the vector that holds how many more parts each part will be divided into more
     * for the parts at the beginning of this coordinate partitioning
     * \param next_future_num_parts_in_parts is the vector that holds how many more parts each part will be divided into more
     * for the parts that will be obtained at the end of this coordinate partitioning.
     * \param concurrent_current_part is the index of the part in the future_num_part_in_parts vector.
     * \param obtained_part_index holds the amount of shift in the next_future_num_parts_in_parts for the output parts.
     */
    void mj_get_initial_cut_coords_target_weights(
        pq_scalar_t min_coord,
        pq_scalar_t max_coord,
        partId_t num_cuts/*p-1*/ ,
        pq_scalar_t global_weight,
        pq_scalar_t *initial_cut_coords /*p - 1 sized, coordinate of each cut line*/,
        pq_scalar_t *target_part_weights /*cumulative weights, at left side of each cut line. p-1 sized*/,

        vector <partId_t> *future_num_part_in_parts, //the vecto
        vector <partId_t> *next_future_num_parts_in_parts,
        partId_t concurrent_current_part,
        partId_t obtained_part_index);

    /*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
     * \param max_coordinate maximum coordinate in the range.
     * \param min_coordinate minimum coordinate in the range.
     *
     * \param concurrent_current_part_index is the index of the part in the inTotalCounts vector.
     * \param coordinate_begin_index holds the beginning of the coordinates in current part.
     * \param coordinate_end_index holds end of the coordinates in current part.
     * \param mj_current_coordinate_permutations is the permutation array, holds the real indices of coordinates on pqCoord array.
     * \param mj_current_dim_coords is the 1D array holding the coordinates.
     * \param mj_part_ids is the array holding the partIds of each coordinate.
     * \param partition_count is the number of parts that the current part will be partitioned into.
     */
    void set_initial_coordinate_parts(
        pq_scalar_t &max_coordinate,
        pq_scalar_t &min_coordinate,
        partId_t &concurrent_current_part_index,
        pq_lno_t coordinate_begin_index,
        pq_lno_t coordinate_end_index,
        pq_lno_t *mj_current_coordinate_permutations,
        pq_scalar_t *mj_current_dim_coords,
        partId_t *mj_part_ids,
        partId_t &partition_count);

    /*! \brief Function that is responsible from 1D partitioning of the given range of coordinates.
     * \param pqJagged_coordinates is 1 dimensional array holding coordinate values.
     * \param imbalanceTolerance is the maximum allowed imbalance ratio.
     * \param current_work_part is the beginning index of concurrentPartCount parts.
     * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
     * \param current_cut_coordinates is the array holding the coordinates of the cut.
     * \param total_incomplete_cut_count is the number of cut lines whose positions should be calculated.
     * \param num_partitioning_in_current_dim is the vector that holds how many parts each part will be divided into.
     *
     */
    void mj_1D_part(
        pq_scalar_t *mj_current_dim_coords,
        pq_scalar_t imbalanceTolerance,
        partId_t current_work_part,
        partId_t current_concurrent_num_parts,
        pq_scalar_t *current_cut_coordinates,
        partId_t total_incomplete_cut_count,
        vector <partId_t> &num_partitioning_in_current_dim);

    /*! \brief Function that calculates the weights of each part according to given part cut coordinates.
     * Function is called inside the parallel region. Thread specific work arrays are provided
     * as function parameter.
     *
     * \param total_part_count is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
     * \param num_cuts is the number of cut lines. P - 1.
     * \param max_coord is the maximum coordinate in the part.
     * \param min_coord is the min coordinate in the part.
     * \param coordinate_begin_index is the index of the first coordinate in current part.
     * \param coordinate_end_index is the index of the last coordinate in current part.
     * \param mj_current_dim_coords is 1 dimensional array holding coordinate values.
     *
     * \param temp_current_cut_coords is the array holding the coordinates of each cut line. Sized P - 1.
     * \param current_cut_status is the boolean array to determine if the correct position for a cut line is found.
     * \param my_current_part_weights is the array holding the part weights for the calling thread.
     * \param my_current_left_closest is the array holding the coordinate of the closest points to the cut lines from left for the calling thread..
     * \param my_current_right_closest is the array holding the coordinate of the closest points to the cut lines from right for the calling thread.
     * \param partIds is the array that holds the part ids of the coordinates
     */
    void mj_1D_part_get_thread_part_weights(
        size_t total_part_count,
        partId_t num_cuts,
        pq_scalar_t max_coord,
        pq_scalar_t min_coord,
        pq_lno_t coordinate_begin_index,
        pq_lno_t coordinate_end_index,
        pq_scalar_t *mj_current_dim_coords,
        pq_scalar_t *temp_current_cut_coords,
        bool *current_cut_status,
        double *my_current_part_weights,
        pq_scalar_t *my_current_left_closest,
        pq_scalar_t *my_current_right_closest);

    /*! \brief Function that reduces the result of multiple threads
     * for left and right closest points and part weights in a single mpi process.
     *
     * \param num_partitioning_in_current_dim is the vector that holds the number of cut lines in current dimension for each part.
     * \param current_work_part holds the index of the first part (important when concurrent parts are used.)
     * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
     */
    void mj_accumulate_thread_results(
        const vector <partId_t> &num_partitioning_in_current_dim,
        partId_t current_work_part,
        partId_t current_concurrent_num_parts);

    /*! \brief Function that calculates the new coordinates for the cut lines.
     * Function is called inside the parallel region. Write the new cut coordinates
     * to new_current_cut_coordinates, and determines if the final position of a cut is found.
     *
     * \param num_total_part is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
     * \param num_cuts is the number of cut lines. P - 1.
     * \param max_coordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
     * \param min_coordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
     * \param global_total_weight is the global total weight in the current range of coordinates.
     * \param used_imbalance_tolerance is the maximum allowed imbalance ratio.
     *
     *
     * \param current_global_part_weights is the array holding the weight of parts. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
     * \param current_local_part_weights is the local totalweight of the processor.
     * \param current_part_target_weights are the desired cumulative part ratios, sized P.
     * \param current_cut_line_determined is the boolean array to determine if the correct position for a cut line is found.
     *
     * \param current_cut_coordinates is the array holding the coordinates of each cut line. Sized P - 1.
     * \param current_cut_upper_bounds is the array holding the upper bound coordinate for each cut line. Sized P - 1.
     * \param current_cut_lower_bounds is the array holding the lower bound coordinate for each cut line. Sized P - 1.
     * \param current_global_left_closest_points is the array holding the closest points to the cut lines from left.
     * \param current_global_right_closest_points is the array holding the closest points to the cut lines from right.
     * \param current_cut_lower_bound_weights is the array holding the weight of the parts at the left of lower bound coordinates.
     * \param current_cut_upper_weights is the array holding the weight of the parts at the left of upper bound coordinates.
     * \param new_current_cut_coordinates is the work array, sized P - 1.
     *
     * \param current_part_cut_line_weight_ratio holds how much weight of the coordinates on the cutline should be put on left side.
     * \param rectelinear_cut_count is the count of cut lines whose balance can be achived via distributing the points in same coordinate to different parts.
     * \param my_num_incomplete_cut is the number of cutlines whose position has not been determined yet. For K > 1 it is the count in a single part (whose cut lines are determined).
     */
    void mj_get_new_cut_coordinates(
        const size_t &num_total_part,
        const partId_t &num_cuts,
        const pq_scalar_t &max_coordinate,
        const pq_scalar_t &min_coordinate,
        const pq_scalar_t &global_total_weight,
        const pq_scalar_t &used_imbalance_tolerance,
        pq_scalar_t * current_global_part_weights,
        const pq_scalar_t * current_local_part_weights,
        const pq_scalar_t *current_part_target_weights,
        bool *current_cut_line_determined,
        pq_scalar_t *current_cut_coordinates,
        pq_scalar_t *current_cut_upper_bounds,
        pq_scalar_t *current_cut_lower_bounds,
        pq_scalar_t *current_global_left_closest_points,
        pq_scalar_t *current_global_right_closest_points,
        pq_scalar_t * current_cut_lower_bound_weights,
        pq_scalar_t * current_cut_upper_weights,
        pq_scalar_t *new_current_cut_coordinates,
        pq_scalar_t *current_part_cut_line_weight_to_put_left,
        partId_t *rectelinear_cut_count,
        partId_t &my_num_incomplete_cut);

    /*! \brief
     * Function that calculates the next pivot position,
     * according to given coordinates of upper bound and lower bound, the weights at upper and lower bounds, and the expected weight.
     * \param cut_upper_bound is the upper bound coordinate of the cut.
     * \param cut_lower_bound is the lower bound coordinate of the cut.
     * \param cut_upper_weight is the weights at the upper bound of the cut.
     * \param cut_lower_weight is the weights at the lower bound of the cut.
     * \param expected_weight is the expected weight that should be placed on the left of the cut line.
     */
    void mj_calculate_new_cut_position (
    	pq_scalar_t cut_upper_bound,
        pq_scalar_t cut_lower_bound,
        pq_scalar_t cut_upper_weight,
        pq_scalar_t cut_lower_weight,
        pq_scalar_t expected_weight,
        pq_scalar_t &new_cut_position);

    /*! \brief Function that determines the permutation indices of the coordinates.
     * \param num_parts is the number of parts.
     * \param mj_current_dim_coords is 1 dimensional array holding the coordinate values.
     * \param current_concurrent_cut_coordinate is 1 dimensional array holding the cut coordinates.
     * \param coordinate_begin is the start index of the given partition on partitionedPointPermutations.
     * \param coordinate_end is the end index of the given partition on partitionedPointPermutations.
     * \param used_local_cut_line_weight_to_left holds how much weight of the coordinates on the cutline should be put on left side.
     * \param used_thread_part_weight_work is the two dimensional array holding the weight of parts for each thread. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
     * \param out_part_xadj is the indices of coordinates calculated for the partition on next dimension. \param pqCoordDim is dimension of the input
     */
    void mj_create_new_partitions(
        partId_t num_parts,
        pq_scalar_t *mj_current_dim_coords,
        pq_scalar_t *current_concurrent_cut_coordinate,
        pq_lno_t coordinate_begin,
        pq_lno_t coordinate_end,
        pq_scalar_t *used_local_cut_line_weight_to_left,
        double **used_thread_part_weight_work,
        pq_lno_t *out_part_xadj);

    /*! \brief Function checks if should do migration or not.
     * It returns true to point that migration should be done when
     * -migration_reduce_all_population are higher than a predetermined value
     * -num_coords_for_last_dim_part that left for the last dimension partitioning is less than a predetermined value
     * -the imbalance of the processors on the parts are higher than given threshold.

     * \param input_num_parts is the number of parts when migration is called.
     * \param output_num_parts is the output number of parts after migration.
     * \param next_future_num_parts_in_parts is the number of total future parts each
     * part is partitioned into. This will be updated when migration is performed.
     * \param output_part_begin_index is the number that will be used as beginning part number
     * when final solution part numbers are assigned.
     * \param migration_reduce_all_population is the estimated total number of reduceall operations
     * multiplied with number of processors to be used for determining migration.
     *
     * \param num_coords_for_last_dim_part is the estimated number of points in each part,
     * when last dimension partitioning is performed.
     * \param iteration is the string that gives information about the dimension for printing purposes.
     * \param input_part_boxes is the array that holds the part boxes after the migration. (swapped)
     * \param output_part_boxes is the array that holds the part boxes before the migration. (swapped)
     *
     */
    bool mj_perform_migration(
        partId_t in_num_parts, //current umb parts
        partId_t &out_num_parts, //output umb parts.
        vector<partId_t> *next_future_num_parts_in_parts,
        partId_t &output_part_begin_index,
        size_t migration_reduce_all_population,
        pq_lno_t num_coords_for_last_dim_part,
        string iteration,
        RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > &input_part_boxes,
        RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > &output_part_boxes);

    /*! \brief Function fills up the num_points_in_all_processor_parts, so that
     * it has the number of coordinates in each processor of each part.
     * to access how many points processor i has on part j, num_points_in_all_processor_parts[i * num_parts + j].
     *
     * \param num_procs is the number of processor attending to migration operation.
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param num_points_in_all_processor_parts is the output array that holds
     * the number of coordinates in each part in each processor.
     */
    void get_processor_num_points_in_parts(
    		partId_t num_procs,
    		partId_t num_parts,
    		pq_gno_t *&num_points_in_all_processor_parts);

    /*! \brief Function checks if should do migration or not.
     * It returns true to point that migration should be done when
     * -migration_reduce_all_population are higher than a predetermined value
     * -num_coords_for_last_dim_part that left for the last dimension partitioning is less than a predetermined value
     * -the imbalance of the processors on the parts are higher than given threshold.
     * \param migration_reduce_all_population is the multiplication of the number of reduceall operations estimated and the number of processors.
     * \param num_coords_for_last_dim_part is the estimated number of coordinates in a part per processor in the last dimension partitioning.
     * \param num_procs is the number of processor attending to migration operation.
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param num_points_in_all_processor_parts is the input array that holds
     * the number of coordinates in each part in each processor.
     */
    bool mj_check_to_migrate(
    		size_t migration_reduce_all_population,
    		pq_lno_t num_coords_for_last_dim_part,
    		partId_t num_procs,
    		partId_t num_parts,
    		pq_gno_t *num_points_in_all_processor_parts);


    /*! \brief Function fills up coordinate_destionations is the output array
     * that holds which part each coordinate should be sent. In addition it calculates
     * the shift amount (output_part_numbering_begin_index) to be done when
     * final numberings of the parts are performed.
     *
     * \param num_points_in_all_processor_parts is the array holding the num points in each part in each proc.
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param num_procs is the number of processor attending to migration operation.

     * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
     * \param processor_ranks_for_subcomm is the ranks of the processors that will be in the subcommunicator with me.
     * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
     * \param out_num_part is the number of parts assigned to the process.
     * \param out_part_indices is the indices of the part to which the processor is assigned.
     * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
     * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
     */
    void mj_migration_part_proc_assignment(
    		pq_gno_t * num_points_in_all_processor_parts,
    		partId_t num_parts,
    		partId_t num_procs,
    		pq_lno_t *send_count_to_each_proc,
    		vector<partId_t> &processor_ranks_for_subcomm,
    		vector<partId_t> *next_future_num_parts_in_parts,
    		partId_t &out_num_part,
    		vector<partId_t> &out_part_indices,
    		partId_t &output_part_numbering_begin_index,
    		int *coordinate_destionations);

    /*! \brief Function that assigned the processors to parts, when there are more processors then parts.
     *	sets the destination of each coordinate in coordinate_destionations, also edits output_part_numbering_begin_index,
     *	and out_part_index, and returns the processor_ranks_for_subcomm which represents the ranks of the processors
     *	that will be used for creating the subcommunicator.
     *
     * \param num_points_in_all_processor_parts is the array holding the num points in each part in each proc.
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param num_procs is the number of processor attending to migration operation.

     * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
     * \param processor_ranks_for_subcomm is the ranks of the processors that will be in the subcommunicator with me.
     * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
     * \param out_part_index is the index of the part to which the processor is assigned.
     * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
     * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
     */
    void mj_assign_proc_to_parts(
    		pq_gno_t * num_points_in_all_processor_parts,
    		partId_t num_parts,
    		partId_t num_procs,
    		pq_lno_t *send_count_to_each_proc,
    		vector<partId_t> &processor_ranks_for_subcomm,
    		vector<partId_t> *next_future_num_parts_in_parts,
    		partId_t &out_part_index,
    		partId_t &output_part_numbering_begin_index,
    		int *coordinate_destionations);

    /*! \brief Function fills up coordinate_destionations is the output array
     * that holds which part each coordinate should be sent.
     *
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param num_procs is the number of processor attending to migration operation.
     * \param part_assignment_proc_begin_indices ([i]) points to the first processor index that part i will be sent to.
     * \param processor_chains_in_parts the array that holds the linked list structure, started from part_assignment_proc_begin_indices ([i]).
     * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
     * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
     */
    void assign_send_destinations(
    		partId_t num_parts,
    		partId_t *part_assignment_proc_begin_indices,
    		partId_t *processor_chains_in_parts,
    		pq_lno_t *send_count_to_each_proc,
    		int *coordinate_destionations);

    /*! \brief Function fills up coordinate_destionations is the output array
     * that holds which part each coordinate should be sent. In addition it calculates
     * the shift amount (output_part_numbering_begin_index) to be done when
     * final numberings of the parts are performed.
     *
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param sort_item_part_to_proc_assignment is the sorted parts with respect to the assigned processors.
     * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
     * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
     * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
     *
     */
    void assign_send_destionations2(
        partId_t num_parts,
        uSortItem<partId_t, partId_t> * sort_item_part_to_proc_assignment, //input sorted wrt processors
        int *coordinate_destionations,
        partId_t &output_part_numbering_begin_index,
        vector<partId_t> *next_future_num_parts_in_parts);

    /*! \brief Function fills up coordinate_destionations is the output array
     * that holds which part each coordinate should be sent. In addition it calculates
     * the shift amount (output_part_numbering_begin_index) to be done when
     * final numberings of the parts are performed.
     *
     * \param num_points_in_all_processor_parts is the array holding the num points in each part in each proc.
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param num_procs is the number of processor attending to migration operation.

     * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
     * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
     * \param out_num_part is the number of parts assigned to the process.
     * \param out_part_indices is the indices of the part to which the processor is assigned.
     * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
     * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
     */
    void mj_assign_parts_to_procs(
        pq_gno_t * num_points_in_all_processor_parts,
        partId_t num_parts,
        partId_t num_procs,
        pq_lno_t *send_count_to_each_proc, //output: sized nprocs, show the number of send point counts to each proc.
        vector<partId_t> *next_future_num_parts_in_parts,//input how many more partitions the part will be partitioned into.
        partId_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
        vector<partId_t> &out_part_indices, //output: the part indices which the processor is assigned to.
        partId_t &output_part_numbering_begin_index, //output: how much the part number should be shifted when setting the solution
        int *coordinate_destionations);

    /*! \brief Function fills up coordinate_destionations is the output array
     * that holds which part each coordinate should be sent. In addition it calculates
     * the shift amount (output_part_numbering_begin_index) to be done when
     * final numberings of the parts are performed.
     *
     *
     * \param num_procs is the number of processor attending to migration operation.
     * \param num_new_local_points is the output to represent the new number of local points.
     * \param iteration is the string for the current iteration.
     * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
     * \param num_parts is the number of parts that exist in the current partitioning.
     */
    void mj_migrate_coords(
        partId_t num_procs,
        pq_gno_t &num_new_local_points,
        string iteration,
        int *coordinate_destionations,
        partId_t num_parts);

    /*! \brief Function creates the new subcomminicator for the processors
     * given in processor_ranks_for_subcomm.
     *
     * \param processor_ranks_for_subcomm is the vector that has the ranks of
     * the processors that will be in the same group.
     */
    void create_sub_communicator(vector<partId_t> &processor_ranks_for_subcomm);


    /*! \brief Function writes the new permutation arrays after the migration.
     *
     * \param output_num_parts is the number of parts that is assigned to the processor.
     * \param num_parts is the number of parts right before migration.
     */
    void fill_permutation_array(
        partId_t output_num_parts,
        partId_t num_parts);

    /*! \brief Function checks if should do migration or not.
     * \param current_num_parts is the number of parts in the process.
     * \param output_part_begin_index is the number that will be used as beginning part number
     * \param output_part_boxes is the array that holds the part boxes
     * \param is_data_ever_migrated is the boolean value which is true
     * if the data is ever migrated during the partitioning.
     *
     */
    void set_final_parts(
    		partId_t current_num_parts,
    		partId_t output_part_begin_index,
    		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > output_part_boxes,
    		bool is_data_ever_migrated);
    /*! \brief Function frees all allocated work memory.
     */
    void free_work_memory();
    /*! \brief Function creates consistent chunks for task partitioning. Used only in the case of
     * sequential task partitioning, where consistent handle of the points on the cuts are required.
     *
     * \param num_parts is the number of parts.
     * \param mj_current_dim_coords is 1 dimensional array holding the coordinate values.
     * \param current_concurrent_cut_coordinate is 1 dimensional array holding the cut coordinates.
     * \param coordinate_begin is the start index of the given partition on partitionedPointPermutations.
     * \param coordinate_end is the end index of the given partition on partitionedPointPermutations.
     * \param used_local_cut_line_weight_to_left holds how much weight of the coordinates on the cutline should be put on left side.
     *
     * \param out_part_xadj is the indices of begginning and end of the parts in the output partition.
     * \param coordInd is the index according to which the partitioning is done.
     */
    void create_consistent_chunks(
        partId_t num_parts,
        pq_scalar_t *mj_current_dim_coords,
        pq_scalar_t *current_concurrent_cut_coordinate,
        pq_lno_t coordinate_begin,
        pq_lno_t coordinate_end,
        pq_scalar_t *used_local_cut_line_weight_to_left,
        pq_lno_t *out_part_xadj,
        int coordInd);
public:
    AlgMJ();

    /*! \brief Multi Jagged  coordinate partitioning algorithm.
     *
     *  \param env   library configuration and problem parameters
     *  \param problemComm the communicator for the problem
     *  \param imbalance_tolerance : the input provided imbalance tolerance.
     *  \param num_global_parts: number of target global parts.
     *  \param part_no_array: part no array, if provided this will be used for partitioning.
     *  \param recursion_depth: if part no array is provided, it is the length of part no array,
     *  						if part no is not provided than it is the number of steps that algorithm will divide into num_global_parts parts.
     *
     *  \param coord_dim: coordinate dimension
     *  \param num_local_coords: number of local coordinates
     *  \param num_global_coords: number of global coordinates
     *  \param initial_mj_gnos: the list of initial global id's
     *  \param mj_coordinates: the two dimensional coordinate array.
     *
     *  \param weight_dim: weight dimension
     *  \param mj_uniform_weights: if weight dimension [i] has uniform weight or not.
     *  \param mj_weights: the two dimensional array for weights in every weight dimension
     *  \param mj_uniform_parts: if the target partitioning aims uniform parts
     *  \param mj_part_sizes: if the target partitioning does not aim uniform parts, then weight of each part.
     *
     *  \param result_assigned_part_ids: Output - 1D pointer, should be provided as null.
     *  			the result partids corresponding to the coordinates given in result_mj_gnos.
     *  \param result_mj_gnos: Output - 1D pointer, should be provided as null.
     *  			the result coordinate global id's corresponding to the part_ids array.
     *
     */
    void multi_jagged_part(
    		const RCP<const Environment> &env,
        	RCP<Comm<int> > &problemComm,

        	double imbalance_tolerance,
        	size_t num_global_parts,
        	partId_t *part_no_array,
        	int recursion_depth,

        	int coord_dim,
        	pq_lno_t num_local_coords,
        	pq_gno_t num_global_coords,
        	const pq_gno_t *initial_mj_gnos,
        	pq_scalar_t **mj_coordinates,

        	int weight_dim,
        	bool *mj_uniform_weights,
        	pq_scalar_t **mj_weights,
        	bool *mj_uniform_parts,
        	pq_scalar_t **mj_part_sizes,

        	partId_t *&result_assigned_part_ids,
        	pq_gno_t *&result_mj_gnos

    		);
    /*! \brief Multi Jagged  coordinate partitioning algorithm.
     *
     *  \param distribute_points_on_cut_lines_ :  if partitioning can distribute points on same coordinate to different parts.
     *  \param max_concurrent_part_calculation_ : how many parts we can calculate concurrently.
     *  \param check_migrate_avoid_migration_option_ : whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
     *  \param minimum_migration_imbalance_  : when MJ decides whether to migrate, the minimum imbalance for migration.
     */
    void set_partitioning_parameters(
    		bool distribute_points_on_cut_lines_,
    		int max_concurrent_part_calculation_,
    		int check_migrate_avoid_migration_option_,
    		pq_scalar_t minimum_migration_imbalance_);
    /*! \brief Function call, if the part boxes are intended to be kept.
     *
     */
    void set_to_keep_part_boxes();
    /*! \brief Function returns the part boxes stored.
     * returns null if boxes are not stored, and prints warning message.
     * User needs to call set_to_keep_part_boxes() before partitioning if part boxes are needed.
     */
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > get_part_boxes();

    /*! \brief Special function for partitioning for task mapping.
     * Runs sequential, and performs deterministic partitioning for the
     * partitioning the points along a cutline.
     *
     *  \param env library configuration and problem parameters
     *  \param num_total_coords number of total coordinates
     *  \param num_selected_coords : the number of selected coordinates. This is to set,
     *  							if there are n processors, but only m<n processors
     *  							are selected for mapping.
     *
     *  \param num_target_part: number of target global parts.
     *  \param coord_dim_: coordinate dimension for coordinates
     *  \param mj_coordinates_: the coordinates
     *
     *  \param inital_adjList_output_adjlist: Array allocated by caller, in the size of num_total_coords,
     *  						first num_selected_coords elements should list the indices of the selected processors.
     *  						This is output for output permutation array.
     *  \param output_xadj: The output part xadj array, pointing beginning and end of each part on
     *  	output permutation array (inital_adjList_output_adjlist).
     *
     *  \param rd: recursion depth
     *  \param part_no_array_: possibly null part_no_array, specifying how many parts each should be divided during partitioning.
     */
    void sequential_task_partitioning(
        const RCP<const Environment> &env,
        pq_lno_t num_total_coords,
        pq_lno_t num_selected_coords,
        size_t num_target_part,
        int coord_dim,
        pq_scalar_t **mj_coordinates,
        pq_lno_t *initial_selected_coords_output_permutation,
        pq_lno_t *output_xadj,
        int recursion_depth,
        const partId_t *part_no_array);

};

/*! \brief Special function for partitioning for task mapping.
 * Runs sequential, and performs deterministic partitioning for the
 * partitioning the points along a cutline.
 *
 *  \param env library configuration and problem parameters
 *  \param num_total_coords number of total coordinates
 *  \param num_selected_coords : the number of selected coordinates. This is to set,
 *  							if there are n processors, but only m<n processors
 *  							are selected for mapping.
 *
 *  \param num_target_part: number of target global parts.
 *  \param coord_dim_: coordinate dimension for coordinates
 *  \param mj_coordinates_: the coordinates
 *
 *  \param inital_adjList_output_adjlist: Array allocated by caller, in the size of num_total_coords,
 *  						first num_selected_coords elements should list the indices of the selected processors.
 *  						This is output for output permutation array.
 *  \param output_xadj: The output part xadj array, pointing beginning and end of each part on
 *  	output permutation array (inital_adjList_output_adjlist).
 *
 *  \param rd: recursion depth
 *  \param part_no_array_: possibly null part_no_array, specifying how many parts each should be divided during partitioning.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::sequential_task_partitioning(
    const RCP<const Environment> &env,
    pq_lno_t num_total_coords,
    pq_lno_t num_selected_coords,
    size_t num_target_part,
    int coord_dim_,
    pq_scalar_t **mj_coordinates_,
    pq_lno_t *inital_adjList_output_adjlist,
    pq_lno_t *output_xadj,
    int rd,
    const partId_t *part_no_array_
){

	this->mj_env = env;
	const RCP<Comm<int> > commN;
	this->comm = this->mj_problemComm =  Teuchos::rcp_const_cast<Comm<int> >
	(Teuchos::DefaultComm<int>::getDefaultSerialComm(commN));
	this->myActualRank = this->myRank = 1;


#ifdef HAVE_ZOLTAN2_OMP
	int actual_num_threads = omp_get_num_threads();
	omp_set_num_threads(1);
#endif

    cout << "new mapper" << endl;

    //weights are uniform for task mapping

    //parts are uniform for task mapping
    //as input indices.


    this->imbalance_tolerance = 0;
    this->num_global_parts = num_target_part;
    this->part_no_array = (partId_t *)part_no_array_;
    this->recursion_depth = rd;


    this->coord_dim = coord_dim_;
    this->num_local_coords = num_total_coords;
    this->num_global_coords = num_total_coords;
    this->mj_coordinates = mj_coordinates_;  //will copy the memory to this->mj_coordinates.


    ////temporary memory. It is not used for this, but the functions reuqire these to be allocated.
    ////will copy the memory to this->current_mj_gnos[j].
    this->initial_mj_gnos = allocMemory<pq_gno_t>(this->num_local_coords);


    this->weight_dim = 0;
    bool *tmp_mj_uniform_weights = new bool[1];
    this->mj_uniform_weights = tmp_mj_uniform_weights ;
    this->mj_uniform_weights[0] = true;


    pq_scalar_t **tmp_mj_weights = new pq_scalar_t *[1];
    this->mj_weights = tmp_mj_weights; //will copy the memory to this->mj_weights

    bool *tmp_mj_uniform_parts = new bool[1];
    this->mj_uniform_parts = tmp_mj_uniform_parts;
    this->mj_uniform_parts[0] = true;

    pq_scalar_t **tmp_mj_part_sizes = new pq_scalar_t * [1];
    this->mj_part_sizes = tmp_mj_part_sizes;
    this->mj_part_sizes[0] = NULL;

    this->num_threads = 1;
    this->set_part_specifications();

    this->allocate_set_work_memory();
    //the end of the initial partition is the end of coordinates.
    this->part_xadj[0] = static_cast<pq_lno_t>(num_selected_coords);
    for(size_t i = 0; i < static_cast<size_t>(num_total_coords); ++i){
        this->coordinate_permutations[i] = inital_adjList_output_adjlist[i];
    }

    partId_t current_num_parts = 1;

    pq_scalar_t *current_cut_coordinates =  this->all_cut_coordinates;




    partId_t future_num_parts = this->total_num_cut;

    vector<partId_t> *future_num_part_in_parts = new vector<partId_t> ();
    vector<partId_t> *next_future_num_parts_in_parts = new vector<partId_t> ();
    next_future_num_parts_in_parts->push_back(this->num_global_parts);
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > t1;
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > t2;


    for (int i = 0; i < this->recursion_depth; ++i){

        //partitioning array. size will be as the number of current partitions and this
        //holds how many parts that each part will be in the current dimension partitioning.
        vector <partId_t> num_partitioning_in_current_dim;


        //number of parts that will be obtained at the end of this partitioning.
        //future_num_part_in_parts is as the size of current number of parts.
        //holds how many more parts each should be divided in the further
        //iterations. this will be used to calculate num_partitioning_in_current_dim,
        //as the number of parts that the part will be partitioned
        //in the current dimension partitioning.

        //next_future_num_parts_in_parts will be as the size of outnumParts,
        //and this will hold how many more parts that each output part
        //should be divided. this array will also be used to determine the weight ratios
        //of the parts.
        //swap the arrays to use iteratively..
        vector<partId_t> *tmpPartVect= future_num_part_in_parts;
        future_num_part_in_parts = next_future_num_parts_in_parts;
        next_future_num_parts_in_parts = tmpPartVect;

        //clear next_future_num_parts_in_parts array as
        //getPartitionArrays expects it to be empty.
        //it also expects num_partitioning_in_current_dim to be empty as well.
        next_future_num_parts_in_parts->clear();


        //returns the total number of output parts for this dimension partitioning.
        partId_t output_part_count_in_dimension =
        		this->update_part_num_arrays(
        				num_partitioning_in_current_dim,
        				future_num_part_in_parts,
        				next_future_num_parts_in_parts,
        				future_num_parts,
        				current_num_parts,
        				i,
        				t1,
        				t2);

        //if the number of obtained parts equal to current number of parts,
        //skip this dimension. For example, this happens when 1 is given in the input
        //part array is given. P=4,5,1,2
        if(output_part_count_in_dimension == current_num_parts) {
            tmpPartVect= future_num_part_in_parts;
            future_num_part_in_parts = next_future_num_parts_in_parts;
            next_future_num_parts_in_parts = tmpPartVect;
            continue;
        }

        //get the coordinate axis along which the partitioning will be done.
        int coordInd = i % this->coord_dim;
        pq_scalar_t * mj_current_dim_coords = this->mj_coordinates[coordInd];
        //convert i to string to be used for debugging purposes.
        string istring = toString<int>(i);

        //alloc Memory to point the indices
        //of the parts in the permutation array.
        this->new_part_xadj = allocMemory<pq_lno_t>(output_part_count_in_dimension);

        //the index where in the outtotalCounts will be written.
        partId_t output_part_index = 0;
        //whatever is written to outTotalCounts will be added with previousEnd
        //so that the points will be shifted.
        partId_t output_coordinate_end_index = 0;

        partId_t current_work_part = 0;
        partId_t current_concurrent_num_parts = min(current_num_parts - current_work_part,
                                         this->max_concurrent_part_calculation);

        partId_t obtained_part_index = 0;

        //run for all available parts.
        for (; current_work_part < current_num_parts;
                     current_work_part += current_concurrent_num_parts){

            current_concurrent_num_parts = min(current_num_parts - current_work_part,
            this->max_concurrent_part_calculation);

            partId_t actual_work_part_count = 0;
            //initialization for 1D partitioning.
            //get the min and max coordinates of each part
            //together with the part weights of each part.
            for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                partId_t current_work_part_in_concurrent_parts = current_work_part + kk;

                //if this part wont be partitioned any further
                //dont do any work for this part.
                if (num_partitioning_in_current_dim[current_work_part_in_concurrent_parts] == 1){
                    continue;
                }
                ++actual_work_part_count;
                pq_lno_t coordinate_end_index= this->part_xadj[current_work_part_in_concurrent_parts];
                pq_lno_t coordinate_begin_index = current_work_part_in_concurrent_parts==0 ? 0: this->part_xadj[current_work_part_in_concurrent_parts -1];


                this->mj_get_local_min_max_coord_totW(
                		coordinate_begin_index,
                		coordinate_end_index,
                		this->coordinate_permutations,
                		mj_current_dim_coords,
                		this->process_local_min_max_coord_total_weight[kk], //min coordinate
                        this->process_local_min_max_coord_total_weight[kk + current_concurrent_num_parts], //max coordinate
                        this->process_local_min_max_coord_total_weight[kk + 2*current_concurrent_num_parts] //total weight);
                );
            }


            //1D partitioning
            if (actual_work_part_count > 0){
                //obtain global Min max of the part.
            	this->mj_get_global_min_max_coord_totW(
                		current_concurrent_num_parts,
                		this->process_local_min_max_coord_total_weight,
                		this->global_min_max_coord_total_weight);

                //represents the total number of cutlines
                //whose coordinate should be determined.
                partId_t total_incomplete_cut_count = 0;

                //Compute weight ratios for parts & cuts:
                //e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
                //part0  cut0  part1 cut1 part2 cut2 part3
                partId_t concurrent_part_cut_shift = 0;
                partId_t concurrent_part_part_shift = 0;
                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    pq_scalar_t min_coordinate = this->global_min_max_coord_total_weight[kk];
                    pq_scalar_t max_coordinate = this->global_min_max_coord_total_weight[kk +
                                                     current_concurrent_num_parts];
                    pq_scalar_t global_total_weight =
                    					this->global_min_max_coord_total_weight[kk +
                                                     2 * current_concurrent_num_parts];

                    partId_t concurrent_current_part_index = current_work_part + kk;

                    partId_t partition_count = num_partitioning_in_current_dim[concurrent_current_part_index];

                    pq_scalar_t *usedCutCoordinate = current_cut_coordinates + concurrent_part_cut_shift;
                    pq_scalar_t *current_target_part_weights = this->target_part_weights +
                                                                     concurrent_part_part_shift;
                    //shift the usedCutCoordinate array as noCuts.
                    concurrent_part_cut_shift += partition_count - 1;
                    //shift the partRatio array as noParts.
                    concurrent_part_part_shift += partition_count;

                    //calculate only if part is not empty,
                    //and part will be further partitioend.
                    if(partition_count > 1 && min_coordinate <= max_coordinate){

                        //increase allDone by the number of cuts of the current
                        //part's cut line number.
                        total_incomplete_cut_count += partition_count - 1;
                        //set the number of cut lines that should be determined
                        //for this part.
                        this->my_incomplete_cut_count[kk] = partition_count - 1;

                        //get the target weights of the parts.
                        this->mj_get_initial_cut_coords_target_weights(
                            min_coordinate,
                            max_coordinate,
                            partition_count - 1,
                            global_total_weight,
                            usedCutCoordinate,
                            current_target_part_weights,
                            future_num_part_in_parts,
                            next_future_num_parts_in_parts,
                            concurrent_current_part_index,
                            obtained_part_index);

                        pq_lno_t coordinate_end_index= this->part_xadj[concurrent_current_part_index];
                        pq_lno_t coordinate_begin_index = concurrent_current_part_index==0 ? 0: this->part_xadj[concurrent_current_part_index -1];

                        //get the initial estimated part assignments of the coordinates.
                        this->set_initial_coordinate_parts(
                            max_coordinate,
                            min_coordinate,
                            concurrent_current_part_index,
                            coordinate_begin_index, coordinate_end_index,
                            this->coordinate_permutations,
                            mj_current_dim_coords,
                            this->assigned_part_ids,
                            partition_count);
                    }
                    else {
                        // e.g., if have fewer coordinates than parts, don't need to do next dim.
                        this->my_incomplete_cut_count[kk] = 0;
                    }
                    obtained_part_index += partition_count;
                }

                //used imbalance, it is always 0, as it is difficult to estimate a range.
                pq_scalar_t used_imbalance = 0;

                // Determine cut lines for k parts here.
                this->mj_1D_part(
                    mj_current_dim_coords,
                    used_imbalance,
                    current_work_part,
                    current_concurrent_num_parts,
                    current_cut_coordinates,
                    total_incomplete_cut_count,
                    num_partitioning_in_current_dim);
            }

            //create part chunks
            {

                partId_t output_array_shift = 0;
                partId_t cut_shift = 0;
                size_t tlr_shift = 0;
                size_t partweight_array_shift = 0;

                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    partId_t current_concurrent_work_part = current_work_part + kk;
                    partId_t num_parts = num_partitioning_in_current_dim[current_concurrent_work_part];

                    //if the part is empty, skip the part.
                    if((num_parts != 1  ) && this->global_min_max_coord_total_weight[kk] >
                             this->global_min_max_coord_total_weight[kk + current_concurrent_num_parts]) {

                        for(partId_t jj = 0; jj < num_parts; ++jj){
                            this->new_part_xadj[output_part_index + output_array_shift + jj] = 0;
                        }
                        cut_shift += num_parts - 1;
                        tlr_shift += (4 *(num_parts - 1) + 1);
                        output_array_shift += num_parts;
                        partweight_array_shift += (2 * (num_parts - 1) + 1);
                        continue;
                    }

                    pq_lno_t coordinate_end = this->part_xadj[current_concurrent_work_part];
                    pq_lno_t coordinate_begin = current_concurrent_work_part==0 ? 0: this->part_xadj[current_concurrent_work_part
                                                             -1];
                    pq_scalar_t *current_concurrent_cut_coordinate = current_cut_coordinates + cut_shift;
                    pq_scalar_t *used_local_cut_line_weight_to_left = this->process_cut_line_weight_to_put_left +
                                                         cut_shift;

                    for(int ii = 0; ii < this->num_threads; ++ii){
                        this->thread_part_weight_work[ii] = this->thread_part_weights[ii] +  partweight_array_shift;
                    }

                    if(num_parts > 1){
                        // Rewrite the indices based on the computed cuts.
                    	this->create_consistent_chunks(
                            num_parts,
                            mj_current_dim_coords,
                            current_concurrent_cut_coordinate,
                            coordinate_begin,
                            coordinate_end,
                            used_local_cut_line_weight_to_left,
                            this->new_part_xadj + output_part_index + output_array_shift,
                            coordInd );
                    }
                    else {
                        //if this part is partitioned into 1 then just copy
                        //the old values.
                        pq_lno_t part_size = coordinate_end - coordinate_begin;
                        *(this->new_part_xadj + output_part_index + output_array_shift) = part_size;
                        memcpy(this->new_coordinate_permutations + coordinate_begin,
                        this->coordinate_permutations + coordinate_begin,
                        part_size * sizeof(pq_lno_t));
                    }
                    cut_shift += num_parts - 1;
                    tlr_shift += (4 *(num_parts - 1) + 1);
                    output_array_shift += num_parts;
                    partweight_array_shift += (2 * (num_parts - 1) + 1);
                }

                //shift cut coordinates so that all cut coordinates are stored.
                //current_cut_coordinates += cutShift;

                //getChunks from coordinates partitioned the parts and
                //wrote the indices as if there were a single part.
                //now we need to shift the beginning indices.
                for(partId_t kk = 0; kk < current_concurrent_num_parts; ++kk){
                    partId_t num_parts = num_partitioning_in_current_dim[ current_work_part + kk];
                    for (partId_t ii = 0;ii < num_parts ; ++ii){
                        //shift it by previousCount
                        this->new_part_xadj[output_part_index+ii] += output_coordinate_end_index;
                    }
                    //increase the previous count by current end.
                    output_coordinate_end_index = this->new_part_xadj[output_part_index + num_parts - 1];
                    //increase the current out.
                    output_part_index += num_parts ;
                }
            }
        }
        // end of this partitioning dimension

        //set the current num parts for next dim partitioning
        current_num_parts = output_part_count_in_dimension;

        //swap the coordinate permutations for the next dimension.
        pq_lno_t * tmp = this->coordinate_permutations;
        this->coordinate_permutations = this->new_coordinate_permutations;
        this->new_coordinate_permutations = tmp;

        freeArray<pq_lno_t>(this->part_xadj);
        this->part_xadj = this->new_part_xadj;
    }

    for(pq_lno_t i = 0; i < num_total_coords; ++i){
    	inital_adjList_output_adjlist[i] = this->coordinate_permutations[i];
    }

    for(size_t i = 0; i < this->num_global_parts ; ++i){
        output_xadj[i] = this->part_xadj[i];
    }


    delete future_num_part_in_parts;
    delete next_future_num_parts_in_parts;

    //free the extra memory that we allocated.
    freeArray<partId_t>(this->assigned_part_ids);
    freeArray<pq_gno_t>(this->initial_mj_gnos);
    freeArray<pq_gno_t>(this->current_mj_gnos);
    freeArray<bool>(tmp_mj_uniform_weights);
    freeArray<bool>(tmp_mj_uniform_parts);
    freeArray<pq_scalar_t *>(tmp_mj_weights);
    freeArray<pq_scalar_t *>(tmp_mj_part_sizes);

    this->free_work_memory();

    /*
    freeArray<partId_t>(this->assigned_part_ids);
    freeArray<pq_lno_t>(this->coordinate_permutations);
    freeArray<pq_lno_t>(this->new_coordinate_permutations);
    freeArray<pq_lno_t>(this->part_xadj);

    freeArray<pq_scalar_t>(this->all_cut_coordinates);
    freeArray<pq_scalar_t>(this->cut_coordinates_work_array);
    freeArray<pq_scalar_t>(this->max_min_coords);

    if(this->distribute_points_on_cut_lines){
        freeArray<float>(this->process_cut_line_weight_to_put_left);
        for(int i = 0; i < this->num_threads; ++i){
            freeArray<float>(this->thread_cut_line_weight_to_put_left[i]);
        }
        freeArray<float *>(this->thread_cut_line_weight_to_put_left);
    }
    freeArray<pq_scalar_t>(this->target_part_weights);

    freeArray<pq_scalar_t>(this->cut_upper_bound_coordinates);
    freeArray<pq_scalar_t>(this->cut_lower_bound_coordinates);
    freeArray<pq_scalar_t>(this->cut_lower_bound_weights);
    freeArray<pq_scalar_t>(this->cut_upper_bound_weights);
    freeArray<pq_scalar_t> (this->process_local_min_max_coord_total_weight);
    freeArray<pq_scalar_t> (this->global_min_max_coord_total_weight);
    freeArray<bool>(this->is_cut_line_determined);
    freeArray<partId_t>(this->my_incomplete_cut_count);

    for(int i = 0; i < this->num_threads; ++i){
        freeArray<double>(this->thread_part_weights[i]);
        freeArray<pq_scalar_t>(this->thread_cut_right_closest_point[i]);
        freeArray<pq_scalar_t>(this->thread_cut_left_closest_point[i]);
    }
    freeArray<double *>(this->thread_part_weights);
    freeArray<pq_scalar_t *>(this->thread_cut_left_closest_point);
    freeArray<pq_scalar_t *>(this->thread_cut_right_closest_point);

    freeArray<double *> (this->thread_part_weight_work);

    for(int i = 0; i < this->num_threads; ++i){
        freeArray<pq_lno_t>(this->thread_point_counts[i]);
    }
    freeArray<pq_lno_t *>(this->thread_point_counts);
    freeArray<pq_scalar_t>(this->process_rectelinear_cut_weight);

    freeArray<pq_scalar_t>(this->global_rectelinear_cut_weight);
    freeArray<pq_scalar_t>(this->total_part_weight_left_right_closests);
    freeArray<pq_scalar_t>(this->global_total_part_weight_left_right_closests);
    //env->timerStop(MACRO_TIMERS, "PQJagged - " +partitioningName+"-Problem_Partitioning");
    */

#ifdef HAVE_ZOLTAN2_OMP
    omp_set_num_threads(actual_num_threads);
#endif
}

/*! \brief Multi Jagged  coordinate partitioning algorithm default constructor.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::AlgMJ(){
	this->fEpsilon = numeric_limits<float>::epsilon();
	this->sEpsilon = numeric_limits<pq_scalar_t>::epsilon() * 100;

    this->maxScalar_t = numeric_limits<pq_scalar_t>::max();
    this->minScalar_t = -numeric_limits<pq_scalar_t>::max();

	this->distribute_points_on_cut_lines = true;
	this->max_concurrent_part_calculation = 1;
	this->check_migrate_avoid_migration_option = 0;
	this->minimum_migration_imbalance = 0.35;
	this->mj_keep_part_boxes = 0;
	this->num_threads = 1;

}

/*! \brief Function returns the part boxes stored
 * returns null if boxes are not stored, and prints warning mesage.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::get_part_boxes(){

	if (this->mj_keep_part_boxes){
		return this->kept_boxes;
	} else {
		cerr << "Warning: part boxes are not stored - "
				<< "too store part boxes call set_to_keep_part_boxes() before partitioning" << endl;
		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > a;
		return a;
	}
}
/*! \brief Function call, if the part boxes are intended to be kept.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::set_to_keep_part_boxes(){
	this->mj_keep_part_boxes = 1;
}





/* \brief Either the mj array (part_no_array) or num_global_parts should be provided in
 * the input. part_no_array takes
 * precedence if both are provided.
 * Depending on these parameters, total cut/part number,
 * maximum part/cut number along a dimension, estimated number of reduceAlls,
 * and the number of parts before the last dimension is calculated.
 * */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::set_part_specifications(){

	this->total_num_cut = 0; //how many cuts will be totally
	this->total_num_part = 1;    //how many parts will be totally
	this->max_num_part_along_dim = 0;         //maximum part count along a dimension.
	this->total_dim_num_reduce_all = 0;    //estimate on #reduceAlls can be done.
	this->last_dim_num_part = 1; //max no of parts that might occur
	//during the partition before the
	//last partitioning dimension.
	this->max_num_cut_along_dim = 0;
	this->max_num_total_part_along_dim = 0;


	if (this->part_no_array){
		//if user provided part array, traverse the array and set variables.
		for (int i = 0; i < this->recursion_depth; ++i){
			this->total_dim_num_reduce_all += this->total_num_part;
			this->total_num_part *= this->part_no_array[i];
			if(this->part_no_array[i] > this->max_num_part_along_dim) {
				this->max_num_part_along_dim = this->part_no_array[i];
			}
		}
		this->last_dim_num_part = this->total_num_part / this->part_no_array[recursion_depth-1];
		this->num_global_parts = this->total_num_part;
	} else {
		partId_t future_num_parts = this->num_global_parts;

		//we need to calculate the part numbers now, to determine the maximum along the dimensions.
		for (int i = 0; i < this->recursion_depth; ++i){

			partId_t maxNoPartAlongI = this->get_part_count(
					future_num_parts, 1.0f / (this->recursion_depth - i));

			if (maxNoPartAlongI > this->max_num_part_along_dim){
				this->max_num_part_along_dim = maxNoPartAlongI;
			}

			partId_t nfutureNumParts = future_num_parts / maxNoPartAlongI;
			if (future_num_parts % maxNoPartAlongI){
				++nfutureNumParts;
			}
			future_num_parts = nfutureNumParts;
		}
		this->total_num_part = this->num_global_parts;
		//estimate reduceAll Count here.
		//we find the upperbound instead.
		partId_t p = 1;
		for (int i = 0; i < this->recursion_depth; ++i){
			this->total_dim_num_reduce_all += p;
			p *= this->max_num_part_along_dim;
		}

		this->last_dim_num_part  = p / this->max_num_part_along_dim;
	}

	this->total_num_cut = this->total_num_part - 1;
	this->max_num_cut_along_dim = this->max_num_part_along_dim - 1;
	this->max_num_total_part_along_dim = this->max_num_part_along_dim + size_t(this->max_num_cut_along_dim);
	//maxPartNo is P, maxCutNo = P-1, matTotalPartcount = 2P-1



	//refine the concurrent part count, if it is given bigger than the maximum possible part count.
    if(this->max_concurrent_part_calculation > this->last_dim_num_part){
        if(this->mj_problemComm->getRank() == 0){
            cerr << "Warning: Concurrent part count ("<< this->max_concurrent_part_calculation <<
            ") has been set bigger than maximum amount that can be used." <<
            " Setting to:" << this->last_dim_num_part << "." << endl;
        }
        this->max_concurrent_part_calculation = this->last_dim_num_part;
    }

}
/* \brief Tries to determine the part number for current dimension,
 * by trying to make the partitioning as square as possible.
 * \param num_total_future how many more partitionings are required.
 * \param root how many more recursion depth is left.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
inline partId_t AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::get_part_count(
		partId_t num_total_future,
		double root){
	double fp = pow(num_total_future, root);
	partId_t ip = partId_t (fp);
	if (fp - ip < this->fEpsilon * 100){
		return ip;
	}
	else {
		return ip  + 1;
	}
}

/* \brief Function returns how many parts that will be obtained after this dimension partitioning.
 * It sets how many parts each current part will be partitioned into in this dimension to num_partitioning_in_current_dim vector,
 * sets how many total future parts each obtained part will be partitioned into in next_future_num_parts_in_parts vector,
 * If part boxes are kept, then sets initializes the output_part_boxes as its ancestor.
 *
 *  \param num_partitioning_in_current_dim: output. How many parts each current part will be partitioned into.
 *  \param future_num_part_in_parts: input, how many future parts each current part will be partitioned into.
 *  \param next_future_num_parts_in_parts: output, how many future parts each obtained part will be partitioned into.
 *  \param future_num_parts: output, max number of future parts that will be obtained from a single
 *  \param current_num_parts: input, how many parts are there currently.
 *  \param current_iteration: input, current dimension iteration number.
 *  \param input_part_boxes: input, if boxes are kept, current boxes.
 *  \param output_part_boxes: output, if boxes are kept, the initial box boundaries for obtained parts.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
partId_t AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::update_part_num_arrays(

	vector <partId_t> &num_partitioning_in_current_dim, //assumes this vector is empty.
    vector<partId_t> *future_num_part_in_parts,
    vector<partId_t> *next_future_num_parts_in_parts, //assumes this vector is empty.
    partId_t &future_num_parts,
    partId_t current_num_parts,

    int current_iteration,

    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > input_part_boxes,
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > output_part_boxes
){
	//how many parts that will be obtained after this dimension.
    partId_t output_num_parts = 0;
    if(this->part_no_array){
        //when the partNo array is provided as input,
        //each current partition will be partition to the same number of parts.
        //we dont need to use the future_num_part_in_parts vector in this case.

        partId_t p = this->part_no_array[current_iteration];
        if (p < 1){
            cout << "i:" << current_iteration << " p is given as:" << p << endl;
            exit(1);
        }
        if (p == 1){
            return current_num_parts;
        }

        for (partId_t ii = 0; ii < current_num_parts; ++ii){
            num_partitioning_in_current_dim.push_back(p);

        }
        //cout << "me:" << this->myRank << " current_iteration" << current_iteration << " current_num_parts:" << current_num_parts << endl;
        //cout << "num_partitioning_in_current_dim[0]:" << num_partitioning_in_current_dim[0] << endl;
        //set the new value of future_num_parts.
        future_num_parts /= num_partitioning_in_current_dim[0];
        output_num_parts = current_num_parts * num_partitioning_in_current_dim[0];

        if (this->mj_keep_part_boxes){
            for (partId_t k = 0; k < current_num_parts; ++k){
            	//initialized the output boxes as its ancestor.
            	for (partId_t j = 0; j < num_partitioning_in_current_dim[0]; ++j){
                    output_part_boxes->push_back((*input_part_boxes)[k]);
                }
            }
        }

        //set the how many more parts each part will be divided.
        //this is obvious when partNo array is provided as input.
        //however, fill this so that weights will be calculated according to this array.
        for (partId_t ii = 0; ii < output_num_parts; ++ii){
            next_future_num_parts_in_parts->push_back(future_num_parts);
        }
    }
    else {
        //if partNo array is not provided as input,
        //future_num_part_in_parts  holds how many parts each part should be divided.
        //initially it holds a single number equal to the total number of global parts.

        //calculate the future_num_parts from beginning,
        //since each part might be divided into different number of parts.
        future_num_parts = 1;

        //cout << "i:" << i << endl;

        for (partId_t ii = 0; ii < current_num_parts; ++ii){
            //get how many parts a part should be divided.
            partId_t future_num_parts_of_part_ii = (*future_num_part_in_parts)[ii];

            //get the ideal number of parts that is close to the
            //(recursion_depth - i) root of the future_num_parts_of_part_ii.
            partId_t num_partitions_in_current_dim =
            					this->get_part_count(
            							future_num_parts_of_part_ii,
            							1.0 / (this->recursion_depth - current_iteration)
                                	);

            if (num_partitions_in_current_dim > this->max_num_part_along_dim){
                cerr << "ERROR: maxPartNo calculation is wrong." << endl;
                exit(1);
            }
            //add this number to num_partitioning_in_current_dim vector.
            num_partitioning_in_current_dim.push_back(num_partitions_in_current_dim);


            //increase the output number of parts.
            output_num_parts += num_partitions_in_current_dim;

            //ideal number of future partitions for each part.
            partId_t ideal_num_future_parts_in_part = future_num_parts_of_part_ii / num_partitions_in_current_dim;
            for (partId_t iii = 0; iii < num_partitions_in_current_dim; ++iii){
                partId_t num_future_parts_for_part_iii = ideal_num_future_parts_in_part;

                //if there is a remainder in the part increase the part weight.
                if (iii < future_num_parts_of_part_ii % num_partitions_in_current_dim){
                    //if not uniform, add 1 for the extra parts.
                    ++num_future_parts_for_part_iii;
                }
                next_future_num_parts_in_parts->push_back(num_future_parts_for_part_iii);

                //if part boxes are stored, initialize the box of the parts as the ancestor.
                if (this->mj_keep_part_boxes){
                    output_part_boxes->push_back((*input_part_boxes)[ii]);
                }

                //set num future_num_parts to maximum in this part.
                if (num_future_parts_for_part_iii > future_num_parts) future_num_parts = num_future_parts_for_part_iii;
            }
        }
    }
    return output_num_parts;
}


/* \brief Allocates and initializes the work memory that will be used by MJ.
 *
 * */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::allocate_set_work_memory(){

	//points to process that initially owns the coordinate.
	this->owner_of_coordinate  = NULL;

	//Throughout the partitioning execution,
	//instead of the moving the coordinates, hold a permutation array for parts.
	//coordinate_permutations holds the current permutation.
	this->coordinate_permutations =  allocMemory< pq_lno_t>(this->num_local_coords);
	//initial configuration, set each pointer-i to i.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
	for(pq_lno_t i = 0; i < this->num_local_coords; ++i){
		this->coordinate_permutations[i] = i;
	}

	//new_coordinate_permutations holds the current permutation.
	this->new_coordinate_permutations = allocMemory< pq_lno_t>(this->num_local_coords);



	this->assigned_part_ids = NULL;
	if(this->num_local_coords > 0){
		this->assigned_part_ids = allocMemory<partId_t>(this->num_local_coords);
	}



	//single partition starts at index-0, and ends at numLocalCoords
	//inTotalCounts array holds the end points in coordinate_permutations array
	//for each partition. Initially sized 1, and single element is set to numLocalCoords.
	this->part_xadj = allocMemory<pq_lno_t>(1);
	this->part_xadj[0] = static_cast<pq_lno_t>(this->num_local_coords);//the end of the initial partition is the end of coordinates.
	//the ends points of the output, this is allocated later.
	this->new_part_xadj = NULL;

	// only store this much if cuts are needed to be stored.
	//this->all_cut_coordinates = allocMemory< pq_scalar_t>(this->total_num_cut);


	this->all_cut_coordinates  = allocMemory< pq_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

	this->max_min_coords =  allocMemory< pq_scalar_t>(this->num_threads * 2);

	this->process_cut_line_weight_to_put_left = NULL; //how much weight percentage should a MPI put left side of the each cutline
	this->thread_cut_line_weight_to_put_left = NULL; //how much weight percentage should each thread in MPI put left side of the each outline
	//distribute_points_on_cut_lines = false;
	if(this->distribute_points_on_cut_lines){
		this->process_cut_line_weight_to_put_left = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
		this->thread_cut_line_weight_to_put_left = allocMemory<pq_scalar_t *>(this->num_threads);
		for(int i = 0; i < this->num_threads; ++i){
			this->thread_cut_line_weight_to_put_left[i] = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim);
		}
	    this->process_rectelinear_cut_weight = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim);
	    this->global_rectelinear_cut_weight = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim);
	}


	// work array to manipulate coordinate of cutlines in different iterations.
	//necessary because previous cut line information is used for determining
	//the next cutline information. therefore, cannot update the cut work array
	//until all cutlines are determined.
	this->cut_coordinates_work_array = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim *
			this->max_concurrent_part_calculation);


	//cumulative part weight array.
	this->target_part_weights = allocMemory<pq_scalar_t>(
					this->max_num_part_along_dim * this->max_concurrent_part_calculation);
	// the weight from left to write.

    this->cut_upper_bound_coordinates = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);  //upper bound coordinate of a cut line
    this->cut_lower_bound_coordinates = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim* this->max_concurrent_part_calculation);  //lower bound coordinate of a cut line
    this->cut_lower_bound_weights = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim* this->max_concurrent_part_calculation);  //lower bound weight of a cut line
    this->cut_upper_bound_weights = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim* this->max_concurrent_part_calculation);  //upper bound weight of a cut line

    this->process_local_min_max_coord_total_weight = allocMemory<pq_scalar_t>(3 * this->max_concurrent_part_calculation); //combined array to exchange the min and max coordinate, and total weight of part.
    this->global_min_max_coord_total_weight = allocMemory<pq_scalar_t>(3 * this->max_concurrent_part_calculation);//global combined array with the results for min, max and total weight.

    //is_cut_line_determined is used to determine if a cutline is determined already.
    //If a cut line is already determined, the next iterations will skip this cut line.
    this->is_cut_line_determined = allocMemory<bool>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
    //my_incomplete_cut_count count holds the number of cutlines that have not been finalized for each part
    //when concurrentPartCount>1, using this information, if my_incomplete_cut_count[x]==0, then no work is done for this part.
    this->my_incomplete_cut_count =  allocMemory<partId_t>(this->max_concurrent_part_calculation);
    //local part weights of each thread.
    this->thread_part_weights = allocMemory<double *>(this->num_threads);
    //the work manupulation array for partweights.
    this->thread_part_weight_work = allocMemory<double *>(this->num_threads);

    //thread_cut_left_closest_point to hold the closest coordinate to a cutline from left (for each thread).
    this->thread_cut_left_closest_point = allocMemory<pq_scalar_t *>(this->num_threads);
    //thread_cut_right_closest_point to hold the closest coordinate to a cutline from right (for each thread)
    this->thread_cut_right_closest_point = allocMemory<pq_scalar_t *>(this->num_threads);

    //to store how many points in each part a thread has.
    this->thread_point_counts = allocMemory<pq_lno_t *>(this->num_threads);

    for(int i = 0; i < this->num_threads; ++i){
        //partWeights[i] = allocMemory<pq_scalar_t>(maxTotalPartCount);
        this->thread_part_weights[i] = allocMemory < double >(this->max_num_total_part_along_dim * this->max_concurrent_part_calculation);
        this->thread_cut_right_closest_point[i] = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
        this->thread_cut_left_closest_point[i] = allocMemory<pq_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
        this->thread_point_counts[i] =  allocMemory<pq_lno_t>(this->max_num_part_along_dim);
    }
    //for faster communication, concatanation of
    //totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
    //leftClosest distances sized P-1, since P-1 cut lines
    //rightClosest distances size P-1, since P-1 cut lines.
    this->total_part_weight_left_right_closests = allocMemory<pq_scalar_t>((this->max_num_total_part_along_dim + this->max_num_cut_along_dim * 2) * this->max_concurrent_part_calculation);
    this->global_total_part_weight_left_right_closests = allocMemory<pq_scalar_t>((this->max_num_total_part_along_dim + this->max_num_cut_along_dim * 2) * this->max_concurrent_part_calculation);


    pq_scalar_t **coord = allocMemory<pq_scalar_t *>(this->coord_dim);
    for (int i=0; i < this->coord_dim; i++){
    	coord[i] = allocMemory<pq_scalar_t>(this->num_local_coords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    	for (pq_lno_t j=0; j < this->num_local_coords; j++)
    		coord[i][j] = this->mj_coordinates[i][j];
    }
    this->mj_coordinates = coord;


    int criteria_dim = (this->weight_dim ? this->weight_dim : 1);
    pq_scalar_t **weights = allocMemory<pq_scalar_t *>(criteria_dim);

    for (int i=0; i < criteria_dim; i++){
    	weights[i] = NULL;
    }
    for (int i=0; i < this->weight_dim; i++){
    	weights[i] = allocMemory<pq_scalar_t>(this->num_local_coords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    	for (pq_lno_t j=0; j < this->num_local_coords; j++)
    		weights[i][j] = this->mj_weights[i][j];

    }
	this->mj_weights = weights;
    this->current_mj_gnos = allocMemory<pq_gno_t>(this->num_local_coords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for (pq_lno_t j=0; j < this->num_local_coords; j++)
    	this->current_mj_gnos[j] = this->initial_mj_gnos[j];

    this->owner_of_coordinate = allocMemory<int>(this->num_local_coords);

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for (pq_lno_t j=0; j < this->num_local_coords; j++)
    	this->owner_of_coordinate[j] = this->myActualRank;
}

/* \brief for part communication we keep track of the box boundaries.
 * This is performed when either asked specifically, or when geometric mapping is performed afterwards.
 * This function initializes a single box with all global min and max coordinates.
 * \param initial_partitioning_boxes the input and output vector for boxes.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::init_part_boxes(
		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > & initial_partitioning_boxes
		){

	//local min coords
    pq_scalar_t *mins = allocMemory<pq_scalar_t>(this->coord_dim);
    //global min coords
    pq_scalar_t *gmins = allocMemory<pq_scalar_t>(this->coord_dim);
    //local max coords
    pq_scalar_t *maxs = allocMemory<pq_scalar_t>(this->coord_dim);
    //global max coords
    pq_scalar_t *gmaxs = allocMemory<pq_scalar_t>(this->coord_dim);

    for (int i = 0; i < this->coord_dim; ++i){
        pq_scalar_t localMin = this->mj_coordinates[i][0];
        pq_scalar_t localMax = this->mj_coordinates[i][0];
        for (pq_lno_t j = 1; j < this->num_local_coords; ++j){
            //cout << " pqJagged_coordinates[i][i]:" <<
            //pqJagged_coordinates[i][j] << endl;
            if (this->mj_coordinates[i][j] < localMin){
                localMin = this->mj_coordinates[i][j];
            }
            if (this->mj_coordinates[i][j] > localMax){
                localMax = this->mj_coordinates[i][j];
            }
        }
        //cout << " localMin:" << localMin << endl;
        //cout << " localMax:" << localMax << endl;
        mins[i] = localMin;
        maxs[i] = localMax;
    }
    reduceAll<int, pq_scalar_t>(*this->comm, Teuchos::REDUCE_MIN,
            this->coord_dim, mins, gmins
    );

    reduceAll<int, pq_scalar_t>(*this->comm, Teuchos::REDUCE_MAX,
            this->coord_dim, maxs, gmaxs
    );

    //create single box with all areas.
    coordinateModelPartBox <pq_scalar_t, partId_t> tmpBox (0, this->coord_dim,
                                                        gmins, gmaxs);
    //coordinateModelPartBox <pq_scalar_t, partId_t> tmpBox (0, coordDim);
    freeArray<pq_scalar_t>(mins);
    freeArray<pq_scalar_t>(gmins);
    freeArray<pq_scalar_t>(maxs);
    freeArray<pq_scalar_t>(gmaxs);
    initial_partitioning_boxes->push_back(tmpBox);
}




/*! \brief Function to determine the local minimum and maximum coordinate, and local total weight
 *  in the given set of local points.
 * \param coordinate_begin_index is the start index of the given partition on partitionedPointPermutations.
 * \param coordinate_end_index is the end index of the given partition on partitionedPointPermutations.
 * \param mj_current_coordinate_permutations is the permutation array that point to the actual coordinate index. Sized as numLocalCoords.
 * \param mj_current_dim_coords float-like array representing the coordinates in a single dimension. Sized as numLocalCoords.
 * \param min_coordinate is the output to represent the local minimumCoordinate in  given range of coordinates.
 * \param max_coordinate is the output to represent the local maximum coordinate in the given range of coordinates.
 * \param total_weight is the output to represent the local total weight in the coordinate in the given range of coordinates.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_get_local_min_max_coord_totW(
		pq_lno_t coordinate_begin_index,
		pq_lno_t coordinate_end_index,
		pq_lno_t *mj_current_coordinate_permutations,
		pq_scalar_t *mj_current_dim_coords,
		pq_scalar_t &min_coordinate,
		pq_scalar_t &max_coordinate,
		pq_scalar_t &total_weight){

    //if the part is empty.
    //set the min and max coordinates as reverse.
    if(coordinate_begin_index >= coordinate_end_index)
    {
        min_coordinate = this->maxScalar_t;
        max_coordinate = this->minScalar_t;
        total_weight = 0;
    }
    else {
        pq_scalar_t my_total_weight = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
        {
            //if uniform weights are used, then weight is equal to count.
            if (this->mj_uniform_weights[0]) {
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
                {
                    my_total_weight = coordinate_end_index - coordinate_begin_index;
                }

            }
            else {
                //if not uniform, then weights are reducted from threads.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for reduction(+:my_total_weight)
#endif
                for (pq_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){
                    int i = mj_current_coordinate_permutations[ii];
                    my_total_weight += this->mj_weights[0][i];
                }
            }

            int my_thread_id = 0;
#ifdef HAVE_ZOLTAN2_OMP
            my_thread_id = omp_get_thread_num();
#endif
            pq_scalar_t my_thread_min_coord, my_thread_max_coord;
            my_thread_min_coord=my_thread_max_coord
                =mj_current_dim_coords[mj_current_coordinate_permutations[coordinate_begin_index]];


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
            for(pq_lno_t j = coordinate_begin_index + 1; j < coordinate_end_index; ++j){
                int i = mj_current_coordinate_permutations[j];
                if(mj_current_dim_coords[i] > my_thread_max_coord)
                    my_thread_max_coord = mj_current_dim_coords[i];
                if(mj_current_dim_coords[i] < my_thread_min_coord)
                    my_thread_min_coord = mj_current_dim_coords[i];
            }
            this->max_min_coords[my_thread_id] = my_thread_min_coord;
            this->max_min_coords[my_thread_id + this->num_threads] = my_thread_max_coord;

#ifdef HAVE_ZOLTAN2_OMP
//we need a barrier here, because max_min_array might not be filled by some of the threads.
#pragma omp barrier
#pragma omp single nowait
#endif
            {
                min_coordinate = this->max_min_coords[0];
                for(int i = 1; i < this->num_threads; ++i){
                    if(this->max_min_coords[i] < min_coordinate)
                        min_coordinate = this->max_min_coords[i];
                }
            }

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single nowait
#endif
            {
                max_coordinate = this->max_min_coords[this->num_threads];
                for(int i = this->num_threads + 1; i < this->num_threads * 2; ++i){
                    if(this->max_min_coords[i] > max_coordinate)
                        max_coordinate = this->max_min_coords[i];
                }
            }
        }
        total_weight = my_total_weight;
    }
}


/*! \brief Function that reduces global minimum and maximum coordinates with global total weight from given local arrays.
 * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
 * \param local_min_max_total is the array holding local min and max coordinate values with local total weight.
 * First concurrentPartCount entries are minimums of the parts, next concurrentPartCount entries are max, and then the total weights.
 * \param global_min_max_total is the output array holding global min and global coordinate values with global total weight.
 * The structure is same as localMinMaxTotal.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_get_global_min_max_coord_totW(
    partId_t current_concurrent_num_parts,
    pq_scalar_t *local_min_max_total,
    pq_scalar_t *global_min_max_total){

	//reduce min for first current_concurrent_num_parts elements, reduce max for next
	//concurrentPartCount elements,
	//reduce sum for the last concurrentPartCount elements.
	if(this->comm->getSize()  > 1){
		Teuchos::PQJaggedCombinedMinMaxTotalReductionOp<int, pq_scalar_t>
			reductionOp(
					current_concurrent_num_parts,
					current_concurrent_num_parts,
					current_concurrent_num_parts);
		try{
			reduceAll<int, pq_scalar_t>(
					*(this->comm),
					reductionOp,
					3 * current_concurrent_num_parts,
					local_min_max_total,
					global_min_max_total);
		}
		Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))
	}
	else {
		partId_t s = 3 * current_concurrent_num_parts;
		for (partId_t i = 0; i < s; ++i){
			global_min_max_total[i] = local_min_max_total[i];
		}
	}
}



/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 * \param min_coord minimum coordinate in the range.
 * \param max_coord maximum coordinate in the range.
 *
 * \param num_cuts holds the number of cuts in the current partitioning dimension.
 * \param global_weight holds the global total weight in the current part.
 *
 * \param initial_cut_coords is the output array for the initial cut lines.
 * \param target_part_weights is the output array holding the cumulative ratios of parts in current partitioning.
 * For partitioning to 4 uniformly, target_part_weights will be (0.25 * globalTotalWeight, 0.5 *globalTotalWeight , 0.75 * globalTotalWeight, globalTotalWeight).
 *
 * \param future_num_part_in_parts is the vector that holds how many more parts each part will be divided into more
 * for the parts at the beginning of this coordinate partitioning
 * \param next_future_num_parts_in_parts is the vector that holds how many more parts each part will be divided into more
 * for the parts that will be obtained at the end of this coordinate partitioning.
 * \param concurrent_current_part is the index of the part in the future_num_part_in_parts vector.
 * \param obtained_part_index holds the amount of shift in the next_future_num_parts_in_parts for the output parts.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_get_initial_cut_coords_target_weights(
    pq_scalar_t min_coord,
    pq_scalar_t max_coord,
    partId_t num_cuts/*p-1*/ ,
    pq_scalar_t global_weight,
    pq_scalar_t *initial_cut_coords /*p - 1 sized, coordinate of each cut line*/,
    pq_scalar_t *current_target_part_weights /*cumulative weights, at left side of each cut line. p-1 sized*/,

    vector <partId_t> *future_num_part_in_parts, //the vecto
    vector <partId_t> *next_future_num_parts_in_parts,
    partId_t concurrent_current_part,
    partId_t obtained_part_index
){

    pq_scalar_t coord_range = max_coord - min_coord;
    if(this->mj_uniform_parts[0]){
        {
            partId_t cumulative = 0;
            //how many total future parts the part will be partitioned into.
            pq_scalar_t total_future_part_count_in_part = pq_scalar_t((*future_num_part_in_parts)[concurrent_current_part]);
            //how much each part should weigh in ideal case.
            pq_scalar_t unit_part_weight = global_weight / total_future_part_count_in_part;
            for(partId_t i = 0; i < num_cuts; ++i){
                cumulative += (*next_future_num_parts_in_parts)[i + obtained_part_index];
                //set target part weight.
                current_target_part_weights[i] = cumulative * unit_part_weight;
                //set initial cut coordinate.
                initial_cut_coords[i] = min_coord + (coord_range *
                                         cumulative) / total_future_part_count_in_part;
            }
            current_target_part_weights[num_cuts] = 1;
        }

        //round the target part weights.
        if (this->mj_uniform_weights[0]){
            for(partId_t i = 0; i < num_cuts + 1; ++i){
                current_target_part_weights[i] = long(current_target_part_weights[i] + 0.5);
            }
        }
    }
    else {
    	cerr << "MJ does not support non uniform part weights" << endl;
    	exit(1);
    }
}


/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 * \param max_coordinate maximum coordinate in the range.
 * \param min_coordinate minimum coordinate in the range.
 *
 * \param concurrent_current_part_index is the index of the part in the inTotalCounts vector.
 * \param coordinate_begin_index holds the beginning of the coordinates in current part.
 * \param coordinate_end_index holds end of the coordinates in current part.
 * \param mj_current_coordinate_permutations is the permutation array, holds the real indices of coordinates on pqCoord array.
 * \param mj_current_dim_coords is the 1D array holding the coordinates.
 * \param mj_part_ids is the array holding the partIds of each coordinate.
 * \param partition_count is the number of parts that the current part will be partitioned into.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::set_initial_coordinate_parts(
    pq_scalar_t &max_coordinate,
    pq_scalar_t &min_coordinate,
    partId_t &concurrent_current_part_index,
    pq_lno_t coordinate_begin_index,
    pq_lno_t coordinate_end_index,
    pq_lno_t *mj_current_coordinate_permutations,
    pq_scalar_t *mj_current_dim_coords,
    partId_t *mj_part_ids,
    partId_t &partition_count
){
    pq_scalar_t coordinate_range = max_coordinate - min_coordinate;

    //if there is single point, or if all points are along a line.
    //set initial part to 0 for all.
    if(ABS(coordinate_range) < this->sEpsilon ){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(pq_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){
        	mj_part_ids[mj_current_coordinate_permutations[ii]] = 0;
        }
    }
    else{

        //otherwise estimate an initial part for each coordinate.
        //assuming uniform distribution of points.
        pq_scalar_t slice = coordinate_range / partition_count;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(pq_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){

            pq_lno_t iii = mj_current_coordinate_permutations[ii];
            partId_t pp = partId_t((mj_current_dim_coords[iii] - min_coordinate) / slice);
            mj_part_ids[iii] = 2 * pp;
        }
    }
}


/*! \brief Function that is responsible from 1D partitioning of the given range of coordinates.
 * \param pqJagged_coordinates is 1 dimensional array holding coordinate values.
 * \param imbalanceTolerance is the maximum allowed imbalance ratio.
 * \param current_work_part is the beginning index of concurrentPartCount parts.
 * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
 * \param current_cut_coordinates is the array holding the coordinates of the cut.
 * \param total_incomplete_cut_count is the number of cut lines whose positions should be calculated.
 * \param num_partitioning_in_current_dim is the vector that holds how many parts each part will be divided into.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_1D_part(
    pq_scalar_t *mj_current_dim_coords,
    pq_scalar_t used_imbalance_tolerance,
    partId_t current_work_part,
    partId_t current_concurrent_num_parts,
    pq_scalar_t *current_cut_coordinates,
    partId_t total_incomplete_cut_count,
    vector <partId_t> &num_partitioning_in_current_dim
){


    partId_t rectelinear_cut_count = 0;
    pq_scalar_t *temp_cut_coords = current_cut_coordinates;

    Teuchos::PQJaggedCombinedReductionOp<partId_t, pq_scalar_t>
                 *reductionOp = NULL;
    reductionOp = new Teuchos::PQJaggedCombinedReductionOp
                     <partId_t, pq_scalar_t>(
                    		 &num_partitioning_in_current_dim ,
                    		 current_work_part ,
                    		 current_concurrent_num_parts);

    size_t total_reduction_size = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(total_incomplete_cut_count,  rectelinear_cut_count)
#endif
    {
        int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
        me = omp_get_thread_num();
#endif
        double *my_thread_part_weights = this->thread_part_weights[me];
        pq_scalar_t *my_thread_left_closest = this->thread_cut_left_closest_point[me];
        pq_scalar_t *my_thread_right_closest = this->thread_cut_right_closest_point[me];

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            {
                //initialize the lower and upper bounds of the cuts.
                partId_t next = 0;
                for(partId_t i = 0; i < current_concurrent_num_parts; ++i){

                    partId_t num_part_in_dim =  num_partitioning_in_current_dim[current_work_part + i];
                    partId_t num_cut_in_dim = num_part_in_dim - 1;
                    total_reduction_size += (4 * num_cut_in_dim + 1);

                    for(partId_t ii = 0; ii < num_cut_in_dim; ++ii){
                        this->is_cut_line_determined[next] = false;
                        this->cut_lower_bound_coordinates[next] = global_min_max_coord_total_weight[i]; //min coordinate
                        this->cut_upper_bound_coordinates[next] = global_min_max_coord_total_weight[i + current_concurrent_num_parts]; //max coordinate

                        this->cut_upper_bound_weights[next] = global_min_max_coord_total_weight[i + 2 * current_concurrent_num_parts]; //total weight
                        this->cut_lower_bound_weights[next] = 0;

                        if(this->distribute_points_on_cut_lines){
                            this->process_cut_line_weight_to_put_left[next] = 0;
                        }
                        ++next;
                    }
                }
            }

        //no need to have barrier here.
        //pragma omp single have implicit barrier.

        int iteration = 0;
        while (total_incomplete_cut_count != 0){
            iteration += 1;
            //cout << "\niteration:" << iteration  << " ";
            partId_t concurrent_cut_shifts = 0;
            size_t total_part_shift = 0;

            for (partId_t kk = 0; kk < current_concurrent_num_parts; ++kk){
                partId_t num_parts =  -1;
                num_parts =  num_partitioning_in_current_dim[current_work_part + kk];

                partId_t num_cuts = num_parts - 1;
                size_t total_part_count = num_parts + size_t (num_cuts) ;
                if (this->my_incomplete_cut_count[kk] > 0){

                    //although isDone shared, currentDone is private and same for all.
                    bool *current_cut_status = this->is_cut_line_determined + concurrent_cut_shifts;
                    double *my_current_part_weights = my_thread_part_weights + total_part_shift;
                    pq_scalar_t *my_current_left_closest = my_thread_left_closest + concurrent_cut_shifts;
                    pq_scalar_t *my_current_right_closest = my_thread_right_closest + concurrent_cut_shifts;

                    partId_t conccurent_current_part = current_work_part + kk;
                    pq_lno_t coordinate_begin_index = conccurent_current_part ==0 ? 0: this->part_xadj[conccurent_current_part -1];
                    pq_lno_t coordinate_end_index = this->part_xadj[conccurent_current_part];
                    pq_scalar_t *temp_current_cut_coords = temp_cut_coords + concurrent_cut_shifts;

                    pq_scalar_t min_coord = global_min_max_coord_total_weight[kk];
                    pq_scalar_t max_coord = global_min_max_coord_total_weight[kk + current_concurrent_num_parts];

                    // compute part weights using existing cuts
                    this->mj_1D_part_get_thread_part_weights(
                        total_part_count,
                        num_cuts,
                        max_coord,//globalMinMaxTotal[kk + concurrentPartCount],//maxScalar,
                        min_coord,//globalMinMaxTotal[kk]//minScalar,
                        coordinate_begin_index,
                        coordinate_end_index,
                        mj_current_dim_coords,
                        temp_current_cut_coords,
                        current_cut_status,
                        my_current_part_weights,
                        my_current_left_closest,
                        my_current_right_closest);

                }

                concurrent_cut_shifts += num_cuts;
                total_part_shift += total_part_count;
            }

            //sum up the results of threads
            this->mj_accumulate_thread_results(
                num_partitioning_in_current_dim,
                current_work_part,
                current_concurrent_num_parts);

            //now sum up the results of mpi processors.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            {
                if(this->comm->getSize() > 1){
                    try{


                        reduceAll<int, pq_scalar_t>( *(this->comm), *reductionOp,
                        		total_reduction_size,
                        		this->total_part_weight_left_right_closests,
                        		this->global_total_part_weight_left_right_closests);

                    }
                    Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))
                }
                else {
                        memcpy(
                        	this->global_total_part_weight_left_right_closests,
                            this->total_part_weight_left_right_closests,
                            total_reduction_size * sizeof(pq_scalar_t));
                }
            }

            //how much cut will be shifted for the next part in the concurrent part calculation.
            partId_t cut_shift = 0;

            //how much the concantaneted array will be shifted for the next part in concurrent part calculation.
            size_t tlr_shift = 0;
            for (partId_t kk = 0; kk < current_concurrent_num_parts; ++kk){
            	partId_t num_parts =  num_partitioning_in_current_dim[current_work_part + kk];
            	partId_t num_cuts = num_parts - 1;
            	size_t num_total_part = num_parts + size_t (num_cuts) ;

            	//if the cuts of this cut has already been completed.
            	//nothing to do for this part.
            	//just update the shift amount and proceed.
            	if (this->my_incomplete_cut_count[kk] == 0) {
            		cut_shift += num_cuts;
            		tlr_shift += (num_total_part + 2 * num_cuts);
            		continue;
            	}

            	pq_scalar_t *current_local_part_weights = this->total_part_weight_left_right_closests  + tlr_shift ;
            	pq_scalar_t *current_global_tlr = this->global_total_part_weight_left_right_closests + tlr_shift;
            	pq_scalar_t *current_global_left_closest_points = current_global_tlr + num_total_part; //left closest points
            	pq_scalar_t *current_global_right_closest_points = current_global_tlr + num_total_part + num_cuts; //right closest points
            	pq_scalar_t *current_global_part_weights = current_global_tlr;
            	bool *current_cut_line_determined = this->is_cut_line_determined + cut_shift;

            	pq_scalar_t *current_part_target_weights = this->target_part_weights + cut_shift + kk;
            	pq_scalar_t *current_part_cut_line_weight_to_put_left = this->process_cut_line_weight_to_put_left + cut_shift;

            	pq_scalar_t min_coordinate = global_min_max_coord_total_weight[kk];
            	pq_scalar_t max_coordinate = global_min_max_coord_total_weight[kk + current_concurrent_num_parts];
            	pq_scalar_t global_total_weight = global_min_max_coord_total_weight[kk + current_concurrent_num_parts * 2];
            	pq_scalar_t *current_cut_lower_bound_weights = this->cut_lower_bound_weights + cut_shift;
            	pq_scalar_t *current_cut_upper_weights = this->cut_upper_bound_weights + cut_shift;
            	pq_scalar_t *current_cut_upper_bounds = this->cut_upper_bound_coordinates + cut_shift;
            	pq_scalar_t *current_cut_lower_bounds = this->cut_lower_bound_coordinates + cut_shift;

            	partId_t initial_incomplete_cut_count = this->my_incomplete_cut_count[kk];

            	// Now compute the new cut coordinates.
            	this->mj_get_new_cut_coordinates(
            			num_total_part,
            			num_cuts,
            			max_coordinate,
            			min_coordinate,
            			global_total_weight,
            			used_imbalance_tolerance,
            			current_global_part_weights,
            			current_local_part_weights,
            			current_part_target_weights,
            			current_cut_line_determined,
            			temp_cut_coords + cut_shift,
            			current_cut_upper_bounds,
            			current_cut_lower_bounds,
            			current_global_left_closest_points,
            			current_global_right_closest_points,
            			current_cut_lower_bound_weights,
            			current_cut_upper_weights,
            			this->cut_coordinates_work_array +cut_shift, //new cut coordinates
            			current_part_cut_line_weight_to_put_left,
            			&rectelinear_cut_count,
            			this->my_incomplete_cut_count[kk]);

            	cut_shift += num_cuts;
            	tlr_shift += (num_total_part + 2 * num_cuts);
            	partId_t iteration_complete_cut_count = initial_incomplete_cut_count - this->my_incomplete_cut_count[kk];
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            	{
            		total_incomplete_cut_count -= iteration_complete_cut_count;
            	}

            }
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp barrier
#pragma omp single
#endif
            {
            	//swap the cut coordinates for next iteration.
                pq_scalar_t *t = temp_cut_coords;
                temp_cut_coords = this->cut_coordinates_work_array;
                this->cut_coordinates_work_array = t;
            }
        }

        // Needed only if keep_cuts; otherwise can simply swap array pointers
        // cutCoordinates and cutCoordinatesWork.
        // (at first iteration, cutCoordinates == cutCoorindates_tmp).
        // computed cuts must be in cutCoordinates.
        if (current_cut_coordinates != temp_cut_coords){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
        	{
        		partId_t next = 0;
        		for(partId_t i = 0; i < current_concurrent_num_parts; ++i){
        			partId_t num_parts = -1;
        			num_parts = num_partitioning_in_current_dim[current_work_part + i];
        			partId_t num_cuts = num_parts - 1;

        			for(partId_t ii = 0; ii < num_cuts; ++ii){
        				current_cut_coordinates[next + ii] = temp_cut_coords[next + ii];
        			}
        			next += num_cuts;
        		}
        	}

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            {
                this->cut_coordinates_work_array = temp_cut_coords;
            }
        }
    }
    delete reductionOp;
}


/*! \brief Function that calculates the weights of each part according to given part cut coordinates.
 * Function is called inside the parallel region. Thread specific work arrays are provided
 * as function parameter.
 *
 * \param total_part_count is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
 * \param num_cuts is the number of cut lines. P - 1.
 * \param max_coord is the maximum coordinate in the part.
 * \param min_coord is the min coordinate in the part.
 * \param coordinate_begin_index is the index of the first coordinate in current part.
 * \param coordinate_end_index is the index of the last coordinate in current part.
 * \param mj_current_dim_coords is 1 dimensional array holding coordinate values.
 *
 * \param temp_current_cut_coords is the array holding the coordinates of each cut line. Sized P - 1.
 * \param current_cut_status is the boolean array to determine if the correct position for a cut line is found.
 * \param my_current_part_weights is the array holding the part weights for the calling thread.
 * \param my_current_left_closest is the array holding the coordinate of the closest points to the cut lines from left for the calling thread..
 * \param my_current_right_closest is the array holding the coordinate of the closest points to the cut lines from right for the calling thread.
 * \param partIds is the array that holds the part ids of the coordinates
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_1D_part_get_thread_part_weights(
    size_t total_part_count,
    partId_t num_cuts,
    pq_scalar_t max_coord,
    pq_scalar_t min_coord,
    pq_lno_t coordinate_begin_index,
    pq_lno_t coordinate_end_index,
    pq_scalar_t *mj_current_dim_coords,
    pq_scalar_t *temp_current_cut_coords,
    bool *current_cut_status,
    double *my_current_part_weights,
    pq_scalar_t *my_current_left_closest,
    pq_scalar_t *my_current_right_closest){

	// initializations for part weights, left/right closest
	for (size_t i = 0; i < total_part_count; ++i){
		my_current_part_weights[i] = 0;
	}

	//initialize the left and right closest coordinates
	//to their max value.
	for(partId_t i = 0; i < num_cuts; ++i){
		my_current_left_closest[i] = min_coord - 1;
		my_current_right_closest[i] = max_coord + 1;
	}
	//pq_lno_t comparison_count = 0;
	pq_scalar_t minus_EPSILON = -this->sEpsilon;
#ifdef HAVE_ZOLTAN2_OMP
	//no need for the barrier as all threads uses their local memories.
	//dont change the static scheduling here, as it is assumed when the new
	//partitions are created later.
#pragma omp for
#endif
	for (pq_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){
		int i = this->coordinate_permutations[ii];

		//the accesses to assigned_part_ids are thread safe
		//since each coordinate is assigned to only a single thread.
		partId_t j = this->assigned_part_ids[i] / 2;

		if(j >= num_cuts){
			j = num_cuts - 1;
		}

		partId_t lower_cut_index = 0;
		partId_t upper_cut_index = num_cuts - 1;

		pq_scalar_t w = this->mj_uniform_weights[0]? 1:this->mj_weights[0][i];
		bool is_inserted = false;
		bool is_on_left_of_cut = false;
		bool is_on_right_of_cut = false;
		partId_t last_compared_part = -1;

		pq_scalar_t coord = mj_current_dim_coords[i];

		while(upper_cut_index >= lower_cut_index)
		{
			//comparison_count++;
			last_compared_part = -1;
			is_on_left_of_cut = false;
			is_on_right_of_cut = false;
			pq_scalar_t cut = temp_current_cut_coords[j];
			pq_scalar_t distance_to_cut = coord - cut;
			pq_scalar_t abs_distance_to_cut = ABS(distance_to_cut);

			//if it is on the line.
			if(abs_distance_to_cut < this->sEpsilon){

				my_current_part_weights[j * 2 + 1] += w;
				this->assigned_part_ids[i] = j * 2 + 1;

				//assign left and right closest point to cut as the point is on the cut.
				my_current_left_closest[j] = coord;
				my_current_right_closest[j] = coord;
				//now we need to check if there are other cuts on the same cut coordinate.
				//if there are, then we add the weight of the cut to all cuts in the same coordinate.
				partId_t kk = j + 1;
				while(kk < num_cuts){
					// Needed when cuts shared the same position
					distance_to_cut =ABS(temp_current_cut_coords[kk] - cut);
					if(distance_to_cut < this->sEpsilon){
						my_current_part_weights[2 * kk + 1] += w;
						my_current_left_closest[kk] = coord;
						my_current_right_closest[kk] = coord;
						kk++;
					}
					else{
						//cut is far away.
						//just check the left closest point for the next cut.
						if(coord - my_current_left_closest[kk] > this->sEpsilon){
							my_current_left_closest[kk] = coord;
						}
						break;
					}
				}


				kk = j - 1;
				//continue checking for the cuts on the left if they share the same coordinate.
				while(kk >= 0){
					distance_to_cut =ABS(temp_current_cut_coords[kk] - cut);
					if(distance_to_cut < this->sEpsilon){
						my_current_part_weights[2 * kk + 1] += w;
						//try to write the partId as the leftmost cut.
						this->assigned_part_ids[i] = kk * 2 + 1;
						my_current_left_closest[kk] = coord;
						my_current_right_closest[kk] = coord;
						kk--;
					}
					else{
						//if cut is far away on the left of the point.
						//then just compare for right closest point.
						if(my_current_right_closest[kk] - coord > this->sEpsilon){
							my_current_right_closest[kk] = coord;
						}
						break;
					}
				}

				is_inserted = true;
				break;
			}
			else {
				//if point is on the left of the cut.
				if (distance_to_cut < 0) {
					bool _break = false;
					if(j > 0){
						//check distance to the cut on the left the current cut compared.
						//if point is on the right, then we find the part of the point.
						pq_scalar_t distance_to_next_cut = coord - temp_current_cut_coords[j - 1];
						if(distance_to_next_cut > this->sEpsilon){
							_break = true;
						}
					}
					//if point is not on the right of the next cut, then
					//set the upper bound to this cut.
					upper_cut_index = j - 1;
					//set the last part, and mark it as on the left of the last part.
					is_on_left_of_cut = true;
					last_compared_part = j;
					if(_break) break;
				}
				else {
					//if point is on the right of the cut.
					bool _break = false;
					if(j < num_cuts - 1){
						//check distance to the cut on the left the current cut compared.
						//if point is on the right, then we find the part of the point.
						pq_scalar_t distance_to_next_cut = coord - temp_current_cut_coords[j + 1];
						if(distance_to_next_cut < minus_EPSILON){
                    	 _break = true;
                     }
					}

					//if point is not on the left of the next cut, then
					//set the upper bound to this cut.
					lower_cut_index = j + 1;
					//set the last part, and mark it as on the right of the last part.
					is_on_right_of_cut = true;
					last_compared_part = j;
					if(_break) break;
				}
			}

			j = (upper_cut_index + lower_cut_index) / 2;
		}
		if(!is_inserted){
			if(is_on_right_of_cut){

				//add it to the right of the last compared part.
				my_current_part_weights[2 * last_compared_part + 2] += w;
				this->assigned_part_ids[i] = 2 * last_compared_part + 2;

				//update the right closest point of last compared cut.
				if(my_current_right_closest[last_compared_part] - coord > this->sEpsilon){
					my_current_right_closest[last_compared_part] = coord;
				}
				//update the left closest point of the cut on the right of the last compared cut.
				if(last_compared_part+1 < num_cuts){

					if(coord - my_current_left_closest[last_compared_part + 1] > this->sEpsilon){
						my_current_left_closest[last_compared_part + 1] = coord;
					}
				}

			}
			else if(is_on_left_of_cut){

				//add it to the left of the last compared part.
				my_current_part_weights[2 * last_compared_part] += w;
				this->assigned_part_ids[i] = 2 * last_compared_part;


				//update the left closest point of last compared cut.
				if(coord - my_current_left_closest[last_compared_part] > this->sEpsilon){
					my_current_left_closest[last_compared_part] = coord;
				}

				//update the right closest point of the cut on the left of the last compared cut.
				if(last_compared_part-1 >= 0){
					if(my_current_right_closest[last_compared_part -1] - coord > this->sEpsilon){
						my_current_right_closest[last_compared_part -1] = coord;
					}
				}
			}
		}
	}

	// prefix sum computation.
	//we need prefix sum for each part to determine cut positions.
	for (size_t i = 1; i < total_part_count; ++i){
		// check for cuts sharing the same position; all cuts sharing a position
		// have the same weight == total weight for all cuts sharing the position.
		// don't want to accumulate that total weight more than once.
		if(i % 2 == 0 && i > 1 && i < total_part_count - 1 &&
				ABS(temp_current_cut_coords[i / 2] - temp_current_cut_coords[i /2 - 1])
		< this->sEpsilon){
			//i % 2 = 0 when part i represents the cut coordinate.
			//if it is a cut, and if the next cut also have the same coordinate, then
			//dont addup.
			my_current_part_weights[i] = my_current_part_weights[i-2];
			continue;
		}
		//otherwise do the prefix sum.
		my_current_part_weights[i] += my_current_part_weights[i-1];
	}
}


/*! \brief Function that reduces the result of multiple threads
 * for left and right closest points and part weights in a single mpi process.
 *
 * \param num_partitioning_in_current_dim is the vector that holds the number of cut lines in current dimension for each part.
 * \param current_work_part holds the index of the first part (important when concurrent parts are used.)
 * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_accumulate_thread_results(
    const vector <partId_t> &num_partitioning_in_current_dim,
    partId_t current_work_part,
    partId_t current_concurrent_num_parts){

#ifdef HAVE_ZOLTAN2_OMP
	//needs barrier here, as it requires all threads to finish mj_1D_part_get_thread_part_weights
	//using parallel region here reduces the performance because of the cache invalidates.
#pragma omp barrier
#pragma omp single
#endif
	{
		size_t tlr_array_shift = 0;
		partId_t cut_shift = 0;

		//iterate for all concurrent parts to find the left and right closest points in the process.
		for(partId_t i = 0; i < current_concurrent_num_parts; ++i){

			partId_t num_parts_in_part =  num_partitioning_in_current_dim[current_work_part + i];
			partId_t num_cuts_in_part = num_parts_in_part - 1;
			size_t num_total_part_in_part = num_parts_in_part + size_t (num_cuts_in_part) ;

			//iterate for cuts in a single part.
			for(partId_t ii = 0; ii < num_cuts_in_part ; ++ii){
				partId_t next = tlr_array_shift + ii;
				partId_t cut_index = cut_shift + ii;
				if(this->is_cut_line_determined[cut_index]) continue;
				pq_scalar_t left_closest_in_process = this->thread_cut_left_closest_point[0][cut_index],
						right_closest_in_process = this->thread_cut_right_closest_point[0][cut_index];

				//find the closest points from left and right for the cut in the process.
				for (int j = 1; j < this->num_threads; ++j){
					if (this->thread_cut_right_closest_point[j][cut_index] < right_closest_in_process ){
						right_closest_in_process = this->thread_cut_right_closest_point[j][cut_index];
					}
					if (this->thread_cut_left_closest_point[j][cut_index] > left_closest_in_process ){
						left_closest_in_process = this->thread_cut_left_closest_point[j][cut_index];
					}
				}
				//store the left and right closes points.
				this->total_part_weight_left_right_closests[num_total_part_in_part +
				                                            next] = left_closest_in_process;
				this->total_part_weight_left_right_closests[num_total_part_in_part +
				                                            num_cuts_in_part + next] = right_closest_in_process;
			}
			//set the shift position in the arrays
			tlr_array_shift += (num_total_part_in_part + 2 * num_cuts_in_part);
			cut_shift += num_cuts_in_part;
		}

		tlr_array_shift = 0;
		cut_shift = 0;
		size_t total_part_array_shift = 0;

		//iterate for all concurrent parts to find the total weight in the process.
		for(partId_t i = 0; i < current_concurrent_num_parts; ++i){

			partId_t num_parts_in_part =  num_partitioning_in_current_dim[current_work_part + i];
			partId_t num_cuts_in_part = num_parts_in_part - 1;
			size_t num_total_part_in_part = num_parts_in_part + size_t (num_cuts_in_part) ;

			for(size_t j = 0; j < num_total_part_in_part; ++j){

				partId_t cut_ind = j / 2 + cut_shift;

				//need to check j !=  num_total_part_in_part - 1
						// which is same as j/2 != num_cuts_in_part.
				//we cannot check it using cut_ind, because of the concurrent part concantanetion.
				if(j !=  num_total_part_in_part - 1 && this->is_cut_line_determined[cut_ind]) continue;
				double pwj = 0;
				for (int k = 0; k < this->num_threads; ++k){
					pwj += this->thread_part_weights[k][total_part_array_shift + j];
				}
				//size_t jshift = j % total_part_count + i * (total_part_count + 2 * noCuts);
				this->total_part_weight_left_right_closests[tlr_array_shift + j] = pwj;
			}
			cut_shift += num_cuts_in_part;
			tlr_array_shift += num_total_part_in_part + 2 * num_cuts_in_part;
			total_part_array_shift += num_total_part_in_part;
		}
	}
	//the other threads needs to wait here.
	//but we don't need a pragma omp barrier.
	//as omp single has already have implicit barrier.
}


/*! \brief
 * Function that calculates the next pivot position,
 * according to given coordinates of upper bound and lower bound, the weights at upper and lower bounds, and the expected weight.
 * \param cut_upper_bound is the upper bound coordinate of the cut.
 * \param cut_lower_bound is the lower bound coordinate of the cut.
 * \param cut_upper_weight is the weights at the upper bound of the cut.
 * \param cut_lower_weight is the weights at the lower bound of the cut.
 * \param expected_weight is the expected weight that should be placed on the left of the cut line.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_calculate_new_cut_position (
	pq_scalar_t cut_upper_bound,
    pq_scalar_t cut_lower_bound,
    pq_scalar_t cut_upper_weight,
    pq_scalar_t cut_lower_weight,
    pq_scalar_t expected_weight,
    pq_scalar_t &new_cut_position){

    if(ABS(cut_upper_bound - cut_lower_bound) < this->sEpsilon){
    	new_cut_position = cut_upper_bound; //or lower bound does not matter.
    }


    if(ABS(cut_upper_weight - cut_lower_weight) < this->sEpsilon){
    	new_cut_position = cut_lower_bound;
    }

    pq_scalar_t coordinate_range = (cut_upper_bound - cut_lower_bound);
    pq_scalar_t weight_range = (cut_upper_weight - cut_lower_weight);
    pq_scalar_t my_weight_diff = (expected_weight - cut_lower_weight);

    pq_scalar_t required_shift = (my_weight_diff / weight_range);
    int scale_constant = 20;
    int shiftint= int (required_shift * scale_constant);
    if (shiftint == 0) shiftint = 1;
    required_shift = pq_scalar_t (shiftint) / scale_constant;
    new_cut_position = coordinate_range * required_shift + cut_lower_bound;
}


/*! \brief Function that determines the permutation indices of the coordinates.
 * \param num_parts is the number of parts.
 * \param mj_current_dim_coords is 1 dimensional array holding the coordinate values.
 * \param current_concurrent_cut_coordinate is 1 dimensional array holding the cut coordinates.
 * \param coordinate_begin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinate_end is the end index of the given partition on partitionedPointPermutations.
 * \param used_local_cut_line_weight_to_left holds how much weight of the coordinates on the cutline should be put on left side.
 * \param used_thread_part_weight_work is the two dimensional array holding the weight of parts for each thread. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param out_part_xadj is the indices of coordinates calculated for the partition on next dimension. \param pqCoordDim is dimension of the input
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_create_new_partitions(
    partId_t num_parts,
    pq_scalar_t *mj_current_dim_coords,
    pq_scalar_t *current_concurrent_cut_coordinate,
    pq_lno_t coordinate_begin,
    pq_lno_t coordinate_end,
    pq_scalar_t *used_local_cut_line_weight_to_left,
    double **used_thread_part_weight_work,
    pq_lno_t *out_part_xadj){

	partId_t num_cuts = num_parts - 1;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
	{
		int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
		me = omp_get_thread_num();
#endif

		pq_lno_t *thread_num_points_in_parts = this->thread_point_counts[me];
		pq_scalar_t *my_local_thread_cut_weights_to_put_left = NULL;

		//now if the rectelinear partitioning is allowed we decide how
		//much weight each thread should put to left and right.
		if (this->distribute_points_on_cut_lines){
			my_local_thread_cut_weights_to_put_left = this->thread_cut_line_weight_to_put_left[me];
			// this for assumes the static scheduling in mj_1D_part calculation.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
			for (partId_t i = 0; i < num_cuts; ++i){
				//the left to be put on the left of the cut.
				pq_scalar_t left_weight = used_local_cut_line_weight_to_left[i];
				for(int ii = 0; ii < this->num_threads; ++ii){
					if(left_weight > this->sEpsilon){
						//the weight of thread ii on cut.
						pq_scalar_t thread_ii_weight_on_cut = used_thread_part_weight_work[ii][i * 2 + 1] - used_thread_part_weight_work[ii][i * 2 ];
						if(thread_ii_weight_on_cut < left_weight){
							//if left weight is bigger than threads weight on cut.
							this->thread_cut_line_weight_to_put_left[ii][i] = thread_ii_weight_on_cut;
						}
						else {
							//if thread's weight is bigger than space, then put only a portion.
							this->thread_cut_line_weight_to_put_left[ii][i] = left_weight ;
						}
						left_weight -= thread_ii_weight_on_cut;
					}
					else {
						this->thread_cut_line_weight_to_put_left[ii][i] = 0;
					}
				}
			}

			if(num_cuts > 0){
				//this is a special case. If cutlines share the same coordinate, their weights are equal.
				//we need to adjust the ratio for that.
				for (partId_t i = num_cuts - 1; i > 0 ; --i){
					if(ABS(current_concurrent_cut_coordinate[i] - current_concurrent_cut_coordinate[i -1]) < this->sEpsilon){
						my_local_thread_cut_weights_to_put_left[i] -= my_local_thread_cut_weights_to_put_left[i - 1] ;
					}
					my_local_thread_cut_weights_to_put_left[i] = int ((my_local_thread_cut_weights_to_put_left[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL)
                    						/ pq_scalar_t(SIGNIFICANCE_MUL);
				}
			}
		}

		for(partId_t ii = 0; ii < num_parts; ++ii){
			thread_num_points_in_parts[ii] = 0;
		}


#ifdef HAVE_ZOLTAN2_OMP
		//dont change static scheduler. the static partitioner used later as well.
#pragma omp for
#endif
		for (pq_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){

			pq_lno_t coordinate_index = this->coordinate_permutations[ii];
			pq_scalar_t coordinate_weight = this->mj_uniform_weights[0]? 1:this->mj_weights[0][coordinate_index];
			partId_t coordinate_assigned_place = this->assigned_part_ids[coordinate_index];
			partId_t coordinate_assigned_part = coordinate_assigned_place / 2;
			if(coordinate_assigned_place % 2 == 1){
				//if it is on the cut.
				if(this->distribute_points_on_cut_lines
						&& my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] > this->sEpsilon * EPS_SCALE){
					//if the rectelinear partitioning is allowed,
					//and the thread has still space to put on the left of the cut
					//then thread puts the vertex to left.
					my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] -= coordinate_weight;
					//if putting the vertex to left increased the weight more than expected.
					//and if the next cut is on the same coordinate,
					//then we need to adjust how much weight next cut puts to its left as well,
					//in order to take care of the imbalance.
					if(my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] < 0
							&& coordinate_assigned_part < num_cuts - 1
							&& ABS(current_concurrent_cut_coordinate[coordinate_assigned_part+1] -
									current_concurrent_cut_coordinate[coordinate_assigned_part]) < this->sEpsilon){
						my_local_thread_cut_weights_to_put_left[coordinate_assigned_part + 1] += my_local_thread_cut_weights_to_put_left[coordinate_assigned_part];
					}
					++thread_num_points_in_parts[coordinate_assigned_part];
					this->assigned_part_ids[coordinate_index] = coordinate_assigned_part;
				}
				else{
					//if there is no more space on the left, put the coordinate to the right of the cut.
					++coordinate_assigned_part;
					//this while loop is necessary when a line is partitioned into more than 2 parts.
					while(this->distribute_points_on_cut_lines &&
							coordinate_assigned_part < num_cuts){
						//traverse all the cut lines having the same partitiong
						if(ABS(current_concurrent_cut_coordinate[coordinate_assigned_part] -
								current_concurrent_cut_coordinate[coordinate_assigned_part - 1])
								< this->sEpsilon){
							//if line has enough space on left, put it there.
							if(my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] >
							this->sEpsilon &&
							my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] >=
							ABS(my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] - coordinate_weight)){
								my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] -= coordinate_weight;
								//Again if it put too much on left of the cut,
								//update how much the next cut sharing the same coordinate will put to its left.
								if(my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] < 0 &&
										coordinate_assigned_part < num_cuts - 1 &&
										ABS(current_concurrent_cut_coordinate[coordinate_assigned_part+1] -
												current_concurrent_cut_coordinate[coordinate_assigned_part]) < this->sEpsilon){
									my_local_thread_cut_weights_to_put_left[coordinate_assigned_part + 1] += my_local_thread_cut_weights_to_put_left[coordinate_assigned_part];
								}
								break;
							}
						}
						else {
							break;
						}
						++coordinate_assigned_part;
					}
					++thread_num_points_in_parts[coordinate_assigned_part];
					this->assigned_part_ids[coordinate_index] = coordinate_assigned_part;
				}
			}
			else {
				//if it is already assigned to a part, then just put it to the corresponding part.
				++thread_num_points_in_parts[coordinate_assigned_part];
				this->assigned_part_ids[coordinate_index] = coordinate_assigned_part;
			}
		}



		//now we calculate where each thread will write in new_coordinate_permutations array.
		//first we find the out_part_xadj, by marking the begin and end points of each part found.
		//the below loop find the number of points in each part, and writes it to out_part_xadj
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
		for(partId_t j = 0; j < num_parts; ++j){
			pq_lno_t num_points_in_part_j_upto_thread_i = 0;
			for (int i = 0; i < this->num_threads; ++i){
				pq_lno_t thread_num_points_in_part_j = this->thread_point_counts[i][j];
				//prefix sum to thread point counts, so that each will have private space to write.
				this->thread_point_counts[i][j] = num_points_in_part_j_upto_thread_i;
				num_points_in_part_j_upto_thread_i += thread_num_points_in_part_j;

			}
			out_part_xadj[j] = num_points_in_part_j_upto_thread_i;// + prev2; //+ coordinateBegin;
		}

		//now we need to do a prefix sum to out_part_xadj[j], to point begin and end of each part.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
		{
			//perform prefix sum for num_points in parts.
			for(partId_t j = 1; j < num_parts; ++j){
				out_part_xadj[j] += out_part_xadj[j - 1];
			}
		}

		//shift the num points in threads thread to obtain the
		//beginning index of each thread's private space.
		for(partId_t j = 1; j < num_parts; ++j){
			thread_num_points_in_parts[j] += out_part_xadj[j - 1] ;
		}


		//now thread gets the coordinate and writes the index of coordinate to the permutation array
		//using the part index we calculated.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
		for (pq_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){
			pq_lno_t i = this->coordinate_permutations[ii];
			partId_t p =  this->assigned_part_ids[i];
			this->new_coordinate_permutations[coordinate_begin +
			                                  thread_num_points_in_parts[p]++] = i;
		}
	}
}



/*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
 *
 * \param num_total_part is the sum of number of cutlines and number of parts. Simply it is 2*P - 1.
 * \param num_cuts is the number of cut lines. P - 1.
 * \param max_coordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
 * \param min_coordinate is the maximum coordinate in the current range of coordinates and in the current dimension.
 * \param global_total_weight is the global total weight in the current range of coordinates.
 * \param used_imbalance_tolerance is the maximum allowed imbalance ratio.
 *
 *
 * \param current_global_part_weights is the array holding the weight of parts. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param current_local_part_weights is the local totalweight of the processor.
 * \param current_part_target_weights are the desired cumulative part ratios, sized P.
 * \param current_cut_line_determined is the boolean array to determine if the correct position for a cut line is found.
 *
 * \param current_cut_coordinates is the array holding the coordinates of each cut line. Sized P - 1.
 * \param current_cut_upper_bounds is the array holding the upper bound coordinate for each cut line. Sized P - 1.
 * \param current_cut_lower_bounds is the array holding the lower bound coordinate for each cut line. Sized P - 1.
 * \param current_global_left_closest_points is the array holding the closest points to the cut lines from left.
 * \param current_global_right_closest_points is the array holding the closest points to the cut lines from right.
 * \param current_cut_lower_bound_weights is the array holding the weight of the parts at the left of lower bound coordinates.
 * \param current_cut_upper_weights is the array holding the weight of the parts at the left of upper bound coordinates.
 * \param new_current_cut_coordinates is the work array, sized P - 1.
 *
 * \param current_part_cut_line_weight_ratio holds how much weight of the coordinates on the cutline should be put on left side.
 * \param rectelinear_cut_count is the count of cut lines whose balance can be achived via distributing the points in same coordinate to different parts.
 * \param my_num_incomplete_cut is the number of cutlines whose position has not been determined yet. For K > 1 it is the count in a single part (whose cut lines are determined).
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_get_new_cut_coordinates(
		const size_t &num_total_part,
		const partId_t &num_cuts,
		const pq_scalar_t &max_coordinate,
		const pq_scalar_t &min_coordinate,
		const pq_scalar_t &global_total_weight,
		const pq_scalar_t &used_imbalance_tolerance,
		pq_scalar_t * current_global_part_weights,
		const pq_scalar_t * current_local_part_weights,
		const pq_scalar_t *current_part_target_weights,
		bool *current_cut_line_determined,
		pq_scalar_t *current_cut_coordinates,
		pq_scalar_t *current_cut_upper_bounds,
		pq_scalar_t *current_cut_lower_bounds,
		pq_scalar_t *current_global_left_closest_points,
		pq_scalar_t *current_global_right_closest_points,
		pq_scalar_t * current_cut_lower_bound_weights,
		pq_scalar_t * current_cut_upper_weights,
		pq_scalar_t *new_current_cut_coordinates,
		pq_scalar_t *current_part_cut_line_weight_to_put_left,
		partId_t *rectelinear_cut_count,
		partId_t &my_num_incomplete_cut){

	//seen weight in the part
	pq_scalar_t seen_weight_in_part = 0;
	//expected weight for part.
	pq_scalar_t expected_weight_in_part = 0;
	//imbalance for the left and right side of the cut.
	pq_scalar_t imbalance_on_left = 0, imbalance_on_right = 0;


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
	for (partId_t i = 0; i < num_cuts; i++){
		//if left and right closest points are not set yet,
		//set it to the cut itself.
		if(min_coordinate - current_global_left_closest_points[i] > this->sEpsilon)
			current_global_left_closest_points[i] = current_cut_coordinates[i];
		if(current_global_right_closest_points[i] - max_coordinate > this->sEpsilon)
			current_global_right_closest_points[i] = current_cut_coordinates[i];

	}

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
	for (partId_t i = 0; i < num_cuts; i++){

		if(this->distribute_points_on_cut_lines){
			//init the weight on the cut.
			this->global_rectelinear_cut_weight[i] = 0;
			this->process_rectelinear_cut_weight[i] = 0;
		}
		//if already determined at previous iterations,
		//then just write the coordinate to new array, and proceed.
		if(current_cut_line_determined[i]) {
			new_current_cut_coordinates[i] = current_cut_coordinates[i];
			continue;
		}

		//current weight of the part at the left of the cut line.
		seen_weight_in_part = current_global_part_weights[i * 2];

		/*
		cout << "seen_weight_in_part:" << i << " is "<< seen_weight_in_part << endl;
		cout << "\tcut:" << current_cut_coordinates[i]
		       << " current_cut_lower_bounds:" << current_cut_lower_bounds[i]
               << " current_cut_upper_bounds:" << current_cut_upper_bounds[i] << endl;
               */
		//expected ratio
		expected_weight_in_part = current_part_target_weights[i];
		//leftImbalance = imbalanceOf(seenW, globalTotalWeight, expected);
		imbalance_on_left = imbalanceOf2(seen_weight_in_part, expected_weight_in_part);
		//rightImbalance = imbalanceOf(globalTotalWeight - seenW, globalTotalWeight, 1 - expected);
		imbalance_on_right = imbalanceOf2(global_total_weight - seen_weight_in_part, global_total_weight - expected_weight_in_part);

		bool is_left_imbalance_valid = ABS(imbalance_on_left) - used_imbalance_tolerance < this->sEpsilon ;
		bool is_right_imbalance_valid = ABS(imbalance_on_right) - used_imbalance_tolerance < this->sEpsilon;

		//if the cut line reaches to desired imbalance.
		if(is_left_imbalance_valid && is_right_imbalance_valid){
			current_cut_line_determined[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
			my_num_incomplete_cut -= 1;
			new_current_cut_coordinates [i] = current_cut_coordinates[i];
			continue;
		}
		else if(imbalance_on_left < 0){
			//if left imbalance < 0 then we need to move the cut to right.

			if(this->distribute_points_on_cut_lines){
				//if it is okay to distribute the coordinate on
				//the same coordinate to left and right.
				//then check if we can reach to the target weight by including the
				//coordinates in the part.
				if (current_global_part_weights[i * 2 + 1] == expected_weight_in_part){
					//if it is we are done.
					current_cut_line_determined[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
					my_num_incomplete_cut -= 1;

					//then assign everything on the cut to the left of the cut.
					new_current_cut_coordinates [i] = current_cut_coordinates[i];

					//for this cut all the weight on cut will be put to left.

					current_part_cut_line_weight_to_put_left[i] = current_local_part_weights[i * 2 + 1] - current_local_part_weights[i * 2];
					continue;
				}
				else if (current_global_part_weights[i * 2 + 1] > expected_weight_in_part){

					//if the weight is larger than the expected weight,
					//then we need to distribute some points to left, some to right.
					current_cut_line_determined[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
					*rectelinear_cut_count += 1;
					//increase the num cuts to be determined with rectelinear partitioning.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
					my_num_incomplete_cut -= 1;
					new_current_cut_coordinates [i] = current_cut_coordinates[i];
					this->process_rectelinear_cut_weight[i] = current_local_part_weights[i * 2 + 1] -
							current_local_part_weights[i * 2];
					continue;
				}
			}
			//we need to move further right,so set lower bound to current line, and shift it to the closes point from right.
			current_cut_lower_bounds[i] = current_global_right_closest_points[i];
			//set the lower bound weight to the weight we have seen.
			current_cut_lower_bound_weights[i] = seen_weight_in_part;

			//compare the upper bound with what has been found in the last iteration.
			//we try to make more strict bounds for the cut here.
			for (partId_t ii = i + 1; ii < num_cuts ; ++ii){
				pq_scalar_t p_weight = current_global_part_weights[ii * 2];
				pq_scalar_t line_weight = current_global_part_weights[ii * 2 + 1];

				if(p_weight >= expected_weight_in_part){
					//if a cut on the right has the expected weight, then we found
					//our cut position. Set up and low coordiantes to this new cut coordinate.
					//but we need one more iteration to finalize the cut position,
					//as wee need to update the part ids.
					if(p_weight == expected_weight_in_part){
						current_cut_upper_bounds[i] = current_cut_coordinates[ii];
						current_cut_upper_weights[i] = p_weight;
						current_cut_lower_bounds[i] = current_cut_coordinates[ii];
						current_cut_lower_bound_weights[i] = p_weight;
					} else if (p_weight < current_cut_upper_weights[i]){
						//if a part weight is larger then my expected weight,
						//but lower than my upper bound weight, update upper bound.
						current_cut_upper_bounds[i] = current_global_left_closest_points[ii];
						current_cut_upper_weights[i] = p_weight;
					}
					break;
				}
				//if comes here then pw < ew
				//then compare the weight against line weight.
				if(line_weight >= expected_weight_in_part){
					//if the line is larger than the expected weight,
					//then we need to reach to the balance by distributing coordinates on this line.
					current_cut_upper_bounds[i] = current_cut_coordinates[ii];
					current_cut_upper_weights[i] = line_weight;
					current_cut_lower_bounds[i] = current_cut_coordinates[ii];
					current_cut_lower_bound_weights[i] = p_weight;
					break;
				}
				//if a stricter lower bound is found,
				//update the lower bound.
				if (p_weight <= expected_weight_in_part && p_weight >= current_cut_lower_bound_weights[i]){
					current_cut_lower_bounds[i] = current_global_right_closest_points[ii] ;
					current_cut_lower_bound_weights[i] = p_weight;
				}
			}


			pq_scalar_t new_cut_position = 0;
			this->mj_calculate_new_cut_position(
					current_cut_upper_bounds[i],
					current_cut_lower_bounds[i],
					current_cut_upper_weights[i],
					current_cut_lower_bound_weights[i],
					expected_weight_in_part, new_cut_position);

			//if cut line does not move significantly.
			//then finalize the search.
			if (ABS(current_cut_coordinates[i] - new_cut_position) < this->sEpsilon
				/*|| current_cut_lower_bounds[i] - current_cut_upper_bounds[i] > this->sEpsilon*/
				){
				current_cut_line_determined[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
				my_num_incomplete_cut -= 1;

				//set the cut coordinate and proceed.
				new_current_cut_coordinates [i] = current_cut_coordinates[i];
			} else {
				new_current_cut_coordinates [i] = new_cut_position;
			}
		} else {

			//need to move the cut line to left.
			//set upper bound to current line.
			current_cut_upper_bounds[i] = current_global_left_closest_points[i];
			current_cut_upper_weights[i] = seen_weight_in_part;

			// compare the current cut line weights with previous upper and lower bounds.
			for (int ii = i - 1; ii >= 0; --ii){
				pq_scalar_t p_weight = current_global_part_weights[ii * 2];
				pq_scalar_t line_weight = current_global_part_weights[ii * 2 + 1];
				if(p_weight <= expected_weight_in_part){
					if(p_weight == expected_weight_in_part){
						//if the weight of the part is my expected weight
						//then we find the solution.
						current_cut_upper_bounds[i] = current_cut_coordinates[ii];
						current_cut_upper_weights[i] = p_weight;
						current_cut_lower_bounds[i] = current_cut_coordinates[ii];
						current_cut_lower_bound_weights[i] = p_weight;
					}
					else if (p_weight > current_cut_lower_bound_weights[i]){
						//if found weight is bigger than the lower bound
						//then update the lower bound.
						current_cut_lower_bounds[i] = current_global_right_closest_points[ii];
						current_cut_lower_bound_weights[i] = p_weight;

						//at the same time, if weight of line is bigger than the
						//expected weight, then update the upper bound as well.
						//in this case the balance will be obtained by distributing weightss
						//on this cut position.
						if(line_weight > expected_weight_in_part){
							current_cut_upper_bounds[i] = current_global_right_closest_points[ii];
							current_cut_upper_weights[i] = line_weight;
						}
					}
					break;
				}
				//if the weight of the cut on the left is still bigger than my weight,
				//and also if the weight is smaller than the current upper weight,
				//or if the weight is equal to current upper weight, but on the left of
				// the upper weight, then update upper bound.
				if (p_weight >= expected_weight_in_part &&
						(p_weight < current_cut_upper_weights[i] ||
								(p_weight == current_cut_upper_weights[i] &&
										current_cut_upper_bounds[i] > current_global_left_closest_points[ii]
								)
						)
					){
					current_cut_upper_bounds[i] = current_global_left_closest_points[ii] ;
					current_cut_upper_weights[i] = p_weight;
				}
			}
			pq_scalar_t new_cut_position = 0;
			this->mj_calculate_new_cut_position(
					current_cut_upper_bounds[i],
					current_cut_lower_bounds[i],
					current_cut_upper_weights[i],
					current_cut_lower_bound_weights[i],
					expected_weight_in_part,
					new_cut_position);

			//if cut line does not move significantly.
			if (ABS(current_cut_coordinates[i] - new_cut_position) < this->sEpsilon
					/*|| current_cut_lower_bounds[i] - current_cut_upper_bounds[i] > this->sEpsilon*/ ){
				current_cut_line_determined[i] = true;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
				my_num_incomplete_cut -= 1;
				//set the cut coordinate and proceed.
				new_current_cut_coordinates [ i] = current_cut_coordinates[i];
			} else {
				new_current_cut_coordinates [ i] = new_cut_position;
			}
		}
	}

	//communication to determine the ratios of processors for the distribution
	//of coordinates on the cut lines.
#ifdef HAVE_ZOLTAN2_OMP
	//no need barrier here as it is implicit.
#pragma omp single
#endif
	{
		if(*rectelinear_cut_count > 0){

			try{
				Teuchos::scan<int,pq_scalar_t>(
						*comm, Teuchos::REDUCE_SUM,
						num_cuts,
						this->process_rectelinear_cut_weight,
						this->global_rectelinear_cut_weight
				);
			}
			Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))

			for (partId_t i = 0; i < num_cuts; ++i){
				//if cut line weight to be distributed.
				if(this->global_rectelinear_cut_weight[i] > 0) {
					//expected weight to go to left of the cut.
					pq_scalar_t expected_part_weight = current_part_target_weights[i];
					//the weight that should be put to left of the cut.
					pq_scalar_t necessary_weight_on_line_for_left = expected_part_weight - current_global_part_weights[i * 2];
					//the weight of the cut in the process
					pq_scalar_t my_weight_on_line = this->process_rectelinear_cut_weight[i];
					//the sum of the cut weights upto this process, including the weight of this process.
					pq_scalar_t weight_on_line_upto_process_inclusive = this->global_rectelinear_cut_weight[i];
					//the space on the left side of the cut after all processes before this process (including this process)
					//puts their weights on cut to left.
					pq_scalar_t space_to_put_left = necessary_weight_on_line_for_left - weight_on_line_upto_process_inclusive;
					//add my weight to this space to find out how much space is left to me.
					pq_scalar_t space_left_to_me = space_to_put_left + my_weight_on_line;
					if(space_left_to_me < 0){
						//space_left_to_me is negative and i dont need to put anything to left.
						current_part_cut_line_weight_to_put_left[i] = 0;
					}
					else if(space_left_to_me >= my_weight_on_line){
						//space left to me is bigger than the weight of the processor on cut.
						//so put everything to left.
						current_part_cut_line_weight_to_put_left[i] = my_weight_on_line;
					}
					else {
						//put only the weight as much as the space.
						current_part_cut_line_weight_to_put_left[i] = space_left_to_me ;
					}

				}
			}
			*rectelinear_cut_count = 0;
		}
	}
}

/*! \brief Function fills up the num_points_in_all_processor_parts, so that
 * it has the number of coordinates in each processor of each part.
 * to access how many points processor i has on part j, num_points_in_all_processor_parts[i * num_parts + j].
 *
 * \param num_procs is the number of processor attending to migration operation.
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param num_points_in_all_processor_parts is the output array that holds
 * the number of coordinates in each part in each processor.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::get_processor_num_points_in_parts(
    partId_t num_procs,
    partId_t num_parts,
    pq_gno_t *&num_points_in_all_processor_parts){

	//initially allocation_size is num_parts
	size_t allocation_size = num_parts * (num_procs + 1);

	//this will be output
	//holds how many each processor has in each part.
	//last portion is the sum of all processor points in each part.

	//allocate memory for the local num coordinates in each part.
	pq_gno_t *num_local_points_in_each_part_to_reduce_sum = allocMemory<pq_gno_t>(allocation_size);


	//this is the portion of the memory which will be used
	//at the summation to obtain total number of processors' points in each part.
	pq_gno_t *my_local_points_to_reduce_sum = num_local_points_in_each_part_to_reduce_sum + num_procs * num_parts;
	//this is the portion of the memory where each stores its local number.
	//this information is needed by other processors.
	pq_gno_t *my_local_point_counts_in_each_art = num_local_points_in_each_part_to_reduce_sum + this->myRank * num_parts;

	//initialize the array with 0's.
	memset(num_local_points_in_each_part_to_reduce_sum, 0, sizeof(pq_gno_t)*allocation_size);

	//write the number of coordinates in each part.
	for (partId_t i = 0; i < num_parts; ++i){
		pq_lno_t part_begin_index = 0;
		if (i > 0){
			part_begin_index = this->new_part_xadj[i - 1];
		}
		pq_lno_t part_end_index = this->new_part_xadj[i];
		my_local_points_to_reduce_sum[i] = part_end_index - part_begin_index;
	}

	//copy the local num parts to the last portion of array,
	//so that this portion will represent the global num points in each part after the reduction.
	memcpy (my_local_point_counts_in_each_art,
			my_local_points_to_reduce_sum,
			sizeof(pq_gno_t) * (num_parts) );


	//reduceAll operation.
	//the portion that belongs to a processor with index p
	//will start from myRank * num_parts.
	//the global number of points will be held at the index
	try{
		reduceAll<int, pq_gno_t>(
				*(this->comm),
				Teuchos::REDUCE_SUM,
				allocation_size,
				num_local_points_in_each_part_to_reduce_sum,
				num_points_in_all_processor_parts);
	}
	Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))
	freeArray<pq_gno_t>(num_local_points_in_each_part_to_reduce_sum);
}



/*! \brief Function checks if should do migration or not.
 * It returns true to point that migration should be done when
 * -migration_reduce_all_population are higher than a predetermined value
 * -num_coords_for_last_dim_part that left for the last dimension partitioning is less than a predetermined value
 * -the imbalance of the processors on the parts are higher than given threshold.
 * \param migration_reduce_all_population is the multiplication of the number of reduceall operations estimated and the number of processors.
 * \param num_coords_for_last_dim_part is the estimated number of coordinates in a part per processor in the last dimension partitioning.
 * \param num_procs is the number of processor attending to migration operation.
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param num_points_in_all_processor_parts is the input array that holds
 * the number of coordinates in each part in each processor.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
bool AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_check_to_migrate(
    size_t migration_reduce_all_population,
    pq_lno_t num_coords_for_last_dim_part,
    partId_t num_procs,
    partId_t num_parts,
    pq_gno_t *num_points_in_all_processor_parts){

	//if reduce all count and population in the last dim is too high
    if (migration_reduce_all_population > FUTURE_REDUCEALL_CUTOFF) return true;
    //if the work in a part per processor in the last dim is too low.
    if (num_coords_for_last_dim_part < MIN_WORK_LAST_DIM) return true;

	//if migration is to be checked and the imbalance is too high
    if (this->check_migrate_avoid_migration_option == 0){
    	double global_imbalance = 0;
    	//global shift to reach the sum of coordiante count in each part.
    	size_t global_shift = num_procs * num_parts;

    	for (partId_t ii = 0; ii < num_procs; ++ii){
    		for (partId_t i = 0; i < num_parts; ++i){
    			double ideal_num = num_points_in_all_processor_parts[global_shift + i]
    								/ double(num_procs);

    			global_imbalance += ABS(ideal_num -
    					num_points_in_all_processor_parts[ii * num_parts + i]) /  (ideal_num);
    		}
    	}
    	global_imbalance /= num_parts;
    	global_imbalance /= num_procs;

		/*
    	if (this->myRank == 0) {
    		cout << "imbalance for next iteration:" << global_imbalance << endl;
    	}
    	*/

    	if(global_imbalance <= this->minimum_migration_imbalance){
    		return false;
    	}
    	else {
    		return true;
    	}
    }
    else {
    	//if migration is forced
        return true;
    }
}


/*! \brief Function fills up coordinate_destionations is the output array
 * that holds which part each coordinate should be sent.
 *
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param part_assignment_proc_begin_indices ([i]) points to the first processor index that part i will be sent to.
 * \param processor_chains_in_parts the array that holds the linked list structure, started from part_assignment_proc_begin_indices ([i]).
 * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
 * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::assign_send_destinations(
    partId_t num_parts,
    partId_t *part_assignment_proc_begin_indices,
    partId_t *processor_chains_in_parts,
    pq_lno_t *send_count_to_each_proc,
    int *coordinate_destionations){

    for (partId_t p = 0; p < num_parts; ++p){
        pq_lno_t part_begin = 0;
        if (p > 0) part_begin = this->new_part_xadj[p - 1];
        pq_lno_t part_end = this->new_part_xadj[p];

        //get the first part that current processor will send its part-p.
        partId_t proc_to_sent = part_assignment_proc_begin_indices[p];
        //initialize how many point I sent to this processor.
        pq_lno_t num_total_send = 0;
        for (pq_lno_t j=part_begin; j < part_end; j++){
            pq_lno_t local_ind = this->new_coordinate_permutations[j];
            while (num_total_send >= send_count_to_each_proc[proc_to_sent]){
            	//then get the next processor to send the points in part p.
                num_total_send = 0;
                //assign new processor to part_assign_begin[p]
                part_assignment_proc_begin_indices[p] = processor_chains_in_parts[proc_to_sent];
                //remove the previous processor
                processor_chains_in_parts[proc_to_sent] = -1;
                //choose the next processor as the next one to send.
                proc_to_sent = part_assignment_proc_begin_indices[p];
            }
            //write the gno index to corresponding position in sendBuf.
            coordinate_destionations[local_ind] = proc_to_sent;
            ++num_total_send;
        }
    }
}

/*! \brief Function fills up coordinate_destionations is the output array
 * that holds which part each coordinate should be sent.
 *
 * \param num_points_in_all_processor_parts is the array holding the num points in each part in each proc.
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param num_procs is the number of processor attending to migration operation.

 * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
 * \param processor_ranks_for_subcomm is the ranks of the processors that will be in the subcommunicator with me.
 * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
 * \param out_part_index is the index of the part to which the processor is assigned.
 * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
 * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_assign_proc_to_parts(
		pq_gno_t * num_points_in_all_processor_parts,
		partId_t num_parts,
		partId_t num_procs,
		pq_lno_t *send_count_to_each_proc,
		vector<partId_t> &processor_ranks_for_subcomm,
		vector<partId_t> *next_future_num_parts_in_parts,
		partId_t &out_part_index,
		partId_t &output_part_numbering_begin_index,
		int *coordinate_destionations){


    pq_gno_t *global_num_points_in_parts = num_points_in_all_processor_parts + num_procs * num_parts;
    partId_t *num_procs_assigned_to_each_part = allocMemory<partId_t>(num_parts);

    //boolean variable if the process finds its part to be assigned.
    bool did_i_find_my_group = false;

    partId_t num_free_procs = num_procs;
    partId_t minimum_num_procs_required_for_rest_of_parts = num_parts - 1;

    double max_imbalance_difference = 0;
    partId_t max_differing_part = 0;

    //find how many processor each part requires.
    for (partId_t i=0; i < num_parts; i++){

        //scalar portion of the required processors
        double scalar_required_proc = num_procs *
                (double (global_num_points_in_parts[i]) / double (this->num_global_coords));

        //round it to closest integer.
        partId_t required_proc = static_cast<partId_t> (0.5 + scalar_required_proc);

        //if assigning the required num procs, creates problems for the rest of the parts.
        //then only assign {num_free_procs - (minimum_num_procs_required_for_rest_of_parts)} procs to this part.
        if (num_free_procs - required_proc < minimum_num_procs_required_for_rest_of_parts){
            required_proc = num_free_procs - (minimum_num_procs_required_for_rest_of_parts);
        }

        //reduce the free processor count
        num_free_procs -= required_proc;
        //reduce the free minimum processor count required for the rest of the part by 1.
        --minimum_num_procs_required_for_rest_of_parts;

        //part (i) is assigned to (required_proc) processors.
        num_procs_assigned_to_each_part[i] = required_proc;

        //because of the roundings some processors might be left as unassigned.
        //we want to assign those processors to the part with most imbalance.
        //find the part with the maximum imbalance here.
        double imbalance_wrt_ideal = (scalar_required_proc - required_proc) /  required_proc;
        if (imbalance_wrt_ideal > max_imbalance_difference){
            max_imbalance_difference = imbalance_wrt_ideal;
            max_differing_part = i;
        }
    }

    //assign extra processors to the part with maximum imbalance than the ideal.
    if (num_free_procs > 0){
        num_procs_assigned_to_each_part[max_differing_part] +=  num_free_procs;
    }

    //now find what are the best processors with least migration for each part.

    //part_assignment_proc_begin_indices ([i]) is the array that holds the beginning
    //index of a processor that processor sends its data for part - i
    partId_t *part_assignment_proc_begin_indices = allocMemory<partId_t>(num_parts);
    //the next processor send is found in processor_chains_in_parts, in linked list manner.
    partId_t *processor_chains_in_parts = allocMemory<partId_t>(num_procs);
    partId_t *processor_part_assignments = allocMemory<partId_t>(num_procs);

    //initialize the assignment of each processor.
    //this has a linked list implementation.
    //the beginning of processors assigned
    //to each part is hold at  part_assignment_proc_begin_indices[part].
    //then the next processor assigned to that part is located at
    //proc_part_assignments[part_assign_begins[part]], this is a chain
    //until the value of -1 is reached.
    for (int i = 0; i < num_procs; ++i ){
        processor_part_assignments[i] = -1;
        processor_chains_in_parts[i] = -1;
    }
    for (int i = 0; i < num_parts; ++i ){
        part_assignment_proc_begin_indices[i] = -1;
    }


    //Allocate memory for sorting data structure.
    uSortItem<partId_t, pq_lno_t> * sort_item_num_part_points_in_procs = allocMemory <uSortItem<partId_t, partId_t> > (num_procs);
    for(partId_t i = 0; i < num_parts; ++i){
        //the algorithm tries to minimize the cost of migration,
    	//by assigning the processors with highest number of coordinates on that part.
    	//here we might want to implement a maximum weighted bipartite matching algorithm.
    	for(partId_t ii = 0; ii < num_procs; ++ii){
    		sort_item_num_part_points_in_procs[ii].id = ii;
    		//if processor is not assigned yet.
    		//add its num points to the sort data structure.
    		if (processor_part_assignments[ii] == -1){
    			sort_item_num_part_points_in_procs[ii].val =
    					num_points_in_all_processor_parts[ii * num_parts + i];
    		}
    		else {
    			//if processor is already assigned, insert -nLocal - 1 so that it won't be selected again.
    			//would be same if we simply set it to -1,
    			//but more information with no extra cost (which is used later) is provided.
    			sort_item_num_part_points_in_procs[ii].val = -num_points_in_all_processor_parts[ii * num_parts + i] - 1;
    		}
    	}
        //sort the processors in the part.
        uqsort<partId_t, partId_t>(num_procs, sort_item_num_part_points_in_procs);

        partId_t required_proc_count =  num_procs_assigned_to_each_part[i];
        pq_gno_t total_num_points_in_part = global_num_points_in_parts[i];
        pq_gno_t ideal_num_points_in_a_proc =
                ceil (total_num_points_in_part / double (required_proc_count));

        //starts sending to least heaviest part.
        partId_t next_proc_to_send_index = num_procs - required_proc_count;
        partId_t next_proc_to_send_id = sort_item_num_part_points_in_procs[next_proc_to_send_index].id;
        pq_lno_t space_left_in_sent_proc =  ideal_num_points_in_a_proc - sort_item_num_part_points_in_procs[next_proc_to_send_index].val;

        //find the processors that will be assigned to this part, which are the heaviest
        //non assigned processors.
        for(partId_t ii = num_procs - 1; ii >= num_procs - required_proc_count; --ii){
            partId_t proc_id = sort_item_num_part_points_in_procs[ii].id;
            //assign processor to part - i.
            processor_part_assignments[proc_id] = i;
        }

        bool did_change_sign = false;
        //if processor has a minus count, reverse it.
        for(partId_t ii = 0; ii < num_procs; ++ii){
            if (sort_item_num_part_points_in_procs[ii].val < 0){
                did_change_sign = true;
                sort_item_num_part_points_in_procs[ii].val = -sort_item_num_part_points_in_procs[ii].val - 1;
            }
            else {
                break;
            }
        }
        if(did_change_sign){
            //resort the processors in the part for the rest of the processors that is not assigned.
            uqsort<partId_t, partId_t>(num_procs - required_proc_count, sort_item_num_part_points_in_procs);
        }

        //check if this processors is one of the procs assigned to this part.
        //if it is, then get the group.
        if (!did_i_find_my_group){
            for(partId_t ii = num_procs - 1; ii >= num_procs - required_proc_count; --ii){

            	partId_t proc_id_to_assign = sort_item_num_part_points_in_procs[ii].id;
            	//add the proc to the group.
                processor_ranks_for_subcomm.push_back(proc_id_to_assign);

                if(proc_id_to_assign == this->myRank){
                	//if the assigned process is me, then I find my group.
                    did_i_find_my_group = true;
                    //set the beginning of part i to my rank.
                    part_assignment_proc_begin_indices[i] = this->myRank;
                    processor_chains_in_parts[this->myRank] = -1;

                    //set send count to myself to the number of points that I have in part i.
                    send_count_to_each_proc[this->myRank] = sort_item_num_part_points_in_procs[ii].val;

                    //calculate the shift required for the output_part_numbering_begin_index
                    for (partId_t in = 0; in < i; ++in){
                        output_part_numbering_begin_index += (*next_future_num_parts_in_parts)[in];
                    }
                    out_part_index = i;
                }
            }
            //if these was not my group,
            //clear the subcomminicator processor array.
            if (!did_i_find_my_group){
                processor_ranks_for_subcomm.clear();
            }
        }

        //send points of the nonassigned coordinates to the assigned coordinates.
        //starts from the heaviest nonassigned processor.
        //TODO we might want to play with this part, that allows more computational imbalance
        //but having better communication balance.
        for(partId_t ii = num_procs - required_proc_count - 1; ii >= 0; --ii){
            partId_t nonassigned_proc_id = sort_item_num_part_points_in_procs[ii].id;
            pq_lno_t num_points_to_sent = sort_item_num_part_points_in_procs[ii].val;

            //we set number of points to -to_sent - 1 for the assigned processors.
            //we reverse it here. This should not happen, as we have already reversed them above.
#ifdef MJ_DEBUG
            if (num_points_to_sent < 0) {
            	cout << "Migration - processor assignments - for part:" << i << "from proc:" << nonassigned_proc_id << " num_points_to_sent:" << num_points_to_sent << endl;
            	exit(1);
            }
#endif

            //now sends the points to the assigned processors.
            while (num_points_to_sent > 0){
                //if the processor has enough space.
                if (num_points_to_sent <= space_left_in_sent_proc){
                	//reduce the space left in the processor.
                	space_left_in_sent_proc -= num_points_to_sent;
                	//if my rank is the one that is sending the coordinates.
                    if (this->myRank == nonassigned_proc_id){
                        //set my sent count to the sent processor.
                        send_count_to_each_proc[next_proc_to_send_id] = num_points_to_sent;
                        //save the processor in the list (processor_chains_in_parts and part_assignment_proc_begin_indices)
                        //that the processor will send its point in part-i.
                        partId_t prev_begin = part_assignment_proc_begin_indices[i];
                        part_assignment_proc_begin_indices[i] = next_proc_to_send_id;
                        processor_chains_in_parts[next_proc_to_send_id] = prev_begin;
                    }
                    num_points_to_sent = 0;
                }
                else {
                    //there might be no space left in the processor.
                    if(space_left_in_sent_proc > 0){
                        num_points_to_sent -= space_left_in_sent_proc;

                        //send as the space left in the processor.
                        if (this->myRank == nonassigned_proc_id){
                        	//send as much as the space in this case.
                            send_count_to_each_proc[next_proc_to_send_id] = space_left_in_sent_proc;
                            partId_t prev_begin = part_assignment_proc_begin_indices[i];
                            part_assignment_proc_begin_indices[i] = next_proc_to_send_id;
                            processor_chains_in_parts[next_proc_to_send_id] = prev_begin;

                        }
                    }
                    //change the sent part
                    ++next_proc_to_send_index;

#ifdef MJ_DEBUG
                    if(next_part_to_send_index <  nprocs - required_proc_count ){
                    	cout << "Migration - processor assignments - for part:"
                    			<< i
                    			<<  " next_part_to_send :" << next_part_to_send_index
                    			<< " nprocs:" << nprocs
                    			<< " required_proc_count:" << required_proc_count
                    			<< " Error: next_part_to_send_index <  nprocs - required_proc_count" << endl;
                    	exit(1)l

                    }
#endif
                    //send the new id.
                    next_proc_to_send_id =  sort_item_num_part_points_in_procs[next_proc_to_send_index].id;
                    //set the new space in the processor.
                    space_left_in_sent_proc = ideal_num_points_in_a_proc - sort_item_num_part_points_in_procs[next_proc_to_send_index].val;
                }
            }
        }
    }



    this->assign_send_destinations(
    		num_parts,
            part_assignment_proc_begin_indices,
            processor_chains_in_parts,
            send_count_to_each_proc,
            coordinate_destionations);

    freeArray<partId_t>(part_assignment_proc_begin_indices);
    freeArray<partId_t>(processor_chains_in_parts);
    freeArray<partId_t>(processor_part_assignments);
    freeArray<uSortItem<partId_t, partId_t> > (sort_item_num_part_points_in_procs);
    freeArray<partId_t > (num_procs_assigned_to_each_part);

}


/*! \brief Function fills up coordinate_destionations is the output array
 * that holds which part each coordinate should be sent. In addition it calculates
 * the shift amount (output_part_numbering_begin_index) to be done when
 * final numberings of the parts are performed.
 *
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param sort_item_part_to_proc_assignment is the sorted parts with respect to the assigned processors.
 * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
 * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
 * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::assign_send_destionations2(
    partId_t num_parts,
    uSortItem<partId_t, partId_t> * sort_item_part_to_proc_assignment, //input sorted wrt processors
    int *coordinate_destionations,
    partId_t &output_part_numbering_begin_index,
    vector<partId_t> *next_future_num_parts_in_parts){

    partId_t part_shift_amount = output_part_numbering_begin_index;
    partId_t previous_processor = -1;
    for(partId_t i = 0; i < num_parts; ++i){
        partId_t p = sort_item_part_to_proc_assignment[i].id;
        //assigned processors are sorted.
        pq_lno_t part_begin_index = 0;
        if (p > 0) part_begin_index = this->new_part_xadj[p - 1];
        pq_lno_t part_end_index = this->new_part_xadj[p];

        partId_t assigned_proc = sort_item_part_to_proc_assignment[i].val;
        if (this->myRank == assigned_proc && previous_processor != assigned_proc){
            output_part_numbering_begin_index =  part_shift_amount;
        }
        previous_processor = assigned_proc;
        part_shift_amount += (*next_future_num_parts_in_parts)[p];

        for (pq_lno_t j=part_begin_index; j < part_end_index; j++){
            pq_lno_t localInd = this->new_coordinate_permutations[j];
            coordinate_destionations[localInd] = assigned_proc;
        }
    }
}


/*! \brief Function fills up coordinate_destionations is the output array
 * that holds which part each coordinate should be sent. In addition it calculates
 * the shift amount (output_part_numbering_begin_index) to be done when
 * final numberings of the parts are performed.
 *
 * \param num_points_in_all_processor_parts is the array holding the num points in each part in each proc.
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param num_procs is the number of processor attending to migration operation.

 * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
 * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
 * \param out_num_part is the number of parts assigned to the process.
 * \param out_part_indices is the indices of the part to which the processor is assigned.
 * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
 * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_assign_parts_to_procs(
    pq_gno_t * num_points_in_all_processor_parts,
    partId_t num_parts,
    partId_t num_procs,
    pq_lno_t *send_count_to_each_proc, //output: sized nprocs, show the number of send point counts to each proc.
    vector<partId_t> *next_future_num_parts_in_parts,//input how many more partitions the part will be partitioned into.
    partId_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
    vector<partId_t> &out_part_indices, //output: the part indices which the processor is assigned to.
    partId_t &output_part_numbering_begin_index, //output: how much the part number should be shifted when setting the solution
    int *coordinate_destionations){
    out_num_part = 0;

    pq_gno_t *global_num_points_in_parts = num_points_in_all_processor_parts + num_procs * num_parts;
    out_part_indices.clear();

    //to sort the parts that is assigned to the processors.
    //id is the part number, sort value is the assigned processor id.
    uSortItem<partId_t, partId_t> * sort_item_part_to_proc_assignment  = allocMemory <uSortItem<partId_t, partId_t> >(num_parts);
    uSortItem<partId_t, pq_gno_t> * sort_item_num_points_of_proc_in_part_i = allocMemory <uSortItem<partId_t, pq_gno_t> >(num_procs);


    //calculate the optimal number of coordinates that should be assigned to each processor.
    pq_lno_t work_each = pq_lno_t (this->num_global_coords / (double (num_procs)) + 0.5f);
    //to hold the left space as the number of coordinates to the optimal number in each proc.
    pq_lno_t *space_in_each_processor = allocMemory <pq_lno_t>(num_procs);
    //initialize left space in each.
    for (partId_t i = 0; i < num_procs; ++i){
        space_in_each_processor[i] = work_each;
    }

    //we keep track of how many parts each processor is assigned to.
    //because in some weird inputs, it might be possible that some
    //processors is not assigned to any part. Using these variables,
    //we force each processor to have at least one part.
    partId_t *num_parts_proc_assigned = allocMemory <partId_t>(num_procs);
    memset(num_parts_proc_assigned, 0, sizeof(partId_t) * num_procs);
    int empty_proc_count = num_procs;

    //to sort the parts with decreasing order of their coordiantes.
    //id are the part numbers, sort value is the number of points in each.
    uSortItem<partId_t, pq_gno_t> * sort_item_point_counts_in_parts  = allocMemory <uSortItem<partId_t, pq_gno_t> >(num_parts);

    //initially we will sort the parts according to the number of coordinates they have.
    //so that we will start assigning with the part that has the most number of coordinates.
    for (partId_t i = 0; i < num_parts; ++i){
        sort_item_point_counts_in_parts[i].id = i;
        sort_item_point_counts_in_parts[i].val = global_num_points_in_parts[i];
    }
    //sort parts with increasing order of loads.
    uqsort<partId_t, pq_gno_t>(num_parts, sort_item_point_counts_in_parts);


    //assigning parts to the processors
    //traverse the part win decreasing order of load.
    //first assign the heaviest part.
    for (partId_t j = 0; j < num_parts; ++j){
        //sorted with increasing order, traverse inverse.
        partId_t i = sort_item_point_counts_in_parts[num_parts - 1 - j].id;
        //load of the part
        pq_gno_t load = global_num_points_in_parts[i];

        //assigned processors
        partId_t assigned_proc = -1;
        //if not fit best processor.
        partId_t best_proc_to_assign = 0;


        //sort processors with increasing number of points in this part.
        for (partId_t ii = 0; ii < num_procs; ++ii){
            sort_item_num_points_of_proc_in_part_i[ii].id = ii;

            //if there are still enough parts to fill empty processors, than proceed normally.
            //but if empty processor count is equal to the number of part, then
            //we force to part assignments only to empty processors.
            if (empty_proc_count < num_parts - j || num_parts_proc_assigned[ii] == 0){
            	//how many points processor ii has in part i?
            	sort_item_num_points_of_proc_in_part_i[ii].val =  num_points_in_all_processor_parts[ii * num_parts + i];
            }
            else {
            	sort_item_num_points_of_proc_in_part_i[ii].val  = -1;
            }
        }
        uqsort<partId_t, pq_gno_t>(num_procs, sort_item_num_points_of_proc_in_part_i);

        //traverse all processors with decreasing load.
        for (partId_t iii = num_procs - 1; iii >= 0; --iii){
            partId_t ii = sort_item_num_points_of_proc_in_part_i[iii].id;
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

        //if none had enough space, then assign it to best part.
        if (assigned_proc == -1){
            assigned_proc = best_proc_to_assign;
        }

        if (num_parts_proc_assigned[assigned_proc]++ == 0){
        	--empty_proc_count;
        }
        space_in_each_processor[assigned_proc] -= load;
        //to sort later, part-i is assigned to the proccessor - assignment.
        sort_item_part_to_proc_assignment[j].id = i; //part i
        sort_item_part_to_proc_assignment[j].val = assigned_proc; //assigned to processor - assignment.


        //if assigned processor is me, increase the number.
        if (assigned_proc == this->myRank){
            out_num_part++;//assigned_part_count;
            out_part_indices.push_back(i);
        }
        //increase the send to that processor by the number of points in that part.
        //as everyone send their coordiantes in this part to the processor assigned to this part.
        send_count_to_each_proc[assigned_proc] += num_points_in_all_processor_parts[this->myRank * num_parts + i];
    }
    freeArray<partId_t>(num_parts_proc_assigned);
    freeArray< uSortItem<partId_t, pq_gno_t> > (sort_item_num_points_of_proc_in_part_i);
    freeArray<uSortItem<partId_t, pq_gno_t> >(sort_item_point_counts_in_parts);
    freeArray<pq_lno_t >(space_in_each_processor);


    //sort assignments with respect to the assigned processors.
    uqsort<partId_t, partId_t>(num_parts, sort_item_part_to_proc_assignment);
    //fill sendBuf.


    this->assign_send_destionations2(
            num_parts,
            sort_item_part_to_proc_assignment,
            coordinate_destionations,
            output_part_numbering_begin_index,
            next_future_num_parts_in_parts);

    freeArray<uSortItem<partId_t, partId_t> >(sort_item_part_to_proc_assignment);
}


/*! \brief Function fills up coordinate_destionations is the output array
 * that holds which part each coordinate should be sent. In addition it calculates
 * the shift amount (output_part_numbering_begin_index) to be done when
 * final numberings of the parts are performed.
 *
 * \param num_points_in_all_processor_parts is the array holding the num points in each part in each proc.
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param num_procs is the number of processor attending to migration operation.

 * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
 * \param processor_ranks_for_subcomm is the ranks of the processors that will be in the subcommunicator with me.
 * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
 * \param out_num_part is the number of parts assigned to the process.
 * \param out_part_indices is the indices of the part to which the processor is assigned.
 * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
 * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_migration_part_proc_assignment(
    pq_gno_t * num_points_in_all_processor_parts,
    partId_t num_parts,
    partId_t num_procs,
    pq_lno_t *send_count_to_each_proc,
    vector<partId_t> &processor_ranks_for_subcomm,
    vector<partId_t> *next_future_num_parts_in_parts,
    partId_t &out_num_part,
    vector<partId_t> &out_part_indices,
    partId_t &output_part_numbering_begin_index,
    int *coordinate_destionations){



	processor_ranks_for_subcomm.clear();
	// if (this->num_local_coords > 0)
	if (num_procs > num_parts){
		//if there are more processors than the number of current part
		//then processors share the existing parts.
		//at the end each processor will have a single part,
		//but a part will be shared by a group of processors.
		partId_t out_part_index = 0;
		this->mj_assign_proc_to_parts(
				num_points_in_all_processor_parts,
				num_parts,
				num_procs,
				send_count_to_each_proc,
				processor_ranks_for_subcomm,
				next_future_num_parts_in_parts,
				out_part_index,
				output_part_numbering_begin_index,
				coordinate_destionations
		);

		out_num_part = 1;
		out_part_indices.clear();
		out_part_indices.push_back(out_part_index);
	}
	else {

		//there are more parts than the processors.
		//therefore a processor will be assigned multiple parts,
		//the subcommunicators will only have a single processor.
		processor_ranks_for_subcomm.push_back(this->myRank);

		//since there are more parts then procs,
		//assign multiple parts to processors.
		this->mj_assign_parts_to_procs(
				num_points_in_all_processor_parts,
				num_parts,
				num_procs,
				send_count_to_each_proc,
				next_future_num_parts_in_parts,
				out_num_part,
				out_part_indices,
				output_part_numbering_begin_index,
				coordinate_destionations);
	}
}

/*! \brief Function fills up coordinate_destionations is the output array
 * that holds which part each coordinate should be sent. In addition it calculates
 * the shift amount (output_part_numbering_begin_index) to be done when
 * final numberings of the parts are performed.
 *
 *
 * \param num_procs is the number of processor attending to migration operation.
 * \param num_new_local_points is the output to represent the new number of local points.
 * \param iteration is the string for the current iteration.
 * \param coordinate_destionations is the output array that holds which part each coordinate should be sent.
 * \param num_parts is the number of parts that exist in the current partitioning.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_migrate_coords(
    partId_t num_procs,
    pq_gno_t &num_new_local_points,
    string iteration,
    int *coordinate_destionations,
    partId_t num_parts){

	ZOLTAN_COMM_OBJ *plan = NULL;
	MPI_Comm mpi_comm = Teuchos2MPI (this->comm);
	pq_lno_t num_incoming_gnos = 0;
	int message_tag = 7859;

	this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Migration Z1PlanCreating-" + iteration);
	int ierr = Zoltan_Comm_Create(
			&plan,
			this->num_local_coords,
			coordinate_destionations,
			mpi_comm,
			message_tag,
			&num_incoming_gnos);
	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
	this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Migration Z1PlanCreating-" + iteration);

	this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Migration Z1Migration-" + iteration);
	pq_gno_t *incoming_gnos = allocMemory< pq_gno_t>(num_incoming_gnos);

	//migrate gnos.
	message_tag++;
	ierr = Zoltan_Comm_Do(
			plan,
			message_tag,
			(char *) this->current_mj_gnos,
			sizeof(pq_gno_t),
			(char *) incoming_gnos);
	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

	freeArray<pq_gno_t>(this->current_mj_gnos);
	this->current_mj_gnos = incoming_gnos;


	//migrate coordinates
	for (int i = 0; i < this->coord_dim; ++i){
		message_tag++;
		pq_scalar_t *coord = this->mj_coordinates[i];

		this->mj_coordinates[i] = allocMemory<pq_scalar_t>(num_incoming_gnos);
		ierr = Zoltan_Comm_Do(
				plan,
				message_tag,
				(char *) coord,
				sizeof(pq_scalar_t),
				(char *) this->mj_coordinates[i]);
		Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
		freeArray<pq_scalar_t>(coord);
	}

	//migrate weights.
	for (int i = 0; i < this->weight_dim; ++i){
		message_tag++;
		pq_scalar_t *weight = this->mj_weights[i];

		this->mj_weights[i] = allocMemory<pq_scalar_t>(num_incoming_gnos);
		ierr = Zoltan_Comm_Do(
				plan,
				message_tag,
				(char *) weight,
				sizeof(pq_scalar_t),
				(char *) this->mj_weights[i]);
		Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
		freeArray<pq_scalar_t>(weight);
	}


	//migrate owners.
	int *coord_own = allocMemory<int>(num_incoming_gnos);
	message_tag++;
	ierr = Zoltan_Comm_Do(
			plan,
			message_tag,
			(char *) this->owner_of_coordinate,
			sizeof(int), (char *) coord_own);
	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
	freeArray<int>(this->owner_of_coordinate);
	this->owner_of_coordinate = coord_own;


	//if num procs is less than num parts,
	//we need the part assigment arrays as well, since
	//there will be multiple parts in processor.
	partId_t *new_parts = allocMemory<int>(num_incoming_gnos);
	if(num_procs < num_parts){
		message_tag++;
		ierr = Zoltan_Comm_Do(
				plan,
				message_tag,
				(char *) this->assigned_part_ids,
				sizeof(partId_t),
				(char *) new_parts);
		Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
	}
	freeArray<partId_t>(this->assigned_part_ids);
	this->assigned_part_ids = new_parts;

	ierr = Zoltan_Comm_Destroy(&plan);
	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
	num_new_local_points = num_incoming_gnos;
	this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Migration Z1Migration-" + iteration);
}

/*! \brief Function creates the new subcomminicator for the processors
 * given in processor_ranks_for_subcomm.
 *
 * \param processor_ranks_for_subcomm is the vector that has the ranks of
 * the processors that will be in the same group.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::create_sub_communicator(vector<partId_t> &processor_ranks_for_subcomm){
    partId_t group_size = processor_ranks_for_subcomm.size();
    partId_t *ids = allocMemory<partId_t>(group_size);
    for(partId_t i = 0; i < group_size; ++i) {
        ids[i] = processor_ranks_for_subcomm[i];
    }
    ArrayView<const partId_t> idView(ids, group_size);
    this->comm = this->comm->createSubcommunicator(idView);
    freeArray<partId_t>(ids);
}


/*! \brief Function writes the new permutation arrays after the migration.
 *
 * \param output_num_parts is the number of parts that is assigned to the processor.
 * \param num_parts is the number of parts right before migration.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::fill_permutation_array(
    partId_t output_num_parts,
    partId_t num_parts){
	//if there is single output part, then simply fill the permutation array.
    if (output_num_parts == 1){
        for(pq_lno_t i = 0; i < this->num_local_coords; ++i){
            this->new_coordinate_permutations[i] = i;
        }
        this->new_part_xadj[0] = this->num_local_coords;
    }
    else {

    	//otherwise we need to count how many points are there in each part.
    	//we allocate here as num_parts, because the sent partids are up to num_parts,
    	//although there are outout_num_parts different part.
        pq_lno_t *num_points_in_parts = allocMemory<pq_lno_t>(num_parts);
        //part shift holds the which part number an old part number corresponds to.
        partId_t *part_shifts = allocMemory<partId_t>(num_parts);

        memset(num_points_in_parts, 0, sizeof(pq_lno_t) * num_parts);

        for(pq_lno_t i = 0; i < this->num_local_coords; ++i){
            partId_t ii = this->assigned_part_ids[i];
            ++num_points_in_parts[ii];
        }

        //write the end points of the parts.
        partId_t p = 0;
        pq_lno_t prev_index = 0;
        for(partId_t i = 0; i < num_parts; ++i){
            if(num_points_in_parts[i] > 0)  {
                this->new_part_xadj[p] =  prev_index + num_points_in_parts[i];
                prev_index += num_points_in_parts[i];
                part_shifts[i] = p++;
            }
        }

        //for the rest of the parts write the end index as end point.
        partId_t assigned_num_parts = p - 1;
        for (;p < num_parts; ++p){
            this->new_part_xadj[p] =  this->new_part_xadj[assigned_num_parts];
        }
        for(partId_t i = 0; i < output_num_parts; ++i){
            num_points_in_parts[i] = this->new_part_xadj[i];
        }

        //write the permutation array here.
        //get the part of the coordinate i, shift it to obtain the new part number.
        //assign it to the end of the new part numbers pointer.
        for(pq_lno_t i = this->num_local_coords - 1; i >= 0; --i){
            partId_t part = part_shifts[partId_t(this->assigned_part_ids[i])];
            this->new_coordinate_permutations[--num_points_in_parts[part]] = i;
        }

        freeArray<pq_lno_t>(num_points_in_parts);
        freeArray<partId_t>(part_shifts);
    }
}


/*! \brief Function checks if should do migration or not.
 * It returns true to point that migration should be done when
 * -migration_reduce_all_population are higher than a predetermined value
 * -num_coords_for_last_dim_part that left for the last dimension partitioning is less than a predetermined value
 * -the imbalance of the processors on the parts are higher than given threshold.

 * \param input_num_parts is the number of parts when migration is called.
 * \param output_num_parts is the output number of parts after migration.
 * \param next_future_num_parts_in_parts is the number of total future parts each
 * part is partitioned into. This will be updated when migration is performed.
 * \param output_part_begin_index is the number that will be used as beginning part number
 * when final solution part numbers are assigned.
 * \param migration_reduce_all_population is the estimated total number of reduceall operations
 * multiplied with number of processors to be used for determining migration.
 *
 * \param num_coords_for_last_dim_part is the estimated number of points in each part,
 * when last dimension partitioning is performed.
 * \param iteration is the string that gives information about the dimension for printing purposes.
 * \param input_part_boxes is the array that holds the part boxes after the migration. (swapped)
 * \param output_part_boxes is the array that holds the part boxes before the migration. (swapped)
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
bool AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::mj_perform_migration(
    partId_t input_num_parts, //current umb parts
    partId_t &output_num_parts, //output umb parts.
    vector<partId_t> *next_future_num_parts_in_parts,
    partId_t &output_part_begin_index,
    size_t migration_reduce_all_population,
    pq_lno_t num_coords_for_last_dim_part,
    string iteration,
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > &input_part_boxes,
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > &output_part_boxes){

	partId_t num_procs = this->comm->getSize();
	this->myRank = this->comm->getRank();


	//this array holds how many points each processor has in each part.
	//to access how many points processor i has on part j,
	//num_points_in_all_processor_parts[i * num_parts + j]
	pq_gno_t *num_points_in_all_processor_parts = allocMemory<pq_gno_t>(input_num_parts * (num_procs + 1));

	//get the number of coordinates in each part in each processor.
	this->get_processor_num_points_in_parts(
			num_procs,
			input_num_parts,
			num_points_in_all_processor_parts);


	//check if migration will be performed or not.
	if (!this->mj_check_to_migrate(
			migration_reduce_all_population,
			num_coords_for_last_dim_part,
			num_procs,
			input_num_parts,
			num_points_in_all_processor_parts)){
		freeArray<pq_gno_t>(num_points_in_all_processor_parts);
		return false;
	}


	pq_lno_t *send_count_to_each_proc = NULL;
	int *coordinate_destionations = allocMemory<int>(this->num_local_coords);
	send_count_to_each_proc = allocMemory<pq_lno_t>(num_procs);
	for (int i = 0; i < num_procs; ++i) send_count_to_each_proc[i] = 0;

	vector<partId_t> processor_ranks_for_subcomm;
	vector<partId_t> out_part_indices;

	//determine which processors are assigned to which parts
	this->mj_migration_part_proc_assignment(
			num_points_in_all_processor_parts,
			input_num_parts,
			num_procs,
			send_count_to_each_proc,
			processor_ranks_for_subcomm,
			next_future_num_parts_in_parts,
			output_num_parts,
			out_part_indices,
			output_part_begin_index,
			coordinate_destionations);




	freeArray<pq_lno_t>(send_count_to_each_proc);
	vector <partId_t> tmpv;

	std::sort (out_part_indices.begin(), out_part_indices.end());
	partId_t outP = out_part_indices.size();

	pq_gno_t new_global_num_points = 0;
	pq_gno_t *global_num_points_in_parts = num_points_in_all_processor_parts + num_procs * input_num_parts;

	if (this->mj_keep_part_boxes){
		input_part_boxes->clear();
	}

	//now we calculate the new values for next_future_num_parts_in_parts.
	//same for the part boxes.
	for (partId_t i = 0; i < outP; ++i){
		partId_t ind = out_part_indices[i];
		new_global_num_points += global_num_points_in_parts[ind];
		tmpv.push_back((*next_future_num_parts_in_parts)[ind]);
		if (this->mj_keep_part_boxes){
			input_part_boxes->push_back((*output_part_boxes)[ind]);
		}
	}
	//swap the input and output part boxes.
	if (this->mj_keep_part_boxes){
		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > tmpPartBoxes = input_part_boxes;
		input_part_boxes = output_part_boxes;
		output_part_boxes = tmpPartBoxes;
	}
	next_future_num_parts_in_parts->clear();
	for (partId_t i = 0; i < outP; ++i){
		partId_t p = tmpv[i];
		next_future_num_parts_in_parts->push_back(p);
	}

	freeArray<pq_gno_t>(num_points_in_all_processor_parts);

	pq_gno_t num_new_local_points = 0;


	//perform the actual migration operation here.
	this->mj_migrate_coords(
			num_procs,
			num_new_local_points,
			iteration,
			coordinate_destionations,
			input_num_parts);


	freeArray<int>(coordinate_destionations);

	if(this->num_local_coords != num_new_local_points){
		freeArray<pq_lno_t>(this->new_coordinate_permutations);
		freeArray<pq_lno_t>(this->coordinate_permutations);

		this->new_coordinate_permutations = allocMemory<pq_lno_t>(num_new_local_points);
		this->coordinate_permutations = allocMemory<pq_lno_t>(num_new_local_points);
	}
	this->num_local_coords = num_new_local_points;
	this->num_global_coords = new_global_num_points;



	//create subcommunicator.
	this->create_sub_communicator(processor_ranks_for_subcomm);
	processor_ranks_for_subcomm.clear();

	//fill the new permutation arrays.
	this->fill_permutation_array(
			output_num_parts,
			input_num_parts);
	return true;
}


/*! \brief Function creates consistent chunks for task partitioning. Used only in the case of
 * sequential task partitioning, where consistent handle of the points on the cuts are required.
 *
 * \param num_parts is the number of parts.
 * \param mj_current_dim_coords is 1 dimensional array holding the coordinate values.
 * \param current_concurrent_cut_coordinate is 1 dimensional array holding the cut coordinates.
 * \param coordinate_begin is the start index of the given partition on partitionedPointPermutations.
 * \param coordinate_end is the end index of the given partition on partitionedPointPermutations.
 * \param used_local_cut_line_weight_to_left holds how much weight of the coordinates on the cutline should be put on left side.
 *
 * \param out_part_xadj is the indices of begginning and end of the parts in the output partition.
 * \param coordInd is the index according to which the partitioning is done.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::create_consistent_chunks(
    partId_t num_parts,
    pq_scalar_t *mj_current_dim_coords,
    pq_scalar_t *current_concurrent_cut_coordinate,
    pq_lno_t coordinate_begin,
    pq_lno_t coordinate_end,
    pq_scalar_t *used_local_cut_line_weight_to_left,
    pq_lno_t *out_part_xadj,
    int coordInd){

	//pq_lno_t numCoordsInPart =  coordinateEnd - coordinateBegin;
	partId_t no_cuts = num_parts - 1;



	int me = 0;
	pq_lno_t *thread_num_points_in_parts = this->thread_point_counts[me];
	float *my_local_thread_cut_weights_to_put_left = NULL;


	//now if the rectelinear partitioning is allowed we decide how
	//much weight each thread should put to left and right.
	if (this->distribute_points_on_cut_lines){
		my_local_thread_cut_weights_to_put_left = this->thread_cut_line_weight_to_put_left[me];
		for (partId_t i = 0; i < no_cuts; ++i){
			//the left to be put on the left of the cut.
			pq_scalar_t left_weight = used_local_cut_line_weight_to_left[i];
			for(int ii = 0; ii < this->num_threads; ++ii){
				if(left_weight > this->sEpsilon){
					//the weight of thread ii on cut.
					pq_scalar_t thread_ii_weight_on_cut = this->thread_part_weight_work[ii][i * 2 + 1] - this->thread_part_weight_work[ii][i * 2 ];
					if(thread_ii_weight_on_cut < left_weight){
						this->thread_cut_line_weight_to_put_left[ii][i] = thread_ii_weight_on_cut;
					}
					else {
						this->thread_cut_line_weight_to_put_left[ii][i] = left_weight ;
					}
					left_weight -= thread_ii_weight_on_cut;
				}
				else {
					this->thread_cut_line_weight_to_put_left[ii][i] = 0;
				}
			}
		}

		if(no_cuts > 0){
			//this is a special case. If cutlines share the same coordinate, their weights are equal.
			//we need to adjust the ratio for that.
			for (partId_t i = no_cuts - 1; i > 0 ; --i){
				if(ABS(current_concurrent_cut_coordinate[i] - current_concurrent_cut_coordinate[i -1]) < this->sEpsilon){
					my_local_thread_cut_weights_to_put_left[i] -= my_local_thread_cut_weights_to_put_left[i - 1] ;
				}
				my_local_thread_cut_weights_to_put_left[i] = int ((my_local_thread_cut_weights_to_put_left[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL)
                    										/ pq_scalar_t(SIGNIFICANCE_MUL);
			}
		}
	}

	for(partId_t ii = 0; ii < num_parts; ++ii){
		thread_num_points_in_parts[ii] = 0;
	}

	//for this specific case we dont want to distribute the points along the cut position
	//randomly, as we need a specific ordering of them. Instead,
	//we put the coordinates into a sort item, where we sort those
	//using the coordinates of points on other dimensions and the index.


	//some of the cuts might share the same position.
	//in this case, if cut i and cut j share the same position
	//cut_map[i] = cut_map[j] = sort item index.
	partId_t *cut_map = allocMemory<partId_t> (no_cuts);


	typedef uMultiSortItem<pq_lno_t, int, pq_scalar_t> multiSItem;
	typedef std::vector< multiSItem > multiSVector;
	typedef std::vector<multiSVector> multiS2Vector;

	//to keep track of the memory allocated.
	std::vector<pq_scalar_t *>allocated_memory;

	//vector for which the coordinates will be sorted.
	multiS2Vector sort_vector_points_on_cut;

	//the number of cuts that have different coordinates.
	partId_t different_cut_count = 1;
	cut_map[0] = 0;

	//now we insert 1 sort vector for all cuts on the different
	//positins.if multiple cuts are on the same position, they share sort vectors.
	multiSVector tmpMultiSVector;
	sort_vector_points_on_cut.push_back(tmpMultiSVector);

	for (partId_t i = 1; i < no_cuts ; ++i){
		//if cuts share the same cut coordinates
		//set the cutmap accordingly.
		if(ABS(current_concurrent_cut_coordinate[i] - current_concurrent_cut_coordinate[i -1]) < this->sEpsilon){
			cut_map[i] = cut_map[i-1];
		}
		else {
			cut_map[i] = different_cut_count++;
			multiSVector tmp2MultiSVector;
			sort_vector_points_on_cut.push_back(tmp2MultiSVector);
		}
	}


	//now the actual part assigment.
	for (pq_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){

		pq_lno_t i = this->coordinate_permutations[ii];

		partId_t pp = this->assigned_part_ids[i];
		partId_t p = pp / 2;
		//if the coordinate is on a cut.
		if(pp % 2 == 1 ){
			pq_scalar_t *vals = allocMemory<pq_scalar_t>(this->coord_dim -1);
			allocated_memory.push_back(vals);

			//we insert the coordinates to the sort item here.
			int val_ind = 0;
			for(int dim = coordInd + 1; dim < this->coord_dim; ++dim){
				vals[val_ind++] = this->mj_coordinates[dim][i];
			}
			for(int dim = 0; dim < coordInd; ++dim){
				vals[val_ind++] = this->mj_coordinates[dim][i];
			}
			multiSItem tempSortItem(i, this->coord_dim -1, vals);
			//inser the point to the sort vector pointed by the cut_map[p].
			partId_t cmap = cut_map[p];
			sort_vector_points_on_cut[cmap].push_back(tempSortItem);
		}
		else {
			//if it is not on the cut, simple sorting.
			++thread_num_points_in_parts[p];
			this->assigned_part_ids[i] = p;
		}
	}

	//sort all the sort vectors.
	for (partId_t i = 0; i < different_cut_count; ++i){
		std::sort (sort_vector_points_on_cut[i].begin(), sort_vector_points_on_cut[i].end());
	}

	//we do the part assignment for the points on cuts here.
	partId_t previous_cut_map = cut_map[0];

	//this is how much previous part owns the weight of the current part.
	//when target part weight is 1.6, and the part on the left is given 2,
	//the left has an extra 0.4, while the right has missing 0.4 from the previous cut.
	//this parameter is used to balance this issues.
	//in the above example weight_stolen_from_previous_part will be 0.4.
	//if the left part target is 2.2 but it is given 2,
	//then weight_stolen_from_previous_part will be -0.2.
	pq_scalar_t weight_stolen_from_previous_part = 0;
	for (partId_t p = 0; p < no_cuts; ++p){

		partId_t mapped_cut = cut_map[p];

		//if previous cut map is done, and it does not have the same index,
		//then assign all points left on that cut to its right.
		if (previous_cut_map != mapped_cut){
			pq_lno_t sort_vector_end = (pq_lno_t)sort_vector_points_on_cut[previous_cut_map].size() - 1;
			for (; sort_vector_end >= 0; --sort_vector_end){
				multiSItem t = sort_vector_points_on_cut[previous_cut_map][sort_vector_end];
				pq_lno_t i = t.index;
				++thread_num_points_in_parts[p];
				this->assigned_part_ids[i] = p;
			}
			sort_vector_points_on_cut[previous_cut_map].clear();
		}
		pq_lno_t sort_vector_end = (pq_lno_t)sort_vector_points_on_cut[mapped_cut].size() - 1;

		for (; sort_vector_end >= 0; --sort_vector_end){
			multiSItem t = sort_vector_points_on_cut[mapped_cut][sort_vector_end];
			pq_lno_t i = t.index;
			pq_scalar_t w = this->mj_uniform_weights[0]? 1:this->mj_weights[0][i];


			//part p has enough space for point i, then put it to point i.
			if(	my_local_thread_cut_weights_to_put_left[p] + weight_stolen_from_previous_part> this->sEpsilon * EPS_SCALE &&
				my_local_thread_cut_weights_to_put_left[p] + weight_stolen_from_previous_part - ABS(my_local_thread_cut_weights_to_put_left[p] + weight_stolen_from_previous_part - w)
					> this->sEpsilon){

				my_local_thread_cut_weights_to_put_left[p] -= w;
				sort_vector_points_on_cut[mapped_cut].pop_back();
				++thread_num_points_in_parts[p];
				this->assigned_part_ids[i] = p;
				//if putting this weight to left overweights the left cut, then
				//increase the space for the next cut using weight_stolen_from_previous_part.
				if(p < no_cuts - 1 && my_local_thread_cut_weights_to_put_left[p] < this->sEpsilon){
					if(mapped_cut == cut_map[p + 1] ){
						//if the cut before the cut indexed at p was also at the same position
						//special case, as we handle the weight differently here.
						if (previous_cut_map != mapped_cut){
							weight_stolen_from_previous_part = my_local_thread_cut_weights_to_put_left[p];
						}
						else {
							//if the cut before the cut indexed at p was also at the same position
							//we assign extra weights cumulatively in this case.
							weight_stolen_from_previous_part += my_local_thread_cut_weights_to_put_left[p];
						}
					}
					else{
						weight_stolen_from_previous_part = -my_local_thread_cut_weights_to_put_left[p];
					}
					//end assignment for part p
					break;
				}
			} else {
				//if part p does not have enough space for this point
				//and if there is another cut sharing the same positon,
				//again increase the space for the next
				if(p < no_cuts - 1 && mapped_cut == cut_map[p + 1]){
					if (previous_cut_map != mapped_cut){
						weight_stolen_from_previous_part = my_local_thread_cut_weights_to_put_left[p];
					}
					else {
						weight_stolen_from_previous_part += my_local_thread_cut_weights_to_put_left[p];
					}
				}
				else{
					weight_stolen_from_previous_part = -my_local_thread_cut_weights_to_put_left[p];
				}
				//end assignment for part p
				break;
			}
		}
		previous_cut_map = mapped_cut;
	}
	//put everything left on the last cut to the last part.
	pq_lno_t sort_vector_end = (pq_lno_t)sort_vector_points_on_cut[previous_cut_map].size() - 1;
	for (; sort_vector_end >= 0; --sort_vector_end){
		multiSItem t = sort_vector_points_on_cut[previous_cut_map][sort_vector_end];
		pq_lno_t i = t.index;
		++thread_num_points_in_parts[no_cuts];
		this->assigned_part_ids[i] = no_cuts;
	}
	sort_vector_points_on_cut[previous_cut_map].clear();
	freeArray<partId_t> (cut_map);

	//free the memory allocated for vertex sort items .
	pq_lno_t vSize = (pq_lno_t) allocated_memory.size();
	for(pq_lno_t i = 0; i < vSize; ++i){
		freeArray<pq_scalar_t> (allocated_memory[i]);
	}

	//creation of part_xadj as in usual case.
	for(partId_t j = 0; j < num_parts; ++j){
		pq_lno_t num_points_in_part_j_upto_thread_i = 0;
		for (int i = 0; i < this->num_threads; ++i){
			pq_lno_t thread_num_points_in_part_j = this->thread_point_counts[i][j];
			this->thread_point_counts[i][j] = num_points_in_part_j_upto_thread_i;
			num_points_in_part_j_upto_thread_i += thread_num_points_in_part_j;

		}
		out_part_xadj[j] = num_points_in_part_j_upto_thread_i;// + prev2; //+ coordinateBegin;
	}

	//perform prefix sum for num_points in parts.
	for(partId_t j = 1; j < num_parts; ++j){
		out_part_xadj[j] += out_part_xadj[j - 1];
	}


	//shift the num points in threads thread to obtain the
	//beginning index of each thread's private space.
	for(partId_t j = 1; j < num_parts; ++j){
		thread_num_points_in_parts[j] += out_part_xadj[j - 1] ;
	}

	//now thread gets the coordinate and writes the index of coordinate to the permutation array
	//using the part index we calculated.
	for (pq_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){
		pq_lno_t i = this->coordinate_permutations[ii];
		partId_t p =  this->assigned_part_ids[i];
		this->new_coordinate_permutations[coordinate_begin +
		                                  thread_num_points_in_parts[p]++] = i;
	}
}



/*! \brief Function sends the found partids to the owner of the coordinates,
 * if the data is ever migrated. otherwise, it seets the part numbers and returns.
 * \param current_num_parts is the number of parts in the process.
 * \param output_part_begin_index is the number that will be used as beginning part number
 * \param output_part_boxes is the array that holds the part boxes
 * \param is_data_ever_migrated is the boolean value which is true
 * if the data is ever migrated during the partitioning.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::set_final_parts(
		partId_t current_num_parts,
		partId_t output_part_begin_index,
		RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > output_part_boxes,
		bool is_data_ever_migrated){
    this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Part_Assignment");

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(partId_t i = 0; i < current_num_parts;++i){

    	pq_lno_t begin = 0;
    	pq_lno_t end = this->part_xadj[i];

    	if(i > 0) begin = this->part_xadj[i -1];
    	partId_t part_to_set_index = i + output_part_begin_index;
    	if (this->mj_keep_part_boxes){
    		(*output_part_boxes)[i].setpId(part_to_set_index);
    	}
    	for (pq_lno_t ii = begin; ii < end; ++ii){
    		pq_lno_t k = this->coordinate_permutations[ii];
    		this->assigned_part_ids[k] = part_to_set_index;
    	}
    }

    //ArrayRCP<const pq_gno_t> gnoList;
    if(!is_data_ever_migrated){
    	//freeArray<pq_gno_t>(this->current_mj_gnos);
        //if(this->num_local_coords > 0){
        //    gnoList = arcpFromArrayView(this->mj_gnos);
        //}
    }
    else {
    	//if data is migrated, then send part numbers to the original owners.
    	ZOLTAN_COMM_OBJ *plan = NULL;
    	MPI_Comm mpi_comm = Teuchos2MPI (this->mj_problemComm);

    	pq_lno_t incoming = 0;
    	int message_tag = 7856;

    	this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Final Z1PlanCreating");
    	int ierr = Zoltan_Comm_Create( &plan, this->num_local_coords,
    			this->owner_of_coordinate, mpi_comm, message_tag,
    			&incoming);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
    	this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Final Z1PlanCreating" );

    	pq_gno_t *incoming_gnos = allocMemory< pq_gno_t>(incoming);

    	message_tag++;
    	this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Final Z1PlanComm");
    	ierr = Zoltan_Comm_Do( plan, message_tag, (char *) this->current_mj_gnos,
    			sizeof(pq_gno_t), (char *) incoming_gnos);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

    	freeArray<pq_gno_t>(this->current_mj_gnos);
    	this->current_mj_gnos = incoming_gnos;

    	partId_t *incoming_partIds = allocMemory< partId_t>(incoming);

    	message_tag++;
    	ierr = Zoltan_Comm_Do( plan, message_tag, (char *) this->assigned_part_ids,
    			sizeof(partId_t), (char *) incoming_partIds);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
    	freeArray<partId_t>(this->assigned_part_ids);
    	this->assigned_part_ids = incoming_partIds;

    	this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Final Z1PlanComm");
    	ierr = Zoltan_Comm_Destroy(&plan);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

    	this->num_local_coords = incoming;
    	//gnoList = arcp(this->current_mj_gnos, 0, this->num_local_coords, true);
    }

    this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Part_Assignment");

    this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Solution_Part_Assignment");

    //ArrayRCP<partId_t> partId;
    //partId = arcp(this->assigned_part_ids, 0, this->num_local_coords, true);


    if (this->mj_keep_part_boxes){
    	this->kept_boxes = output_part_boxes;
        //this->mj_solution->setPartBoxes(output_part_boxes);
    }

    //this->mj_solution->setParts(gnoList, partId, true);

    this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Solution_Part_Assignment");
}

/*! \brief Function frees all allocated work memory.
*/
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::free_work_memory(){
	this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Free");

	for (int i=0; i < this->coord_dim; i++){

		freeArray<pq_scalar_t>(this->mj_coordinates[i]);
	}
	freeArray<pq_scalar_t *>(this->mj_coordinates);

	for (int i=0; i < this->weight_dim; i++){
		freeArray<pq_scalar_t>(this->mj_weights[i]);
	}
	freeArray<pq_scalar_t *>(this->mj_weights);

	freeArray<int>(this->owner_of_coordinate);

	for(int i = 0; i < this->num_threads; ++i){
		freeArray<pq_lno_t>(this->thread_point_counts[i]);
	}

	freeArray<pq_lno_t *>(this->thread_point_counts);
	freeArray<double *> (this->thread_part_weight_work);

	if(this->distribute_points_on_cut_lines){
		freeArray<pq_scalar_t>(this->process_cut_line_weight_to_put_left);
		for(int i = 0; i < this->num_threads; ++i){
			freeArray<pq_scalar_t>(this->thread_cut_line_weight_to_put_left[i]);
		}
		freeArray<pq_scalar_t *>(this->thread_cut_line_weight_to_put_left);
		freeArray<pq_scalar_t>(this->process_rectelinear_cut_weight);
		freeArray<pq_scalar_t>(this->global_rectelinear_cut_weight);
	}

	freeArray<partId_t>(this->my_incomplete_cut_count);


	freeArray<pq_scalar_t>(this->max_min_coords);

	freeArray<pq_lno_t>(this->new_part_xadj);

	freeArray<pq_lno_t>(this->coordinate_permutations);

	freeArray<pq_lno_t>(this->new_coordinate_permutations);

	freeArray<pq_scalar_t>(this->all_cut_coordinates);



	freeArray<pq_scalar_t> (this->process_local_min_max_coord_total_weight);

	freeArray<pq_scalar_t> (this->global_min_max_coord_total_weight);





	freeArray<pq_scalar_t>(this->cut_coordinates_work_array);

	freeArray<pq_scalar_t>(this->target_part_weights);

	freeArray<pq_scalar_t>(this->cut_upper_bound_coordinates);

	freeArray<pq_scalar_t>(this->cut_lower_bound_coordinates);


	freeArray<pq_scalar_t>(this->cut_lower_bound_weights);
	freeArray<pq_scalar_t>(this->cut_upper_bound_weights);
	freeArray<bool>(this->is_cut_line_determined);
	freeArray<pq_scalar_t>(this->total_part_weight_left_right_closests);
	freeArray<pq_scalar_t>(this->global_total_part_weight_left_right_closests);

	for(int i = 0; i < this->num_threads; ++i){
		freeArray<double>(this->thread_part_weights[i]);
		freeArray<pq_scalar_t>(this->thread_cut_right_closest_point[i]);
		freeArray<pq_scalar_t>(this->thread_cut_left_closest_point[i]);
	}

	freeArray<double *>(this->thread_part_weights);
	freeArray<pq_scalar_t *>(this->thread_cut_left_closest_point);
	freeArray<pq_scalar_t *>(this->thread_cut_right_closest_point);

	this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Free");
}

/*! \brief Multi Jagged  coordinate partitioning algorithm.
 *
 *  \param distribute_points_on_cut_lines_ :  if partitioning can distribute points on same coordinate to different parts.
 *  \param max_concurrent_part_calculation_ : how many parts we can calculate concurrently.
 *  \param check_migrate_avoid_migration_option_ : whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
 *  \param minimum_migration_imbalance_  : when MJ decides whether to migrate, the minimum imbalance for migration.
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::set_partitioning_parameters(
		bool distribute_points_on_cut_lines_,
		int max_concurrent_part_calculation_,
		int check_migrate_avoid_migration_option_,
		pq_scalar_t minimum_migration_imbalance_){
	this->distribute_points_on_cut_lines = distribute_points_on_cut_lines_;
	this->max_concurrent_part_calculation = max_concurrent_part_calculation_;
	this->check_migrate_avoid_migration_option = check_migrate_avoid_migration_option_;
	this->minimum_migration_imbalance = minimum_migration_imbalance_;

}


/*! \brief Multi Jagged  coordinate partitioning algorithm.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm the communicator for the problem
 *  \param imbalance_tolerance : the input provided imbalance tolerance.
 *  \param num_global_parts: number of target global parts.
 *  \param part_no_array: part no array, if provided this will be used for partitioning.
 *  \param recursion_depth: if part no array is provided, it is the length of part no array,
 *  						if part no is not provided than it is the number of steps that algorithm will divide into num_global_parts parts.
 *
 *  \param coord_dim: coordinate dimension
 *  \param num_local_coords: number of local coordinates
 *  \param num_global_coords: number of global coordinates
 *  \param initial_mj_gnos: the list of initial global id's
 *  \param mj_coordinates: the two dimensional coordinate array.
 *
 *  \param weight_dim: weight dimension
 *  \param mj_uniform_weights: if weight dimension [i] has uniform weight or not.
 *  \param mj_weights: the two dimensional array for weights in every weight dimension
 *  \param mj_uniform_parts: if the target partitioning aims uniform parts
 *  \param mj_part_sizes: if the target partitioning does not aim uniform parts, then weight of each part.
 *
 *  \param result_assigned_part_ids: Output - 1D pointer, should be provided as null. Memory is given in the function.
 *  			the result partids corresponding to the coordinates given in result_mj_gnos.
 *  \param result_mj_gnos: Output - 1D pointer, should be provided as null. Memory is given in the function.
 *  			the result coordinate global id's corresponding to the part_ids array.
 *
 */
template <typename pq_scalar_t, typename pq_lno_t, typename pq_gno_t>
void AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t>::multi_jagged_part(

		const RCP<const Environment> &env,
    	RCP<Comm<int> > &problemComm,

    	double imbalance_tolerance_,
    	size_t num_global_parts_,
    	partId_t *part_no_array_,
    	int recursion_depth_,

    	int coord_dim_,
    	pq_lno_t num_local_coords_,
    	pq_gno_t num_global_coords_,
    	const pq_gno_t *initial_mj_gnos_,
    	pq_scalar_t **mj_coordinates_,

    	int weight_dim_,
    	bool *mj_uniform_weights_,
    	pq_scalar_t **mj_weights_,
    	bool *mj_uniform_parts_,
    	pq_scalar_t **mj_part_sizes_,

    	partId_t *&result_assigned_part_ids_,
    	pq_gno_t *&result_mj_gnos_

		){

#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL

    Z2_THROW_EXPERIMENTAL("Zoltan2 PQJagged is strictly experimental software "
            "while it is being developed and tested.")

#else

#ifdef print_debug
    if(comm->getRank() == 0){
    	cout << "size of gno:" << sizeof(pq_gno_t) << endl;
    	cout << "size of lno:" << sizeof(pq_lno_t) << endl;
    	cout << "size of pq_scalar_t:" << sizeof(pq_scalar_t) << endl;
    }
#endif

	cout << "multi jagged refactoring" << endl;


    this->mj_env = env;
    this->mj_problemComm = problemComm;
    this->myActualRank = this->myRank = this->mj_problemComm->getRank();
    //this->mj_coords = coords;
    //this->mj_solution = solution;

    this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Total");
    this->mj_env->debug(3, "In PQ Jagged");


    //this->set_input_parameters(this->mj_env->getParameters());
    {
    	this->imbalance_tolerance = imbalance_tolerance_;
    	this->num_global_parts = num_global_parts_;
    	this->part_no_array =  part_no_array_;
    	this->recursion_depth = recursion_depth_;


    	this->coord_dim = coord_dim_;
    	this->num_local_coords = num_local_coords_;
    	this->num_global_coords = num_global_coords_;
    	this->mj_coordinates = mj_coordinates_; //will copy the memory to this->mj_coordinates.
    	this->initial_mj_gnos = (pq_gno_t *) initial_mj_gnos_; //will copy the memory to this->current_mj_gnos[j].


    	this->weight_dim = weight_dim_;
    	this->mj_uniform_weights = mj_uniform_weights_;
    	this->mj_weights = mj_weights_; //will copy the memory to this->mj_weights
    	this->mj_uniform_parts = mj_uniform_parts_;
    	this->mj_part_sizes = mj_part_sizes_;

    	this->num_threads = 1;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel

    	{
    		this->num_threads = omp_get_num_threads();
    	}
#endif

    }
    //this->set_input_data();
    this->set_part_specifications();

    this->allocate_set_work_memory();



    //We duplicate the comm as we create subcommunicators during migration.
    //We keep the problemComm as it is, while comm changes after each migration.
    this->comm = this->mj_problemComm->duplicate();

    //initially there is a single partition
    partId_t current_num_parts = 1;
    pq_scalar_t *current_cut_coordinates =  this->all_cut_coordinates;


    this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Partitioning");

    partId_t output_part_begin_index = 0;
    partId_t future_num_parts = this->total_num_part;
    bool is_data_ever_migrated = false;

    vector<partId_t> *future_num_part_in_parts = new vector<partId_t> ();
    vector<partId_t> *next_future_num_parts_in_parts = new vector<partId_t> ();
    next_future_num_parts_in_parts->push_back(this->num_global_parts);

    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
     input_part_boxes(new vector <coordinateModelPartBox <pq_scalar_t, partId_t> >
      (), true) ;
    RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
     output_part_boxes(new vector <coordinateModelPartBox <pq_scalar_t, partId_t> >
      (), true);

    if(this->mj_keep_part_boxes){
    	this->init_part_boxes(output_part_boxes);
    }

    for (int i = 0; i < this->recursion_depth; ++i){
        //partitioning array. size will be as the number of current partitions and this
        //holds how many parts that each part will be in the current dimension partitioning.
        vector <partId_t> num_partitioning_in_current_dim;

        //number of parts that will be obtained at the end of this partitioning.
        //future_num_part_in_parts is as the size of current number of parts.
        //holds how many more parts each should be divided in the further
        //iterations. this will be used to calculate num_partitioning_in_current_dim,
        //as the number of parts that the part will be partitioned
        //in the current dimension partitioning.

        //next_future_num_parts_in_parts will be as the size of outnumParts,
        //and this will hold how many more parts that each output part
        //should be divided. this array will also be used to determine the weight ratios
        //of the parts.
        //swap the arrays to use iteratively..
        vector<partId_t> *tmpPartVect= future_num_part_in_parts;
        future_num_part_in_parts = next_future_num_parts_in_parts;
        next_future_num_parts_in_parts = tmpPartVect;

        //clear next_future_num_parts_in_parts array as
        //getPartitionArrays expects it to be empty.
        //it also expects num_partitioning_in_current_dim to be empty as well.
        next_future_num_parts_in_parts->clear();

        if(this->mj_keep_part_boxes){
            RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
                                                 tmpPartBoxes = input_part_boxes;
            input_part_boxes = output_part_boxes;
            output_part_boxes = tmpPartBoxes;
            output_part_boxes->clear();
        }

        //returns the total no. of output parts for this dimension partitioning.
        partId_t output_part_count_in_dimension =
        		this->update_part_num_arrays(
        				num_partitioning_in_current_dim,
        				future_num_part_in_parts,
        				next_future_num_parts_in_parts,
        				future_num_parts,
        				current_num_parts,
        				i,
        				input_part_boxes,
        				output_part_boxes);

        //if the number of obtained parts equal to current number of parts,
        //skip this dimension. For example, this happens when 1 is given in the input
        //part array is given. P=4,5,1,2
        if(output_part_count_in_dimension == current_num_parts) {
        	//still need to swap the input output arrays.
            tmpPartVect= future_num_part_in_parts;
            future_num_part_in_parts = next_future_num_parts_in_parts;
            next_future_num_parts_in_parts = tmpPartVect;

            if(this->mj_keep_part_boxes){
                RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > >
                                             tmpPartBoxes = input_part_boxes;
                input_part_boxes = output_part_boxes;
                output_part_boxes = tmpPartBoxes;
            }
            continue;
        }


        //get the coordinate axis along which the partitioning will be done.
        int coordInd = i % this->coord_dim;
        pq_scalar_t * mj_current_dim_coords = this->mj_coordinates[coordInd];

        //convert i to string to be used for debugging purposes.
        string istring = toString<int>(i);
        this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Partitioning_" + istring);

        //alloc Memory to point the indices
        //of the parts in the permutation array.
        this->new_part_xadj = allocMemory<pq_lno_t>(output_part_count_in_dimension);

        //the index where in the new_part_xadj will be written.
        partId_t output_part_index = 0;
        //whatever is written to output_part_index will be added with putput_coordinate_end_index
        //so that the points will be shifted.
        partId_t output_coordinate_end_index = 0;

        partId_t current_work_part = 0;
        partId_t current_concurrent_num_parts =
        		min(current_num_parts - current_work_part, this->max_concurrent_part_calculation);

        partId_t obtained_part_index = 0;

        //run for all available parts.
        for (; current_work_part < current_num_parts;
                 current_work_part += current_concurrent_num_parts){

            current_concurrent_num_parts = min(current_num_parts - current_work_part,
                                 this->max_concurrent_part_calculation);

            partId_t actual_work_part_count = 0;
            //initialization for 1D partitioning.
            //get the min and max coordinates of each part
            //together with the part weights of each part.
            for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                partId_t current_work_part_in_concurrent_parts = current_work_part + kk;

                //if this part wont be partitioned any further
                //dont do any work for this part.
                if (num_partitioning_in_current_dim[current_work_part_in_concurrent_parts] == 1){
                    continue;
                }
                ++actual_work_part_count;
                pq_lno_t coordinate_end_index= this->part_xadj[current_work_part_in_concurrent_parts];
                pq_lno_t coordinate_begin_index = current_work_part_in_concurrent_parts==0 ? 0: this->part_xadj[current_work_part_in_concurrent_parts -1];


                //cout << "i:" << i << " j:" << current_work_part + kk
                //		<< " coordinate_begin_index:" << coordinate_begin_index
                //		<< " coordinate_end_index:" << coordinate_end_index << endl;

                this->mj_get_local_min_max_coord_totW(
                    		coordinate_begin_index,
                    		coordinate_end_index,
                    		this->coordinate_permutations,
                    		mj_current_dim_coords,
                            this->process_local_min_max_coord_total_weight[kk], //min_coordinate
                            this->process_local_min_max_coord_total_weight[kk + current_concurrent_num_parts], //max_coordinate
                            this->process_local_min_max_coord_total_weight[kk + 2*current_concurrent_num_parts]); //total_weight

            }

            //1D partitioning
            if (actual_work_part_count > 0){
                //obtain global Min max of the part.
                this->mj_get_global_min_max_coord_totW(
                		current_concurrent_num_parts,
                		this->process_local_min_max_coord_total_weight,
                		this->global_min_max_coord_total_weight);

                //represents the total number of cutlines
                //whose coordinate should be determined.
                partId_t total_incomplete_cut_count = 0;

                //Compute weight ratios for parts & cuts:
                //e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
                //part0  cut0  part1 cut1 part2 cut2 part3
                partId_t concurrent_part_cut_shift = 0;
                partId_t concurrent_part_part_shift = 0;
                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    pq_scalar_t min_coordinate = this->global_min_max_coord_total_weight[kk];
                    pq_scalar_t max_coordinate = this->global_min_max_coord_total_weight[kk +
                                                     current_concurrent_num_parts];

                    pq_scalar_t global_total_weight = this->global_min_max_coord_total_weight[kk +
                                                        2 * current_concurrent_num_parts];

                    partId_t concurrent_current_part_index = current_work_part + kk;

                    partId_t partition_count = num_partitioning_in_current_dim[concurrent_current_part_index];

                    pq_scalar_t *usedCutCoordinate = current_cut_coordinates + concurrent_part_cut_shift;
                    pq_scalar_t *current_target_part_weights = this->target_part_weights +
                                                        concurrent_part_part_shift;
                    //shift the usedCutCoordinate array as noCuts.
                    concurrent_part_cut_shift += partition_count - 1;
                    //shift the partRatio array as noParts.
                    concurrent_part_part_shift += partition_count;


                    //calculate only if part is not empty,
                    //and part will be further partitioned.
                    if(partition_count > 1 && min_coordinate <= max_coordinate){

                        //increase num_cuts_do_be_determined by the number of cuts of the current
                        //part's cut line number.
                        total_incomplete_cut_count += partition_count - 1;
                        //set the number of cut lines that should be determined
                        //for this part.
                        this->my_incomplete_cut_count[kk] = partition_count - 1;

                        //get the target weights of the parts.
                        this->mj_get_initial_cut_coords_target_weights(
                        		min_coordinate,
                        		max_coordinate,
                        		partition_count - 1,
                        		global_total_weight,
                        		usedCutCoordinate,
                        		current_target_part_weights,
                        		future_num_part_in_parts,
                        		next_future_num_parts_in_parts,
                        		concurrent_current_part_index,
                        		obtained_part_index);

                        pq_lno_t coordinate_end_index= this->part_xadj[concurrent_current_part_index];
                        pq_lno_t coordinate_begin_index = concurrent_current_part_index==0 ? 0: this->part_xadj[concurrent_current_part_index -1];

                        //get the initial estimated part assignments of the
                        //coordinates.
                        this->set_initial_coordinate_parts(
                            max_coordinate,
                            min_coordinate,
                            concurrent_current_part_index,
                            coordinate_begin_index, coordinate_end_index,
                            this->coordinate_permutations,
                            mj_current_dim_coords,
                            this->assigned_part_ids,
                            partition_count);
                    }
                    else {
                        // e.g., if have fewer coordinates than parts, don't need to do next dim.
                        this->my_incomplete_cut_count[kk] = 0;
                    }
                    obtained_part_index += partition_count;
                }



                //used imbalance, it is always 0, as it is difficult to
                //estimate a range.
                pq_scalar_t used_imbalance = 0;


                // Determine cut lines for all concurrent parts parts here.
                this->mj_1D_part(
                    mj_current_dim_coords,
                    used_imbalance,
                    current_work_part,
                    current_concurrent_num_parts,
                    current_cut_coordinates,
                    total_incomplete_cut_count,
                    num_partitioning_in_current_dim);
            }

            //create new part chunks
            {
                partId_t output_array_shift = 0;
                partId_t cut_shift = 0;
                size_t tlr_shift = 0;
                size_t partweight_array_shift = 0;

                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    partId_t current_concurrent_work_part = current_work_part + kk;
                    partId_t num_parts = num_partitioning_in_current_dim[current_concurrent_work_part];

                    //if the part is empty, skip the part.
                    if((num_parts != 1  )
                    		&&
                    		this->global_min_max_coord_total_weight[kk] >
                             this->global_min_max_coord_total_weight[kk + current_concurrent_num_parts]) {

                    	//we still need to write the begin and end point of the
                    	//empty part. simply set it zero, the array indices will be shifted later.
                    	for(partId_t jj = 0; jj < num_parts; ++jj){
                    		this->new_part_xadj[output_part_index + output_array_shift + jj] = 0;
                    	}
                        cut_shift += num_parts - 1;
                        tlr_shift += (4 *(num_parts - 1) + 1);
                        output_array_shift += num_parts;
                        partweight_array_shift += (2 * (num_parts - 1) + 1);
                        continue;
                    }

                    pq_lno_t coordinate_end= this->part_xadj[current_concurrent_work_part];
                    pq_lno_t coordinate_begin = current_concurrent_work_part==0 ? 0: this->part_xadj[
                                                                current_concurrent_work_part -1];
                    pq_scalar_t *current_concurrent_cut_coordinate = current_cut_coordinates + cut_shift;
                    pq_scalar_t *used_local_cut_line_weight_to_left = this->process_cut_line_weight_to_put_left +
                                                            cut_shift;

                    //pq_scalar_t *used_tlr_array =  this->total_part_weight_left_right_closests + tlr_shift;

                    for(int ii = 0; ii < this->num_threads; ++ii){
                        this->thread_part_weight_work[ii] = this->thread_part_weights[ii] +  partweight_array_shift;
                    }

                    if(num_parts > 1){
                        if(this->mj_keep_part_boxes){
                        	//if part boxes are to be stored update the boundaries.
                            for (partId_t j = 0; j < num_parts - 1; ++j){
                                (*output_part_boxes)[output_array_shift + output_part_index +
                                 j].updateMinMax(current_concurrent_cut_coordinate[j], 1
                                  /*update max*/, coordInd);

                                (*output_part_boxes)[output_array_shift + output_part_index + j +
                                 1].updateMinMax(current_concurrent_cut_coordinate[j], 0
                                  /*update min*/, coordInd);
                            }
                        }

                        // Rewrite the indices based on the computed cuts.
                        this->mj_create_new_partitions(
                            num_parts,
                            mj_current_dim_coords,
                            current_concurrent_cut_coordinate,
                            coordinate_begin,
                            coordinate_end,
                            used_local_cut_line_weight_to_left,
                            this->thread_part_weight_work,
                            this->new_part_xadj + output_part_index + output_array_shift
                            );

                    }
                    else {
                        //if this part is partitioned into 1 then just copy
                        //the old values.
                        pq_lno_t part_size = coordinate_end - coordinate_begin;
                        *(this->new_part_xadj + output_part_index + output_array_shift) = part_size;
                        memcpy(
                        	this->new_coordinate_permutations + coordinate_begin,
                            this->coordinate_permutations + coordinate_begin,
                            part_size * sizeof(pq_lno_t));
                    }
                    cut_shift += num_parts - 1;
                    tlr_shift += (4 *(num_parts - 1) + 1);
                    output_array_shift += num_parts;
                    partweight_array_shift += (2 * (num_parts - 1) + 1);
                }

                //shift cut coordinates so that all cut coordinates are stored.
                //no shift now because we dont keep the cuts.
                //current_cut_coordinates += cut_shift;

                //mj_create_new_partitions from coordinates partitioned the parts and
                //write the indices as if there were a single part.
                //now we need to shift the beginning indices.
                for(partId_t kk = 0; kk < current_concurrent_num_parts; ++kk){
                    partId_t num_parts = num_partitioning_in_current_dim[ current_work_part + kk];
                    for (partId_t ii = 0;ii < num_parts ; ++ii){
                        //shift it by previousCount
                        this->new_part_xadj[output_part_index+ii] += output_coordinate_end_index;
                    }
                    //increase the previous count by current end.
                    output_coordinate_end_index = this->new_part_xadj[output_part_index + num_parts - 1];
                    //increase the current out.
                    output_part_index += num_parts ;
                }
            }
        }
        // end of this partitioning dimension


        int current_world_size = this->comm->getSize();
        long migration_reduce_all_population = this->total_dim_num_reduce_all * current_world_size;


        bool is_migrated_in_current_dimension = false;

        //we migrate if there are more partitionings to be done after this step
        //and if the migration is not forced to be avoided.
        //and the operation is not sequential.
        if (future_num_parts > 1 &&
            this->check_migrate_avoid_migration_option >= 0 &&
            current_world_size > 1){

        	this->mj_env->timerStart(MACRO_TIMERS, "PQJagged - Problem_Migration-" + istring);
        	partId_t num_parts = output_part_count_in_dimension;
        	if ( this->mj_perform_migration(
        					num_parts,
        					current_num_parts, //output
        					next_future_num_parts_in_parts, //output
        					output_part_begin_index,
        					migration_reduce_all_population,
        					this->num_local_coords / (future_num_parts * current_num_parts),
        					istring,
        					input_part_boxes, output_part_boxes) ) {
        		is_migrated_in_current_dimension = true;
        		is_data_ever_migrated = true;
        		this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Migration-" +
        				istring);
        		//since data is migrated, we reduce the number of reduceAll operations for the last part.
        		this->total_dim_num_reduce_all /= num_parts;
        	}
        	else {
        		is_migrated_in_current_dimension = false;
        		this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Migration-" + istring);
        	}
        }

        //swap the coordinate permutations for the next dimension.
        pq_lno_t * tmp = this->coordinate_permutations;
        this->coordinate_permutations = this->new_coordinate_permutations;
        this->new_coordinate_permutations = tmp;

        if(!is_migrated_in_current_dimension){
            this->total_dim_num_reduce_all -= current_num_parts;
            current_num_parts = output_part_count_in_dimension;
        }
        freeArray<pq_lno_t>(this->part_xadj);
        this->part_xadj = this->new_part_xadj;

        this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Partitioning_" + istring);
    }

    // Partitioning is done
    delete future_num_part_in_parts;
    delete next_future_num_parts_in_parts;

    this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Problem_Partitioning");
    /////////////////////////////End of the partitioning////////////////////////


    //get the final parts of each initial coordinate
    //the results will be written to
    //this->assigned_part_ids for gnos given in this->current_mj_gnos
    this->set_final_parts(
    		current_num_parts,
    		output_part_begin_index,
    		output_part_boxes,
    		is_data_ever_migrated);

	result_assigned_part_ids_ = this->assigned_part_ids;
	result_mj_gnos_ = this->current_mj_gnos;

    this->free_work_memory();
    this->mj_env->timerStop(MACRO_TIMERS, "PQJagged - Total");
    this->mj_env->debug(3, "Out of PQ Jagged");
#endif // INCLUDE_ZOLTAN2_EXPERIMENTAL

}


/*! \brief Multi Jagged coordinate partitioning algorithm.
 *
 */
template <typename Adapter>
class Zoltan2_AlgMJ{
private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
	typedef typename Adapter::scalar_t pq_scalar_t;
	typedef typename Adapter::gno_t pq_gno_t;
	typedef typename Adapter::lno_t pq_lno_t;
	typedef typename Adapter::node_t pq_node_t;
#endif
	AlgMJ<pq_scalar_t, pq_lno_t, pq_gno_t> mj_partitioner;

	RCP<const Environment> mj_env; //the environment object
	RCP<const CoordinateModel<typename Adapter::base_adapter_t> > mj_coords; //coordinate adapter
	RCP<PartitioningSolution<Adapter> > mj_solution; //solution object

	//PARAMETERS
	double imbalance_tolerance; //input imbalance tolerance.
	size_t num_global_parts; //the targeted number of parts
	partId_t *part_no_array; //input part array specifying num part to divide along each dim.
	int recursion_depth; //the number of steps that partitioning will be solved in.

	int coord_dim; // coordinate dimension.
    pq_lno_t num_local_coords; //number of local coords.
    pq_gno_t num_global_coords; //number of global coords.
    const pq_gno_t *initial_mj_gnos; //initial global ids of the coordinates.
    pq_scalar_t **mj_coordinates; //two dimension coordinate array

    int weight_dim; //weight dimension.
    bool *mj_uniform_weights; //if the coordinates have uniform weights.
    pq_scalar_t **mj_weights; //two dimension weight array
    bool *mj_uniform_parts; //if the target parts are uniform
    pq_scalar_t **mj_part_sizes; //target part weight sizes.

    bool distribute_points_on_cut_lines; //if partitioning can distribute points on same coordiante to different parts.
    int max_concurrent_part_calculation; // how many parts we can calculate concurrently.
    int check_migrate_avoid_migration_option; //whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
    pq_scalar_t minimum_migration_imbalance; //when MJ decides whether to migrate, the minimum imbalance for migration.
    int mj_keep_part_boxes; //if the boxes need to be kept.

    int num_threads;

    int mj_run_as_rcb; //if this is set, then recursion depth is adjusted to its maximum value.

	void set_input_parameters(const Teuchos::ParameterList &p);

	void free_work_memory();
public:
	Zoltan2_AlgMJ(){}
	~Zoltan2_AlgMJ(){}

	/*! \brief Multi Jagged  coordinate partitioning algorithm.
	 *
	 *  \param env   library configuration and problem parameters
	 *  \param problemComm the communicator for the problem
	 *  \param coords    a CoordinateModel with user data
	 *  \param solution  a PartitioningSolution, on input it
	 *      contains part information, on return it also contains
	 *      the solution and quality metrics.
	 */
	void multi_jagged_part(
			const RCP<const Environment> &env,
			RCP<Comm<int> > &problemComm,
			const RCP<const CoordinateModel<typename Adapter::base_adapter_t> > &mj_coords,
			RCP<PartitioningSolution<Adapter> > &solution);
};


/*! \brief Multi Jagged  coordinate partitioning algorithm.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm the communicator for the problem
 *  \param coords    a CoordinateModel with user data
 *  \param solution  a PartitioningSolution, on input it
 *      contains part information, on return it also contains
 *      the solution and quality metrics.
 */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::multi_jagged_part(
		const RCP<const Environment> &env,
		RCP<Comm<int> > &problemComm,
		const RCP<const CoordinateModel<typename Adapter::base_adapter_t> > &coords,
		RCP<PartitioningSolution<Adapter> > &solution){
    this->mj_env = env;
    this->mj_coords = coords;
    this->mj_solution = solution;
    this->set_input_parameters(this->mj_env->getParameters());
    if (this->mj_keep_part_boxes){
    	this->mj_partitioner.set_to_keep_part_boxes();
    }
    this->mj_partitioner.set_partitioning_parameters(
    		this->distribute_points_on_cut_lines,
    		this->max_concurrent_part_calculation,
    		this->check_migrate_avoid_migration_option,
    		this->minimum_migration_imbalance);

	partId_t *result_assigned_part_ids = NULL;
	pq_gno_t *result_mj_gnos = NULL;
    this->mj_partitioner.multi_jagged_part(
    		env,
    		problemComm,

    		this->imbalance_tolerance,
    		this->num_global_parts,
    		this->part_no_array,
    		this->recursion_depth,

    		this->coord_dim,
    		this->num_local_coords,
    		this->num_global_coords,
    		this->initial_mj_gnos,
    		this->mj_coordinates,

    		this->weight_dim,
    		this->mj_uniform_weights,
    		this->mj_weights,
    		this->mj_uniform_parts,
    		this->mj_part_sizes,

    		result_assigned_part_ids,
        	result_mj_gnos
    		);

    ArrayRCP<const pq_gno_t> gnoList = arcp(result_mj_gnos, 0, this->num_local_coords, true);
    ArrayRCP<partId_t> partId = arcp(result_assigned_part_ids, 0, this->num_local_coords, true);
    this->mj_solution->setParts(gnoList, partId, true);
    if (this->mj_keep_part_boxes){
    	RCP < vector <coordinateModelPartBox <pq_scalar_t, partId_t> > > output_part_boxes = this->mj_partitioner.get_part_boxes();
        this->mj_solution->setPartBoxes(output_part_boxes);
    }
    this->free_work_memory();
}
/* \brief Freeing the memory allocated.
 * */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::free_work_memory(){
	freeArray<pq_scalar_t *>(this->mj_coordinates);
	freeArray<pq_scalar_t *>(this->mj_weights);
	freeArray<bool>(this->mj_uniform_parts);
	freeArray<pq_scalar_t *>(this->mj_part_sizes);
	freeArray<bool>(this->mj_uniform_weights);

}

/* \brief Sets the partitioning data, and the partitioning parameters.
 * \param pl: is the parameter list provided to zoltan2 call
 * */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::set_input_parameters(const Teuchos::ParameterList &pl){

	this->coord_dim = this->mj_coords->getCoordinateDim();
	this->weight_dim = this->mj_coords->getCoordinateWeightDim();
	this->num_local_coords = this->mj_coords->getLocalNumCoordinates();
	this->num_global_coords = this->mj_coords->getGlobalNumCoordinates();
	int criteria_dim = (this->weight_dim ? this->weight_dim : 1);

	// From the Solution we get part information.
	// If the part sizes for a given criteria are not uniform,
	// then they are values that sum to 1.0.
	this->num_global_parts = this->mj_solution->getTargetGlobalNumberOfParts();
	//allocate only two dimensional pointer.
	//raw pointer addresess will be obtained from multivector.
	this->mj_coordinates = allocMemory<pq_scalar_t *>(this->coord_dim);
	this->mj_weights = allocMemory<pq_scalar_t *>(criteria_dim);

	//if the partitioning results are to be uniform.
	this->mj_uniform_parts = allocMemory< bool >(criteria_dim);
	//if in a criteria dimension, uniform part is false this shows ratios of
	//the target part weights.
	this->mj_part_sizes =  allocMemory<pq_scalar_t *>(criteria_dim);
	//if the weights of coordinates are uniform in a criteria dimension.
	this->mj_uniform_weights = allocMemory< bool >(criteria_dim);


	typedef StridedData<pq_lno_t, pq_scalar_t> input_t;
	ArrayView<const pq_gno_t> gnos;
	ArrayView<input_t> xyz;
	ArrayView<input_t> wgts;

	this->mj_coords->getCoordinates(gnos, xyz, wgts);
	//obtain global ids.
	ArrayView<const pq_gno_t> mj_gnos = gnos;
	this->initial_mj_gnos = mj_gnos.getRawPtr();

	//extract coordinates from multivector.
	for (int dim=0; dim < this->coord_dim; dim++){
		ArrayRCP<const pq_scalar_t> ar;
		xyz[dim].getInputArray(ar);
		//pqJagged coordinate values assignment
		this->mj_coordinates[dim] =  (pq_scalar_t *)ar.getRawPtr();
	}

	//if no weights are provided set uniform weight.
	if (this->weight_dim == 0){
		this->mj_uniform_weights[0] = true;
		this->mj_weights[0] = NULL;
	}
	else{
		//if weights are provided get weights for all weight dimensions.
		for (int wdim = 0; wdim < this->weight_dim; wdim++){
			if (wgts[wdim].size() == 0){
				//if the weight dimension has uniform weights.
				this->mj_uniform_weights[wdim] = true;
				this->mj_weights[wdim] = NULL;
			}
			else{
				//if the weight dimension doesn't have uniform weights, set weights of the objects.
				ArrayRCP<const pq_scalar_t> ar;
				wgts[wdim].getInputArray(ar);
				this->mj_uniform_weights[wdim] = false;
				this->mj_weights[wdim] = (pq_scalar_t *) ar.getRawPtr();
			}
		}
	}


	for (int wdim = 0; wdim < criteria_dim; wdim++){
		if (this->mj_solution->criteriaHasUniformPartSizes(wdim)){
			this->mj_uniform_parts[wdim] = true;
			this->mj_part_sizes[wdim] = NULL;
		}
		else{
			cerr << "MJ does not support non uniform target part weights" << endl;
			exit(1);
		}
	}



	const Teuchos::ParameterEntry *pe = pl.getEntryPtr("imbalance_tolerance");
	if (pe){
		double tol;
		tol = pe->getValue(&tol);
		this->imbalance_tolerance = tol - 1.0;
	}

	if (this->imbalance_tolerance <= 0)
		this->imbalance_tolerance= 10e-4;


	//if an input partitioning array is provided.
	this->part_no_array = NULL;
	//the length of the input partitioning array.
	this->recursion_depth = 0;

	if (pl.getPtr<Array <partId_t> >("pqParts")){
		this->part_no_array = (partId_t *) pl.getPtr<Array <partId_t> >("pqParts")->getRawPtr();
		this->recursion_depth = pl.getPtr<Array <partId_t> >("pqParts")->size() - 1;
		this->mj_env->debug(2, "PQparts provided by user");
	}

	//get mj specific parameters.
	this->distribute_points_on_cut_lines = true;
	this->max_concurrent_part_calculation = 1;

	this->mj_run_as_rcb = 0;
	int mj_user_recursion_depth = -1;
	this->mj_keep_part_boxes = 0;

	this->check_migrate_avoid_migration_option = 0;
	this->minimum_migration_imbalance = 0.35;

	// TODO: Append all the parameters with mj_
	pe = pl.getEntryPtr("migration_imbalance_cut_off");
	if (pe){
		double imb;
		imb = pe->getValue(&imb);
		this->minimum_migration_imbalance = imb - 1.0;
	}

	pe = pl.getEntryPtr("migration_check_option");
	if (pe){
		this->check_migrate_avoid_migration_option = pe->getValue(&this->check_migrate_avoid_migration_option);
	}else {
		this->check_migrate_avoid_migration_option = 0;
	}
	if (this->check_migrate_avoid_migration_option > 1) this->check_migrate_avoid_migration_option = -1;


	pe = pl.getEntryPtr("parallel_part_calculation_count");
	if (pe){
		this->max_concurrent_part_calculation = pe->getValue(&this->max_concurrent_part_calculation);
	}else {
		this->max_concurrent_part_calculation = 1; // Set to 1 if not provided.
	}

	pe = pl.getEntryPtr("keep_part_boxes");
	if (pe){
		this->mj_keep_part_boxes = pe->getValue(&this->mj_keep_part_boxes);
	}else {
		this->mj_keep_part_boxes = 0; // Set to invalid value
	}

	//need to keep part boxes if mapping type is geometric.
	if (this->mj_keep_part_boxes == 0){
		pe = pl.getEntryPtr("mapping_type");
		if (pe){
			int mapping_type = -1;
			mapping_type = pe->getValue(&mapping_type);
			if (mapping_type == 0){
				mj_keep_part_boxes  = 1;
			}
		}
	}

	//need to keep part boxes if mapping type is geometric.
	pe = pl.getEntryPtr("mj_enable_rcb");
	if (pe){
		this->mj_run_as_rcb = pe->getValue(&this->mj_run_as_rcb);
	}else {
		this->mj_run_as_rcb = 0; // Set to invalid value
	}

	pe = pl.getEntryPtr("recursion_depth");
	if (pe){
		mj_user_recursion_depth = pe->getValue(&mj_user_recursion_depth);
	}else {
		mj_user_recursion_depth = -1; // Set to invalid value
	}

	int val = 0;
	pe = pl.getEntryPtr("rectilinear_blocks");
	if (pe)
		val = pe->getValue(&val);

	if (val == 1){
		this->distribute_points_on_cut_lines = false;
	} else {
		this->distribute_points_on_cut_lines = true;
	}



	if (this->mj_run_as_rcb){
		mj_user_recursion_depth = (int)(ceil(log ((this->num_global_parts)) / log (2.0)));
	}
	if (this->recursion_depth < 1){
		if (mj_user_recursion_depth > 0){
			this->recursion_depth = mj_user_recursion_depth;
		}
		else {
			this->recursion_depth = this->coord_dim;
		}
	}

	this->num_threads = 1;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
	{
		this->num_threads = omp_get_num_threads();
	}

#endif

}


} // namespace Zoltan2

#endif
