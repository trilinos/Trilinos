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

/*! \file Zoltan2_AlgMultiJagged.hpp
  \brief Contains the Multi-jagged algorthm.
 */

#ifndef _ZOLTAN2_ALGMultiJagged_HPP_
#define _ZOLTAN2_ALGMultiJagged_HPP_

#include <Zoltan2_MultiJagged_ReductionOps.hpp>
#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_Algorithm.hpp>

#include <Tpetra_Distributor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Zoltan2_CoordinatePartitioningGraph.hpp>
#include <new>          // ::operator new[]
#include <algorithm>    // std::sort
#include <Zoltan2_Util.hpp>
#include <vector>
#include <unordered_map>

#ifdef HAVE_ZOLTAN2_ZOLTAN
#ifdef HAVE_ZOLTAN2_MPI
#define ENABLE_ZOLTAN_MIGRATION
#include "zoltan_comm_cpp.h"
#include "zoltan_types.h" // for error codes
#endif
#endif


#ifdef HAVE_ZOLTAN2_OMP
#include <omp.h>
#endif

#define LEAST_SIGNIFICANCE 0.0001
#define SIGNIFICANCE_MUL 1000

//if the (last dimension reduce all count) x the mpi world size
//estimated to be bigger than this number then migration will be forced
//in earlier iterations.
#define FUTURE_REDUCEALL_CUTOFF 1500000
//if parts right before last dimension are estimated to have less than
//MIN_WORK_LAST_DIM many coords, migration will be forced in earlier iterations.
#define MIN_WORK_LAST_DIM 1000




#define ABS(x) ((x) >= 0 ? (x) : -(x))
//imbalance calculation. Wreal / Wexpected - 1
#define imbalanceOf(Wachieved, totalW, expectedRatio) \
        (Wachieved) / ((totalW) * (expectedRatio)) - 1
#define imbalanceOf2(Wachieved, wExpected) \
        (Wachieved) / (wExpected) - 1



using std::vector;

namespace Teuchos{

/*! \brief Zoltan2_BoxBoundaries is a reduction operation
 * to all reduce the all box boundaries.
*/

template <typename Ordinal, typename T>
class Zoltan2_BoxBoundaries  : public ValueTypeReductionOp<Ordinal,T>
{
private:
    Ordinal size;
    T _EPSILON;

public:
    /*! \brief Default Constructor
     */
    Zoltan2_BoxBoundaries ():size(0), _EPSILON (std::numeric_limits<T>::epsilon()){}

    /*! \brief Constructor
     *   \param nsum  the count of how many sums will be computed at the
     *             start of the list.
     *   \param nmin  following the sums, this many minimums will be computed.
     *   \param nmax  following the minimums, this many maximums will be computed.
     */
    Zoltan2_BoxBoundaries (Ordinal s_):
        size(s_), _EPSILON (std::numeric_limits<T>::epsilon()){}

    /*! \brief Implement Teuchos::ValueTypeReductionOp interface
     */
    void reduce( const Ordinal count, const T inBuffer[], T inoutBuffer[]) const
    {
        for (Ordinal i=0; i < count; i++){
            if (Z2_ABS(inBuffer[i]) >  _EPSILON){
                inoutBuffer[i] = inBuffer[i];
            }
        }
    }
};
} // namespace Teuchos

namespace Zoltan2{

/*! \brief Allocates memory for the given size.
 *
 */
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

/*! \brief Frees the given array.
 *
 */
template <typename T>
void freeArray(T *&array){
    if(array != NULL){
        delete [] array;
        array = NULL;
    }
}

/*! \brief Converts the given object to string.
 *
 */
template <typename tt>
std::string toString(tt obj){
    std::stringstream ss (std::stringstream::in |std::stringstream::out);
    ss << obj;
    std::string tmp = "";
    ss >> tmp;
    return tmp;
}

/*! \brief Class for sorting items with multiple values.
 * First sorting with respect to val[0], then val[1] then ... val[count-1].
 * The last tie breaking is done with index values.
 * Used for task mapping partitioning where the points on a cut line needs to be
 * distributed consistently.
 *
 */
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
        this->_EPSILON = std::numeric_limits<WT>::epsilon() * 100;
    }


    uMultiSortItem(IT index_ ,CT count_, WT *vals_){
        this->index = index_;
        this->count = count_;
        this->val = vals_;
        this->_EPSILON = std::numeric_limits<WT>::epsilon() * 100;
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

/*! \brief Sort items for quick sort function.
 *
 */
template <class IT, class WT>
struct uSortItem
{
    IT id;
    //unsigned int val;
    WT val;
};// uSortItem;

/*! \brief Quick sort function.
 *	Sorts the arr of uSortItems, with respect to increasing vals.
 */
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
                std::cout << "uqsort: NSTACK too small in sort." << std::endl;
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



/*! \brief Multi Jagged coordinate partitioning algorithm.
 *
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
class AlgMJ
{
private:
    typedef coordinateModelPartBox<mj_scalar_t, mj_part_t> mj_partBox_t;
    typedef std::vector<mj_partBox_t> mj_partBoxVector_t;

    RCP<const Environment> mj_env; //the environment object
    RCP<Comm<int> > mj_problemComm; //initial comm object

    double imbalance_tolerance; //input imbalance tolerance.
    mj_part_t *part_no_array; //input part array specifying num part to divide along each dim.
    int recursion_depth; //the number of steps that partitioning will be solved in.
    int coord_dim, num_weights_per_coord; //coordinate dim and # of weights per coord

    size_t initial_num_loc_coords; //initial num local coords.
    global_size_t initial_num_glob_coords; //initial num global coords.

    mj_lno_t num_local_coords; //number of local coords.
    mj_gno_t num_global_coords; //number of global coords.

    mj_scalar_t **mj_coordinates; //two dimension coordinate array
    mj_scalar_t **mj_weights; //two dimension weight array
    bool *mj_uniform_parts; //if the target parts are uniform
    mj_scalar_t **mj_part_sizes; //target part weight sizes.
    bool *mj_uniform_weights; //if the coordinates have uniform weights.

    ArrayView<const mj_gno_t> mj_gnos; //global ids of the coordinates, comes from the input
    size_t num_global_parts; //the targeted number of parts

    mj_gno_t *initial_mj_gnos; //initial global ids of the coordinates.
    mj_gno_t *current_mj_gnos; //current global ids of the coordinates, might change during migration.
    int *owner_of_coordinate; //the actual processor owner of the coordinate, to track after migrations.

    mj_lno_t *coordinate_permutations; //permutation of coordinates, for partitioning.
    mj_lno_t *new_coordinate_permutations; //permutation work array.
    mj_part_t *assigned_part_ids; //the part ids assigned to coordinates.

    mj_lno_t *part_xadj; //beginning and end of each part.
    mj_lno_t *new_part_xadj; // work array for beginning and end of each part.

    //get mj specific parameters.
    bool distribute_points_on_cut_lines; //if partitioning can distribute points on same coordiante to different parts.
    mj_part_t max_concurrent_part_calculation; // how many parts we can calculate concurrently.

    int mj_run_as_rcb; //if this is set, then recursion depth is adjusted to its maximum value.
    int mj_user_recursion_depth; //the recursion depth value provided by user.
    int mj_keep_part_boxes; //if the boxes need to be kept.

    int check_migrate_avoid_migration_option; //whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
    mj_scalar_t minimum_migration_imbalance; //when MJ decides whether to migrate, the minimum imbalance for migration.
    int num_threads; //num threads

    mj_part_t total_num_cut ; //how many cuts will be totally
    mj_part_t total_num_part;    //how many parts will be totally

    mj_part_t max_num_part_along_dim ;         //maximum part count along a dimension.
    mj_part_t max_num_cut_along_dim; //maximum cut count along a dimension.
    size_t max_num_total_part_along_dim; //maximum part+cut count along a dimension.

    mj_part_t total_dim_num_reduce_all;    //estimate on #reduceAlls can be done.
    mj_part_t last_dim_num_part; //max no of parts that might occur
                                //during the partition before the
                                //last partitioning dimension.

    RCP<Comm<int> > comm; //comm object than can be altered during execution
    float fEpsilon; //epsilon for float
    mj_scalar_t sEpsilon; //epsilon for mj_scalar_t

    mj_scalar_t maxScalar_t; //max possible scalar
    mj_scalar_t minScalar_t; //min scalar

    mj_scalar_t *all_cut_coordinates;
    mj_scalar_t *max_min_coords;
    mj_scalar_t *process_cut_line_weight_to_put_left; //how much weight should a MPI put left side of the each cutline
    mj_scalar_t **thread_cut_line_weight_to_put_left; //how much weight percentage should each thread in MPI put left side of the each outline

    // work array to manipulate coordinate of cutlines in different iterations.
    //necessary because previous cut line information is used for determining
    //the next cutline information. therefore, cannot update the cut work array
    //until all cutlines are determined.
    mj_scalar_t *cut_coordinates_work_array;

    //cumulative part weight array.
    mj_scalar_t *target_part_weights;

    mj_scalar_t *cut_upper_bound_coordinates ;  //upper bound coordinate of a cut line
    mj_scalar_t *cut_lower_bound_coordinates ;  //lower bound coordinate of a cut line
    mj_scalar_t *cut_lower_bound_weights ;  //lower bound weight of a cut line
    mj_scalar_t *cut_upper_bound_weights ;  //upper bound weight of a cut line

    mj_scalar_t *process_local_min_max_coord_total_weight ; //combined array to exchange the min and max coordinate, and total weight of part.
    mj_scalar_t *global_min_max_coord_total_weight ;//global combined array with the results for min, max and total weight.

    //isDone is used to determine if a cutline is determined already.
    //If a cut line is already determined, the next iterations will skip this cut line.
    bool *is_cut_line_determined;
    //my_incomplete_cut_count count holds the number of cutlines that have not been finalized for each part
    //when concurrentPartCount>1, using this information, if my_incomplete_cut_count[x]==0, then no work is done for this part.
    mj_part_t *my_incomplete_cut_count;
    //local part weights of each thread.
    double **thread_part_weights;
    //the work manupulation array for partweights.
    double **thread_part_weight_work;

    //thread_cut_left_closest_point to hold the closest coordinate to a cutline from left (for each thread).
    mj_scalar_t **thread_cut_left_closest_point;
    //thread_cut_right_closest_point to hold the closest coordinate to a cutline from right (for each thread)
    mj_scalar_t **thread_cut_right_closest_point;

    //to store how many points in each part a thread has.
    mj_lno_t **thread_point_counts;

    mj_scalar_t *process_rectilinear_cut_weight;
    mj_scalar_t *global_rectilinear_cut_weight;

    //for faster communication, concatanation of
    //totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
    //leftClosest distances sized P-1, since P-1 cut lines
    //rightClosest distances size P-1, since P-1 cut lines.
    mj_scalar_t *total_part_weight_left_right_closests ;
    mj_scalar_t *global_total_part_weight_left_right_closests;

    RCP<mj_partBoxVector_t> kept_boxes;  // vector of all boxes for all parts;
                                         // constructed only if 
                                         // mj_keep_part_boxes == true
    RCP<mj_partBox_t> global_box;
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
    inline mj_part_t get_part_count(
    		mj_part_t num_total_future,
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
    void init_part_boxes(RCP<mj_partBoxVector_t> & outPartBoxes);

    /* \brief compute global bounding box:  min/max coords of global domain */
    void compute_global_box();

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
    mj_part_t update_part_num_arrays(
    		std::vector<mj_part_t> &num_partitioning_in_current_dim, //assumes this vector is empty.
    		std::vector<mj_part_t> *future_num_part_in_parts,
    		std::vector<mj_part_t> *next_future_num_parts_in_parts, //assumes this vector is empty.
    		mj_part_t &future_num_parts,
    		mj_part_t current_num_parts,
    		int current_iteration,
    		RCP<mj_partBoxVector_t> input_part_boxes,
    		RCP<mj_partBoxVector_t> output_part_boxes);

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
    		mj_lno_t coordinate_begin_index,
    		mj_lno_t coordinate_end_index,
    		mj_lno_t *mj_current_coordinate_permutations,
    		mj_scalar_t *mj_current_dim_coords,
    		mj_scalar_t &min_coordinate,
    		mj_scalar_t &max_coordinate,
    		mj_scalar_t &total_weight);

    /*! \brief Function that reduces global minimum and maximum coordinates with global total weight from given local arrays.
     * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
     * \param local_min_max_total is the array holding local min and max coordinate values with local total weight.
     * First current_concurrent_num_parts entries are minimums of the parts, next current_concurrent_num_parts entries are max, and then the total weights.
     * \param global_min_max_total is the output array holding global min and global coordinate values with global total weight.
     * The structure is same as local_min_max_total.
     */
    void mj_get_global_min_max_coord_totW(
        mj_part_t current_concurrent_num_parts,
        mj_scalar_t *local_min_max_total,
        mj_scalar_t *global_min_max_total);

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
        mj_scalar_t min_coord,
        mj_scalar_t max_coord,
        mj_part_t num_cuts/*p-1*/ ,
        mj_scalar_t global_weight,
        mj_scalar_t *initial_cut_coords /*p - 1 sized, coordinate of each cut line*/,
        mj_scalar_t *target_part_weights /*cumulative weights, at left side of each cut line. p-1 sized*/,

        std::vector <mj_part_t> *future_num_part_in_parts, //the vecto
        std::vector <mj_part_t> *next_future_num_parts_in_parts,
        mj_part_t concurrent_current_part,
        mj_part_t obtained_part_index);

    /*! \brief Function that calculates the new coordinates for the cut lines. Function is called inside the parallel region.
     * \param max_coordinate maximum coordinate in the range.
     * \param min_coordinate minimum coordinate in the range.
     *
     * \param concurrent_current_part_index is the index of the part in the inTotalCounts vector.
     * \param coordinate_begin_index holds the beginning of the coordinates in current part.
     * \param coordinate_end_index holds end of the coordinates in current part.
     * \param mj_current_coordinate_permutations is the permutation array, holds the real indices of coordinates on mj_current_dim_coords array.
     * \param mj_current_dim_coords is the 1D array holding the coordinates.
     * \param mj_part_ids is the array holding the partIds of each coordinate.
     * \param partition_count is the number of parts that the current part will be partitioned into.
     */
    void set_initial_coordinate_parts(
        mj_scalar_t &max_coordinate,
        mj_scalar_t &min_coordinate,
        mj_part_t &concurrent_current_part_index,
        mj_lno_t coordinate_begin_index,
        mj_lno_t coordinate_end_index,
        mj_lno_t *mj_current_coordinate_permutations,
        mj_scalar_t *mj_current_dim_coords,
        mj_part_t *mj_part_ids,
        mj_part_t &partition_count);

    /*! \brief Function that is responsible from 1D partitioning of the given range of coordinates.
     * \param mj_current_dim_coords is 1 dimensional array holding coordinate values.
     * \param imbalanceTolerance is the maximum allowed imbalance ratio.
     * \param current_work_part is the beginning index of concurrentPartCount parts.
     * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
     * \param current_cut_coordinates is the array holding the coordinates of the cut.
     * \param total_incomplete_cut_count is the number of cut lines whose positions should be calculated.
     * \param num_partitioning_in_current_dim is the vector that holds how many parts each part will be divided into.
     *
     */
    void mj_1D_part(
        mj_scalar_t *mj_current_dim_coords,
        mj_scalar_t imbalanceTolerance,
        mj_part_t current_work_part,
        mj_part_t current_concurrent_num_parts,
        mj_scalar_t *current_cut_coordinates,
        mj_part_t total_incomplete_cut_count,
        std::vector <mj_part_t> &num_partitioning_in_current_dim);

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
        mj_part_t num_cuts,
        mj_scalar_t max_coord,
        mj_scalar_t min_coord,
        mj_lno_t coordinate_begin_index,
        mj_lno_t coordinate_end_index,
        mj_scalar_t *mj_current_dim_coords,
        mj_scalar_t *temp_current_cut_coords,
        bool *current_cut_status,
        double *my_current_part_weights,
        mj_scalar_t *my_current_left_closest,
        mj_scalar_t *my_current_right_closest);

    /*! \brief Function that reduces the result of multiple threads
     * for left and right closest points and part weights in a single mpi process.
     *
     * \param num_partitioning_in_current_dim is the vector that holds the number of cut lines in current dimension for each part.
     * \param current_work_part holds the index of the first part (important when concurrent parts are used.)
     * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
     */
    void mj_accumulate_thread_results(
        const std::vector <mj_part_t> &num_partitioning_in_current_dim,
        mj_part_t current_work_part,
        mj_part_t current_concurrent_num_parts);

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
     * \param rectilinear_cut_count is the count of cut lines whose balance can be achived via distributing the points in same coordinate to different parts.
     * \param my_num_incomplete_cut is the number of cutlines whose position has not been determined yet. For K > 1 it is the count in a single part (whose cut lines are determined).
     */
    void mj_get_new_cut_coordinates(
        const size_t &num_total_part,
        const mj_part_t &num_cuts,
        const mj_scalar_t &max_coordinate,
        const mj_scalar_t &min_coordinate,
        const mj_scalar_t &global_total_weight,
        const mj_scalar_t &used_imbalance_tolerance,
        mj_scalar_t * current_global_part_weights,
        const mj_scalar_t * current_local_part_weights,
        const mj_scalar_t *current_part_target_weights,
        bool *current_cut_line_determined,
        mj_scalar_t *current_cut_coordinates,
        mj_scalar_t *current_cut_upper_bounds,
        mj_scalar_t *current_cut_lower_bounds,
        mj_scalar_t *current_global_left_closest_points,
        mj_scalar_t *current_global_right_closest_points,
        mj_scalar_t * current_cut_lower_bound_weights,
        mj_scalar_t * current_cut_upper_weights,
        mj_scalar_t *new_current_cut_coordinates,
        mj_scalar_t *current_part_cut_line_weight_to_put_left,
        mj_part_t *rectilinear_cut_count,
        mj_part_t &my_num_incomplete_cut);

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
    	mj_scalar_t cut_upper_bound,
        mj_scalar_t cut_lower_bound,
        mj_scalar_t cut_upper_weight,
        mj_scalar_t cut_lower_weight,
        mj_scalar_t expected_weight,
        mj_scalar_t &new_cut_position);

    /*! \brief Function that determines the permutation indices of the coordinates.
     * \param num_parts is the number of parts.
     * \param mj_current_dim_coords is 1 dimensional array holding the coordinate values.
     * \param current_concurrent_cut_coordinate is 1 dimensional array holding the cut coordinates.
     * \param coordinate_begin is the start index of the given partition on partitionedPointPermutations.
     * \param coordinate_end is the end index of the given partition on partitionedPointPermutations.
     * \param used_local_cut_line_weight_to_left holds how much weight of the coordinates on the cutline should be put on left side.
     * \param used_thread_part_weight_work is the two dimensional array holding the weight of parts for each thread. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
     * \param out_part_xadj is the indices of coordinates calculated for the partition on next dimension.
     */
    void mj_create_new_partitions(
        mj_part_t num_parts,
        mj_scalar_t *mj_current_dim_coords,
        mj_scalar_t *current_concurrent_cut_coordinate,
        mj_lno_t coordinate_begin,
        mj_lno_t coordinate_end,
        mj_scalar_t *used_local_cut_line_weight_to_left,
        double **used_thread_part_weight_work,
        mj_lno_t *out_part_xadj);

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
        mj_part_t in_num_parts, //current umb parts
        mj_part_t &out_num_parts, //output umb parts.
        std::vector<mj_part_t> *next_future_num_parts_in_parts,
        mj_part_t &output_part_begin_index,
        size_t migration_reduce_all_population,
        mj_lno_t num_coords_for_last_dim_part,
        std::string iteration,
        RCP<mj_partBoxVector_t> &input_part_boxes,
        RCP<mj_partBoxVector_t> &output_part_boxes);

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
    		mj_part_t num_procs,
    		mj_part_t num_parts,
    		mj_gno_t *&num_points_in_all_processor_parts);

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
    		mj_lno_t num_coords_for_last_dim_part,
    		mj_part_t num_procs,
    		mj_part_t num_parts,
    		mj_gno_t *num_points_in_all_processor_parts);


    /*! \brief Function fills up coordinate_destinations is the output array
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
     * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
     */
    void mj_migration_part_proc_assignment(
    		mj_gno_t * num_points_in_all_processor_parts,
    		mj_part_t num_parts,
    		mj_part_t num_procs,
    		mj_lno_t *send_count_to_each_proc,
    		std::vector<mj_part_t> &processor_ranks_for_subcomm,
    		std::vector<mj_part_t> *next_future_num_parts_in_parts,
    		mj_part_t &out_num_part,
    		std::vector<mj_part_t> &out_part_indices,
    		mj_part_t &output_part_numbering_begin_index,
    		int *coordinate_destinations);

    /*! \brief Function that assigned the processors to parts, when there are more processors then parts.
     *	sets the destination of each coordinate in coordinate_destinations, also edits output_part_numbering_begin_index,
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
     * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
     */
    void mj_assign_proc_to_parts(
    		mj_gno_t * num_points_in_all_processor_parts,
    		mj_part_t num_parts,
    		mj_part_t num_procs,
    		mj_lno_t *send_count_to_each_proc,
    		std::vector<mj_part_t> &processor_ranks_for_subcomm,
    		std::vector<mj_part_t> *next_future_num_parts_in_parts,
    		mj_part_t &out_part_index,
    		mj_part_t &output_part_numbering_begin_index,
    		int *coordinate_destinations);

    /*! \brief Function fills up coordinate_destinations is the output array
     * that holds which part each coordinate should be sent.
     *
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param num_procs is the number of processor attending to migration operation.
     * \param part_assignment_proc_begin_indices ([i]) points to the first processor index that part i will be sent to.
     * \param processor_chains_in_parts the array that holds the linked list structure, started from part_assignment_proc_begin_indices ([i]).
     * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
     * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
     */
    void assign_send_destinations(
    		mj_part_t num_parts,
    		mj_part_t *part_assignment_proc_begin_indices,
    		mj_part_t *processor_chains_in_parts,
    		mj_lno_t *send_count_to_each_proc,
    		int *coordinate_destinations);

    /*! \brief Function fills up coordinate_destinations is the output array
     * that holds which part each coordinate should be sent. In addition it calculates
     * the shift amount (output_part_numbering_begin_index) to be done when
     * final numberings of the parts are performed.
     *
     * \param num_parts is the number of parts that exist in the current partitioning.
     * \param sort_item_part_to_proc_assignment is the sorted parts with respect to the assigned processors.
     * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
     * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
     * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
     *
     */
    void assign_send_destinations2(
        mj_part_t num_parts,
        uSortItem<mj_part_t, mj_part_t> * sort_item_part_to_proc_assignment, //input sorted wrt processors
        int *coordinate_destinations,
        mj_part_t &output_part_numbering_begin_index,
        std::vector<mj_part_t> *next_future_num_parts_in_parts);

    /*! \brief Function fills up coordinate_destinations is the output array
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
     * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
     */
    void mj_assign_parts_to_procs(
        mj_gno_t * num_points_in_all_processor_parts,
        mj_part_t num_parts,
        mj_part_t num_procs,
        mj_lno_t *send_count_to_each_proc, //output: sized nprocs, show the number of send point counts to each proc.
        std::vector<mj_part_t> *next_future_num_parts_in_parts,//input how many more partitions the part will be partitioned into.
        mj_part_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
        std::vector<mj_part_t> &out_part_indices, //output: the part indices which the processor is assigned to.
        mj_part_t &output_part_numbering_begin_index, //output: how much the part number should be shifted when setting the solution
        int *coordinate_destinations);

    /*! \brief Function fills up coordinate_destinations is the output array
     * that holds which part each coordinate should be sent. In addition it calculates
     * the shift amount (output_part_numbering_begin_index) to be done when
     * final numberings of the parts are performed.
     *
     *
     * \param num_procs is the number of processor attending to migration operation.
     * \param num_new_local_points is the output to represent the new number of local points.
     * \param iteration is the string for the current iteration.
     * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
     * \param num_parts is the number of parts that exist in the current partitioning.
     */
    void mj_migrate_coords(
        mj_part_t num_procs,
        mj_lno_t &num_new_local_points,
        std::string iteration,
        int *coordinate_destinations,
        mj_part_t num_parts);

    /*! \brief Function creates the new subcomminicator for the processors
     * given in processor_ranks_for_subcomm.
     *
     * \param processor_ranks_for_subcomm is the vector that has the ranks of
     * the processors that will be in the same group.
     */
    void create_sub_communicator(vector<mj_part_t> &processor_ranks_for_subcomm);


    /*! \brief Function writes the new permutation arrays after the migration.
     *
     * \param output_num_parts is the number of parts that is assigned to the processor.
     * \param num_parts is the number of parts right before migration.
     */
    void fill_permutation_array(
        mj_part_t output_num_parts,
        mj_part_t num_parts);

    /*! \brief Function checks if should do migration or not.
     * \param current_num_parts is the number of parts in the process.
     * \param output_part_begin_index is the number that will be used as beginning part number
     * \param output_part_boxes is the array that holds the part boxes
     * \param is_data_ever_migrated is the boolean value which is true
     * if the data is ever migrated during the partitioning.
     *
     */
    void set_final_parts(
    		mj_part_t current_num_parts,
    		mj_part_t output_part_begin_index,
    		RCP<mj_partBoxVector_t> &output_part_boxes,
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
        mj_part_t num_parts,
        mj_scalar_t *mj_current_dim_coords,
        mj_scalar_t *current_concurrent_cut_coordinate,
        mj_lno_t coordinate_begin,
        mj_lno_t coordinate_end,
        mj_scalar_t *used_local_cut_line_weight_to_left,
        mj_lno_t *out_part_xadj,
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
     *  \param num_weights_per_coord: number of weights per coordinate
     *  \param mj_uniform_weights: if weight index [i] has uniform weight or not.
     *  \param mj_weights: the two dimensional array for weights
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
        	mj_part_t *part_no_array,
        	int recursion_depth,

        	int coord_dim,
        	mj_lno_t num_local_coords,
        	mj_gno_t num_global_coords,
        	const mj_gno_t *initial_mj_gnos,
        	mj_scalar_t **mj_coordinates,

        	int num_weights_per_coord,
        	bool *mj_uniform_weights,
        	mj_scalar_t **mj_weights,
        	bool *mj_uniform_parts,
        	mj_scalar_t **mj_part_sizes,

        	mj_part_t *&result_assigned_part_ids,
        	mj_gno_t *&result_mj_gnos

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
    		mj_scalar_t minimum_migration_imbalance_);
    /*! \brief Function call, if the part boxes are intended to be kept.
     *
     */
    void set_to_keep_part_boxes();

    /*! \brief Return the global bounding box: min/max coords of global domain
     */
    RCP<mj_partBox_t> get_global_box() const;

    RCP<mj_partBoxVector_t> get_kept_boxes() const;
    
    RCP<mj_partBoxVector_t> compute_global_box_boundaries(
        RCP<mj_partBoxVector_t> &localPartBoxes) const;

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
        mj_lno_t num_total_coords,
        mj_lno_t num_selected_coords,
        size_t num_target_part,
        int coord_dim,
        mj_scalar_t **mj_coordinates,
        mj_lno_t *initial_selected_coords_output_permutation,
        mj_lno_t *output_xadj,
        int recursion_depth,
        const mj_part_t *part_no_array);

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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::sequential_task_partitioning(
    const RCP<const Environment> &env,
    mj_lno_t num_total_coords,
    mj_lno_t num_selected_coords,
    size_t num_target_part,
    int coord_dim_,
    mj_scalar_t **mj_coordinates_,
    mj_lno_t *inital_adjList_output_adjlist,
    mj_lno_t *output_xadj,
    int rd,
    const mj_part_t *part_no_array_
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

    //weights are uniform for task mapping
 
    //parts are uniform for task mapping
    //as input indices.

    this->imbalance_tolerance = 0;
    this->num_global_parts = num_target_part;
    this->part_no_array = (mj_part_t *)part_no_array_;
    this->recursion_depth = rd;

    this->coord_dim = coord_dim_;
    this->num_local_coords = num_total_coords;
    this->num_global_coords = num_total_coords;
    this->mj_coordinates = mj_coordinates_;  //will copy the memory to this->mj_coordinates.

    ////temporary memory. It is not used here, but the functions require these to be allocated.
    ////will copy the memory to this->current_mj_gnos[j].
    this->initial_mj_gnos = allocMemory<mj_gno_t>(this->num_local_coords);

    this->num_weights_per_coord = 0;
    bool *tmp_mj_uniform_weights = new bool[1];
    this->mj_uniform_weights = tmp_mj_uniform_weights ;
    this->mj_uniform_weights[0] = true;

    mj_scalar_t **tmp_mj_weights = new mj_scalar_t *[1];
    this->mj_weights = tmp_mj_weights; //will copy the memory to this->mj_weights

    bool *tmp_mj_uniform_parts = new bool[1];
    this->mj_uniform_parts = tmp_mj_uniform_parts;
    this->mj_uniform_parts[0] = true;

    mj_scalar_t **tmp_mj_part_sizes = new mj_scalar_t * [1];
    this->mj_part_sizes = tmp_mj_part_sizes;
    this->mj_part_sizes[0] = NULL;

    this->num_threads = 1;
    this->set_part_specifications();

    this->allocate_set_work_memory();
    //the end of the initial partition is the end of coordinates.
    this->part_xadj[0] = static_cast<mj_lno_t>(num_selected_coords);
    for(size_t i = 0; i < static_cast<size_t>(num_total_coords); ++i){
        this->coordinate_permutations[i] = inital_adjList_output_adjlist[i];
    }

    mj_part_t current_num_parts = 1;

    mj_scalar_t *current_cut_coordinates =  this->all_cut_coordinates;

    mj_part_t future_num_parts = this->total_num_part;

    std::vector<mj_part_t> *future_num_part_in_parts = new std::vector<mj_part_t> ();
    std::vector<mj_part_t> *next_future_num_parts_in_parts = new std::vector<mj_part_t> ();
    next_future_num_parts_in_parts->push_back(this->num_global_parts);
    RCP<mj_partBoxVector_t> t1;
    RCP<mj_partBoxVector_t> t2;

    for (int i = 0; i < this->recursion_depth; ++i){

        //partitioning array. size will be as the number of current partitions and this
        //holds how many parts that each part will be in the current dimension partitioning.
        std::vector <mj_part_t> num_partitioning_in_current_dim;

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
        std::vector<mj_part_t> *tmpPartVect= future_num_part_in_parts;
        future_num_part_in_parts = next_future_num_parts_in_parts;
        next_future_num_parts_in_parts = tmpPartVect;

        //clear next_future_num_parts_in_parts array as
        //getPartitionArrays expects it to be empty.
        //it also expects num_partitioning_in_current_dim to be empty as well.
        next_future_num_parts_in_parts->clear();


        //returns the total number of output parts for this dimension partitioning.
        mj_part_t output_part_count_in_dimension =
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
        mj_scalar_t * mj_current_dim_coords = this->mj_coordinates[coordInd];
        //convert i to string to be used for debugging purposes.
        std::string istring = toString<int>(i);

        //alloc Memory to point the indices
        //of the parts in the permutation array.
        this->new_part_xadj = allocMemory<mj_lno_t>(output_part_count_in_dimension);

        //the index where in the outtotalCounts will be written.
        mj_part_t output_part_index = 0;
        //whatever is written to outTotalCounts will be added with previousEnd
        //so that the points will be shifted.
        mj_part_t output_coordinate_end_index = 0;

        mj_part_t current_work_part = 0;
        mj_part_t current_concurrent_num_parts = std::min(current_num_parts - current_work_part,
                                         this->max_concurrent_part_calculation);

        mj_part_t obtained_part_index = 0;

        //run for all available parts.
        for (; current_work_part < current_num_parts;
                     current_work_part += current_concurrent_num_parts){

            current_concurrent_num_parts = std::min(current_num_parts - current_work_part,
            this->max_concurrent_part_calculation);

            mj_part_t actual_work_part_count = 0;
            //initialization for 1D partitioning.
            //get the min and max coordinates of each part
            //together with the part weights of each part.
            for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                mj_part_t current_work_part_in_concurrent_parts = current_work_part + kk;

                //if this part wont be partitioned any further
                //dont do any work for this part.
                if (num_partitioning_in_current_dim[current_work_part_in_concurrent_parts] == 1){
                    continue;
                }
                ++actual_work_part_count;
                mj_lno_t coordinate_end_index= this->part_xadj[current_work_part_in_concurrent_parts];
                mj_lno_t coordinate_begin_index = current_work_part_in_concurrent_parts==0 ? 0: this->part_xadj[current_work_part_in_concurrent_parts -1];

                /*
                std::cout << "i:" << i << " j:" << current_work_part + kk
                		<< " coordinate_begin_index:" << coordinate_begin_index
                		<< " coordinate_end_index:" << coordinate_end_index
                		<< " total:" << coordinate_end_index - coordinate_begin_index<< std::endl;
                		*/
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
                mj_part_t total_incomplete_cut_count = 0;

                //Compute weight ratios for parts & cuts:
                //e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
                //part0  cut0  part1 cut1 part2 cut2 part3
                mj_part_t concurrent_part_cut_shift = 0;
                mj_part_t concurrent_part_part_shift = 0;
                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    mj_scalar_t min_coordinate = this->global_min_max_coord_total_weight[kk];
                    mj_scalar_t max_coordinate = this->global_min_max_coord_total_weight[kk +
                                                     current_concurrent_num_parts];
                    mj_scalar_t global_total_weight =
                    					this->global_min_max_coord_total_weight[kk +
                                                     2 * current_concurrent_num_parts];

                    mj_part_t concurrent_current_part_index = current_work_part + kk;

                    mj_part_t partition_count = num_partitioning_in_current_dim[concurrent_current_part_index];

                    mj_scalar_t *usedCutCoordinate = current_cut_coordinates + concurrent_part_cut_shift;
                    mj_scalar_t *current_target_part_weights = this->target_part_weights +
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

                        mj_lno_t coordinate_end_index= this->part_xadj[concurrent_current_part_index];
                        mj_lno_t coordinate_begin_index = concurrent_current_part_index==0 ? 0: this->part_xadj[concurrent_current_part_index -1];

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
                mj_scalar_t used_imbalance = 0;

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

                mj_part_t output_array_shift = 0;
                mj_part_t cut_shift = 0;
                size_t tlr_shift = 0;
                size_t partweight_array_shift = 0;

                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    mj_part_t current_concurrent_work_part = current_work_part + kk;
                    mj_part_t num_parts = num_partitioning_in_current_dim[current_concurrent_work_part];

                    //if the part is empty, skip the part.
                    if((num_parts != 1  ) && this->global_min_max_coord_total_weight[kk] >
                             this->global_min_max_coord_total_weight[kk + current_concurrent_num_parts]) {

                        for(mj_part_t jj = 0; jj < num_parts; ++jj){
                            this->new_part_xadj[output_part_index + output_array_shift + jj] = 0;
                        }
                        cut_shift += num_parts - 1;
                        tlr_shift += (4 *(num_parts - 1) + 1);
                        output_array_shift += num_parts;
                        partweight_array_shift += (2 * (num_parts - 1) + 1);
                        continue;
                    }

                    mj_lno_t coordinate_end = this->part_xadj[current_concurrent_work_part];
                    mj_lno_t coordinate_begin = current_concurrent_work_part==0 ? 0: this->part_xadj[current_concurrent_work_part
                                                             -1];
                    mj_scalar_t *current_concurrent_cut_coordinate = current_cut_coordinates + cut_shift;
                    mj_scalar_t *used_local_cut_line_weight_to_left = this->process_cut_line_weight_to_put_left +
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
                        mj_lno_t part_size = coordinate_end - coordinate_begin;
                        *(this->new_part_xadj + output_part_index + output_array_shift) = part_size;
                        memcpy(this->new_coordinate_permutations + coordinate_begin,
                        this->coordinate_permutations + coordinate_begin,
                        part_size * sizeof(mj_lno_t));
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
                for(mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk){
                    mj_part_t num_parts = num_partitioning_in_current_dim[ current_work_part + kk];
                    for (mj_part_t ii = 0;ii < num_parts ; ++ii){
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
        mj_lno_t * tmp = this->coordinate_permutations;
        this->coordinate_permutations = this->new_coordinate_permutations;
        this->new_coordinate_permutations = tmp;

        freeArray<mj_lno_t>(this->part_xadj);
        this->part_xadj = this->new_part_xadj;
    }

    for(mj_lno_t i = 0; i < num_total_coords; ++i){
    	inital_adjList_output_adjlist[i] = this->coordinate_permutations[i];
    }

    for(size_t i = 0; i < this->num_global_parts ; ++i){
        output_xadj[i] = this->part_xadj[i];
    }

    delete future_num_part_in_parts;
    delete next_future_num_parts_in_parts;

    //free the extra memory that we allocated.
    freeArray<mj_part_t>(this->assigned_part_ids);
    freeArray<mj_gno_t>(this->initial_mj_gnos);
    freeArray<mj_gno_t>(this->current_mj_gnos);
    freeArray<bool>(tmp_mj_uniform_weights);
    freeArray<bool>(tmp_mj_uniform_parts);
    freeArray<mj_scalar_t *>(tmp_mj_weights);
    freeArray<mj_scalar_t *>(tmp_mj_part_sizes);

    this->free_work_memory();

#ifdef HAVE_ZOLTAN2_OMP
    omp_set_num_threads(actual_num_threads);
#endif
}

/*! \brief Multi Jagged  coordinate partitioning algorithm default constructor.
 *
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::AlgMJ():
	mj_env(), mj_problemComm(), imbalance_tolerance(0),
	part_no_array(NULL), recursion_depth(0), coord_dim(0),
	num_weights_per_coord(0), initial_num_loc_coords(0),
        initial_num_glob_coords(0),
	num_local_coords(0), num_global_coords(0), mj_coordinates(NULL),
	mj_weights(NULL), mj_uniform_parts(NULL), mj_part_sizes(NULL),
	mj_uniform_weights(NULL), mj_gnos(), num_global_parts(1),
	initial_mj_gnos(NULL), current_mj_gnos(NULL), owner_of_coordinate(NULL),
	coordinate_permutations(NULL), new_coordinate_permutations(NULL),
	assigned_part_ids(NULL), part_xadj(NULL), new_part_xadj(NULL),
	distribute_points_on_cut_lines(true), max_concurrent_part_calculation(1),
	mj_run_as_rcb(0), mj_user_recursion_depth(0), mj_keep_part_boxes(0),
	check_migrate_avoid_migration_option(0), minimum_migration_imbalance(0.30),
	num_threads(1), total_num_cut(0), total_num_part(0), max_num_part_along_dim(0),
	max_num_cut_along_dim(0), max_num_total_part_along_dim(0), total_dim_num_reduce_all(0),
	last_dim_num_part(0), comm(), fEpsilon(0), sEpsilon(0), maxScalar_t(0), minScalar_t(0),
	all_cut_coordinates(NULL), max_min_coords(NULL), process_cut_line_weight_to_put_left(NULL),
	thread_cut_line_weight_to_put_left(NULL), cut_coordinates_work_array(NULL),
	target_part_weights(NULL), cut_upper_bound_coordinates(NULL), cut_lower_bound_coordinates(NULL),
	cut_lower_bound_weights(NULL), cut_upper_bound_weights(NULL),
	process_local_min_max_coord_total_weight(NULL), global_min_max_coord_total_weight(NULL),
	is_cut_line_determined(NULL), my_incomplete_cut_count(NULL),
	thread_part_weights(NULL), thread_part_weight_work(NULL),
	thread_cut_left_closest_point(NULL), thread_cut_right_closest_point(NULL),
	thread_point_counts(NULL), process_rectilinear_cut_weight(NULL),
	global_rectilinear_cut_weight(NULL),total_part_weight_left_right_closests(NULL),
	global_total_part_weight_left_right_closests(NULL),
        kept_boxes(),global_box(),
	myRank(0), myActualRank(0)
{
    this->fEpsilon = std::numeric_limits<float>::epsilon();
    this->sEpsilon = std::numeric_limits<mj_scalar_t>::epsilon() * 100;

    this->maxScalar_t = std::numeric_limits<mj_scalar_t>::max();
    this->minScalar_t = -std::numeric_limits<mj_scalar_t>::max();

}


/*! \brief Function returns the part boxes stored
 * returns null if boxes are not stored, and prints warning mesage.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
RCP<typename AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t>::mj_partBox_t>
AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t>::get_global_box() const 
{
  return this->global_box;
}

/*! \brief Function call, if the part boxes are intended to be kept.
 *
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::set_to_keep_part_boxes(){
  this->mj_keep_part_boxes = 1;
}


/* \brief Either the mj array (part_no_array) or num_global_parts should be provided in
 * the input. part_no_array takes
 * precedence if both are provided.
 * Depending on these parameters, total cut/part number,
 * maximum part/cut number along a dimension, estimated number of reduceAlls,
 * and the number of parts before the last dimension is calculated.
 * */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::set_part_specifications(){

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
		mj_part_t future_num_parts = this->num_global_parts;

		//we need to calculate the part numbers now, to determine the maximum along the dimensions.
		for (int i = 0; i < this->recursion_depth; ++i){

			mj_part_t maxNoPartAlongI = this->get_part_count(
					future_num_parts, 1.0f / (this->recursion_depth - i));

			if (maxNoPartAlongI > this->max_num_part_along_dim){
				this->max_num_part_along_dim = maxNoPartAlongI;
			}

			mj_part_t nfutureNumParts = future_num_parts / maxNoPartAlongI;
			if (future_num_parts % maxNoPartAlongI){
				++nfutureNumParts;
			}
			future_num_parts = nfutureNumParts;
		}
		this->total_num_part = this->num_global_parts;
		//estimate reduceAll Count here.
		//we find the upperbound instead.
		mj_part_t p = 1;
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
            std::cerr << "Warning: Concurrent part count ("<< this->max_concurrent_part_calculation <<
            ") has been set bigger than maximum amount that can be used." <<
            " Setting to:" << this->last_dim_num_part << "." << std::endl;
        }
        this->max_concurrent_part_calculation = this->last_dim_num_part;
    }

}
/* \brief Tries to determine the part number for current dimension,
 * by trying to make the partitioning as square as possible.
 * \param num_total_future how many more partitionings are required.
 * \param root how many more recursion depth is left.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
inline mj_part_t AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::get_part_count(
		mj_part_t num_total_future,
		double root)
{
	double fp = pow(num_total_future, root);
	mj_part_t ip = mj_part_t (fp);
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
mj_part_t AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::update_part_num_arrays(
	std::vector <mj_part_t> &num_partitioning_in_current_dim, //assumes this vector is empty.
    std::vector<mj_part_t> *future_num_part_in_parts,
    std::vector<mj_part_t> *next_future_num_parts_in_parts, //assumes this vector is empty.
    mj_part_t &future_num_parts,
    mj_part_t current_num_parts,
    int current_iteration,
    RCP<mj_partBoxVector_t> input_part_boxes,
    RCP<mj_partBoxVector_t> output_part_boxes
){
	//how many parts that will be obtained after this dimension.
    mj_part_t output_num_parts = 0;
    if(this->part_no_array){
        //when the partNo array is provided as input,
        //each current partition will be partition to the same number of parts.
        //we dont need to use the future_num_part_in_parts vector in this case.

        mj_part_t p = this->part_no_array[current_iteration];
        if (p < 1){
            std::cout << "i:" << current_iteration << " p is given as:" << p << std::endl;
            exit(1);
        }
        if (p == 1){
            return current_num_parts;
        }

        for (mj_part_t ii = 0; ii < current_num_parts; ++ii){
            num_partitioning_in_current_dim.push_back(p);

        }
        //cout << "me:" << this->myRank << " current_iteration" << current_iteration <<
        //" current_num_parts:" << current_num_parts << std::endl;
        //cout << "num_partitioning_in_current_dim[0]:" << num_partitioning_in_current_dim[0] << std::endl;
        //set the new value of future_num_parts.

        /*
        cout << "\tfuture_num_parts:" << future_num_parts
        		<< " num_partitioning_in_current_dim[0]:" << num_partitioning_in_current_dim[0]
        		<< future_num_parts/ num_partitioning_in_current_dim[0] << std::endl;
        */

        future_num_parts /= num_partitioning_in_current_dim[0];
        output_num_parts = current_num_parts * num_partitioning_in_current_dim[0];

        if (this->mj_keep_part_boxes){
            for (mj_part_t k = 0; k < current_num_parts; ++k){
            	//initialized the output boxes as its ancestor.
            	for (mj_part_t j = 0; j < num_partitioning_in_current_dim[0]; ++j){
                    output_part_boxes->push_back((*input_part_boxes)[k]);
                }
            }
        }

        //set the how many more parts each part will be divided.
        //this is obvious when partNo array is provided as input.
        //however, fill this so that weights will be calculated according to this array.
        for (mj_part_t ii = 0; ii < output_num_parts; ++ii){
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

        //cout << "i:" << i << std::endl;

        for (mj_part_t ii = 0; ii < current_num_parts; ++ii){
            //get how many parts a part should be divided.
            mj_part_t future_num_parts_of_part_ii = (*future_num_part_in_parts)[ii];

            //get the ideal number of parts that is close to the
            //(recursion_depth - i) root of the future_num_parts_of_part_ii.
            mj_part_t num_partitions_in_current_dim =
            					this->get_part_count(
            							future_num_parts_of_part_ii,
            							1.0 / (this->recursion_depth - current_iteration)
                                	);

            if (num_partitions_in_current_dim > this->max_num_part_along_dim){
                std::cerr << "ERROR: maxPartNo calculation is wrong." << std::endl;
                exit(1);
            }
            //add this number to num_partitioning_in_current_dim vector.
            num_partitioning_in_current_dim.push_back(num_partitions_in_current_dim);


            //increase the output number of parts.
            output_num_parts += num_partitions_in_current_dim;

            //ideal number of future partitions for each part.
            mj_part_t ideal_num_future_parts_in_part = future_num_parts_of_part_ii / num_partitions_in_current_dim;
            for (mj_part_t iii = 0; iii < num_partitions_in_current_dim; ++iii){
                mj_part_t num_future_parts_for_part_iii = ideal_num_future_parts_in_part;

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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::allocate_set_work_memory(){

	//points to process that initially owns the coordinate.
	this->owner_of_coordinate  = NULL;

	//Throughout the partitioning execution,
	//instead of the moving the coordinates, hold a permutation array for parts.
	//coordinate_permutations holds the current permutation.
	this->coordinate_permutations =  allocMemory< mj_lno_t>(this->num_local_coords);
	//initial configuration, set each pointer-i to i.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
	for(mj_lno_t i = 0; i < this->num_local_coords; ++i){
		this->coordinate_permutations[i] = i;
	}

	//new_coordinate_permutations holds the current permutation.
	this->new_coordinate_permutations = allocMemory< mj_lno_t>(this->num_local_coords);

	this->assigned_part_ids = NULL;
	if(this->num_local_coords > 0){
		this->assigned_part_ids = allocMemory<mj_part_t>(this->num_local_coords);
	}

	//single partition starts at index-0, and ends at numLocalCoords
	//inTotalCounts array holds the end points in coordinate_permutations array
	//for each partition. Initially sized 1, and single element is set to numLocalCoords.
	this->part_xadj = allocMemory<mj_lno_t>(1);
	this->part_xadj[0] = static_cast<mj_lno_t>(this->num_local_coords);//the end of the initial partition is the end of coordinates.
	//the ends points of the output, this is allocated later.
	this->new_part_xadj = NULL;

	// only store this much if cuts are needed to be stored.
	//this->all_cut_coordinates = allocMemory< mj_scalar_t>(this->total_num_cut);


	this->all_cut_coordinates  = allocMemory< mj_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

	this->max_min_coords =  allocMemory< mj_scalar_t>(this->num_threads * 2);

	this->process_cut_line_weight_to_put_left = NULL; //how much weight percentage should a MPI put left side of the each cutline
	this->thread_cut_line_weight_to_put_left = NULL; //how much weight percentage should each thread in MPI put left side of the each outline
	//distribute_points_on_cut_lines = false;
	if(this->distribute_points_on_cut_lines){
		this->process_cut_line_weight_to_put_left = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
		this->thread_cut_line_weight_to_put_left = allocMemory<mj_scalar_t *>(this->num_threads);
		for(int i = 0; i < this->num_threads; ++i){
			this->thread_cut_line_weight_to_put_left[i] = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim);
		}
	    this->process_rectilinear_cut_weight = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim);
	    this->global_rectilinear_cut_weight = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim);
	}


	// work array to manipulate coordinate of cutlines in different iterations.
	//necessary because previous cut line information is used for determining
	//the next cutline information. therefore, cannot update the cut work array
	//until all cutlines are determined.
	this->cut_coordinates_work_array = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim *
			this->max_concurrent_part_calculation);


	//cumulative part weight array.
	this->target_part_weights = allocMemory<mj_scalar_t>(
					this->max_num_part_along_dim * this->max_concurrent_part_calculation);
	// the weight from left to write.

    this->cut_upper_bound_coordinates = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);  //upper bound coordinate of a cut line
    this->cut_lower_bound_coordinates = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim* this->max_concurrent_part_calculation);  //lower bound coordinate of a cut line
    this->cut_lower_bound_weights = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim* this->max_concurrent_part_calculation);  //lower bound weight of a cut line
    this->cut_upper_bound_weights = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim* this->max_concurrent_part_calculation);  //upper bound weight of a cut line

    this->process_local_min_max_coord_total_weight = allocMemory<mj_scalar_t>(3 * this->max_concurrent_part_calculation); //combined array to exchange the min and max coordinate, and total weight of part.
    this->global_min_max_coord_total_weight = allocMemory<mj_scalar_t>(3 * this->max_concurrent_part_calculation);//global combined array with the results for min, max and total weight.

    //is_cut_line_determined is used to determine if a cutline is determined already.
    //If a cut line is already determined, the next iterations will skip this cut line.
    this->is_cut_line_determined = allocMemory<bool>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
    //my_incomplete_cut_count count holds the number of cutlines that have not been finalized for each part
    //when concurrentPartCount>1, using this information, if my_incomplete_cut_count[x]==0, then no work is done for this part.
    this->my_incomplete_cut_count =  allocMemory<mj_part_t>(this->max_concurrent_part_calculation);
    //local part weights of each thread.
    this->thread_part_weights = allocMemory<double *>(this->num_threads);
    //the work manupulation array for partweights.
    this->thread_part_weight_work = allocMemory<double *>(this->num_threads);

    //thread_cut_left_closest_point to hold the closest coordinate to a cutline from left (for each thread).
    this->thread_cut_left_closest_point = allocMemory<mj_scalar_t *>(this->num_threads);
    //thread_cut_right_closest_point to hold the closest coordinate to a cutline from right (for each thread)
    this->thread_cut_right_closest_point = allocMemory<mj_scalar_t *>(this->num_threads);

    //to store how many points in each part a thread has.
    this->thread_point_counts = allocMemory<mj_lno_t *>(this->num_threads);

    for(int i = 0; i < this->num_threads; ++i){
        //partWeights[i] = allocMemory<mj_scalar_t>(maxTotalPartCount);
        this->thread_part_weights[i] = allocMemory < double >(this->max_num_total_part_along_dim * this->max_concurrent_part_calculation);
        this->thread_cut_right_closest_point[i] = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
        this->thread_cut_left_closest_point[i] = allocMemory<mj_scalar_t>(this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
        this->thread_point_counts[i] =  allocMemory<mj_lno_t>(this->max_num_part_along_dim);
    }
    //for faster communication, concatanation of
    //totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
    //leftClosest distances sized P-1, since P-1 cut lines
    //rightClosest distances size P-1, since P-1 cut lines.
    this->total_part_weight_left_right_closests = allocMemory<mj_scalar_t>((this->max_num_total_part_along_dim + this->max_num_cut_along_dim * 2) * this->max_concurrent_part_calculation);
    this->global_total_part_weight_left_right_closests = allocMemory<mj_scalar_t>((this->max_num_total_part_along_dim + this->max_num_cut_along_dim * 2) * this->max_concurrent_part_calculation);


    mj_scalar_t **coord = allocMemory<mj_scalar_t *>(this->coord_dim);
    for (int i=0; i < this->coord_dim; i++){
    	coord[i] = allocMemory<mj_scalar_t>(this->num_local_coords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    	for (mj_lno_t j=0; j < this->num_local_coords; j++)
    		coord[i][j] = this->mj_coordinates[i][j];
    }
    this->mj_coordinates = coord;


    int criteria_dim = (this->num_weights_per_coord ? this->num_weights_per_coord : 1);
    mj_scalar_t **weights = allocMemory<mj_scalar_t *>(criteria_dim);

    for (int i=0; i < criteria_dim; i++){
    	weights[i] = NULL;
    }
    for (int i=0; i < this->num_weights_per_coord; i++){
    	weights[i] = allocMemory<mj_scalar_t>(this->num_local_coords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    	for (mj_lno_t j=0; j < this->num_local_coords; j++)
    		weights[i][j] = this->mj_weights[i][j];

    }
	this->mj_weights = weights;
    this->current_mj_gnos = allocMemory<mj_gno_t>(this->num_local_coords);
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for (mj_lno_t j=0; j < this->num_local_coords; j++)
    	this->current_mj_gnos[j] = this->initial_mj_gnos[j];

    this->owner_of_coordinate = allocMemory<int>(this->num_local_coords);

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for (mj_lno_t j=0; j < this->num_local_coords; j++)
    	this->owner_of_coordinate[j] = this->myActualRank;
}

/* \brief compute the global bounding box
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t>::compute_global_box()
{
    //local min coords
    mj_scalar_t *mins = allocMemory<mj_scalar_t>(this->coord_dim);
    //global min coords
    mj_scalar_t *gmins = allocMemory<mj_scalar_t>(this->coord_dim);
    //local max coords
    mj_scalar_t *maxs = allocMemory<mj_scalar_t>(this->coord_dim);
    //global max coords
    mj_scalar_t *gmaxs = allocMemory<mj_scalar_t>(this->coord_dim);

    for (int i = 0; i < this->coord_dim; ++i){
        mj_scalar_t localMin = this->mj_coordinates[i][0];
        mj_scalar_t localMax = this->mj_coordinates[i][0];
        for (mj_lno_t j = 1; j < this->num_local_coords; ++j){
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
    reduceAll<int, mj_scalar_t>(*this->comm, Teuchos::REDUCE_MIN,
            this->coord_dim, mins, gmins
    );

    reduceAll<int, mj_scalar_t>(*this->comm, Teuchos::REDUCE_MAX,
            this->coord_dim, maxs, gmaxs
    );

    //create single box with all areas.
    global_box = rcp(new mj_partBox_t(0,this->coord_dim,gmins,gmaxs));
    //coordinateModelPartBox <mj_scalar_t, mj_part_t> tmpBox (0, coordDim);
    freeArray<mj_scalar_t>(mins);
    freeArray<mj_scalar_t>(gmins);
    freeArray<mj_scalar_t>(maxs);
    freeArray<mj_scalar_t>(gmaxs);
}

/* \brief for part communication we keep track of the box boundaries.
 * This is performed when either asked specifically, or when geometric mapping is performed afterwards.
 * This function initializes a single box with all global min and max coordinates.
 * \param initial_partitioning_boxes the input and output vector for boxes.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::init_part_boxes(
		RCP<mj_partBoxVector_t> & initial_partitioning_boxes
)
{
    mj_partBox_t tmp_box(*global_box);
    initial_partitioning_boxes->push_back(tmp_box);
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_get_local_min_max_coord_totW(
		mj_lno_t coordinate_begin_index,
		mj_lno_t coordinate_end_index,
		mj_lno_t *mj_current_coordinate_permutations,
		mj_scalar_t *mj_current_dim_coords,
		mj_scalar_t &min_coordinate,
		mj_scalar_t &max_coordinate,
		mj_scalar_t &total_weight){

    //if the part is empty.
    //set the min and max coordinates as reverse.
    if(coordinate_begin_index >= coordinate_end_index)
    {
        min_coordinate = this->maxScalar_t;
        max_coordinate = this->minScalar_t;
        total_weight = 0;
    }
    else {
        mj_scalar_t my_total_weight = 0;
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
                for (mj_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){
                    int i = mj_current_coordinate_permutations[ii];
                    my_total_weight += this->mj_weights[0][i];
                }
            }

            int my_thread_id = 0;
#ifdef HAVE_ZOLTAN2_OMP
            my_thread_id = omp_get_thread_num();
#endif
            mj_scalar_t my_thread_min_coord, my_thread_max_coord;
            my_thread_min_coord=my_thread_max_coord
                =mj_current_dim_coords[mj_current_coordinate_permutations[coordinate_begin_index]];


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
            for(mj_lno_t j = coordinate_begin_index + 1; j < coordinate_end_index; ++j){
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_get_global_min_max_coord_totW(
    mj_part_t current_concurrent_num_parts,
    mj_scalar_t *local_min_max_total,
    mj_scalar_t *global_min_max_total){

	//reduce min for first current_concurrent_num_parts elements, reduce max for next
	//concurrentPartCount elements,
	//reduce sum for the last concurrentPartCount elements.
	if(this->comm->getSize()  > 1){
		Teuchos::MultiJaggedCombinedMinMaxTotalReductionOp<int, mj_scalar_t>
			reductionOp(
					current_concurrent_num_parts,
					current_concurrent_num_parts,
					current_concurrent_num_parts);
		try{
			reduceAll<int, mj_scalar_t>(
					*(this->comm),
					reductionOp,
					3 * current_concurrent_num_parts,
					local_min_max_total,
					global_min_max_total);
		}
		Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))
	}
	else {
		mj_part_t s = 3 * current_concurrent_num_parts;
		for (mj_part_t i = 0; i < s; ++i){
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_get_initial_cut_coords_target_weights(
    mj_scalar_t min_coord,
    mj_scalar_t max_coord,
    mj_part_t num_cuts/*p-1*/ ,
    mj_scalar_t global_weight,
    mj_scalar_t *initial_cut_coords /*p - 1 sized, coordinate of each cut line*/,
    mj_scalar_t *current_target_part_weights /*cumulative weights, at left side of each cut line. p-1 sized*/,

    std::vector <mj_part_t> *future_num_part_in_parts, //the vecto
    std::vector <mj_part_t> *next_future_num_parts_in_parts,
    mj_part_t concurrent_current_part,
    mj_part_t obtained_part_index
){

    mj_scalar_t coord_range = max_coord - min_coord;
    if(this->mj_uniform_parts[0]){
        {
            mj_part_t cumulative = 0;
            //how many total future parts the part will be partitioned into.
            mj_scalar_t total_future_part_count_in_part = mj_scalar_t((*future_num_part_in_parts)[concurrent_current_part]);


            //how much each part should weigh in ideal case.
            mj_scalar_t unit_part_weight = global_weight / total_future_part_count_in_part;
            /*
            cout << "total_future_part_count_in_part:" << total_future_part_count_in_part << endl;
            cout << "global_weight:" << global_weight << endl;
            cout << "unit_part_weight" << unit_part_weight <<endl;
            */
            for(mj_part_t i = 0; i < num_cuts; ++i){
                cumulative += (*next_future_num_parts_in_parts)[i + obtained_part_index];

                /*
                cout << "obtained_part_index:" << obtained_part_index <<
                		" (*next_future_num_parts_in_parts)[i + obtained_part_index]:" << (*next_future_num_parts_in_parts)[i + obtained_part_index] <<
                		" cumulative:" << cumulative << endl;
                */
                //set target part weight.
                current_target_part_weights[i] = cumulative * unit_part_weight;
                //cout <<"i:" << i << " current_target_part_weights:" << current_target_part_weights[i] << endl;
                //set initial cut coordinate.
                initial_cut_coords[i] = min_coord + (coord_range *
                                         cumulative) / total_future_part_count_in_part;
            }
            current_target_part_weights[num_cuts] = 1;
        }

        //round the target part weights.
        if (this->mj_uniform_weights[0]){
        	for(mj_part_t i = 0; i < num_cuts + 1; ++i){
                current_target_part_weights[i] = long(current_target_part_weights[i] + 0.5);
            }
        }
    }
    else {
    	std::cerr << "MJ does not support non uniform part weights" << std::endl;
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
 * \param mj_current_coordinate_permutations is the permutation array, holds the real indices of coordinates on mj_current_dim_coords array.
 * \param mj_current_dim_coords is the 1D array holding the coordinates.
 * \param mj_part_ids is the array holding the partIds of each coordinate.
 * \param partition_count is the number of parts that the current part will be partitioned into.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::set_initial_coordinate_parts(
    mj_scalar_t &max_coordinate,
    mj_scalar_t &min_coordinate,
    mj_part_t &concurrent_current_part_index,
    mj_lno_t coordinate_begin_index,
    mj_lno_t coordinate_end_index,
    mj_lno_t *mj_current_coordinate_permutations,
    mj_scalar_t *mj_current_dim_coords,
    mj_part_t *mj_part_ids,
    mj_part_t &partition_count
){
    mj_scalar_t coordinate_range = max_coordinate - min_coordinate;

    //if there is single point, or if all points are along a line.
    //set initial part to 0 for all.
    if(ABS(coordinate_range) < this->sEpsilon ){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(mj_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){
        	mj_part_ids[mj_current_coordinate_permutations[ii]] = 0;
        }
    }
    else{

        //otherwise estimate an initial part for each coordinate.
        //assuming uniform distribution of points.
        mj_scalar_t slice = coordinate_range / partition_count;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(mj_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){

            mj_lno_t iii = mj_current_coordinate_permutations[ii];
            mj_part_t pp = mj_part_t((mj_current_dim_coords[iii] - min_coordinate) / slice);
            mj_part_ids[iii] = 2 * pp;
        }
    }
}


/*! \brief Function that is responsible from 1D partitioning of the given range of coordinates.
 * \param mj_current_dim_coords is 1 dimensional array holding coordinate values.
 * \param imbalanceTolerance is the maximum allowed imbalance ratio.
 * \param current_work_part is the beginning index of concurrentPartCount parts.
 * \param current_concurrent_num_parts is the number of parts whose cut lines will be calculated concurrently.
 * \param current_cut_coordinates is the array holding the coordinates of the cut.
 * \param total_incomplete_cut_count is the number of cut lines whose positions should be calculated.
 * \param num_partitioning_in_current_dim is the vector that holds how many parts each part will be divided into.
 *
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_1D_part(
    mj_scalar_t *mj_current_dim_coords,
    mj_scalar_t used_imbalance_tolerance,
    mj_part_t current_work_part,
    mj_part_t current_concurrent_num_parts,
    mj_scalar_t *current_cut_coordinates,
    mj_part_t total_incomplete_cut_count,
    std::vector <mj_part_t> &num_partitioning_in_current_dim
){


    mj_part_t rectilinear_cut_count = 0;
    mj_scalar_t *temp_cut_coords = current_cut_coordinates;

    Teuchos::MultiJaggedCombinedReductionOp<mj_part_t, mj_scalar_t>
                 *reductionOp = NULL;
    reductionOp = new Teuchos::MultiJaggedCombinedReductionOp
                     <mj_part_t, mj_scalar_t>(
                    		 &num_partitioning_in_current_dim ,
                    		 current_work_part ,
                    		 current_concurrent_num_parts);

    size_t total_reduction_size = 0;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel shared(total_incomplete_cut_count,  rectilinear_cut_count)
#endif
    {
        int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
        me = omp_get_thread_num();
#endif
        double *my_thread_part_weights = this->thread_part_weights[me];
        mj_scalar_t *my_thread_left_closest = this->thread_cut_left_closest_point[me];
        mj_scalar_t *my_thread_right_closest = this->thread_cut_right_closest_point[me];

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp single
#endif
            {
                //initialize the lower and upper bounds of the cuts.
                mj_part_t next = 0;
                for(mj_part_t i = 0; i < current_concurrent_num_parts; ++i){

                    mj_part_t num_part_in_dim =  num_partitioning_in_current_dim[current_work_part + i];
                    mj_part_t num_cut_in_dim = num_part_in_dim - 1;
                    total_reduction_size += (4 * num_cut_in_dim + 1);

                    for(mj_part_t ii = 0; ii < num_cut_in_dim; ++ii){
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
            mj_part_t concurrent_cut_shifts = 0;
            size_t total_part_shift = 0;

            for (mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk){
                mj_part_t num_parts =  -1;
                num_parts =  num_partitioning_in_current_dim[current_work_part + kk];

                mj_part_t num_cuts = num_parts - 1;
                size_t total_part_count = num_parts + size_t (num_cuts) ;
                if (this->my_incomplete_cut_count[kk] > 0){

                    //although isDone shared, currentDone is private and same for all.
                    bool *current_cut_status = this->is_cut_line_determined + concurrent_cut_shifts;
                    double *my_current_part_weights = my_thread_part_weights + total_part_shift;
                    mj_scalar_t *my_current_left_closest = my_thread_left_closest + concurrent_cut_shifts;
                    mj_scalar_t *my_current_right_closest = my_thread_right_closest + concurrent_cut_shifts;

                    mj_part_t conccurent_current_part = current_work_part + kk;
                    mj_lno_t coordinate_begin_index = conccurent_current_part ==0 ? 0: this->part_xadj[conccurent_current_part -1];
                    mj_lno_t coordinate_end_index = this->part_xadj[conccurent_current_part];
                    mj_scalar_t *temp_current_cut_coords = temp_cut_coords + concurrent_cut_shifts;

                    mj_scalar_t min_coord = global_min_max_coord_total_weight[kk];
                    mj_scalar_t max_coord = global_min_max_coord_total_weight[kk + current_concurrent_num_parts];

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
                        reduceAll<int, mj_scalar_t>( *(this->comm), *reductionOp,
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
                            total_reduction_size * sizeof(mj_scalar_t));
                }
            }

            //how much cut will be shifted for the next part in the concurrent part calculation.
            mj_part_t cut_shift = 0;

            //how much the concantaneted array will be shifted for the next part in concurrent part calculation.
            size_t tlr_shift = 0;
            for (mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk){
            	mj_part_t num_parts =  num_partitioning_in_current_dim[current_work_part + kk];
            	mj_part_t num_cuts = num_parts - 1;
            	size_t num_total_part = num_parts + size_t (num_cuts) ;

            	//if the cuts of this cut has already been completed.
            	//nothing to do for this part.
            	//just update the shift amount and proceed.
            	if (this->my_incomplete_cut_count[kk] == 0) {
            		cut_shift += num_cuts;
            		tlr_shift += (num_total_part + 2 * num_cuts);
            		continue;
            	}

            	mj_scalar_t *current_local_part_weights = this->total_part_weight_left_right_closests  + tlr_shift ;
            	mj_scalar_t *current_global_tlr = this->global_total_part_weight_left_right_closests + tlr_shift;
            	mj_scalar_t *current_global_left_closest_points = current_global_tlr + num_total_part; //left closest points
            	mj_scalar_t *current_global_right_closest_points = current_global_tlr + num_total_part + num_cuts; //right closest points
            	mj_scalar_t *current_global_part_weights = current_global_tlr;
            	bool *current_cut_line_determined = this->is_cut_line_determined + cut_shift;

            	mj_scalar_t *current_part_target_weights = this->target_part_weights + cut_shift + kk;
            	mj_scalar_t *current_part_cut_line_weight_to_put_left = this->process_cut_line_weight_to_put_left + cut_shift;

            	mj_scalar_t min_coordinate = global_min_max_coord_total_weight[kk];
            	mj_scalar_t max_coordinate = global_min_max_coord_total_weight[kk + current_concurrent_num_parts];
            	mj_scalar_t global_total_weight = global_min_max_coord_total_weight[kk + current_concurrent_num_parts * 2];
            	mj_scalar_t *current_cut_lower_bound_weights = this->cut_lower_bound_weights + cut_shift;
            	mj_scalar_t *current_cut_upper_weights = this->cut_upper_bound_weights + cut_shift;
            	mj_scalar_t *current_cut_upper_bounds = this->cut_upper_bound_coordinates + cut_shift;
            	mj_scalar_t *current_cut_lower_bounds = this->cut_lower_bound_coordinates + cut_shift;

            	mj_part_t initial_incomplete_cut_count = this->my_incomplete_cut_count[kk];

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
            			&rectilinear_cut_count,
            			this->my_incomplete_cut_count[kk]);

            	cut_shift += num_cuts;
            	tlr_shift += (num_total_part + 2 * num_cuts);
            	mj_part_t iteration_complete_cut_count = initial_incomplete_cut_count - this->my_incomplete_cut_count[kk];
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
                mj_scalar_t *t = temp_cut_coords;
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
        		mj_part_t next = 0;
        		for(mj_part_t i = 0; i < current_concurrent_num_parts; ++i){
        			mj_part_t num_parts = -1;
        			num_parts = num_partitioning_in_current_dim[current_work_part + i];
        			mj_part_t num_cuts = num_parts - 1;

        			for(mj_part_t ii = 0; ii < num_cuts; ++ii){
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_1D_part_get_thread_part_weights(
    size_t total_part_count,
    mj_part_t num_cuts,
    mj_scalar_t max_coord,
    mj_scalar_t min_coord,
    mj_lno_t coordinate_begin_index,
    mj_lno_t coordinate_end_index,
    mj_scalar_t *mj_current_dim_coords,
    mj_scalar_t *temp_current_cut_coords,
    bool *current_cut_status,
    double *my_current_part_weights,
    mj_scalar_t *my_current_left_closest,
    mj_scalar_t *my_current_right_closest){

	// initializations for part weights, left/right closest
	for (size_t i = 0; i < total_part_count; ++i){
		my_current_part_weights[i] = 0;
	}

	//initialize the left and right closest coordinates
	//to their max value.
	for(mj_part_t i = 0; i < num_cuts; ++i){
		my_current_left_closest[i] = min_coord - 1;
		my_current_right_closest[i] = max_coord + 1;
	}
	//mj_lno_t comparison_count = 0;
	mj_scalar_t minus_EPSILON = -this->sEpsilon;
#ifdef HAVE_ZOLTAN2_OMP
	//no need for the barrier as all threads uses their local memories.
	//dont change the static scheduling here, as it is assumed when the new
	//partitions are created later.
#pragma omp for
#endif
	for (mj_lno_t ii = coordinate_begin_index; ii < coordinate_end_index; ++ii){
		int i = this->coordinate_permutations[ii];

		//the accesses to assigned_part_ids are thread safe
		//since each coordinate is assigned to only a single thread.
		mj_part_t j = this->assigned_part_ids[i] / 2;

		if(j >= num_cuts){
			j = num_cuts - 1;
		}

		mj_part_t lower_cut_index = 0;
		mj_part_t upper_cut_index = num_cuts - 1;

		mj_scalar_t w = this->mj_uniform_weights[0]? 1:this->mj_weights[0][i];
		bool is_inserted = false;
		bool is_on_left_of_cut = false;
		bool is_on_right_of_cut = false;
		mj_part_t last_compared_part = -1;

		mj_scalar_t coord = mj_current_dim_coords[i];

		while(upper_cut_index >= lower_cut_index)
		{
			//comparison_count++;
			last_compared_part = -1;
			is_on_left_of_cut = false;
			is_on_right_of_cut = false;
			mj_scalar_t cut = temp_current_cut_coords[j];
			mj_scalar_t distance_to_cut = coord - cut;
			mj_scalar_t abs_distance_to_cut = ABS(distance_to_cut);

			//if it is on the line.
			if(abs_distance_to_cut < this->sEpsilon){

				my_current_part_weights[j * 2 + 1] += w;
				this->assigned_part_ids[i] = j * 2 + 1;

				//assign left and right closest point to cut as the point is on the cut.
				my_current_left_closest[j] = coord;
				my_current_right_closest[j] = coord;
				//now we need to check if there are other cuts on the same cut coordinate.
				//if there are, then we add the weight of the cut to all cuts in the same coordinate.
				mj_part_t kk = j + 1;
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
						mj_scalar_t distance_to_next_cut = coord - temp_current_cut_coords[j - 1];
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
						mj_scalar_t distance_to_next_cut = coord - temp_current_cut_coords[j + 1];
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_accumulate_thread_results(
    const std::vector <mj_part_t> &num_partitioning_in_current_dim,
    mj_part_t current_work_part,
    mj_part_t current_concurrent_num_parts){

#ifdef HAVE_ZOLTAN2_OMP
	//needs barrier here, as it requires all threads to finish mj_1D_part_get_thread_part_weights
	//using parallel region here reduces the performance because of the cache invalidates.
#pragma omp barrier
#pragma omp single
#endif
	{
		size_t tlr_array_shift = 0;
		mj_part_t cut_shift = 0;

		//iterate for all concurrent parts to find the left and right closest points in the process.
		for(mj_part_t i = 0; i < current_concurrent_num_parts; ++i){

			mj_part_t num_parts_in_part =  num_partitioning_in_current_dim[current_work_part + i];
			mj_part_t num_cuts_in_part = num_parts_in_part - 1;
			size_t num_total_part_in_part = num_parts_in_part + size_t (num_cuts_in_part) ;

			//iterate for cuts in a single part.
			for(mj_part_t ii = 0; ii < num_cuts_in_part ; ++ii){
				mj_part_t next = tlr_array_shift + ii;
				mj_part_t cut_index = cut_shift + ii;
				if(this->is_cut_line_determined[cut_index]) continue;
				mj_scalar_t left_closest_in_process = this->thread_cut_left_closest_point[0][cut_index],
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
		for(mj_part_t i = 0; i < current_concurrent_num_parts; ++i){

			mj_part_t num_parts_in_part =  num_partitioning_in_current_dim[current_work_part + i];
			mj_part_t num_cuts_in_part = num_parts_in_part - 1;
			size_t num_total_part_in_part = num_parts_in_part + size_t (num_cuts_in_part) ;

			for(size_t j = 0; j < num_total_part_in_part; ++j){

				mj_part_t cut_ind = j / 2 + cut_shift;

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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_calculate_new_cut_position (
	mj_scalar_t cut_upper_bound,
    mj_scalar_t cut_lower_bound,
    mj_scalar_t cut_upper_weight,
    mj_scalar_t cut_lower_weight,
    mj_scalar_t expected_weight,
    mj_scalar_t &new_cut_position){

    if(ABS(cut_upper_bound - cut_lower_bound) < this->sEpsilon){
    	new_cut_position = cut_upper_bound; //or lower bound does not matter.
    }


    if(ABS(cut_upper_weight - cut_lower_weight) < this->sEpsilon){
    	new_cut_position = cut_lower_bound;
    }

    mj_scalar_t coordinate_range = (cut_upper_bound - cut_lower_bound);
    mj_scalar_t weight_range = (cut_upper_weight - cut_lower_weight);
    mj_scalar_t my_weight_diff = (expected_weight - cut_lower_weight);

    mj_scalar_t required_shift = (my_weight_diff / weight_range);
    int scale_constant = 20;
    int shiftint= int (required_shift * scale_constant);
    if (shiftint == 0) shiftint = 1;
    required_shift = mj_scalar_t (shiftint) / scale_constant;
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
 * \param out_part_xadj is the indices of coordinates calculated for the partition on next dimension.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_create_new_partitions(
    mj_part_t num_parts,
    mj_scalar_t *mj_current_dim_coords,
    mj_scalar_t *current_concurrent_cut_coordinate,
    mj_lno_t coordinate_begin,
    mj_lno_t coordinate_end,
    mj_scalar_t *used_local_cut_line_weight_to_left,
    double **used_thread_part_weight_work,
    mj_lno_t *out_part_xadj){

	mj_part_t num_cuts = num_parts - 1;

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
	{
		int me = 0;
#ifdef HAVE_ZOLTAN2_OMP
		me = omp_get_thread_num();
#endif

		mj_lno_t *thread_num_points_in_parts = this->thread_point_counts[me];
		mj_scalar_t *my_local_thread_cut_weights_to_put_left = NULL;

		//now if the rectilinear partitioning is allowed we decide how
		//much weight each thread should put to left and right.
		if (this->distribute_points_on_cut_lines){
			my_local_thread_cut_weights_to_put_left = this->thread_cut_line_weight_to_put_left[me];
			// this for assumes the static scheduling in mj_1D_part calculation.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
			for (mj_part_t i = 0; i < num_cuts; ++i){
				//the left to be put on the left of the cut.
				mj_scalar_t left_weight = used_local_cut_line_weight_to_left[i];
				for(int ii = 0; ii < this->num_threads; ++ii){
					if(left_weight > this->sEpsilon){
						//the weight of thread ii on cut.
						mj_scalar_t thread_ii_weight_on_cut = used_thread_part_weight_work[ii][i * 2 + 1] - used_thread_part_weight_work[ii][i * 2 ];
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
				for (mj_part_t i = num_cuts - 1; i > 0 ; --i){
					if(ABS(current_concurrent_cut_coordinate[i] - current_concurrent_cut_coordinate[i -1]) < this->sEpsilon){
						my_local_thread_cut_weights_to_put_left[i] -= my_local_thread_cut_weights_to_put_left[i - 1] ;
					}
					my_local_thread_cut_weights_to_put_left[i] = int ((my_local_thread_cut_weights_to_put_left[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL)
                    						/ mj_scalar_t(SIGNIFICANCE_MUL);
				}
			}
		}

		for(mj_part_t ii = 0; ii < num_parts; ++ii){
			thread_num_points_in_parts[ii] = 0;
		}


#ifdef HAVE_ZOLTAN2_OMP
		//dont change static scheduler. the static partitioner used later as well.
#pragma omp for
#endif
		for (mj_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){

			mj_lno_t coordinate_index = this->coordinate_permutations[ii];
			mj_scalar_t coordinate_weight = this->mj_uniform_weights[0]? 1:this->mj_weights[0][coordinate_index];
			mj_part_t coordinate_assigned_place = this->assigned_part_ids[coordinate_index];
			mj_part_t coordinate_assigned_part = coordinate_assigned_place / 2;
			if(coordinate_assigned_place % 2 == 1){
				//if it is on the cut.
				if(this->distribute_points_on_cut_lines
						&& my_local_thread_cut_weights_to_put_left[coordinate_assigned_part] > this->sEpsilon){
					//if the rectilinear partitioning is allowed,
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
		for(mj_part_t j = 0; j < num_parts; ++j){
			mj_lno_t num_points_in_part_j_upto_thread_i = 0;
			for (int i = 0; i < this->num_threads; ++i){
				mj_lno_t thread_num_points_in_part_j = this->thread_point_counts[i][j];
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
			for(mj_part_t j = 1; j < num_parts; ++j){
				out_part_xadj[j] += out_part_xadj[j - 1];
			}
		}

		//shift the num points in threads thread to obtain the
		//beginning index of each thread's private space.
		for(mj_part_t j = 1; j < num_parts; ++j){
			thread_num_points_in_parts[j] += out_part_xadj[j - 1] ;
		}


		//now thread gets the coordinate and writes the index of coordinate to the permutation array
		//using the part index we calculated.
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
		for (mj_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){
			mj_lno_t i = this->coordinate_permutations[ii];
			mj_part_t p =  this->assigned_part_ids[i];
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
 * \param rectilinear_cut_count is the count of cut lines whose balance can be achived via distributing the points in same coordinate to different parts.
 * \param my_num_incomplete_cut is the number of cutlines whose position has not been determined yet. For K > 1 it is the count in a single part (whose cut lines are determined).
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_get_new_cut_coordinates(
		const size_t &num_total_part,
		const mj_part_t &num_cuts,
		const mj_scalar_t &max_coordinate,
		const mj_scalar_t &min_coordinate,
		const mj_scalar_t &global_total_weight,
		const mj_scalar_t &used_imbalance_tolerance,
		mj_scalar_t * current_global_part_weights,
		const mj_scalar_t * current_local_part_weights,
		const mj_scalar_t *current_part_target_weights,
		bool *current_cut_line_determined,
		mj_scalar_t *current_cut_coordinates,
		mj_scalar_t *current_cut_upper_bounds,
		mj_scalar_t *current_cut_lower_bounds,
		mj_scalar_t *current_global_left_closest_points,
		mj_scalar_t *current_global_right_closest_points,
		mj_scalar_t * current_cut_lower_bound_weights,
		mj_scalar_t * current_cut_upper_weights,
		mj_scalar_t *new_current_cut_coordinates,
		mj_scalar_t *current_part_cut_line_weight_to_put_left,
		mj_part_t *rectilinear_cut_count,
		mj_part_t &my_num_incomplete_cut){

	//seen weight in the part
	mj_scalar_t seen_weight_in_part = 0;
	//expected weight for part.
	mj_scalar_t expected_weight_in_part = 0;
	//imbalance for the left and right side of the cut.
	mj_scalar_t imbalance_on_left = 0, imbalance_on_right = 0;


#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
	for (mj_part_t i = 0; i < num_cuts; i++){
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
	for (mj_part_t i = 0; i < num_cuts; i++){

		if(this->distribute_points_on_cut_lines){
			//init the weight on the cut.
			this->global_rectilinear_cut_weight[i] = 0;
			this->process_rectilinear_cut_weight[i] = 0;
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
					*rectilinear_cut_count += 1;
					//increase the num cuts to be determined with rectilinear partitioning.

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp atomic
#endif
					my_num_incomplete_cut -= 1;
					new_current_cut_coordinates [i] = current_cut_coordinates[i];
					this->process_rectilinear_cut_weight[i] = current_local_part_weights[i * 2 + 1] -
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
			for (mj_part_t ii = i + 1; ii < num_cuts ; ++ii){
				mj_scalar_t p_weight = current_global_part_weights[ii * 2];
				mj_scalar_t line_weight = current_global_part_weights[ii * 2 + 1];

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


			mj_scalar_t new_cut_position = 0;
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
				mj_scalar_t p_weight = current_global_part_weights[ii * 2];
				mj_scalar_t line_weight = current_global_part_weights[ii * 2 + 1];
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
			mj_scalar_t new_cut_position = 0;
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
		if(*rectilinear_cut_count > 0){

			try{
				Teuchos::scan<int,mj_scalar_t>(
						*comm, Teuchos::REDUCE_SUM,
						num_cuts,
						this->process_rectilinear_cut_weight,
						this->global_rectilinear_cut_weight
				);
			}
			Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))

			for (mj_part_t i = 0; i < num_cuts; ++i){
				//if cut line weight to be distributed.
				if(this->global_rectilinear_cut_weight[i] > 0) {
					//expected weight to go to left of the cut.
					mj_scalar_t expected_part_weight = current_part_target_weights[i];
					//the weight that should be put to left of the cut.
					mj_scalar_t necessary_weight_on_line_for_left = expected_part_weight - current_global_part_weights[i * 2];
					//the weight of the cut in the process
					mj_scalar_t my_weight_on_line = this->process_rectilinear_cut_weight[i];
					//the sum of the cut weights upto this process, including the weight of this process.
					mj_scalar_t weight_on_line_upto_process_inclusive = this->global_rectilinear_cut_weight[i];
					//the space on the left side of the cut after all processes before this process (including this process)
					//puts their weights on cut to left.
					mj_scalar_t space_to_put_left = necessary_weight_on_line_for_left - weight_on_line_upto_process_inclusive;
					//add my weight to this space to find out how much space is left to me.
					mj_scalar_t space_left_to_me = space_to_put_left + my_weight_on_line;

					/*
					cout << "expected_part_weight:" << expected_part_weight
							<< " necessary_weight_on_line_for_left:" << necessary_weight_on_line_for_left
							<< " my_weight_on_line" << my_weight_on_line
							<< " weight_on_line_upto_process_inclusive:" << weight_on_line_upto_process_inclusive
							<< " space_to_put_left:" << space_to_put_left
							<< " space_left_to_me" << space_left_to_me << endl;
					 */
					if(space_left_to_me < 0){
						//space_left_to_me is negative and i dont need to put anything to left.
						current_part_cut_line_weight_to_put_left[i] = 0;
					}
					else if(space_left_to_me >= my_weight_on_line){
						//space left to me is bigger than the weight of the processor on cut.
						//so put everything to left.
						current_part_cut_line_weight_to_put_left[i] = my_weight_on_line;
						//cout << "setting current_part_cut_line_weight_to_put_left to my_weight_on_line:" << my_weight_on_line << endl;
					}
					else {
						//put only the weight as much as the space.
						current_part_cut_line_weight_to_put_left[i] = space_left_to_me ;

						//cout << "setting current_part_cut_line_weight_to_put_left to space_left_to_me:" << space_left_to_me << endl;
					}

				}
			}
			*rectilinear_cut_count = 0;
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::get_processor_num_points_in_parts(
    mj_part_t num_procs,
    mj_part_t num_parts,
    mj_gno_t *&num_points_in_all_processor_parts){

	//initially allocation_size is num_parts
	size_t allocation_size = num_parts * (num_procs + 1);

	//this will be output
	//holds how many each processor has in each part.
	//last portion is the sum of all processor points in each part.

	//allocate memory for the local num coordinates in each part.
	mj_gno_t *num_local_points_in_each_part_to_reduce_sum = allocMemory<mj_gno_t>(allocation_size);


	//this is the portion of the memory which will be used
	//at the summation to obtain total number of processors' points in each part.
	mj_gno_t *my_local_points_to_reduce_sum = num_local_points_in_each_part_to_reduce_sum + num_procs * num_parts;
	//this is the portion of the memory where each stores its local number.
	//this information is needed by other processors.
	mj_gno_t *my_local_point_counts_in_each_art = num_local_points_in_each_part_to_reduce_sum + this->myRank * num_parts;

	//initialize the array with 0's.
	memset(num_local_points_in_each_part_to_reduce_sum, 0, sizeof(mj_gno_t)*allocation_size);

	//write the number of coordinates in each part.
	for (mj_part_t i = 0; i < num_parts; ++i){
		mj_lno_t part_begin_index = 0;
		if (i > 0){
			part_begin_index = this->new_part_xadj[i - 1];
		}
		mj_lno_t part_end_index = this->new_part_xadj[i];
		my_local_points_to_reduce_sum[i] = part_end_index - part_begin_index;
	}

	//copy the local num parts to the last portion of array,
	//so that this portion will represent the global num points in each part after the reduction.
	memcpy (my_local_point_counts_in_each_art,
			my_local_points_to_reduce_sum,
			sizeof(mj_gno_t) * (num_parts) );


	//reduceAll operation.
	//the portion that belongs to a processor with index p
	//will start from myRank * num_parts.
	//the global number of points will be held at the index
	try{
		reduceAll<int, mj_gno_t>(
				*(this->comm),
				Teuchos::REDUCE_SUM,
				allocation_size,
				num_local_points_in_each_part_to_reduce_sum,
				num_points_in_all_processor_parts);
	}
	Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))
	freeArray<mj_gno_t>(num_local_points_in_each_part_to_reduce_sum);
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
bool AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_check_to_migrate(
    size_t migration_reduce_all_population,
    mj_lno_t num_coords_for_last_dim_part,
    mj_part_t num_procs,
    mj_part_t num_parts,
    mj_gno_t *num_points_in_all_processor_parts){

	//if reduce all count and population in the last dim is too high
    if (migration_reduce_all_population > FUTURE_REDUCEALL_CUTOFF) return true;
    //if the work in a part per processor in the last dim is too low.
    if (num_coords_for_last_dim_part < MIN_WORK_LAST_DIM) return true;

	//if migration is to be checked and the imbalance is too high
    if (this->check_migrate_avoid_migration_option == 0){
    	double global_imbalance = 0;
    	//global shift to reach the sum of coordiante count in each part.
    	size_t global_shift = num_procs * num_parts;

    	for (mj_part_t ii = 0; ii < num_procs; ++ii){
    		for (mj_part_t i = 0; i < num_parts; ++i){
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


/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent.
 *
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param part_assignment_proc_begin_indices ([i]) points to the first processor index that part i will be sent to.
 * \param processor_chains_in_parts the array that holds the linked list structure, started from part_assignment_proc_begin_indices ([i]).
 * \param send_count_to_each_proc array array storing the number of points to be sent to each part.
 * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::assign_send_destinations(
    mj_part_t num_parts,
    mj_part_t *part_assignment_proc_begin_indices,
    mj_part_t *processor_chains_in_parts,
    mj_lno_t *send_count_to_each_proc,
    int *coordinate_destinations){

    for (mj_part_t p = 0; p < num_parts; ++p){
        mj_lno_t part_begin = 0;
        if (p > 0) part_begin = this->new_part_xadj[p - 1];
        mj_lno_t part_end = this->new_part_xadj[p];

        //get the first part that current processor will send its part-p.
        mj_part_t proc_to_sent = part_assignment_proc_begin_indices[p];
        //initialize how many point I sent to this processor.
        mj_lno_t num_total_send = 0;
        for (mj_lno_t j=part_begin; j < part_end; j++){
            mj_lno_t local_ind = this->new_coordinate_permutations[j];
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
            coordinate_destinations[local_ind] = proc_to_sent;
            ++num_total_send;
        }
    }
}

/*! \brief Function fills up coordinate_destinations is the output array
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
 * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_assign_proc_to_parts(
		mj_gno_t * num_points_in_all_processor_parts,
		mj_part_t num_parts,
		mj_part_t num_procs,
		mj_lno_t *send_count_to_each_proc,
		std::vector<mj_part_t> &processor_ranks_for_subcomm,
		std::vector<mj_part_t> *next_future_num_parts_in_parts,
		mj_part_t &out_part_index,
		mj_part_t &output_part_numbering_begin_index,
		int *coordinate_destinations){


    mj_gno_t *global_num_points_in_parts = num_points_in_all_processor_parts + num_procs * num_parts;
    mj_part_t *num_procs_assigned_to_each_part = allocMemory<mj_part_t>(num_parts);

    //boolean variable if the process finds its part to be assigned.
    bool did_i_find_my_group = false;

    mj_part_t num_free_procs = num_procs;
    mj_part_t minimum_num_procs_required_for_rest_of_parts = num_parts - 1;

    double max_imbalance_difference = 0;
    mj_part_t max_differing_part = 0;

    //find how many processor each part requires.
    for (mj_part_t i=0; i < num_parts; i++){

        //scalar portion of the required processors
        double scalar_required_proc = num_procs *
                (double (global_num_points_in_parts[i]) / double (this->num_global_coords));

        //round it to closest integer.
        mj_part_t required_proc = static_cast<mj_part_t> (0.5 + scalar_required_proc);

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
    mj_part_t *part_assignment_proc_begin_indices = allocMemory<mj_part_t>(num_parts);
    //the next processor send is found in processor_chains_in_parts, in linked list manner.
    mj_part_t *processor_chains_in_parts = allocMemory<mj_part_t>(num_procs);
    mj_part_t *processor_part_assignments = allocMemory<mj_part_t>(num_procs);

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
    uSortItem<mj_part_t, mj_gno_t> * sort_item_num_part_points_in_procs = allocMemory <uSortItem<mj_part_t, mj_gno_t> > (num_procs);
    for(mj_part_t i = 0; i < num_parts; ++i){
        //the algorithm tries to minimize the cost of migration,
    	//by assigning the processors with highest number of coordinates on that part.
    	//here we might want to implement a maximum weighted bipartite matching algorithm.
    	for(mj_part_t ii = 0; ii < num_procs; ++ii){
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
        uqsort<mj_part_t, mj_gno_t>(num_procs, sort_item_num_part_points_in_procs);

        mj_part_t required_proc_count =  num_procs_assigned_to_each_part[i];
        mj_gno_t total_num_points_in_part = global_num_points_in_parts[i];
        mj_gno_t ideal_num_points_in_a_proc =
                Teuchos::as<mj_gno_t>(ceil (total_num_points_in_part / double (required_proc_count)));

        //starts sending to least heaviest part.
        mj_part_t next_proc_to_send_index = num_procs - required_proc_count;
        mj_part_t next_proc_to_send_id = sort_item_num_part_points_in_procs[next_proc_to_send_index].id;
        mj_lno_t space_left_in_sent_proc =  ideal_num_points_in_a_proc - sort_item_num_part_points_in_procs[next_proc_to_send_index].val;

        //find the processors that will be assigned to this part, which are the heaviest
        //non assigned processors.
        for(mj_part_t ii = num_procs - 1; ii >= num_procs - required_proc_count; --ii){
            mj_part_t proc_id = sort_item_num_part_points_in_procs[ii].id;
            //assign processor to part - i.
            processor_part_assignments[proc_id] = i;
        }

        bool did_change_sign = false;
        //if processor has a minus count, reverse it.
        for(mj_part_t ii = 0; ii < num_procs; ++ii){
            // TODO:  THE LINE BELOW PRODUCES A WARNING IF gno_t IS UNSIGNED
            // TODO:  SEE BUG 6194
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
            uqsort<mj_part_t, mj_gno_t>(num_procs - required_proc_count, sort_item_num_part_points_in_procs);
        }

        //check if this processors is one of the procs assigned to this part.
        //if it is, then get the group.
        if (!did_i_find_my_group){
            for(mj_part_t ii = num_procs - 1; ii >= num_procs - required_proc_count; --ii){

            	mj_part_t proc_id_to_assign = sort_item_num_part_points_in_procs[ii].id;
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
                    for (mj_part_t in = 0; in < i; ++in){
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
        for(mj_part_t ii = num_procs - required_proc_count - 1; ii >= 0; --ii){
            mj_part_t nonassigned_proc_id = sort_item_num_part_points_in_procs[ii].id;
            mj_lno_t num_points_to_sent = sort_item_num_part_points_in_procs[ii].val;

            //we set number of points to -to_sent - 1 for the assigned processors.
            //we reverse it here. This should not happen, as we have already reversed them above.
#ifdef MJ_DEBUG
            if (num_points_to_sent < 0) {
            	cout << "Migration - processor assignments - for part:" << i << "from proc:" << nonassigned_proc_id << " num_points_to_sent:" << num_points_to_sent << std::endl;
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
                        mj_part_t prev_begin = part_assignment_proc_begin_indices[i];
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
                            mj_part_t prev_begin = part_assignment_proc_begin_indices[i];
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
                    			<< " Error: next_part_to_send_index <  nprocs - required_proc_count" << std::endl;
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
            coordinate_destinations);

    freeArray<mj_part_t>(part_assignment_proc_begin_indices);
    freeArray<mj_part_t>(processor_chains_in_parts);
    freeArray<mj_part_t>(processor_part_assignments);
    freeArray<uSortItem<mj_part_t, mj_gno_t> > (sort_item_num_part_points_in_procs);
    freeArray<mj_part_t > (num_procs_assigned_to_each_part);

}


/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent. In addition it calculates
 * the shift amount (output_part_numbering_begin_index) to be done when
 * final numberings of the parts are performed.
 *
 * \param num_parts is the number of parts that exist in the current partitioning.
 * \param sort_item_part_to_proc_assignment is the sorted parts with respect to the assigned processors.
 * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
 * \param output_part_numbering_begin_index is how much the numbers should be shifted when numbering the result parts.
 * \param next_future_num_parts_in_parts is the vector, how many more parts each part will be divided into in the future.
 *
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::assign_send_destinations2(
    mj_part_t num_parts,
    uSortItem<mj_part_t, mj_part_t> * sort_item_part_to_proc_assignment, //input sorted wrt processors
    int *coordinate_destinations,
    mj_part_t &output_part_numbering_begin_index,
    std::vector<mj_part_t> *next_future_num_parts_in_parts){

    mj_part_t part_shift_amount = output_part_numbering_begin_index;
    mj_part_t previous_processor = -1;
    for(mj_part_t i = 0; i < num_parts; ++i){
        mj_part_t p = sort_item_part_to_proc_assignment[i].id;
        //assigned processors are sorted.
        mj_lno_t part_begin_index = 0;
        if (p > 0) part_begin_index = this->new_part_xadj[p - 1];
        mj_lno_t part_end_index = this->new_part_xadj[p];

        mj_part_t assigned_proc = sort_item_part_to_proc_assignment[i].val;
        if (this->myRank == assigned_proc && previous_processor != assigned_proc){
            output_part_numbering_begin_index =  part_shift_amount;
        }
        previous_processor = assigned_proc;
        part_shift_amount += (*next_future_num_parts_in_parts)[p];

        for (mj_lno_t j=part_begin_index; j < part_end_index; j++){
            mj_lno_t localInd = this->new_coordinate_permutations[j];
            coordinate_destinations[localInd] = assigned_proc;
        }
    }
}


/*! \brief Function fills up coordinate_destinations is the output array
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
 * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_assign_parts_to_procs(
    mj_gno_t * num_points_in_all_processor_parts,
    mj_part_t num_parts,
    mj_part_t num_procs,
    mj_lno_t *send_count_to_each_proc, //output: sized nprocs, show the number of send point counts to each proc.
    std::vector<mj_part_t> *next_future_num_parts_in_parts,//input how many more partitions the part will be partitioned into.
    mj_part_t &out_num_part, //output, how many parts the processor will have. this is always 1 for this function.
    std::vector<mj_part_t> &out_part_indices, //output: the part indices which the processor is assigned to.
    mj_part_t &output_part_numbering_begin_index, //output: how much the part number should be shifted when setting the solution
    int *coordinate_destinations){
    out_num_part = 0;

    mj_gno_t *global_num_points_in_parts = num_points_in_all_processor_parts + num_procs * num_parts;
    out_part_indices.clear();

    //to sort the parts that is assigned to the processors.
    //id is the part number, sort value is the assigned processor id.
    uSortItem<mj_part_t, mj_part_t> * sort_item_part_to_proc_assignment  = allocMemory <uSortItem<mj_part_t, mj_part_t> >(num_parts);
    uSortItem<mj_part_t, mj_gno_t> * sort_item_num_points_of_proc_in_part_i = allocMemory <uSortItem<mj_part_t, mj_gno_t> >(num_procs);


    //calculate the optimal number of coordinates that should be assigned to each processor.
    mj_lno_t work_each = mj_lno_t (this->num_global_coords / (double (num_procs)) + 0.5f);
    //to hold the left space as the number of coordinates to the optimal number in each proc.
    mj_lno_t *space_in_each_processor = allocMemory <mj_lno_t>(num_procs);
    //initialize left space in each.
    for (mj_part_t i = 0; i < num_procs; ++i){
        space_in_each_processor[i] = work_each;
    }

    //we keep track of how many parts each processor is assigned to.
    //because in some weird inputs, it might be possible that some
    //processors is not assigned to any part. Using these variables,
    //we force each processor to have at least one part.
    mj_part_t *num_parts_proc_assigned = allocMemory <mj_part_t>(num_procs);
    memset(num_parts_proc_assigned, 0, sizeof(mj_part_t) * num_procs);
    int empty_proc_count = num_procs;

    //to sort the parts with decreasing order of their coordiantes.
    //id are the part numbers, sort value is the number of points in each.
    uSortItem<mj_part_t, mj_gno_t> * sort_item_point_counts_in_parts  = allocMemory <uSortItem<mj_part_t, mj_gno_t> >(num_parts);

    //initially we will sort the parts according to the number of coordinates they have.
    //so that we will start assigning with the part that has the most number of coordinates.
    for (mj_part_t i = 0; i < num_parts; ++i){
        sort_item_point_counts_in_parts[i].id = i;
        sort_item_point_counts_in_parts[i].val = global_num_points_in_parts[i];
    }
    //sort parts with increasing order of loads.
    uqsort<mj_part_t, mj_gno_t>(num_parts, sort_item_point_counts_in_parts);


    //assigning parts to the processors
    //traverse the part win decreasing order of load.
    //first assign the heaviest part.
    for (mj_part_t j = 0; j < num_parts; ++j){
        //sorted with increasing order, traverse inverse.
        mj_part_t i = sort_item_point_counts_in_parts[num_parts - 1 - j].id;
        //load of the part
        mj_gno_t load = global_num_points_in_parts[i];

        //assigned processors
        mj_part_t assigned_proc = -1;
        //if not fit best processor.
        mj_part_t best_proc_to_assign = 0;


        //sort processors with increasing number of points in this part.
        for (mj_part_t ii = 0; ii < num_procs; ++ii){
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
        uqsort<mj_part_t, mj_gno_t>(num_procs, sort_item_num_points_of_proc_in_part_i);

        //traverse all processors with decreasing load.
        for (mj_part_t iii = num_procs - 1; iii >= 0; --iii){
            mj_part_t ii = sort_item_num_points_of_proc_in_part_i[iii].id;
            mj_lno_t left_space = space_in_each_processor[ii] - load;
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
    freeArray<mj_part_t>(num_parts_proc_assigned);
    freeArray< uSortItem<mj_part_t, mj_gno_t> > (sort_item_num_points_of_proc_in_part_i);
    freeArray<uSortItem<mj_part_t, mj_gno_t> >(sort_item_point_counts_in_parts);
    freeArray<mj_lno_t >(space_in_each_processor);


    //sort assignments with respect to the assigned processors.
    uqsort<mj_part_t, mj_part_t>(num_parts, sort_item_part_to_proc_assignment);
    //fill sendBuf.


    this->assign_send_destinations2(
            num_parts,
            sort_item_part_to_proc_assignment,
            coordinate_destinations,
            output_part_numbering_begin_index,
            next_future_num_parts_in_parts);

    freeArray<uSortItem<mj_part_t, mj_part_t> >(sort_item_part_to_proc_assignment);
}


/*! \brief Function fills up coordinate_destinations is the output array
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
 * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_migration_part_proc_assignment(
    mj_gno_t * num_points_in_all_processor_parts,
    mj_part_t num_parts,
    mj_part_t num_procs,
    mj_lno_t *send_count_to_each_proc,
    std::vector<mj_part_t> &processor_ranks_for_subcomm,
    std::vector<mj_part_t> *next_future_num_parts_in_parts,
    mj_part_t &out_num_part,
    std::vector<mj_part_t> &out_part_indices,
    mj_part_t &output_part_numbering_begin_index,
    int *coordinate_destinations){



	processor_ranks_for_subcomm.clear();
	// if (this->num_local_coords > 0)
	if (num_procs > num_parts){
		//if there are more processors than the number of current part
		//then processors share the existing parts.
		//at the end each processor will have a single part,
		//but a part will be shared by a group of processors.
		mj_part_t out_part_index = 0;
		this->mj_assign_proc_to_parts(
				num_points_in_all_processor_parts,
				num_parts,
				num_procs,
				send_count_to_each_proc,
				processor_ranks_for_subcomm,
				next_future_num_parts_in_parts,
				out_part_index,
				output_part_numbering_begin_index,
				coordinate_destinations
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
				coordinate_destinations);
	}
}

/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent. In addition it calculates
 * the shift amount (output_part_numbering_begin_index) to be done when
 * final numberings of the parts are performed.
 *
 *
 * \param num_procs is the number of processor attending to migration operation.
 * \param num_new_local_points is the output to represent the new number of local points.
 * \param iteration is the string for the current iteration.
 * \param coordinate_destinations is the output array that holds which part each coordinate should be sent.
 * \param num_parts is the number of parts that exist in the current partitioning.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_migrate_coords(
    mj_part_t num_procs,
    mj_lno_t &num_new_local_points,
    std::string iteration,
    int *coordinate_destinations,
    mj_part_t num_parts)
{
#ifdef ENABLE_ZOLTAN_MIGRATION
    if (sizeof(mj_lno_t) <= sizeof(int)) {

        // Cannot use Zoltan_Comm with local ordinals larger than ints.
        // In Zoltan_Comm_Create, the cast int(this->num_local_coords) 
        // may overflow.

	ZOLTAN_COMM_OBJ *plan = NULL;
	MPI_Comm mpi_comm = Teuchos2MPI (this->comm);
	int num_incoming_gnos = 0;
	int message_tag = 7859;

	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Migration Z1PlanCreating-" + iteration);
	int ierr = Zoltan_Comm_Create(
			&plan,
			int(this->num_local_coords),
			coordinate_destinations,
			mpi_comm,
			message_tag,
			&num_incoming_gnos);
	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Migration Z1PlanCreating-" + iteration);

	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Migration Z1Migration-" + iteration);
	mj_gno_t *incoming_gnos = allocMemory< mj_gno_t>(num_incoming_gnos);

	//migrate gnos.
	message_tag++;
	ierr = Zoltan_Comm_Do(
			plan,
			message_tag,
			(char *) this->current_mj_gnos,
			sizeof(mj_gno_t),
			(char *) incoming_gnos);
	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

	freeArray<mj_gno_t>(this->current_mj_gnos);
	this->current_mj_gnos = incoming_gnos;


	//migrate coordinates
	for (int i = 0; i < this->coord_dim; ++i){
		message_tag++;
		mj_scalar_t *coord = this->mj_coordinates[i];

		this->mj_coordinates[i] = allocMemory<mj_scalar_t>(num_incoming_gnos);
		ierr = Zoltan_Comm_Do(
				plan,
				message_tag,
				(char *) coord,
				sizeof(mj_scalar_t),
				(char *) this->mj_coordinates[i]);
		Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
		freeArray<mj_scalar_t>(coord);
	}

	//migrate weights.
	for (int i = 0; i < this->num_weights_per_coord; ++i){
		message_tag++;
		mj_scalar_t *weight = this->mj_weights[i];

		this->mj_weights[i] = allocMemory<mj_scalar_t>(num_incoming_gnos);
		ierr = Zoltan_Comm_Do(
				plan,
				message_tag,
				(char *) weight,
				sizeof(mj_scalar_t),
				(char *) this->mj_weights[i]);
		Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
		freeArray<mj_scalar_t>(weight);
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
	mj_part_t *new_parts = allocMemory<mj_part_t>(num_incoming_gnos);
	if(num_procs < num_parts){
		message_tag++;
		ierr = Zoltan_Comm_Do(
				plan,
				message_tag,
				(char *) this->assigned_part_ids,
				sizeof(mj_part_t),
				(char *) new_parts);
		Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
	}
	freeArray<mj_part_t>(this->assigned_part_ids);
	this->assigned_part_ids = new_parts;

	ierr = Zoltan_Comm_Destroy(&plan);
	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
	num_new_local_points = num_incoming_gnos;
	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Migration Z1Migration-" + iteration);
    }

    else

#endif  // ENABLE_ZOLTAN_MIGRATION
    {
	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Migration DistributorPlanCreating-" + iteration);
	Tpetra::Distributor distributor(this->comm);
	ArrayView<const mj_part_t> destinations( coordinate_destinations, this->num_local_coords);
	mj_lno_t num_incoming_gnos = distributor.createFromSends(destinations);
	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Migration DistributorPlanCreating-" + iteration);

	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Migration DistributorMigration-" + iteration);
	{
		//migrate gnos.
		ArrayRCP<mj_gno_t> received_gnos(num_incoming_gnos);
		ArrayView<mj_gno_t> sent_gnos(this->current_mj_gnos, this->num_local_coords);
		distributor.doPostsAndWaits<mj_gno_t>(sent_gnos, 1, received_gnos());
		freeArray<mj_gno_t>(this->current_mj_gnos);
		this->current_mj_gnos = allocMemory<mj_gno_t>(num_incoming_gnos);
		memcpy(
				this->current_mj_gnos,
				received_gnos.getRawPtr(),
				num_incoming_gnos * sizeof(mj_gno_t));
	}
	//migrate coordinates
	for (int i = 0; i < this->coord_dim; ++i){

		ArrayView<mj_scalar_t> sent_coord(this->mj_coordinates[i], this->num_local_coords);
		ArrayRCP<mj_scalar_t> received_coord(num_incoming_gnos);
		distributor.doPostsAndWaits<mj_scalar_t>(sent_coord, 1, received_coord());
		freeArray<mj_scalar_t>(this->mj_coordinates[i]);
		this->mj_coordinates[i] = allocMemory<mj_scalar_t>(num_incoming_gnos);
		memcpy(
				this->mj_coordinates[i],
				received_coord.getRawPtr(),
				num_incoming_gnos * sizeof(mj_scalar_t));
	}

	//migrate weights.
	for (int i = 0; i < this->num_weights_per_coord; ++i){

		ArrayView<mj_scalar_t> sent_weight(this->mj_weights[i], this->num_local_coords);
		ArrayRCP<mj_scalar_t> received_weight(num_incoming_gnos);
		distributor.doPostsAndWaits<mj_scalar_t>(sent_weight, 1, received_weight());
		freeArray<mj_scalar_t>(this->mj_weights[i]);
		this->mj_weights[i] = allocMemory<mj_scalar_t>(num_incoming_gnos);
		memcpy(
				this->mj_weights[i],
				received_weight.getRawPtr(),
				num_incoming_gnos * sizeof(mj_scalar_t));
	}

	{
		//migrate the owners of the coordinates
		ArrayView<int> sent_owners(this->owner_of_coordinate, this->num_local_coords);
		ArrayRCP<int> received_owners(num_incoming_gnos);
		distributor.doPostsAndWaits<int>(sent_owners, 1, received_owners());
		freeArray<int>(this->owner_of_coordinate);
		this->owner_of_coordinate = allocMemory<int>(num_incoming_gnos);
		memcpy(
						this->owner_of_coordinate,
						received_owners.getRawPtr(),
						num_incoming_gnos * sizeof(int));
	}

	//if num procs is less than num parts,
	//we need the part assigment arrays as well, since
	//there will be multiple parts in processor.
	if(num_procs < num_parts){
		ArrayView<mj_part_t> sent_partids(this->assigned_part_ids, this->num_local_coords);
		ArrayRCP<mj_part_t> received_partids(num_incoming_gnos);
		distributor.doPostsAndWaits<mj_part_t>(sent_partids, 1, received_partids());
		freeArray<mj_part_t>(this->assigned_part_ids);
		this->assigned_part_ids = allocMemory<mj_part_t>(num_incoming_gnos);
		memcpy(
				this->assigned_part_ids,
				received_partids.getRawPtr(),
				num_incoming_gnos * sizeof(mj_part_t));
	}
	else {
		mj_part_t *new_parts = allocMemory<int>(num_incoming_gnos);
		freeArray<mj_part_t>(this->assigned_part_ids);
		this->assigned_part_ids = new_parts;
	}
	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Migration DistributorMigration-" + iteration);
	num_new_local_points = num_incoming_gnos;

    }
}

/*! \brief Function creates the new subcomminicator for the processors
 * given in processor_ranks_for_subcomm.
 *
 * \param processor_ranks_for_subcomm is the vector that has the ranks of
 * the processors that will be in the same group.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::create_sub_communicator(std::vector<mj_part_t> &processor_ranks_for_subcomm){
    mj_part_t group_size = processor_ranks_for_subcomm.size();
    mj_part_t *ids = allocMemory<mj_part_t>(group_size);
    for(mj_part_t i = 0; i < group_size; ++i) {
        ids[i] = processor_ranks_for_subcomm[i];
    }
    ArrayView<const mj_part_t> idView(ids, group_size);
    this->comm = this->comm->createSubcommunicator(idView);
    freeArray<mj_part_t>(ids);
}


/*! \brief Function writes the new permutation arrays after the migration.
 *
 * \param output_num_parts is the number of parts that is assigned to the processor.
 * \param num_parts is the number of parts right before migration.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::fill_permutation_array(
    mj_part_t output_num_parts,
    mj_part_t num_parts){
	//if there is single output part, then simply fill the permutation array.
    if (output_num_parts == 1){
        for(mj_lno_t i = 0; i < this->num_local_coords; ++i){
            this->new_coordinate_permutations[i] = i;
        }
        this->new_part_xadj[0] = this->num_local_coords;
    }
    else {

    	//otherwise we need to count how many points are there in each part.
    	//we allocate here as num_parts, because the sent partids are up to num_parts,
    	//although there are outout_num_parts different part.
        mj_lno_t *num_points_in_parts = allocMemory<mj_lno_t>(num_parts);
        //part shift holds the which part number an old part number corresponds to.
        mj_part_t *part_shifts = allocMemory<mj_part_t>(num_parts);

        memset(num_points_in_parts, 0, sizeof(mj_lno_t) * num_parts);

        for(mj_lno_t i = 0; i < this->num_local_coords; ++i){
            mj_part_t ii = this->assigned_part_ids[i];
            ++num_points_in_parts[ii];
        }

        //write the end points of the parts.
        mj_part_t p = 0;
        mj_lno_t prev_index = 0;
        for(mj_part_t i = 0; i < num_parts; ++i){
            if(num_points_in_parts[i] > 0)  {
                this->new_part_xadj[p] =  prev_index + num_points_in_parts[i];
                prev_index += num_points_in_parts[i];
                part_shifts[i] = p++;
            }
        }

        //for the rest of the parts write the end index as end point.
        mj_part_t assigned_num_parts = p - 1;
        for (;p < num_parts; ++p){
            this->new_part_xadj[p] =  this->new_part_xadj[assigned_num_parts];
        }
        for(mj_part_t i = 0; i < output_num_parts; ++i){
            num_points_in_parts[i] = this->new_part_xadj[i];
        }

        //write the permutation array here.
        //get the part of the coordinate i, shift it to obtain the new part number.
        //assign it to the end of the new part numbers pointer.
        for(mj_lno_t i = this->num_local_coords - 1; i >= 0; --i){
            mj_part_t part = part_shifts[mj_part_t(this->assigned_part_ids[i])];
            this->new_coordinate_permutations[--num_points_in_parts[part]] = i;
        }

        freeArray<mj_lno_t>(num_points_in_parts);
        freeArray<mj_part_t>(part_shifts);
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
bool AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::mj_perform_migration(
    mj_part_t input_num_parts, //current umb parts
    mj_part_t &output_num_parts, //output umb parts.
    std::vector<mj_part_t> *next_future_num_parts_in_parts,
    mj_part_t &output_part_begin_index,
    size_t migration_reduce_all_population,
    mj_lno_t num_coords_for_last_dim_part,
    std::string iteration,
    RCP<mj_partBoxVector_t> &input_part_boxes,
    RCP<mj_partBoxVector_t> &output_part_boxes
)
{
	mj_part_t num_procs = this->comm->getSize();
	this->myRank = this->comm->getRank();


	//this array holds how many points each processor has in each part.
	//to access how many points processor i has on part j,
	//num_points_in_all_processor_parts[i * num_parts + j]
	mj_gno_t *num_points_in_all_processor_parts = allocMemory<mj_gno_t>(input_num_parts * (num_procs + 1));

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
		freeArray<mj_gno_t>(num_points_in_all_processor_parts);
		return false;
	}


	mj_lno_t *send_count_to_each_proc = NULL;
	int *coordinate_destinations = allocMemory<int>(this->num_local_coords);
	send_count_to_each_proc = allocMemory<mj_lno_t>(num_procs);
	for (int i = 0; i < num_procs; ++i) send_count_to_each_proc[i] = 0;

	std::vector<mj_part_t> processor_ranks_for_subcomm;
	std::vector<mj_part_t> out_part_indices;

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
			coordinate_destinations);




	freeArray<mj_lno_t>(send_count_to_each_proc);
	std::vector <mj_part_t> tmpv;

	std::sort (out_part_indices.begin(), out_part_indices.end());
	mj_part_t outP = out_part_indices.size();

	mj_gno_t new_global_num_points = 0;
	mj_gno_t *global_num_points_in_parts = num_points_in_all_processor_parts + num_procs * input_num_parts;

	if (this->mj_keep_part_boxes){
		input_part_boxes->clear();
	}

	//now we calculate the new values for next_future_num_parts_in_parts.
	//same for the part boxes.
	for (mj_part_t i = 0; i < outP; ++i){
		mj_part_t ind = out_part_indices[i];
		new_global_num_points += global_num_points_in_parts[ind];
		tmpv.push_back((*next_future_num_parts_in_parts)[ind]);
		if (this->mj_keep_part_boxes){
			input_part_boxes->push_back((*output_part_boxes)[ind]);
		}
	}
	//swap the input and output part boxes.
	if (this->mj_keep_part_boxes){
		RCP<mj_partBoxVector_t> tmpPartBoxes = input_part_boxes;
		input_part_boxes = output_part_boxes;
		output_part_boxes = tmpPartBoxes;
	}
	next_future_num_parts_in_parts->clear();
	for (mj_part_t i = 0; i < outP; ++i){
		mj_part_t p = tmpv[i];
		next_future_num_parts_in_parts->push_back(p);
	}

	freeArray<mj_gno_t>(num_points_in_all_processor_parts);

	mj_lno_t num_new_local_points = 0;


	//perform the actual migration operation here.
	this->mj_migrate_coords(
			num_procs,
			num_new_local_points,
			iteration,
			coordinate_destinations,
			input_num_parts);


	freeArray<int>(coordinate_destinations);

	if(this->num_local_coords != num_new_local_points){
		freeArray<mj_lno_t>(this->new_coordinate_permutations);
		freeArray<mj_lno_t>(this->coordinate_permutations);

		this->new_coordinate_permutations = allocMemory<mj_lno_t>(num_new_local_points);
		this->coordinate_permutations = allocMemory<mj_lno_t>(num_new_local_points);
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::create_consistent_chunks(
    mj_part_t num_parts,
    mj_scalar_t *mj_current_dim_coords,
    mj_scalar_t *current_concurrent_cut_coordinate,
    mj_lno_t coordinate_begin,
    mj_lno_t coordinate_end,
    mj_scalar_t *used_local_cut_line_weight_to_left,
    mj_lno_t *out_part_xadj,
    int coordInd){

	//mj_lno_t numCoordsInPart =  coordinateEnd - coordinateBegin;
	mj_part_t no_cuts = num_parts - 1;



	int me = 0;
	mj_lno_t *thread_num_points_in_parts = this->thread_point_counts[me];
	mj_scalar_t *my_local_thread_cut_weights_to_put_left = NULL;


	//now if the rectilinear partitioning is allowed we decide how
	//much weight each thread should put to left and right.
	if (this->distribute_points_on_cut_lines){

		my_local_thread_cut_weights_to_put_left = this->thread_cut_line_weight_to_put_left[me];
		for (mj_part_t i = 0; i < no_cuts; ++i){
			//the left to be put on the left of the cut.
			mj_scalar_t left_weight = used_local_cut_line_weight_to_left[i];
			//cout << "i:" << i << " left_weight:" << left_weight << endl;
			for(int ii = 0; ii < this->num_threads; ++ii){
				if(left_weight > this->sEpsilon){
					//the weight of thread ii on cut.
					mj_scalar_t thread_ii_weight_on_cut = this->thread_part_weight_work[ii][i * 2 + 1] - this->thread_part_weight_work[ii][i * 2 ];
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
			for (mj_part_t i = no_cuts - 1; i > 0 ; --i){
				if(ABS(current_concurrent_cut_coordinate[i] - current_concurrent_cut_coordinate[i -1]) < this->sEpsilon){
					my_local_thread_cut_weights_to_put_left[i] -= my_local_thread_cut_weights_to_put_left[i - 1] ;
				}
				my_local_thread_cut_weights_to_put_left[i] = int ((my_local_thread_cut_weights_to_put_left[i] + LEAST_SIGNIFICANCE) * SIGNIFICANCE_MUL)
                    										/ mj_scalar_t(SIGNIFICANCE_MUL);
			}
		}
	}

	for(mj_part_t ii = 0; ii < num_parts; ++ii){
		thread_num_points_in_parts[ii] = 0;
	}

	//for this specific case we dont want to distribute the points along the cut position
	//randomly, as we need a specific ordering of them. Instead,
	//we put the coordinates into a sort item, where we sort those
	//using the coordinates of points on other dimensions and the index.


	//some of the cuts might share the same position.
	//in this case, if cut i and cut j share the same position
	//cut_map[i] = cut_map[j] = sort item index.
	mj_part_t *cut_map = allocMemory<mj_part_t> (no_cuts);


	typedef uMultiSortItem<mj_lno_t, int, mj_scalar_t> multiSItem;
	typedef std::vector< multiSItem > multiSVector;
	typedef std::vector<multiSVector> multiS2Vector;

	//to keep track of the memory allocated.
	std::vector<mj_scalar_t *>allocated_memory;

	//vector for which the coordinates will be sorted.
	multiS2Vector sort_vector_points_on_cut;

	//the number of cuts that have different coordinates.
	mj_part_t different_cut_count = 1;
	cut_map[0] = 0;

	//now we insert 1 sort vector for all cuts on the different
	//positins.if multiple cuts are on the same position, they share sort vectors.
	multiSVector tmpMultiSVector;
	sort_vector_points_on_cut.push_back(tmpMultiSVector);

	for (mj_part_t i = 1; i < no_cuts ; ++i){
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
	for (mj_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){

		mj_lno_t i = this->coordinate_permutations[ii];

		mj_part_t pp = this->assigned_part_ids[i];
		mj_part_t p = pp / 2;
		//if the coordinate is on a cut.
		if(pp % 2 == 1 ){
			mj_scalar_t *vals = allocMemory<mj_scalar_t>(this->coord_dim -1);
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
			mj_part_t cmap = cut_map[p];
			sort_vector_points_on_cut[cmap].push_back(tempSortItem);
		}
		else {
			//if it is not on the cut, simple sorting.
			++thread_num_points_in_parts[p];
			this->assigned_part_ids[i] = p;
		}
	}

	//sort all the sort vectors.
	for (mj_part_t i = 0; i < different_cut_count; ++i){
		std::sort (sort_vector_points_on_cut[i].begin(), sort_vector_points_on_cut[i].end());
	}

	//we do the part assignment for the points on cuts here.
	mj_part_t previous_cut_map = cut_map[0];

	//this is how much previous part owns the weight of the current part.
	//when target part weight is 1.6, and the part on the left is given 2,
	//the left has an extra 0.4, while the right has missing 0.4 from the previous cut.
	//this parameter is used to balance this issues.
	//in the above example weight_stolen_from_previous_part will be 0.4.
	//if the left part target is 2.2 but it is given 2,
	//then weight_stolen_from_previous_part will be -0.2.
	mj_scalar_t weight_stolen_from_previous_part = 0;
	for (mj_part_t p = 0; p < no_cuts; ++p){

		mj_part_t mapped_cut = cut_map[p];

		//if previous cut map is done, and it does not have the same index,
		//then assign all points left on that cut to its right.
		if (previous_cut_map != mapped_cut){
			mj_lno_t sort_vector_end = (mj_lno_t)sort_vector_points_on_cut[previous_cut_map].size() - 1;
			for (; sort_vector_end >= 0; --sort_vector_end){
				multiSItem t = sort_vector_points_on_cut[previous_cut_map][sort_vector_end];
				mj_lno_t i = t.index;
				++thread_num_points_in_parts[p];
				this->assigned_part_ids[i] = p;
			}
			sort_vector_points_on_cut[previous_cut_map].clear();
		}
		mj_lno_t sort_vector_end = (mj_lno_t)sort_vector_points_on_cut[mapped_cut].size() - 1;

		for (; sort_vector_end >= 0; --sort_vector_end){
			multiSItem t = sort_vector_points_on_cut[mapped_cut][sort_vector_end];
			mj_lno_t i = t.index;
			mj_scalar_t w = this->mj_uniform_weights[0]? 1:this->mj_weights[0][i];


			//part p has enough space for point i, then put it to point i.
			if(	my_local_thread_cut_weights_to_put_left[p] + weight_stolen_from_previous_part> this->sEpsilon &&
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
	mj_lno_t sort_vector_end = (mj_lno_t)sort_vector_points_on_cut[previous_cut_map].size() - 1;
	for (; sort_vector_end >= 0; --sort_vector_end){
		multiSItem t = sort_vector_points_on_cut[previous_cut_map][sort_vector_end];
		mj_lno_t i = t.index;
		++thread_num_points_in_parts[no_cuts];
		this->assigned_part_ids[i] = no_cuts;
	}
	sort_vector_points_on_cut[previous_cut_map].clear();
	freeArray<mj_part_t> (cut_map);

	//free the memory allocated for vertex sort items .
	mj_lno_t vSize = (mj_lno_t) allocated_memory.size();
	for(mj_lno_t i = 0; i < vSize; ++i){
		freeArray<mj_scalar_t> (allocated_memory[i]);
	}

	//creation of part_xadj as in usual case.
	for(mj_part_t j = 0; j < num_parts; ++j){
		mj_lno_t num_points_in_part_j_upto_thread_i = 0;
		for (int i = 0; i < this->num_threads; ++i){
			mj_lno_t thread_num_points_in_part_j = this->thread_point_counts[i][j];
			this->thread_point_counts[i][j] = num_points_in_part_j_upto_thread_i;
			num_points_in_part_j_upto_thread_i += thread_num_points_in_part_j;

		}
		out_part_xadj[j] = num_points_in_part_j_upto_thread_i;// + prev2; //+ coordinateBegin;
	}

	//perform prefix sum for num_points in parts.
	for(mj_part_t j = 1; j < num_parts; ++j){
		out_part_xadj[j] += out_part_xadj[j - 1];
	}


	//shift the num points in threads thread to obtain the
	//beginning index of each thread's private space.
	for(mj_part_t j = 1; j < num_parts; ++j){
		thread_num_points_in_parts[j] += out_part_xadj[j - 1] ;
	}

	//now thread gets the coordinate and writes the index of coordinate to the permutation array
	//using the part index we calculated.
	for (mj_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii){
		mj_lno_t i = this->coordinate_permutations[ii];
		mj_part_t p =  this->assigned_part_ids[i];
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
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::set_final_parts(
		mj_part_t current_num_parts,
		mj_part_t output_part_begin_index,
		RCP<mj_partBoxVector_t> &output_part_boxes,
		bool is_data_ever_migrated)
{
    this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Part_Assignment");

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for(mj_part_t i = 0; i < current_num_parts;++i){

    	mj_lno_t begin = 0;
    	mj_lno_t end = this->part_xadj[i];

    	if(i > 0) begin = this->part_xadj[i -1];
    	mj_part_t part_to_set_index = i + output_part_begin_index;
    	if (this->mj_keep_part_boxes){
    		(*output_part_boxes)[i].setpId(part_to_set_index);
    	}
    	for (mj_lno_t ii = begin; ii < end; ++ii){
    		mj_lno_t k = this->coordinate_permutations[ii];
    		this->assigned_part_ids[k] = part_to_set_index;
    	}
    }

    //ArrayRCP<const mj_gno_t> gnoList;
    if(!is_data_ever_migrated){
    	//freeArray<mj_gno_t>(this->current_mj_gnos);
        //if(this->num_local_coords > 0){
        //    gnoList = arcpFromArrayView(this->mj_gnos);
        //}
    }
    else {
#ifdef ENABLE_ZOLTAN_MIGRATION
      if (sizeof(mj_lno_t) <=  sizeof(int)) {

        // Cannot use Zoltan_Comm with local ordinals larger than ints.
        // In Zoltan_Comm_Create, the cast int(this->num_local_coords) 
        // may overflow.

    	//if data is migrated, then send part numbers to the original owners.
    	ZOLTAN_COMM_OBJ *plan = NULL;
    	MPI_Comm mpi_comm = Teuchos2MPI (this->mj_problemComm);

    	int incoming = 0;
    	int message_tag = 7856;

    	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Final Z1PlanCreating");
    	int ierr = Zoltan_Comm_Create( &plan, int(this->num_local_coords),
    			this->owner_of_coordinate, mpi_comm, message_tag,
    			&incoming);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
    	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Final Z1PlanCreating" );

    	mj_gno_t *incoming_gnos = allocMemory< mj_gno_t>(incoming);

    	message_tag++;
    	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Final Z1PlanComm");
    	ierr = Zoltan_Comm_Do( plan, message_tag, (char *) this->current_mj_gnos,
    			sizeof(mj_gno_t), (char *) incoming_gnos);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

    	freeArray<mj_gno_t>(this->current_mj_gnos);
    	this->current_mj_gnos = incoming_gnos;

    	mj_part_t *incoming_partIds = allocMemory< mj_part_t>(incoming);

    	message_tag++;
    	ierr = Zoltan_Comm_Do( plan, message_tag, (char *) this->assigned_part_ids,
    			sizeof(mj_part_t), (char *) incoming_partIds);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
    	freeArray<mj_part_t>(this->assigned_part_ids);
    	this->assigned_part_ids = incoming_partIds;

    	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Final Z1PlanComm");
    	ierr = Zoltan_Comm_Destroy(&plan);
    	Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

    	this->num_local_coords = incoming;
    	//gnoList = arcp(this->current_mj_gnos, 0, this->num_local_coords, true);
      }
      else

#endif  // !ENABLE_ZOLTAN_MIGRATION
      {
    	//if data is migrated, then send part numbers to the original owners.
    	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Final DistributorPlanCreating");
    	Tpetra::Distributor distributor(this->mj_problemComm);
    	ArrayView<const mj_part_t> owners_of_coords(this->owner_of_coordinate, this->num_local_coords);
    	mj_lno_t incoming = distributor.createFromSends(owners_of_coords);
    	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Final DistributorPlanCreating" );

    	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Final DistributorPlanComm");
	//migrate gnos to actual owners.
	ArrayRCP<mj_gno_t> received_gnos(incoming);
	ArrayView<mj_gno_t> sent_gnos(this->current_mj_gnos, this->num_local_coords);
	distributor.doPostsAndWaits<mj_gno_t>(sent_gnos, 1, received_gnos());
	freeArray<mj_gno_t>(this->current_mj_gnos);
	this->current_mj_gnos = allocMemory<mj_gno_t>(incoming);
	memcpy( this->current_mj_gnos,
		received_gnos.getRawPtr(),
		incoming * sizeof(mj_gno_t));

		//migrate part ids to actual owners.
	ArrayView<mj_part_t> sent_partids(this->assigned_part_ids, this->num_local_coords);
	ArrayRCP<mj_part_t> received_partids(incoming);
	distributor.doPostsAndWaits<mj_part_t>(sent_partids, 1, received_partids());
	freeArray<mj_part_t>(this->assigned_part_ids);
	this->assigned_part_ids = allocMemory<mj_part_t>(incoming);
	memcpy( this->assigned_part_ids,
		received_partids.getRawPtr(),
		incoming * sizeof(mj_part_t));

    	this->num_local_coords = incoming;
    	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Final DistributorPlanComm");

      }
    }

    this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Part_Assignment");

    this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Solution_Part_Assignment");

    //ArrayRCP<mj_part_t> partId;
    //partId = arcp(this->assigned_part_ids, 0, this->num_local_coords, true);

    if (this->mj_keep_part_boxes){
    	this->kept_boxes = compute_global_box_boundaries(output_part_boxes);
                                                      
    }

    this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Solution_Part_Assignment");
}

/*! \brief Function frees all allocated work memory.
*/
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::free_work_memory(){
	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Problem_Free");

	for (int i=0; i < this->coord_dim; i++){
		freeArray<mj_scalar_t>(this->mj_coordinates[i]);
	}
	freeArray<mj_scalar_t *>(this->mj_coordinates);

	for (int i=0; i < this->num_weights_per_coord; i++){
		freeArray<mj_scalar_t>(this->mj_weights[i]);
	}
	freeArray<mj_scalar_t *>(this->mj_weights);

	freeArray<int>(this->owner_of_coordinate);

	for(int i = 0; i < this->num_threads; ++i){
		freeArray<mj_lno_t>(this->thread_point_counts[i]);
	}

	freeArray<mj_lno_t *>(this->thread_point_counts);
	freeArray<double *> (this->thread_part_weight_work);

	if(this->distribute_points_on_cut_lines){
		freeArray<mj_scalar_t>(this->process_cut_line_weight_to_put_left);
		for(int i = 0; i < this->num_threads; ++i){
			freeArray<mj_scalar_t>(this->thread_cut_line_weight_to_put_left[i]);
		}
		freeArray<mj_scalar_t *>(this->thread_cut_line_weight_to_put_left);
		freeArray<mj_scalar_t>(this->process_rectilinear_cut_weight);
		freeArray<mj_scalar_t>(this->global_rectilinear_cut_weight);
	}

	freeArray<mj_part_t>(this->my_incomplete_cut_count);

	freeArray<mj_scalar_t>(this->max_min_coords);

	freeArray<mj_lno_t>(this->new_part_xadj);

	freeArray<mj_lno_t>(this->coordinate_permutations);

	freeArray<mj_lno_t>(this->new_coordinate_permutations);

	freeArray<mj_scalar_t>(this->all_cut_coordinates);

	freeArray<mj_scalar_t> (this->process_local_min_max_coord_total_weight);

	freeArray<mj_scalar_t> (this->global_min_max_coord_total_weight);

	freeArray<mj_scalar_t>(this->cut_coordinates_work_array);

	freeArray<mj_scalar_t>(this->target_part_weights);

	freeArray<mj_scalar_t>(this->cut_upper_bound_coordinates);

	freeArray<mj_scalar_t>(this->cut_lower_bound_coordinates);

	freeArray<mj_scalar_t>(this->cut_lower_bound_weights);
	freeArray<mj_scalar_t>(this->cut_upper_bound_weights);
	freeArray<bool>(this->is_cut_line_determined);
	freeArray<mj_scalar_t>(this->total_part_weight_left_right_closests);
	freeArray<mj_scalar_t>(this->global_total_part_weight_left_right_closests);

	for(int i = 0; i < this->num_threads; ++i){
		freeArray<double>(this->thread_part_weights[i]);
		freeArray<mj_scalar_t>(this->thread_cut_right_closest_point[i]);
		freeArray<mj_scalar_t>(this->thread_cut_left_closest_point[i]);
	}

	freeArray<double *>(this->thread_part_weights);
	freeArray<mj_scalar_t *>(this->thread_cut_left_closest_point);
	freeArray<mj_scalar_t *>(this->thread_cut_right_closest_point);

	this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Problem_Free");
}

/*! \brief Multi Jagged  coordinate partitioning algorithm.
 *
 *  \param distribute_points_on_cut_lines_ :  if partitioning can distribute points on same coordinate to different parts.
 *  \param max_concurrent_part_calculation_ : how many parts we can calculate concurrently.
 *  \param check_migrate_avoid_migration_option_ : whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
 *  \param minimum_migration_imbalance_  : when MJ decides whether to migrate, the minimum imbalance for migration.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::set_partitioning_parameters(
		bool distribute_points_on_cut_lines_,
		int max_concurrent_part_calculation_,
		int check_migrate_avoid_migration_option_,
		mj_scalar_t minimum_migration_imbalance_){
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
 *  \param num_weights_per_coord: number of weights per coordinate
 *  \param mj_uniform_weights: if weight index [i] has uniform weight or not.
 *  \param mj_weights: the two dimensional array for weights
 *  \param mj_uniform_parts: if the target partitioning aims uniform parts
 *  \param mj_part_sizes: if the target partitioning does not aim uniform parts, then weight of each part.
 *
 *  \param result_assigned_part_ids: Output - 1D pointer, should be provided as null. Memory is given in the function.
 *  			the result partids corresponding to the coordinates given in result_mj_gnos.
 *  \param result_mj_gnos: Output - 1D pointer, should be provided as null. Memory is given in the function.
 *  			the result coordinate global id's corresponding to the part_ids array.
 *
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t>::multi_jagged_part(

	const RCP<const Environment> &env,
    	RCP<Comm<int> > &problemComm,

    	double imbalance_tolerance_,
    	size_t num_global_parts_,
    	mj_part_t *part_no_array_,
    	int recursion_depth_,

    	int coord_dim_,
    	mj_lno_t num_local_coords_,
    	mj_gno_t num_global_coords_,
    	const mj_gno_t *initial_mj_gnos_,
    	mj_scalar_t **mj_coordinates_,

    	int num_weights_per_coord_,
    	bool *mj_uniform_weights_,
    	mj_scalar_t **mj_weights_,
    	bool *mj_uniform_parts_,
    	mj_scalar_t **mj_part_sizes_,

    	mj_part_t *&result_assigned_part_ids_,
    	mj_gno_t *&result_mj_gnos_
)
{

#ifdef print_debug
    if(comm->getRank() == 0){
    	std::cout << "size of gno:" << sizeof(mj_gno_t) << std::endl;
    	std::cout << "size of lno:" << sizeof(mj_lno_t) << std::endl;
    	std::cout << "size of mj_scalar_t:" << sizeof(mj_scalar_t) << std::endl;
    }
#endif

    this->mj_env = env;
    this->mj_problemComm = problemComm;
    this->myActualRank = this->myRank = this->mj_problemComm->getRank();

    this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Total");
    this->mj_env->debug(3, "In MultiJagged Jagged");

    {
    	this->imbalance_tolerance = imbalance_tolerance_;
    	this->num_global_parts = num_global_parts_;
    	this->part_no_array =  part_no_array_;
    	this->recursion_depth = recursion_depth_;

    	this->coord_dim = coord_dim_;
    	this->num_local_coords = num_local_coords_;
    	this->num_global_coords = num_global_coords_;
    	this->mj_coordinates = mj_coordinates_; //will copy the memory to this->mj_coordinates.
    	this->initial_mj_gnos = (mj_gno_t *) initial_mj_gnos_; //will copy the memory to this->current_mj_gnos[j].

    	this->num_weights_per_coord = num_weights_per_coord_;
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
    mj_part_t current_num_parts = 1;
    mj_scalar_t *current_cut_coordinates =  this->all_cut_coordinates;

    this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Problem_Partitioning");

    mj_part_t output_part_begin_index = 0;
    mj_part_t future_num_parts = this->total_num_part;
    bool is_data_ever_migrated = false;

    std::vector<mj_part_t> *future_num_part_in_parts = new std::vector<mj_part_t> ();
    std::vector<mj_part_t> *next_future_num_parts_in_parts = new std::vector<mj_part_t> ();
    next_future_num_parts_in_parts->push_back(this->num_global_parts);

    RCP<mj_partBoxVector_t> input_part_boxes(new mj_partBoxVector_t(), true) ;
    RCP<mj_partBoxVector_t> output_part_boxes(new mj_partBoxVector_t(), true);

    compute_global_box();
    if(this->mj_keep_part_boxes){
    	this->init_part_boxes(output_part_boxes);
    }

    for (int i = 0; i < this->recursion_depth; ++i){
        //partitioning array. size will be as the number of current partitions and this
        //holds how many parts that each part will be in the current dimension partitioning.
        std::vector <mj_part_t> num_partitioning_in_current_dim;

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
        std::vector<mj_part_t> *tmpPartVect= future_num_part_in_parts;
        future_num_part_in_parts = next_future_num_parts_in_parts;
        next_future_num_parts_in_parts = tmpPartVect;

        //clear next_future_num_parts_in_parts array as
        //getPartitionArrays expects it to be empty.
        //it also expects num_partitioning_in_current_dim to be empty as well.
        next_future_num_parts_in_parts->clear();

        if(this->mj_keep_part_boxes){
            RCP<mj_partBoxVector_t> tmpPartBoxes = input_part_boxes;
            input_part_boxes = output_part_boxes;
            output_part_boxes = tmpPartBoxes;
            output_part_boxes->clear();
        }

        //returns the total no. of output parts for this dimension partitioning.
        mj_part_t output_part_count_in_dimension =
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
                RCP<mj_partBoxVector_t> tmpPartBoxes = input_part_boxes;
                input_part_boxes = output_part_boxes;
                output_part_boxes = tmpPartBoxes;
            }
            continue;
        }


        //get the coordinate axis along which the partitioning will be done.
        int coordInd = i % this->coord_dim;
        mj_scalar_t * mj_current_dim_coords = this->mj_coordinates[coordInd];

        //convert i to string to be used for debugging purposes.
        std::string istring = toString<int>(i);
        this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Problem_Partitioning_" + istring);

        //alloc Memory to point the indices
        //of the parts in the permutation array.
        this->new_part_xadj = allocMemory<mj_lno_t>(output_part_count_in_dimension);

        //the index where in the new_part_xadj will be written.
        mj_part_t output_part_index = 0;
        //whatever is written to output_part_index will be added with putput_coordinate_end_index
        //so that the points will be shifted.
        mj_part_t output_coordinate_end_index = 0;

        mj_part_t current_work_part = 0;
        mj_part_t current_concurrent_num_parts =
        		std::min(current_num_parts - current_work_part, this->max_concurrent_part_calculation);

        mj_part_t obtained_part_index = 0;

        //run for all available parts.
        for (; current_work_part < current_num_parts;
                 current_work_part += current_concurrent_num_parts){

            current_concurrent_num_parts = std::min(current_num_parts - current_work_part,
                                 this->max_concurrent_part_calculation);

            mj_part_t actual_work_part_count = 0;
            //initialization for 1D partitioning.
            //get the min and max coordinates of each part
            //together with the part weights of each part.
            for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                mj_part_t current_work_part_in_concurrent_parts = current_work_part + kk;

                //if this part wont be partitioned any further
                //dont do any work for this part.
                if (num_partitioning_in_current_dim[current_work_part_in_concurrent_parts] == 1){
                    continue;
                }
                ++actual_work_part_count;
                mj_lno_t coordinate_end_index= this->part_xadj[current_work_part_in_concurrent_parts];
                mj_lno_t coordinate_begin_index = current_work_part_in_concurrent_parts==0 ? 0: this->part_xadj[current_work_part_in_concurrent_parts -1];

/*
                cout << "i:" << i << " j:" << current_work_part + kk
                		<< " coordinate_begin_index:" << coordinate_begin_index
                		<< " coordinate_end_index:" << coordinate_end_index
                		<< " total:" << coordinate_end_index - coordinate_begin_index<< endl;
                		*/
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
                mj_part_t total_incomplete_cut_count = 0;

                //Compute weight ratios for parts & cuts:
                //e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
                //part0  cut0  part1 cut1 part2 cut2 part3
                mj_part_t concurrent_part_cut_shift = 0;
                mj_part_t concurrent_part_part_shift = 0;
                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    mj_scalar_t min_coordinate = this->global_min_max_coord_total_weight[kk];
                    mj_scalar_t max_coordinate = this->global_min_max_coord_total_weight[kk +
                                                     current_concurrent_num_parts];

                    mj_scalar_t global_total_weight = this->global_min_max_coord_total_weight[kk +
                                                        2 * current_concurrent_num_parts];

                    mj_part_t concurrent_current_part_index = current_work_part + kk;

                    mj_part_t partition_count = num_partitioning_in_current_dim[concurrent_current_part_index];

                    mj_scalar_t *usedCutCoordinate = current_cut_coordinates + concurrent_part_cut_shift;
                    mj_scalar_t *current_target_part_weights = this->target_part_weights +
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

                        mj_lno_t coordinate_end_index= this->part_xadj[concurrent_current_part_index];
                        mj_lno_t coordinate_begin_index = concurrent_current_part_index==0 ? 0: this->part_xadj[concurrent_current_part_index -1];

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
                mj_scalar_t used_imbalance = 0;


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
                mj_part_t output_array_shift = 0;
                mj_part_t cut_shift = 0;
                size_t tlr_shift = 0;
                size_t partweight_array_shift = 0;

                for(int kk = 0; kk < current_concurrent_num_parts; ++kk){
                    mj_part_t current_concurrent_work_part = current_work_part + kk;
                    mj_part_t num_parts = num_partitioning_in_current_dim[current_concurrent_work_part];

                    //if the part is empty, skip the part.
                    if((num_parts != 1  )
                    		&&
                    		this->global_min_max_coord_total_weight[kk] >
                             this->global_min_max_coord_total_weight[kk + current_concurrent_num_parts]) {

                    	//we still need to write the begin and end point of the
                    	//empty part. simply set it zero, the array indices will be shifted later.
                    	for(mj_part_t jj = 0; jj < num_parts; ++jj){
                    		this->new_part_xadj[output_part_index + output_array_shift + jj] = 0;
                    	}
                        cut_shift += num_parts - 1;
                        tlr_shift += (4 *(num_parts - 1) + 1);
                        output_array_shift += num_parts;
                        partweight_array_shift += (2 * (num_parts - 1) + 1);
                        continue;
                    }

                    mj_lno_t coordinate_end= this->part_xadj[current_concurrent_work_part];
                    mj_lno_t coordinate_begin = current_concurrent_work_part==0 ? 0: this->part_xadj[
                                                                current_concurrent_work_part -1];
                    mj_scalar_t *current_concurrent_cut_coordinate = current_cut_coordinates + cut_shift;
                    mj_scalar_t *used_local_cut_line_weight_to_left = this->process_cut_line_weight_to_put_left +
                                                            cut_shift;

                    //mj_scalar_t *used_tlr_array =  this->total_part_weight_left_right_closests + tlr_shift;

                    for(int ii = 0; ii < this->num_threads; ++ii){
                        this->thread_part_weight_work[ii] = this->thread_part_weights[ii] +  partweight_array_shift;
                    }

                    if(num_parts > 1){
                        if(this->mj_keep_part_boxes){
                        	//if part boxes are to be stored update the boundaries.
                            for (mj_part_t j = 0; j < num_parts - 1; ++j){
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
                        mj_lno_t part_size = coordinate_end - coordinate_begin;
                        *(this->new_part_xadj + output_part_index + output_array_shift) = part_size;
                        memcpy(
                        	this->new_coordinate_permutations + coordinate_begin,
                            this->coordinate_permutations + coordinate_begin,
                            part_size * sizeof(mj_lno_t));
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
                for(mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk){
                    mj_part_t num_parts = num_partitioning_in_current_dim[ current_work_part + kk];
                    for (mj_part_t ii = 0;ii < num_parts ; ++ii){
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

        	this->mj_env->timerStart(MACRO_TIMERS, "MultiJagged - Problem_Migration-" + istring);
        	mj_part_t num_parts = output_part_count_in_dimension;
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
        		this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Problem_Migration-" +
        				istring);
        		//since data is migrated, we reduce the number of reduceAll operations for the last part.
        		this->total_dim_num_reduce_all /= num_parts;
        	}
        	else {
        		is_migrated_in_current_dimension = false;
        		this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Problem_Migration-" + istring);
        	}
        }

        //swap the coordinate permutations for the next dimension.
        mj_lno_t * tmp = this->coordinate_permutations;
        this->coordinate_permutations = this->new_coordinate_permutations;
        this->new_coordinate_permutations = tmp;

        if(!is_migrated_in_current_dimension){
            this->total_dim_num_reduce_all -= current_num_parts;
            current_num_parts = output_part_count_in_dimension;
        }
        freeArray<mj_lno_t>(this->part_xadj);
        this->part_xadj = this->new_part_xadj;

        this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Problem_Partitioning_" + istring);
    }

    // Partitioning is done
    delete future_num_part_in_parts;
    delete next_future_num_parts_in_parts;

    this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Problem_Partitioning");
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
    this->mj_env->timerStop(MACRO_TIMERS, "MultiJagged - Total");
    this->mj_env->debug(3, "Out of MultiJagged");

}


/*! \brief Multi Jagged coordinate partitioning algorithm.
 *
 */
template <typename Adapter>
class Zoltan2_AlgMJ : public Algorithm<Adapter>
{
private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

    typedef CoordinateModel<typename Adapter::base_adapter_t> coordinateModel_t;
    typedef typename Adapter::scalar_t mj_scalar_t;
    typedef typename Adapter::gno_t mj_gno_t;
    typedef typename Adapter::lno_t mj_lno_t;
    typedef typename Adapter::node_t mj_node_t;
    typedef typename Adapter::part_t mj_part_t;
    typedef coordinateModelPartBox<mj_scalar_t, mj_part_t> mj_partBox_t;
    typedef std::vector<mj_partBox_t> mj_partBoxVector_t;
#endif
    AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t> mj_partitioner;

    RCP<const Environment> mj_env; //the environment object
    RCP<Comm<int> > mj_problemComm; //initial comm object
    RCP<const coordinateModel_t> mj_coords; //coordinate adapter

    //PARAMETERS
    double imbalance_tolerance; //input imbalance tolerance.
    size_t num_global_parts; //the targeted number of parts
    mj_part_t *part_no_array; //input part array specifying num part to divide along each dim.
    int recursion_depth; //the number of steps that partitioning will be solved in.

    int coord_dim; // coordinate dimension.
    mj_lno_t num_local_coords; //number of local coords.
    mj_gno_t num_global_coords; //number of global coords.
    const mj_gno_t *initial_mj_gnos; //initial global ids of the coordinates.
    mj_scalar_t **mj_coordinates; //two dimension coordinate array

    int num_weights_per_coord; // number of weights per coordinate
    bool *mj_uniform_weights; //if the coordinates have uniform weights.
    mj_scalar_t **mj_weights; //two dimensional weight array
    bool *mj_uniform_parts; //if the target parts are uniform
    mj_scalar_t **mj_part_sizes; //target part weight sizes.

    bool distribute_points_on_cut_lines; //if partitioning can distribute points on same coordiante to different parts.
    mj_part_t max_concurrent_part_calculation; // how many parts we can calculate concurrently.
    int check_migrate_avoid_migration_option; //whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
    mj_scalar_t minimum_migration_imbalance; //when MJ decides whether to migrate, the minimum imbalance for migration.
    int mj_keep_part_boxes; //if the boxes need to be kept.

    int num_threads;

    int mj_run_as_rcb; //if this is set, then recursion depth is adjusted to its maximum value.

    ArrayRCP<mj_part_t> comXAdj_; //communication graph xadj
    ArrayRCP<mj_part_t> comAdj_; //communication graph adj.

    void set_up_partitioning_data(
      const RCP<PartitioningSolution<Adapter> >&solution);

    void set_input_parameters(const Teuchos::ParameterList &p);

    void free_work_memory();

    RCP<mj_partBoxVector_t> getGlobalBoxBoundaries() const;

public:

    Zoltan2_AlgMJ(const RCP<const Environment> &env,
                  RCP<Comm<int> > &problemComm,
                  const RCP<const coordinateModel_t> &coords) :
                        mj_partitioner(), mj_env(env),
                        mj_problemComm(problemComm),
			mj_coords(coords),
                        imbalance_tolerance(0),
			num_global_parts(1), part_no_array(NULL),
                        recursion_depth(0),
			coord_dim(0),num_local_coords(0), num_global_coords(0),
			initial_mj_gnos(NULL), mj_coordinates(NULL),
                        num_weights_per_coord(0),
			mj_uniform_weights(NULL), mj_weights(NULL),
                        mj_uniform_parts(NULL),
			mj_part_sizes(NULL),
                        distribute_points_on_cut_lines(true),
			max_concurrent_part_calculation(1),
                        check_migrate_avoid_migration_option(0),
			minimum_migration_imbalance(0.30),
                        mj_keep_part_boxes(0), num_threads(1), mj_run_as_rcb(0),
                        comXAdj_(), comAdj_()
    {}
    ~Zoltan2_AlgMJ(){}

    /*! \brief Multi Jagged  coordinate partitioning algorithm.
     *
     *  \param solution  a PartitioningSolution, on input it
     *      contains part information, on return it also contains
     *      the solution and quality metrics.
     */
    void partition(const RCP<PartitioningSolution<Adapter> > &solution);

    mj_partBoxVector_t &getPartBoxesView() const
    {
      RCP<mj_partBoxVector_t> pBoxes = this->getGlobalBoxBoundaries();
      return *pBoxes;
    }

    mj_part_t pointAssign(int dim, mj_scalar_t *point) const;

    void boxAssign(int dim, mj_scalar_t *lower, mj_scalar_t *upper,
                   size_t &nPartsFound, mj_part_t **partsFound) const;


    /*! \brief returns communication graph resulting from MJ partitioning.
     */
    void getCommunicationGraph(
                         const PartitioningSolution<Adapter> *solution,
                         ArrayRCP<mj_part_t> &comXAdj,
                         ArrayRCP<mj_part_t> &comAdj);
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
void Zoltan2_AlgMJ<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
    this->set_up_partitioning_data(solution);
    this->set_input_parameters(this->mj_env->getParameters());
    if (this->mj_keep_part_boxes){
    	this->mj_partitioner.set_to_keep_part_boxes();
    }
    this->mj_partitioner.set_partitioning_parameters(
    		this->distribute_points_on_cut_lines,
    		this->max_concurrent_part_calculation,
    		this->check_migrate_avoid_migration_option,
    		this->minimum_migration_imbalance);

    mj_part_t *result_assigned_part_ids = NULL;
    mj_gno_t *result_mj_gnos = NULL;
    this->mj_partitioner.multi_jagged_part(
    		this->mj_env,
    		this->mj_problemComm,

    		this->imbalance_tolerance,
    		this->num_global_parts,
    		this->part_no_array,
    		this->recursion_depth,

    		this->coord_dim,
    		this->num_local_coords,
    		this->num_global_coords,
    		this->initial_mj_gnos,
    		this->mj_coordinates,

    		this->num_weights_per_coord,
    		this->mj_uniform_weights,
    		this->mj_weights,
    		this->mj_uniform_parts,
    		this->mj_part_sizes,

    		result_assigned_part_ids,
        	result_mj_gnos
    		);

    // Reorder results so that they match the order of the input
    std::unordered_map<mj_gno_t, mj_lno_t> localGidToLid;
    localGidToLid.reserve(this->num_local_coords);
    for (mj_lno_t i = 0; i < this->num_local_coords; i++)
      localGidToLid[this->initial_mj_gnos[i]] = i;

    ArrayRCP<mj_part_t> partId = arcp(new mj_part_t[this->num_local_coords],
                                      0, this->num_local_coords, true);
 
    for (mj_lno_t i = 0; i < this->num_local_coords; i++) {
      mj_lno_t origLID = localGidToLid[result_mj_gnos[i]];
      partId[origLID] = result_assigned_part_ids[i];
    }
    
    delete [] result_mj_gnos;
    delete [] result_assigned_part_ids;

    solution->setParts(partId);
    this->free_work_memory();
}

/* \brief Freeing the memory allocated.
 * */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::free_work_memory(){
	freeArray<mj_scalar_t *>(this->mj_coordinates);
	freeArray<mj_scalar_t *>(this->mj_weights);
	freeArray<bool>(this->mj_uniform_parts);
	freeArray<mj_scalar_t *>(this->mj_part_sizes);
	freeArray<bool>(this->mj_uniform_weights);

}

/* \brief Sets the partitioning data for multijagged algorithm.
 * */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::set_up_partitioning_data(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
	this->coord_dim = this->mj_coords->getCoordinateDim();
	this->num_weights_per_coord = this->mj_coords->getNumWeightsPerCoordinate();
	this->num_local_coords = this->mj_coords->getLocalNumCoordinates();
	this->num_global_coords = this->mj_coords->getGlobalNumCoordinates();
	int criteria_dim = (this->num_weights_per_coord ? this->num_weights_per_coord : 1);

	// From the Solution we get part information.
	// If the part sizes for a given criteria are not uniform,
	// then they are values that sum to 1.0.
	this->num_global_parts = solution->getTargetGlobalNumberOfParts();
	//allocate only two dimensional pointer.
	//raw pointer addresess will be obtained from multivector.
	this->mj_coordinates = allocMemory<mj_scalar_t *>(this->coord_dim);
	this->mj_weights = allocMemory<mj_scalar_t *>(criteria_dim);

	//if the partitioning results are to be uniform.
	this->mj_uniform_parts = allocMemory< bool >(criteria_dim);
	//if in a criteria dimension, uniform part is false this shows ratios of
	//the target part weights.
	this->mj_part_sizes =  allocMemory<mj_scalar_t *>(criteria_dim);
	//if the weights of coordinates are uniform in a criteria dimension.
	this->mj_uniform_weights = allocMemory< bool >(criteria_dim);

	typedef StridedData<mj_lno_t, mj_scalar_t> input_t;
	ArrayView<const mj_gno_t> gnos;
	ArrayView<input_t> xyz;
	ArrayView<input_t> wgts;

	this->mj_coords->getCoordinates(gnos, xyz, wgts);
	//obtain global ids.
	ArrayView<const mj_gno_t> mj_gnos = gnos;
	this->initial_mj_gnos = mj_gnos.getRawPtr();

	//extract coordinates from multivector.
	for (int dim=0; dim < this->coord_dim; dim++){
		ArrayRCP<const mj_scalar_t> ar;
		xyz[dim].getInputArray(ar);
		//multiJagged coordinate values assignment
		this->mj_coordinates[dim] =  (mj_scalar_t *)ar.getRawPtr();
	}

	//if no weights are provided set uniform weight.
	if (this->num_weights_per_coord == 0){
		this->mj_uniform_weights[0] = true;
		this->mj_weights[0] = NULL;
	}
	else{
		//if weights are provided get weights for all weight indices
		for (int wdim = 0; wdim < this->num_weights_per_coord; wdim++){
			ArrayRCP<const mj_scalar_t> ar;
			wgts[wdim].getInputArray(ar);
			this->mj_uniform_weights[wdim] = false;
			this->mj_weights[wdim] = (mj_scalar_t *) ar.getRawPtr();
		}
	}

	for (int wdim = 0; wdim < criteria_dim; wdim++){
		if (solution->criteriaHasUniformPartSizes(wdim)){
			this->mj_uniform_parts[wdim] = true;
			this->mj_part_sizes[wdim] = NULL;
		}
		else{
			std::cerr << "MJ does not support non uniform target part weights" << std::endl;
			exit(1);
		}
	}
}

/* \brief Sets the partitioning parameters for multijagged algorithm.
 * \param pl: is the parameter list provided to zoltan2 call
 * */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::set_input_parameters(const Teuchos::ParameterList &pl){

	const Teuchos::ParameterEntry *pe = pl.getEntryPtr("imbalance_tolerance");
	if (pe){
		double tol;
		tol = pe->getValue(&tol);
		this->imbalance_tolerance = tol - 1.0;
	}

    // TODO: May be a more relaxed tolerance is needed. RCB uses 10%
	if (this->imbalance_tolerance <= 0)
		this->imbalance_tolerance= 10e-4;

	//if an input partitioning array is provided.
	this->part_no_array = NULL;
	//the length of the input partitioning array.
	this->recursion_depth = 0;

	if (pl.getPtr<Array <mj_part_t> >("mj_parts")){
		this->part_no_array = (mj_part_t *) pl.getPtr<Array <mj_part_t> >("mj_parts")->getRawPtr();
		this->recursion_depth = pl.getPtr<Array <mj_part_t> >("mj_parts")->size() - 1;
		this->mj_env->debug(2, "mj_parts provided by user");
	}

	//get mj specific parameters.
	this->distribute_points_on_cut_lines = true;
	this->max_concurrent_part_calculation = 1;

	this->mj_run_as_rcb = 0;
	int mj_user_recursion_depth = -1;
	this->mj_keep_part_boxes = 0;
	this->check_migrate_avoid_migration_option = 0;
	this->minimum_migration_imbalance = 0.35;

	pe = pl.getEntryPtr("mj_minimum_migration_imbalance");
	if (pe){
		double imb;
		imb = pe->getValue(&imb);
		this->minimum_migration_imbalance = imb - 1.0;
	}

	pe = pl.getEntryPtr("mj_migration_option");
	if (pe){
		this->check_migrate_avoid_migration_option = pe->getValue(&this->check_migrate_avoid_migration_option);
	}else {
		this->check_migrate_avoid_migration_option = 0;
	}
	if (this->check_migrate_avoid_migration_option > 1) this->check_migrate_avoid_migration_option = -1;


	pe = pl.getEntryPtr("mj_concurrent_part_count");
	if (pe){
		this->max_concurrent_part_calculation = pe->getValue(&this->max_concurrent_part_calculation);
	}else {
		this->max_concurrent_part_calculation = 1; // Set to 1 if not provided.
	}

	pe = pl.getEntryPtr("mj_keep_part_boxes");
	if (pe){
		this->mj_keep_part_boxes = pe->getValue(&this->mj_keep_part_boxes);
	}else {
		this->mj_keep_part_boxes = 0; // Set to invalid value
	}

        // For now, need keep_part_boxes to do pointAssign and boxAssign.
	// pe = pl.getEntryPtr("keep_cuts");
	// if (pe){
	// 	int tmp = pe->getValue(&tmp);
	// 	if (tmp) this->mj_keep_part_boxes = 1;
        // }

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

	pe = pl.getEntryPtr("mj_recursion_depth");
	if (pe){
		mj_user_recursion_depth = pe->getValue(&mj_user_recursion_depth);
	}else {
		mj_user_recursion_depth = -1; // Set to invalid value
	}

	int val = 0;
	pe = pl.getEntryPtr("rectilinear");
	if (pe) val = pe->getValue(&val);
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

/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::boxAssign(
  int dim, 
  typename Adapter::scalar_t *lower, 
  typename Adapter::scalar_t *upper,
  size_t &nPartsFound, 
  typename Adapter::part_t **partsFound) const
{
  // TODO:  Implement with cuts rather than boxes to reduce algorithmic
  // TODO:  complexity.  Or at least do a search through the boxes, using
  // TODO:  p x q x r x ... if possible.

  nPartsFound = 0;
  *partsFound = NULL;

  if (this->mj_keep_part_boxes) {

    // Get vector of part boxes
    RCP<mj_partBoxVector_t> partBoxes = this->getGlobalBoxBoundaries();

    size_t nBoxes = (*partBoxes).size();
    if (nBoxes == 0) {
      throw std::logic_error("no part boxes exist");
    }

    // Determine whether the box overlaps the globalBox at all
    RCP<mj_partBox_t> globalBox = this->mj_partitioner.get_global_box();

    if (globalBox->boxesOverlap(dim, lower, upper)) {

      std::vector<typename Adapter::part_t> partlist;

      // box overlaps the global box; find specific overlapping boxes
      for (size_t i = 0; i < nBoxes; i++) {
        try {
          if ((*partBoxes)[i].boxesOverlap(dim, lower, upper)) {
            nPartsFound++;
            partlist.push_back((*partBoxes)[i].getpId());

//            std::cout << "Given box (";
//            for (int j = 0; j < dim; j++)
//              std::cout << lower[j] << " ";
//            std::cout << ") x (";
//            for (int j = 0; j < dim; j++)
//              std::cout << upper[j] << " ";
//            std::cout << ") overlaps PartBox " 
//                      << (*partBoxes)[i].getpId() << " (";
//            for (int j = 0; j < dim; j++)
//              std::cout << (*partBoxes)[i].getlmins()[j] << " ";
//            std::cout << ") x (";
//            for (int j = 0; j < dim; j++)
//              std::cout << (*partBoxes)[i].getlmaxs()[j] << " ";
//            std::cout << ")" << std::endl;
          }
        }
        Z2_FORWARD_EXCEPTIONS;
      }
      if (nPartsFound) {
        *partsFound = new mj_part_t[nPartsFound];
        for (size_t i = 0; i < nPartsFound; i++)
          (*partsFound)[i] = partlist[i];
      }
    }
    else {
      // Box does not overlap the domain at all.  Find the closest part
      // Not sure how to perform this operation for MJ without having the 
      // cuts.  With the RCB cuts, the concept of a part extending to 
      // infinity was natural.  With the boxes, it is much more difficult.
      // TODO:  For now, return information indicating NO OVERLAP.

    }
  }
  else {
    throw std::logic_error("need to use keep_cuts parameter for boxAssign");
  }
}

/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
typename Adapter::part_t Zoltan2_AlgMJ<Adapter>::pointAssign(
  int dim, 
  typename Adapter::scalar_t *point) const
{

  // TODO:  Implement with cuts rather than boxes to reduce algorithmic
  // TODO:  complexity.  Or at least do a search through the boxes, using
  // TODO:  p x q x r x ... if possible.

  if (this->mj_keep_part_boxes) {
    typename Adapter::part_t foundPart = -1;

    // Get vector of part boxes
    RCP<mj_partBoxVector_t> partBoxes = this->getGlobalBoxBoundaries();

    size_t nBoxes = (*partBoxes).size();
    if (nBoxes == 0) {
      throw std::logic_error("no part boxes exist");
    }

    // Determine whether the point is within the global domain
    RCP<mj_partBox_t> globalBox = this->mj_partitioner.get_global_box();

    if (globalBox->pointInBox(dim, point)) {

      // point is in the global domain; determine in which part it is.
      size_t i;
      for (i = 0; i < nBoxes; i++) {
        try {
          if ((*partBoxes)[i].pointInBox(dim, point)) {
            foundPart = (*partBoxes)[i].getpId();
//            std::cout << "Point (";
//            for (int j = 0; j < dim; j++) std::cout << point[j] << " ";
//            std::cout << ") found in box " << i << " part " << foundPart 
//                      << std::endl;
//            (*partBoxes)[i].print();
            break;
          }
        }
        Z2_FORWARD_EXCEPTIONS;
      }

      if (i == nBoxes) {
        // This error should never occur
        std::ostringstream oss;
        oss << "Point (";
        for (int j = 0; j < dim; j++) oss << point[j] << " ";
        oss << ") not found in domain";
        throw std::logic_error(oss.str());
      }
    }
    
    else {
      // Point is outside the global domain.  
      // Determine to which part it is closest.
      // TODO:  with cuts, would not need this special case

      size_t closestBox = 0;
      mj_scalar_t minDistance = std::numeric_limits<mj_scalar_t>::max();
      mj_scalar_t *centroid = new mj_scalar_t[dim];
      for (size_t i = 0; i < nBoxes; i++) {
        (*partBoxes)[i].computeCentroid(centroid);
        mj_scalar_t sum = 0.;
        mj_scalar_t diff;
        for (int j = 0; j < dim; j++) {
          diff = centroid[j] - point[j];
          sum += diff * diff;
        }
        if (sum < minDistance) {
          minDistance = sum;
          closestBox = i;
        }
      }
      foundPart = (*partBoxes)[closestBox].getpId();
      delete [] centroid;
    }

    return foundPart;
  }
  else {
    throw std::logic_error("need to use keep_cuts parameter for pointAssign");
  }
}

template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::getCommunicationGraph(
  const PartitioningSolution<Adapter> *solution,
  ArrayRCP<typename Zoltan2_AlgMJ<Adapter>::mj_part_t> &comXAdj,
  ArrayRCP<typename Zoltan2_AlgMJ<Adapter>::mj_part_t> &comAdj) 
{
  if(comXAdj_.getRawPtr() == NULL && comAdj_.getRawPtr() == NULL){
    RCP<mj_partBoxVector_t> pBoxes = this->getGlobalBoxBoundaries();
    mj_part_t ntasks =  (*pBoxes).size();
    int dim = (*pBoxes)[0].getDim();
    GridHash<mj_scalar_t, mj_part_t> grid(pBoxes, ntasks, dim);
    grid.getAdjArrays(comXAdj_, comAdj_);
  }
  comAdj = comAdj_;
  comXAdj = comXAdj_;
}


template <typename Adapter>
RCP<typename Zoltan2_AlgMJ<Adapter>::mj_partBoxVector_t>
Zoltan2_AlgMJ<Adapter>::getGlobalBoxBoundaries() const
{
  return this->mj_partitioner.get_kept_boxes();
}


template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
RCP<typename AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t>::mj_partBoxVector_t>
AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t>::get_kept_boxes() const
{
  if (this->mj_keep_part_boxes)
    return this->kept_boxes;
  else
    throw std::logic_error("Error: part boxes are not stored.");
}

template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t>
RCP<typename AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t>::mj_partBoxVector_t>
AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t>::compute_global_box_boundaries(
  RCP<mj_partBoxVector_t> &localPartBoxes
) const
{
  mj_part_t ntasks = this->num_global_parts;
  int dim = (*localPartBoxes)[0].getDim();
  mj_scalar_t *localPartBoundaries = new mj_scalar_t[ntasks * 2 *dim];

  memset(localPartBoundaries, 0, sizeof(mj_scalar_t) * ntasks * 2 *dim);

  mj_scalar_t *globalPartBoundaries = new mj_scalar_t[ntasks * 2 *dim];
  memset(globalPartBoundaries, 0, sizeof(mj_scalar_t) * ntasks * 2 *dim);

  mj_scalar_t *localPartMins = localPartBoundaries;
  mj_scalar_t *localPartMaxs = localPartBoundaries + ntasks * dim;

  mj_scalar_t *globalPartMins = globalPartBoundaries;
  mj_scalar_t *globalPartMaxs = globalPartBoundaries + ntasks * dim;

  mj_part_t boxCount = localPartBoxes->size();
  for (mj_part_t i = 0; i < boxCount; ++i){
    mj_part_t pId = (*localPartBoxes)[i].getpId();
      //cout << "me:" << comm->getRank() << " has:" << pId << endl;

    mj_scalar_t *lmins = (*localPartBoxes)[i].getlmins();
    mj_scalar_t *lmaxs = (*localPartBoxes)[i].getlmaxs();

    for (int j = 0; j < dim; ++j){
      localPartMins[dim * pId + j] = lmins[j];
      localPartMaxs[dim * pId + j] = lmaxs[j];
      /*
      cout << "me:" << comm->getRank()  <<
              " dim * pId + j:"<< dim * pId + j <<
              " localMin:" << localPartMins[dim * pId + j] <<
              " localMax:" << localPartMaxs[dim * pId + j] << endl;
      */
    }
  }

  Teuchos::Zoltan2_BoxBoundaries<int, mj_scalar_t> reductionOp(ntasks * 2 *dim);

  reduceAll<int, mj_scalar_t>(*mj_problemComm, reductionOp,
            ntasks * 2 *dim, localPartBoundaries, globalPartBoundaries);
  RCP<mj_partBoxVector_t> pB(new mj_partBoxVector_t(),true);
  for (mj_part_t i = 0; i < ntasks; ++i){
    Zoltan2::coordinateModelPartBox <mj_scalar_t, mj_part_t> tpb(i, dim,
                                               globalPartMins + dim * i,
                                               globalPartMaxs + dim * i);

    /*
    for (int j = 0; j < dim; ++j){
        cout << "me:" << comm->getRank()  <<
                " dim * pId + j:"<< dim * i + j <<
                " globalMin:" << globalPartMins[dim * i + j] <<
                " globalMax:" << globalPartMaxs[dim * i + j] << endl;
    }
    */
    pB->push_back(tpb);
  }
  delete []localPartBoundaries;
  delete []globalPartBoundaries;
  //RCP <mj_partBoxVector_t> tmpRCPBox(pB, true);
  return pB;
}
} // namespace Zoltan2

#endif
