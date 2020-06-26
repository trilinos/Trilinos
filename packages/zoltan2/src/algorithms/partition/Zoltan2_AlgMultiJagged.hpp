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
#include <Zoltan2_IntegerRangeList.hpp>
#include <Zoltan2_CoordinatePartitioningGraph.hpp>
#include <Zoltan2_Util.hpp>
#include <Tpetra_Distributor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Kokkos_Sort.hpp>

#include <algorithm>    // std::sort
#include <vector>
#include <unordered_map>

#ifdef ZOLTAN2_USEZOLTANCOMM
#ifdef HAVE_ZOLTAN2_MPI
#define ZOLTAN2_MJ_ENABLE_ZOLTAN_MIGRATION
#include "zoltan_comm_cpp.h"
#include "zoltan_types.h" // for error codes
#endif
#endif

namespace Teuchos{

/*! \brief Zoltan2_BoxBoundaries is a reduction operation
 *  to all reduce the all box boundaries.
*/
template <typename Ordinal, typename T>
class Zoltan2_BoxBoundaries  : public ValueTypeReductionOp<Ordinal,T>
{
private:
  Ordinal size;
  T epsilon;

public:
  /*! \brief Default Constructor
   */
  Zoltan2_BoxBoundaries() : size(0),
    epsilon(std::numeric_limits<T>::epsilon()) {}

  /*! \brief Constructor
   *  \param Ordinal          DOCWORK: Documentation
   */
  Zoltan2_BoxBoundaries(Ordinal s_):
    size(s_), epsilon(std::numeric_limits<T>::epsilon()) {}

  /*! \brief Implement Teuchos::ValueTypeReductionOp interface
   *  \param count            DOCWORK: Documentation
   *  \param inBuffer         DOCWORK: Documentation
   *  \param inoutBuffer      DOCWORK: Documentation
   */
  void reduce( const Ordinal count, const T inBuffer[], T inoutBuffer[]) const {
    for(Ordinal i = 0; i < count; i++) {
      if(Z2_ABS(inBuffer[i]) > epsilon) {
        inoutBuffer[i] = inBuffer[i];
      }
    }
  }
};

} // namespace Teuchos

namespace Zoltan2{

/*! \brief Class for sorting items with multiple values.
 *  First sorting with respect to val[0], then val[1] then ... val[count-1].
 *  The last tie breaking is done with index values.
 *  Used for task mapping partitioning where the points on a cut line needs to
 *  be distributed consistently.
 */
template <typename IT, typename CT, typename WT>
class uMultiSortItem
{
public:
  // TODO: Why volatile?
  // no idea, another intel compiler failure.
  volatile IT index;
  volatile CT count;
  volatile WT *val;
  volatile WT epsilon;

  uMultiSortItem() {
    this->index = 0;
    this->count = 0;
    this->val = NULL;
    this->epsilon = std::numeric_limits<WT>::epsilon() * 100;
  }

  // TODO: Document these methods?
  uMultiSortItem(IT index_ ,CT count_, WT *vals_) {
    this->index = index_;
    this->count = count_;
    this->val = vals_;
    this->epsilon = std::numeric_limits<WT>::epsilon() * 100;
  }

  ~uMultiSortItem() {
  }

  void set(IT index_ ,CT count_, WT *vals_) {
    this->index = index_;
    this->count = count_;
    this->val = vals_;
  }

  bool operator<(const uMultiSortItem<IT,CT,WT>& other) const {
    assert(this->count == other.count);
    for(CT i = 0; i < this->count; ++i) {
      // if the values are equal go to next one.
      if(std::abs(this->val[i] - other.val[i]) < this->epsilon) {
        continue;
      }
      // if next value is smaller return true;
      if(this->val[i] < other.val[i]) {
        return true;
      }
      // if next value is bigger return false;
      else {
        return false;
      }
    }
    // if they are totally equal.
    return this->index < other.index;
  }
};

/*! \brief Sort items for quick sort function.
 */
template <class IT, class WT>
struct uSortItem
{
  IT id;
  WT val;
};

/*! \brief Quick sort function.
 *  Sorts the arr of uSortItems, with respect to increasing vals.
 *  DOCWORK: Document input params
 */
template <class IT, class WT>
void uqsort(IT n, uSortItem<IT, WT> * arr) {
  int NSTACK = 50;
  int M = 7;
  IT         i, ir=n, j, k, l=1;
  IT         jstack=0, istack[50];
  WT aval;
  uSortItem<IT,WT>    a;

  --arr;
  for(;;) {
    if(ir-l < M) {
      for(j=l+1;j<=ir;j++) {
        a=arr[j];
        aval = a.val;
        for(i=j-1;i>=1;i--) {
          if(arr[i].val <= aval)
            break;
          arr[i+1] = arr[i];
        }
        arr[i+1]=a;
      }
      if(jstack == 0)
          break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else {
      k=(l+ir) >> 1;
      std::swap(arr[k],arr[l+1]);
      if(arr[l+1].val > arr[ir].val) {
        std::swap(arr[l+1],arr[ir]);
      }
      if(arr[l].val > arr[ir].val) {
        std::swap(arr[l],arr[ir]);
      }
      if(arr[l+1].val > arr[l].val) {
        std::swap(arr[l+1],arr[l]);
      }
      i=l+1;
      j=ir;
      a=arr[l];
      aval = a.val;
      for(;;) {
        do i++; while (arr[i].val < aval);
        do j--; while (arr[j].val > aval);
        if(j < i) break;
        std::swap(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      jstack += 2;
      if(jstack > NSTACK) {
        std::cout << "uqsort: NSTACK too small in sort." << std::endl;
        std::terminate();
      }
      if(ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      }
      else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

template <class IT, class WT, class SIGN>
struct uSignedSortItem
{
  IT id;
  WT val;
  SIGN signbit; // 1 means positive, 0 means negative.
  bool operator<(const uSignedSortItem<IT, WT, SIGN>& rhs) const {
    /*if I am negative, the other is positive*/
    if(this->signbit < rhs.signbit) {
      return true;
    }
    /*if both has the same sign*/
    else if(this->signbit == rhs.signbit) {
      if(this->val < rhs.val) {//if my value is smaller,
        return this->signbit;//then if we both are positive return true.
                            //if we both are negative, return false.
      }
      else if(this->val > rhs.val) {//if my value is larger,
        return !this->signbit; //then if we both are positive return false.
                              //if we both are negative, return true.
      }
      else { //if both are equal.
        return false;
      }
    }
    else {
      /*if I am positive, the other is negative*/
      return false;
    }
  }

  bool operator<=(const uSignedSortItem<IT, WT, SIGN>& rhs) {
    return (this->val == rhs.val && this->signbit == rhs.signbit) || (*this < rhs);
  }
};

/*! \brief Quick sort function.
 *  Sorts the arr of uSignedSortItems, with respect to increasing vals.
 */
template <class IT, class WT, class SIGN>
void uqSignsort(IT n, uSignedSortItem<IT, WT, SIGN> * arr) {
  IT NSTACK = 50;
  IT M = 7;
  IT         i, ir=n, j, k, l=1;
  IT         jstack=0, istack[50];
  uSignedSortItem<IT,WT,SIGN>    a;

  --arr;
  for(;;) {
    if(ir < M + l) {
      for(j=l+1;j<=ir;j++) {
        a=arr[j];
        for(i=j-1;i>=1;i--) {
          if(arr[i] <= a) {
              break;
          }
          arr[i+1] = arr[i];
        }
        arr[i+1]=a;
      }
      if(jstack == 0) {
        break;
      }
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else {
      k=(l+ir) >> 1;
      std::swap(arr[k],arr[l+1]);
      if(arr[ir] < arr[l+1]) {
        std::swap(arr[l+1],arr[ir]);
      }
      if(arr[ir] < arr[l] ) {
        std::swap(arr[l],arr[ir]);
      }
      if(arr[l] < arr[l+1]) {
        std::swap(arr[l+1],arr[l]);
      }
      i=l+1;
      j=ir;
      a=arr[l];
      for(;;) {
        do i++; while (arr[i] < a);
        do j--; while (a < arr[j]);
        if(j < i) break;
        std::swap(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      jstack += 2;
      if(jstack > NSTACK) {
        std::cout << "uqsort: NSTACK too small in sort." << std::endl;
        std::terminate();
      }
      if(ir+l+1 >= j+i) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      }
      else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

// This exists only so we can track how many times the MJ algorithm is
// called and put each of those into different timer names.
// Currently the MultiJaggedTest.cpp will actually call it twice.
// First time with data from a Tpetra MultiVector and then a second time using
// a BasicVectorAdapter which allows us to turn UVM off for some tests. The
// results of the two runs are compared which helps to catch a lot of bugs. For
// profiling I'm mostly just interested in the UVM off case and need it to be
// in separate timers. Passing a value through would mess up the API. Possibly
// we could check the Adapter and use that. The statics have to be outside the
// templated class as the two called instances will be different template
// parameters. Another complication is that MultiJagged.cpp will call through
// the Zoltan2_AlgMJ class and we want to time things in both classes. However
// TaskMapper will directly call AlgMJ so I made two counters for the two
// classes to make sure it was always correct. This does not impact any
// behavior and has the sole purpose of generating unique timer names. If you
// run an MJ test you'll see MJ(0) and MJ(1) in the names to distinguish the
// 1st and 2nd run. Right now only MultijaggedTest.cpp cares about this.
struct Zoltan2_AlgMJ_TrackCallsCounter {
  static int get_counter_AlgMJ() {
    static int counter = 0;
    return counter++;
  }
  static int get_counter_Zoltan2_AlgMJ() {
    static int counter = 0;
    return counter++;
  }
};

/*! \brief Multi Jagged coordinate partitioning algorithm.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
class AlgMJ
{
private:
  typedef typename mj_node_t::device_type device_t; // for views
  typedef coordinateModelPartBox mj_partBox_t;
  typedef std::vector<mj_partBox_t> mj_partBoxVector_t;

  //if the (last dimension reduce all count) x the mpi world size
  //estimated to be bigger than this number then migration will be forced
  //in earlier iterations.
  static constexpr size_t future_reduceall_cutoff = 1500000;

  //if parts right before last dimension are estimated to have less than
  //MIN_WORK_LAST_DIM many coords, migration will be forced in earlier iterations.
  static constexpr mj_lno_t min_work_last_dim = 1000;

  static constexpr mj_scalar_t least_signifiance = 0.0001;
  static constexpr int significance_mul = 1000;

  std::string mj_timer_base_string; // for convenience making timer names

  RCP<const Environment> mj_env; // the environment object
  RCP<const Comm<int> > mj_problemComm; // initial comm object
  RCP<Comm<int> > comm; // comm object than can be altered during execution
  double imbalance_tolerance; // input imbalance tolerance.
  int recursion_depth; // number of steps that partitioning will be solved in.
  int coord_dim; // coordinate dim
  int num_weights_per_coord; // # of weights per coord
  size_t initial_num_loc_coords; // initial num local coords.
  global_size_t initial_num_glob_coords; // initial num global coords.
  mj_lno_t num_local_coords; // number of local coords.
  mj_gno_t num_global_coords; // number of global coords.
  mj_scalar_t sEpsilon; // epsilon for mj_scalar_t

  // can distribute points on same coordinant to different parts.
  bool distribute_points_on_cut_lines;

  // how many parts we can calculate concurrently.
  mj_part_t max_concurrent_part_calculation;

  bool mj_run_as_rcb; // means recursion depth is adjusted to maximum value.
  int mj_user_recursion_depth; // the recursion depth value provided by user.
  bool mj_keep_part_boxes; // if the boxes need to be kept.

  // whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
  int check_migrate_avoid_migration_option;

  // when doing the migration, 0 will aim for perfect load-imbalance, 1 - will
  // aim for minimized number of messages with possibly bad load-imbalance
  int migration_type;

  // when MJ decides whether to migrate, the minimum imbalance for migration.
  double minimum_migration_imbalance;

  // Nonuniform first level partitioning
  // (Currently available only for sequential_task_partitioning):
  // Used for Dragonfly task mapping by partitioning Dragonfly RCA
  // machine coordinates and application coordinates.
  // An optimization that completely partitions the most important machine dimension
  // first (i.e. the Dragonfly group coordinate, or RCA's x coordinate). The standard
  // MJ alg follows after the nonuniform first level partitioning.
  //
  // Ex. (first level partitioning): If we have 120 elements,
  // num_first_level_parts = 3, first_level_distribution = [4, 10, 6], then
  // part sizes after first level will be [24, 60, 36]. Standard uniform MJ
  // continues for all subsequent levels.

  // If used, number of parts requested for a nonuniform
  // first level partitioning
  mj_part_t num_first_level_parts;

  // If used, the requested distribution of parts for the
  // nonuniform first level partitioning
  Kokkos::View<mj_part_t*, Kokkos::HostSpace> first_level_distribution;

  mj_part_t total_num_cut ;           // how many cuts will be totally
  mj_part_t total_num_part;           // how many parts will be totally

  mj_part_t max_num_part_along_dim ;  // maximum part count along a dimension.
  mj_part_t max_num_cut_along_dim;    // maximum cut count along a dimension.

  // maximum part+cut count along a dimension.
  size_t max_num_total_part_along_dim;

  mj_part_t total_dim_num_reduce_all;  // estimate on #reduceAlls can be done.

  // max no of parts that might occur during the partition before the last
  // partitioning dimension.
  mj_part_t last_dim_num_part;

  // input part array specifying num part to divide along each dim.
  Kokkos::View<mj_part_t *, Kokkos::HostSpace> part_no_array;

  // two dimension coordinate array
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t **, Kokkos::LayoutLeft, device_t>
    mj_coordinates;

  // two dimension weight array
  Kokkos::View<mj_scalar_t **, device_t> mj_weights;

  // if the target parts are uniform
  Kokkos::View<bool *, Kokkos::HostSpace> mj_uniform_parts;

  // if the coordinates have uniform weights
  Kokkos::View<bool *, Kokkos::HostSpace> mj_uniform_weights;

  int mj_num_teams; // the number of teams

  size_t num_global_parts; // the targeted number of parts

  // vector of all boxes for all parts, constructed if mj_keep_part_boxes true
  RCP<mj_partBoxVector_t> kept_boxes;

  RCP<mj_partBox_t> global_box;

  int myRank;           // processor rank
  int myActualRank;     // initial rank

  bool divide_to_prime_first;

  // initial global ids of the coordinates.
  Kokkos::View<const mj_gno_t*, device_t> initial_mj_gnos;

  // current global ids of the coordinates, might change during migration.
  Kokkos::View<mj_gno_t*, device_t> current_mj_gnos;

  // the actual processor owner of the coordinate, to track after migrations.
  Kokkos::View<int*, Kokkos::HostSpace> owner_of_coordinate;

  // permutation of coordinates, for partitioning.
  Kokkos::View<mj_lno_t*, device_t> coordinate_permutations;

  // permutation work array.
  Kokkos::View<mj_lno_t*, device_t> new_coordinate_permutations;

  // the part ids assigned to coordinates.
  Kokkos::View<mj_part_t*, device_t> assigned_part_ids;

  // beginning and end of each part.
  Kokkos::View<mj_lno_t *, device_t> part_xadj;

  // work array for beginning and end of each part.
  Kokkos::View<mj_lno_t *, device_t> new_part_xadj;

  Kokkos::View<mj_scalar_t *, device_t> all_cut_coordinates;

  // how much weight should a MPI put left side of the each cutline
  Kokkos::View<mj_scalar_t *, device_t>
    process_cut_line_weight_to_put_left;

  // weight percentage each thread in MPI puts left side of the each outline
  Kokkos::View<mj_scalar_t *, device_t>
    thread_cut_line_weight_to_put_left;

  // work array to manipulate coordinate of cutlines in different iterations.
  // necessary because previous cut line information is used for determining
  // the next cutline information. therefore, cannot update the cut work array
  // until all cutlines are determined.
  Kokkos::View<mj_scalar_t *, device_t> cut_coordinates_work_array;

  // Used for swapping above cut_coordinates_work_array
  Kokkos::View<mj_scalar_t *, device_t> temp_cut_coords;

  // cumulative part weight array.
  Kokkos::View<mj_scalar_t *, device_t> target_part_weights;

  // upper bound coordinate of a cut line
  Kokkos::View<mj_scalar_t *, device_t> cut_upper_bound_coordinates;

  // lower bound coordinate of a cut line
  Kokkos::View<mj_scalar_t *, device_t> cut_lower_bound_coordinates;

  // lower bound weight of a cut line
  Kokkos::View<mj_scalar_t *, device_t> cut_lower_bound_weights;

  // upper bound weight of a cut line
  Kokkos::View<mj_scalar_t *, device_t> cut_upper_bound_weights;

  // combined array to exchange the min and max coordinate, and total
  // weight of part.
  Kokkos::View<mj_scalar_t *, device_t>
    process_local_min_max_coord_total_weight;

  // global combined array with the results for min, max and total weight.
  Kokkos::View<mj_scalar_t *, device_t>
    global_min_max_coord_total_weight;

  // isDone is used to determine if a cutline is determined already. If a cut
  // line is already determined, the next iterations will skip this cut line.
  Kokkos::View<bool *, device_t> is_cut_line_determined;

  // incomplete_cut_count count holds the number of cutlines that have not
  // been finalized for each part when concurrentPartCount>1, using this
  // information, if incomplete_cut_count[x]==0, then no work is done
  // for this part.
  Kokkos::View<mj_part_t *, device_t> device_incomplete_cut_count;
  typename decltype(device_incomplete_cut_count)::HostMirror
    incomplete_cut_count;

  // Need a quick accessor for this on host
  typename decltype (part_xadj)::HostMirror host_part_xadj;

  // local part weights of each thread.
  Kokkos::View<double *, device_t>
    thread_part_weights;

  // the work manupulation array for partweights.
  Kokkos::View<double *, device_t>
    thread_part_weight_work;

  // thread_cut_left_closest_point to hold the closest coordinate
  // to a cutline from left (for each thread).
  Kokkos::View<mj_scalar_t *, device_t>
    thread_cut_left_closest_point;

  // thread_cut_right_closest_point to hold the closest coordinate
  // to a cutline from right (for each thread)
  Kokkos::View<mj_scalar_t *, device_t>
    thread_cut_right_closest_point;

  // to store how many points in each part a thread has.
  Kokkos::View<mj_lno_t *, device_t>
    thread_point_counts;

  Kokkos::View<mj_scalar_t *, device_t> process_rectilinear_cut_weight;
  Kokkos::View<mj_scalar_t *, device_t> global_rectilinear_cut_weight;

  // for faster communication, concatanation of
  // totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
  // leftClosest distances sized P-1, since P-1 cut lines
  // rightClosest distances size P-1, since P-1 cut lines.
  Kokkos::View<mj_scalar_t *, device_t>
    total_part_weight_left_right_closests;
  Kokkos::View<mj_scalar_t *, device_t>
    global_total_part_weight_left_right_closests;

  Kokkos::View<mj_part_t*, device_t> device_num_partitioning_in_current_dim;
  typename decltype(device_num_partitioning_in_current_dim)::HostMirror
    host_num_partitioning_in_current_dim; // for quick access on host

  /* \brief helper functio to calculate imbalance.
   * \param achieved balance we achieved.
   * \param expected balance expected.
   */
  KOKKOS_INLINE_FUNCTION
  double calculate_imbalance(mj_scalar_t achieved, mj_scalar_t expected) const {
    return static_cast<double>(achieved) / static_cast<double>(expected) - 1.0;
  }

  /* \brief Either the mj array (part_no_array) or num_global_parts should be
   * provided in the input. part_no_array takes precedence if both are
   * provided. Depending on these parameters, total cut/part number, maximum
   * part/cut number along a dimension, estimated number of reduceAlls,
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

  /* \brief for part communication we keep track of the box boundaries.
   * This is performed when either asked specifically, or when geometric
   * mapping is performed afterwards. This function initializes a single box
   * with all global min and max coordinates.
   * \param initial_partitioning_boxes the input and output vector for boxes.
   */
  void init_part_boxes(RCP<mj_partBoxVector_t> & outPartBoxes);

  /* \brief Function returns how many parts that will be obtained after this
   * dimension partitioning. It sets how many parts each current part will be
   * partitioned into in this dimension to device_num_partitioning_in_current_dim
   * vector, sets how many total future parts each obtained part will be
   * partitioned into in next_future_num_parts_in_parts vector, If part boxes
   * are kept, then sets initializes the output_part_boxes as its ancestor.
   * \param future_num_part_in_parts: input, how many future parts each
   * current part will be partitioned into.
   * \param next_future_num_parts_in_parts: output, how many future parts
   * each obtained part will be partitioned into.
   * \param future_num_parts: output, max number of future parts that will be
   * obtained from a single
   * \param current_num_parts: input, how many parts are there currently.
   * \param current_iteration: input, current dimension iteration number.
   * \param input_part_boxes: input, if boxes are kept, current boxes.
   * \param output_part_boxes: output, if boxes are kept, the initial box
   * boundaries for obtained parts.
   * \param atomic_part_count  // DOCWORK: Documentation
   */
  mj_part_t update_part_num_arrays(
    std::vector<mj_part_t> *future_num_part_in_parts,
    std::vector<mj_part_t> *next_future_num_parts_in_parts,
    mj_part_t &future_num_parts,
    mj_part_t current_num_parts,
    int current_iteration,
    RCP<mj_partBoxVector_t> input_part_boxes,
    RCP<mj_partBoxVector_t> output_part_boxes,
    mj_part_t atomic_part_count);

  /*! \brief Function that calculates the next pivot position,
   * according to given coordinates of upper bound and lower bound, the
   * weights at upper and lower bounds, and the expected weight.
   * \param cut_upper_bound is the upper bound coordinate of the cut.
   * \param cut_lower_bound is the lower bound coordinate of the cut.
   * \param cut_upper_weight is the weights at the upper bound of the cut.
   * \param cut_lower_weight is the weights at the lower bound of the cut.
   * \param expected_weight is the expected weight that should be placed on
   * the left of the cut line.
   * \param new_cut_position DOCWORK: Documentation
   */
  KOKKOS_INLINE_FUNCTION
  void mj_calculate_new_cut_position (
    mj_scalar_t cut_upper_bound,
    mj_scalar_t cut_lower_bound,
    mj_scalar_t cut_upper_weight,
    mj_scalar_t cut_lower_weight,
    mj_scalar_t expected_weight,
    mj_scalar_t &new_cut_position);

  /*! \brief Function checks if should do migration or not.
   * It returns true to point that migration should be done when
   * -migration_reduce_all_population are higher than a predetermined value
   * -num_coords_for_last_dim_part that left for the last dimension
   * partitioning is less than a predetermined value - the imbalance of the
   * processors on the parts are higher than given threshold.
   * \param input_num_parts is the number of parts when migration is called.
   * \param output_num_parts is the output number of parts after migration.
   * \param next_future_num_parts_in_parts is the number of total future parts
   * each part is partitioned into. Updated when migration is performed.
   * \param output_part_begin_index is the number that will be used as
   * beginning part number when final solution part numbers are assigned.
   * \param migration_reduce_all_population is the estimated total number of
   * reduceall operations multiplied with number of processors to be used for
   * determining migration.
   * \param num_coords_for_last_dim_part is the estimated number of points in
   * each part, when last dimension partitioning is performed.
   * \param iteration is the string that gives information about the dimension
   * for printing purposes.
   * \param input_part_boxes is the array that holds the part boxes after the
   * migration. (swapped)
   * \param output_part_boxes is the array that holds the part boxes before
   * the migration. (swapped)
   */
  bool mj_perform_migration(
    mj_part_t in_num_parts, //current number of parts
    mj_part_t &out_num_parts, //output number of parts.
    std::vector<mj_part_t> *next_future_num_parts_in_parts,
    mj_part_t &output_part_begin_index,
    size_t migration_reduce_all_population,
    mj_lno_t num_coords_for_last_dim_part,
    std::string iteration,
    RCP<mj_partBoxVector_t> &input_part_boxes,
    RCP<mj_partBoxVector_t> &output_part_boxes);

  /*! \brief Function checks if should do migration or not.
   * It returns true to point that migration should be done when
   * -migration_reduce_all_population are higher than a predetermined value
   * -num_coords_for_last_dim_part that left for the last dimension
   * partitioning is less than a predetermined value - the imbalance of the
   * processors on the parts are higher than given threshold.
   * \param migration_reduce_all_population is the multiplication of the
   * number of reduceall operations estimated and the number of processors.
   * \param num_coords_for_last_dim_part is the estimated number of
   * coordinates in a part per processor in the last dimension partitioning.
   * \param num_procs is the number of processor attending to migration
   * operation.
   * \param num_parts is the number of parts that exist in the current
   * partitioning.
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
   * that holds which part each coordinate should be sent. In addition it
   * calculates the shift amount (output_part_numbering_begin_index) to be
   * done when final numberings of the parts are performed.
   * \param num_points_in_all_processor_parts is the array holding the num
   * points in each part in each proc.
   * \param num_parts is the number of parts that exist in the current
   * partitioning.
   * \param num_procs is the number of processor attending to migration
   * operation.
   * \param send_count_to_each_proc array array storing the number of points
   * to be sent to each part.
   * \param processor_ranks_for_subcomm is the ranks of the processors that
   * will be in the subcommunicator with me.
   * \param next_future_num_parts_in_parts is the vector, how many more parts
   * each part will be divided into in the future.
   * \param out_num_part is the number of parts assigned to the process.
   * \param out_part_indices is the indices of the part to which the processor
   * is assigned.
   * \param output_part_numbering_begin_index is how much the numbers should
   * be shifted when numbering the result parts.
   * \param coordinate_destinations is the output array that holds which part
   * each coordinate should be sent.
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

  /*! \brief Function that assigned the processors to parts, when there are
   * more processors then parts. Sets the destination of each coordinate in
   * coordinate_destinations, also edits output_part_numbering_begin_index,
   * and out_part_index, and returns the processor_ranks_for_subcomm which
   * represents the ranks of the processors that will be used for creating the
   * subcommunicator.
   * \param num_points_in_all_processor_parts is the array holding the num
   * points in each part in each proc.
   * \param num_parts is the number of parts that exist in the current
   * partitioning.
   * \param num_procs is the number of processor attending to migration
   * operation.
   * \param send_count_to_each_proc array array storing the number of points
   * to be sent to each part.
   * \param processor_ranks_for_subcomm is the ranks of the processors that
   * will be in the subcommunicator with me.
   * \param next_future_num_parts_in_parts is the vector, how many more parts
   * each part will be divided into in the future.
   * \param out_part_index is the index of the part to which the processor
   * is assigned.
   * \param output_part_numbering_begin_index is how much the numbers should
   * be shifted when numbering the result parts.
   * \param coordinate_destinations is the output array that holds which part
   * each coordinate should be sent.
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
   * \param num_parts is the number of parts that exist in the
   * current partitioning.
   * \param num_procs is the number of processors attending to
   * migration operation.
   * \param part_assignment_proc_begin_indices ([i]) points to the first
   * processor index that part i will be sent to.
   * \param processor_chains_in_parts the array that holds the linked list
   * structure, started from part_assignment_proc_begin_indices ([i]).
   * \param send_count_to_each_proc array array storing the number of points to
   * be sent to each part.
   * \param coordinate_destinations is the output array that holds which part
   * each coordinate should be sent.
   */
  void assign_send_destinations(
    mj_part_t num_parts,
    mj_part_t *part_assignment_proc_begin_indices,
    mj_part_t *processor_chains_in_parts,
    mj_lno_t *send_count_to_each_proc,
    int *coordinate_destinations);

  /*! \brief Function fills up coordinate_destinations is the output array
   * that holds which part each coordinate should be sent. In addition it
   * calculates the shift amount (output_part_numbering_begin_index) to be done
   * when final numberings of the parts are performed.
   * \param num_parts is the number of parts in the current partitioning.
   * \param sort_item_part_to_proc_assignment is the sorted parts with respect
   * to the assigned processors.
   * \param coordinate_destinations is the output array that holds which part
   * each coordinate should be sent.
   * \param output_part_numbering_begin_index is how much the numbers should be
   * shifted when numbering the result parts.
   * \param next_future_num_parts_in_parts is the vector, how many more parts
   * each part will be divided into in the future.
   */
  void assign_send_destinations2(
    mj_part_t num_parts,
    uSortItem<mj_part_t, mj_part_t> * sort_item_part_to_proc_assignment,
    int *coordinate_destinations,
    mj_part_t &output_part_numbering_begin_index,
    std::vector<mj_part_t> *next_future_num_parts_in_parts);

  /*! \brief Function fills up coordinate_destinations is the output array
   * that holds which part each coordinate should be sent. In addition it
   * calculates the shift amount (output_part_numbering_begin_index) to be done
   * when final numberings of the parts are performed.
   * \param num_points_in_all_processor_parts is the array holding the num
   * points in each part in each proc.
   * \param num_parts is the number of parts that exist in the current
   * partitioning.
   * \param num_procs is the number of processors attending to
   * migration operation.
   * \param send_count_to_each_proc array array storing the number of points to
   * be sent to each part.
   * \param next_future_num_parts_in_parts is the vector, how many more parts
   * each part will be divided into in the future.
   * \param out_num_part is the number of parts assigned to the process.
   * \param out_part_indices is the indices of the part to which the processor
   * is assigned.
   * \param output_part_numbering_begin_index is how much the numbers should be
   * shifted when numbering the result parts.
   * \param coordinate_destinations is the output array that holds which parta
   * each coordinate should be sent.
   */
  void mj_assign_parts_to_procs(
    mj_gno_t * num_points_in_all_processor_parts,
    mj_part_t num_parts,
    mj_part_t num_procs,
    mj_lno_t *send_count_to_each_proc,
    std::vector<mj_part_t> *next_future_num_parts_in_parts,
    mj_part_t &out_num_part,
    std::vector<mj_part_t> &out_part_indices,
    mj_part_t &output_part_numbering_begin_index,
    int *coordinate_destinations);

  /*! \brief Function fills up coordinate_destinations is the output array
   * that holds which part each coordinate should be sent. In addition it
   * calculates the shift amount (output_part_numbering_begin_index) to be done
   * when final numberings of the parts are performed.
   * \param num_procs is the number of processora attending to
   * migration operation.
   * \param num_new_local_points is the output to represent the new number
   * of local points.
   * \param iteration is the string for the current iteration.
   * \param coordinate_destinations is the output array that holds which part
   * each coordinate should be sent.
   * \param num_parts is the number of parts in the current partitioning.
   */
  void mj_migrate_coords(
    mj_part_t num_procs,
    mj_lno_t &num_new_local_points,
    std::string iteration,
    int *coordinate_destinations,
    mj_part_t num_parts);

  /*! \brief Function creates the new subcomminicator for the processors
   * given in processor_ranks_for_subcomm.
   * \param processor_ranks_for_subcomm is the vector that has the ranks of
   * the processors that will be in the same group.
   */
  void create_sub_communicator(
    std::vector<mj_part_t> &processor_ranks_for_subcomm);

  /*! \brief Function returns the largest prime factor of a given number.
   * input and output are integer-like.
   * \param num_parts  DOCWORK: documentation
   */
  mj_part_t find_largest_prime_factor(mj_part_t num_parts) {
    mj_part_t largest_factor = 1;
    mj_part_t n = num_parts;
    mj_part_t divisor = 2;
    while (n > 1) {
      while (n % divisor == 0) {
        n = n / divisor;
        largest_factor = divisor;
      }
      ++divisor;
      if(divisor * divisor > n) {
        if(n > 1) {
          largest_factor = n;
        }
        break;
      }
    }
    return largest_factor;
  }

public:
  AlgMJ();

  // DOCWORK: Make param documentation use : consistently
  /*! \brief Multi Jagged coordinate partitioning algorithm.
   * \param env: library configuration and problem parameters
   * \param problemComm: the communicator for the problem
   * \param imbalance_tolerance: the input provided imbalance tolerance.
   * \param num_teams: number of teams for CUDA kernels.
   * \param num_global_parts: number of target global parts.
   * \param part_no_array: part no array, if provided this will be used
   * for partitioning.
   * \param recursion_depth: if part no array is provided, it is the length of
   * part no array, if part no is not provided than it is the number of steps
   * that algorithm will divide into num_global_parts parts.
   * \param coord_dim: coordinate dimension
   * \param num_local_coords: number of local coordinates
   * \param num_global_coords: number of global coordinates
   * \param initial_mj_gnos: the list of initial global id's
   * \param mj_coordinates: the two dimensional coordinate array.
   * \param num_weights_per_coord: number of weights per coordinate
   * \param mj_uniform_weights: if weight index [i] has uniform weight or not.
   * \param mj_weights: the two dimensional array for weights
   * \param mj_uniform_parts: if the target partitioning aims uniform parts
   * \param result_assigned_part_ids: Output - the result partids corresponding
   * to the coordinates given im result_mj_gnos.
   * \param result_mj_gnos: Output - the result coordinate global id's
   * corresponding to the part_ids array.
   */
  void multi_jagged_part(
    const RCP<const Environment> &env,
    RCP<const Comm<int> > &problemComm,
    double imbalance_tolerance,
    int num_teams,
    size_t num_global_parts,
    Kokkos::View<mj_part_t*, Kokkos::HostSpace> & part_no_array,
    int recursion_depth,
    int coord_dim,
    mj_lno_t num_local_coords,
    mj_gno_t num_global_coords,
    Kokkos::View<const mj_gno_t*, device_t> & initial_mj_gnos,
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t> & mj_coordinates,
    int num_weights_per_coord,
    Kokkos::View<bool*, Kokkos::HostSpace> & mj_uniform_weights,
    Kokkos::View<mj_scalar_t**, device_t> & mj_weights,
    Kokkos::View<bool*, Kokkos::HostSpace> & mj_uniform_parts,
    Kokkos::View<mj_part_t*, device_t> & result_assigned_part_ids,
    Kokkos::View<mj_gno_t*, device_t> & result_mj_gnos);

  /*! \brief Multi Jagged coordinate partitioning algorithm.
   * \param distribute_points_on_cut_lines_ : if partitioning can distribute
   * points on same coordinate to different parts.
   * \param max_concurrent_part_calculation_ : how many parts we can calculate
   * concurrently.
   * \param check_migrate_avoid_migration_option_ : whether to migrate=1, avoid
   * migrate=2, or leave decision to MJ=0
   * \param minimum_migration_imbalance_ : when MJ decides whether to migrate,
   * the minimum imbalance for migration.
   * \param migration_type_ : when MJ migration whether to migrate for perfect
   * load-imbalance or less messages
   */
  void set_partitioning_parameters(
    bool distribute_points_on_cut_lines_,
    int max_concurrent_part_calculation_,
    int check_migrate_avoid_migration_option_,
    double minimum_migration_imbalance_,
    int migration_type_ = 0);

  /*! \brief Function call, if the part boxes are intended to be kept.
   */
  void set_to_keep_part_boxes();

  /*! \brief DOCWORK: Documentation
   */
  RCP<mj_partBox_t> get_global_box() const;

  /*! \brief DOCWORK: Documentation
   */
  RCP<mj_partBoxVector_t> get_kept_boxes() const;

  /*! \brief DOCWORK: Documentation
   */
  RCP<mj_partBoxVector_t> compute_global_box_boundaries(
    RCP<mj_partBoxVector_t> &localPartBoxes) const;

  /*! \brief Special function for partitioning for task mapping.
   * Runs sequential, and performs deterministic partitioning for the
   * partitioning the points along a cutline.
   * \param env library configuration and problem parameters
   * \param num_total_coords number of total coordinates
   * \param num_selected_coords : the number of selected coordinates. This is
   * to set, if there are n processors, but only m<n processors are selected for
   * mapping.
   * \param num_target_part: number of target global parts.
   * \param coord_dim_: coordinate dimension for coordinates
   * \param mj_coordinates_: the coordinates
   * \param initial_selected_coords_output_permutation: Array allocated by
   * caller, in the size of num_total_coords, first num_selected_coords elements
   * should list the indices of the selected processors. This is output for
   * output permutation array.
   * \param output_xadj: The output part xadj array, pointing beginning and end
   * of each part on output permutation array (inital_adjList_output_adjlist).
   * Returned in CSR format: part i's info in output_xadj[i] : output_xadj[i+1]
   * \param recursion_depth: recursion depth
   * \param part_no_array_: possibly null part_no_array, specifying how many
   * parts each should be divided during partitioning.
   * \param partition_along_longest_dim  DOCWORK: Documentation
   * \param divide_to_prime_first_ DOCWORK: Documentation
   *  Nonuniform first level partitioning (Currently available only for sequential_task_partitioning):
   *  Currently used for Dragonfly task mapping by partitioning Dragonfly RCA
   *  machine coordinates and application coordinates.
   *  An optimization that completely partitions the most important machine dimension
   *  first (i.e. the Dragonfly group coordinate, or RCA's x coordinate). The standard
   *  MJ alg follows after the nonuniform first level partitioning.
   *  \param num_first_level_parts_: the number of parts after the first level of
   *  partitioning (may be nonuniform)
   *  \param first_level_distribution_: a view containing the distribution of
   *  elements in each part after the first cut (used for nonuniform first cuts)
   *  Ex. (first level partitioning): If we have 120 elements,
   *  num_first_level_parts = 3, first_level_distribution = [4, 10, 6], then
   *  part sizes after first level will be [24, 60, 36]. Standard uniform MJ
   *  continues for all subsequent levels.
   */
  void sequential_task_partitioning(
    const RCP<const Environment> &env,
    mj_lno_t num_total_coords,
    mj_lno_t num_selected_coords,
    size_t num_target_part,
    int coord_dim,
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<mj_scalar_t **, Kokkos::LayoutLeft, device_t> & mj_coordinates_,
    Kokkos::View<mj_lno_t *, device_t> &
      initial_selected_coords_output_permutation,
    mj_lno_t *output_xadj,
    int recursion_depth_,
    const Kokkos::View<mj_part_t *, Kokkos::HostSpace> & part_no_array,
    bool partition_along_longest_dim,
    int num_ranks_per_node,
    bool divide_to_prime_first_,
    mj_part_t num_first_level_parts_ = 1,
    const Kokkos::View<mj_part_t *, Kokkos::HostSpace> & first_level_distribution_
      = Kokkos::View<mj_part_t *, Kokkos::HostSpace>());

#ifdef KOKKOS_ENABLE_CUDA
  public:
#else
  private:
#endif

  /* \brief Allocates all required memory for the mj partitioning algorithm.
   */
  void allocate_set_work_memory();

  /* \brief compute global bounding box:  min/max coords of global domain */
  void compute_global_box();

  // DOCWORK: Inconsisent use of ! for descriptive/brief commenting - decide.
  /*! \brief Function to determine the local minimum and maximum coordinate,
   * and local total weight in the given set of local points.
   * \current_work_part  DOCWORK: Documentation
   * \current_concurrent_num_parts  DOCWORK: Documentation
   * \mj_current_dim_coords DOCWORK: Documentation
   */
  void mj_get_local_min_max_coord_totW(
    mj_part_t current_work_part,
    mj_part_t current_concurrent_num_parts,
    Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords);

  /*! \brief Function that reduces global minimum and maximum coordinates with
   * global total weight from given local arrays.
   * \param current_concurrent_num_parts is the number of parts whose cut
   * lines will be calculated concurrently.
   * \param local_min_max_total is the array holding local min and max
   * coordinate values with local total weight.
   * First current_concurrent_num_parts entries are minimums of the parts,
   * next current_concurrent_num_parts entries are max and then total weights.
   * \param global_min_max_total is the output array holding global min and
   * global coordinate values with global total weight.
   * The structure is same as local_min_max_total.
   */
  void mj_get_global_min_max_coord_totW(
    mj_part_t current_concurrent_num_parts,
    Kokkos::View<mj_scalar_t *, device_t> & local_min_max_total,
    Kokkos::View<mj_scalar_t *, device_t> & global_min_max_total);

  /*! \brief Function that calculates the new coordinates for the cut lines.
   * Function is called inside the parallel region.
   * \param min_coord minimum coordinate in the range.
   * \param max_coord maximum coordinate in the range.
   * \param num_cuts holds number of cuts in current partitioning dimension.
   * \param global_weight holds the global total weight in the current part.
   * \param initial_cut_coords is the output array for the initial cut lines.
   * \param target_part_weights is the output array holding the cumulative
   * ratios of parts in current partitioning.
   * For partitioning to 4 uniformly, target_part_weights will be
   * (0.25 * globalTotalWeight, 0.5 *globalTotalWeight , 0.75 *
   * globalTotalWeight, globalTotalWeight).
   * \param future_num_part_in_parts is the vector that holds how many more
   * parts each part will be divided into more
   * for the parts at the beginning of this coordinate partitioning
   * \param next_future_num_parts_in_parts is the vector that holds how many
   * more parts each part will be divided into more for the parts that will be
   * obtained at the end of this coordinate partitioning.
   * \param concurrent_current_part is the index of the part in the
   * future_num_part_in_parts vector.
   * \param obtained_part_index holds the amount of shift in the
   * next_future_num_parts_in_parts for the output parts.
   * Nonuniform first level partitioning:
   * \param num_target_first_level_parts is the number of parts requested
   * after the first level of partitioning (resulting parts may be imbalanced)
   * \param target_first_level_dist is an array requesting the distribution of
   * elements in each part after the first cut (used for nonuniform first cuts)
   * Ex. If we have num_first_level_parts = 3, first_level_dist = [4, 10, 6], then
   * target_part_weights will be [.20, .70, 1.00] * global_weight
   */
  void mj_get_initial_cut_coords_target_weights(
    mj_scalar_t min_coord,
    mj_scalar_t max_coord,
    mj_part_t num_cuts/*p-1*/ ,
    mj_scalar_t global_weight,
    Kokkos::View<mj_scalar_t *, device_t> & initial_cut_coords,
    Kokkos::View<mj_scalar_t *, device_t> & target_part_weights,
    std::vector <mj_part_t> *future_num_part_in_parts,
    std::vector <mj_part_t> *next_future_num_parts_in_parts,
    mj_part_t concurrent_current_part,
    mj_part_t obtained_part_index,
    mj_part_t num_target_first_level_parts = 1,
    const Kokkos::View<mj_part_t *, Kokkos::HostSpace> & target_first_level_dist =
      Kokkos::View<mj_part_t *, Kokkos::HostSpace>());

  /*! \brief Function that calculates the new coordinates for the cut lines.
   * Function is called inside the parallel region.
   * \param max_coordinate maximum coordinate in the range.
   * \param min_coordinate minimum coordinate in the range.
   * \param concurrent_current_part_index is the index of the part in the
   * inTotalCounts vector.
   * \param coordinate_begin_index holds the beginning of the coordinates
   * in current part.
   * \param coordinate_end_index holds end of the coordinates in current part.
   * \param mj_current_coordinate_permutations is the permutation array, holds
   * the real indices of coordinates on mj_current_dim_coords array.
   * \param mj_current_dim_coords is the 1D array holding the coordinates.
   * \param mj_part_ids is the array holding the partIds of each coordinate.
   * \param partition_count is the number of parts that the current part will
   * be partitioned into.
   */
  void set_initial_coordinate_parts(
    mj_scalar_t &max_coordinate,
    mj_scalar_t &min_coordinate,
    mj_lno_t coordinate_begin_index,
    mj_lno_t coordinate_end_index,
    Kokkos::View<mj_lno_t *, device_t> &
      mj_current_coordinate_permutations,
    Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
    Kokkos::View<mj_part_t *, device_t> & mj_part_ids,
    mj_part_t &partition_count);

  /*! \brief Function that is responsible from 1D partitioning of the given
   * range of coordinates.
   * \param mj_current_dim_coords is 1 dimensional array holding coordinate
   * values.
   * \param imbalanceTolerance is the maximum allowed imbalance ratio.
   * \param current_work_part is the beginning index of concurrentPartCount
   * parts.
   * \param current_concurrent_num_parts is the number of parts whose cut
   * lines will be calculated concurrently.
   * \param current_cut_coordinates is the array holding the coordinates of
   * the cut.
   * \param total_incomplete_cut_count is the number of cut lines whose
   * positions should be calculated.
   * \param view_rectilinear_cut_count   DOCWORK: Documentation
   * \param view_total_reduction_size   DOCWORK: Documentation
   */
  void mj_1D_part(
    Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
    double imbalanceTolerance,
    mj_part_t current_work_part,
    mj_part_t current_concurrent_num_parts,
    Kokkos::View<mj_scalar_t *, device_t> & current_cut_coordinates,
    mj_part_t total_incomplete_cut_count,
    Kokkos::View<mj_part_t *, device_t> & view_rectilinear_cut_count,
    Kokkos::View<size_t*, device_t> & view_total_reduction_size);

  /*! \brief Function that calculates the weights of each part according to
   * given part cut coordinates. Function is called inside the parallel
   * region. Thread specific work arrays are provided as function parameter.
   * DOCWORK: Documentation params
   */
  void mj_1D_part_get_part_weights(
    mj_part_t current_concurrent_num_parts,
    mj_part_t current_work_part,
    Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
    int loop_count);

  /*! \brief Function that reduces the result of multiple threads for
   * left and right closest points and part weights in a single mpi process.
   * \param current_work_part holds the index of the first part (important
   * when concurrent parts are used.)
   * \param current_concurrent_num_parts is the number of parts whose cut
   * lines will be calculated concurrently.
   */
  void mj_combine_rightleft_and_weights(
    mj_part_t current_work_part,
    mj_part_t current_concurrent_num_parts);

  /*! \brief Function that determines the permutation indices of coordinates.
   * \param num_parts is the number of parts.
   * \param current_concurrent_work_part DOCWORK: Documentation
   * \param mj_current_dim_coords is 1 dimensional array holding the
   * coordinate values.
   * \param current_concurrent_cut_coordinate is 1 dimensional array holding
   * the cut coordinates.
   * \param used_local_cut_line_weight_to_left holds how much weight of the
   * coordinates on the cutline should be put on left side.
   * \param out_part_xadj is the indices of coordinates calculated for the
   * partition on next dimension.
   */
  void mj_create_new_partitions(
    mj_part_t num_parts,
    mj_part_t current_concurrent_work_part,
    Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
    Kokkos::View<mj_scalar_t *, device_t> & current_concurrent_cut_coordinate,
    Kokkos::View<mj_scalar_t *, device_t> & used_local_cut_line_weight_to_left,
    Kokkos::View<mj_lno_t *, device_t> & out_part_xadj);

  /*! \brief Function that calculates the new coordinates for the cut lines.
   * Function is called inside the parallel region. Write the new cut
   * coordinates to new_current_cut_coordinates, and determines if the final
   * position of a cut is found.
   * \param current_concurrent_num_parts  DOCWORK: Documentation
   * \param kk  DOCWORK: Documentation and pick better name
   * \param num_total_part is the sum of number of cutlines and number of
   * parts. Simply it is 2*P - 1.
   * \param num_cuts is the number of cut lines. P - 1.
   * \param used_imbalance_tolerance is the maximum allowed imbalance ratio.
   * \param current_global_part_weights is the array holding the weight of
   * parts. Assumes there are 2*P - 1 parts (cut lines are seperate parts).
   * \param current_local_part_weights is local totalweight of the processor.
   * \param current_part_target_weights desired cumulative part ratios, size P.
   * \param current_cut_line_determined is the boolean array to determine if
   * the correct position for a cut line is found.
   * \param current_cut_coordinates is the array holding the coordinates of
   * each cut line. Sized P - 1.
   * \param current_cut_upper_bounds is the array holding the upper bound
   * coordinate for each cut line. Sized P - 1.
   * \param current_cut_lower_bounds is the array holding the lower bound
   * coordinate for each cut line. Sized P - 1.
   * \param current_global_left_closest_points is the array holding the
   * closest points to the cut lines from left.
   * \param current_global_right_closest_points is the array holding the
   * closest points to the cut lines from right.
   * \param current_cut_lower_bound_weights is the array holding the weight
   * of the parts at the left of lower bound coordinates.
   * \param current_cut_upper_weights is the array holding the weight of the
   * parts at the left of upper bound coordinates.
   * \param new_current_cut_coordinates is the work array, sized P - 1.
   * \param current_part_cut_line_weight_to_put_left DOCWORK: Documentation
   * \param view_rectilinear_cut_count is the count of cut lines whose balance
   * is acheived via distributing points in same coordinate to different parts.
   */
  void mj_get_new_cut_coordinates(
    mj_part_t current_concurrent_num_parts,
    mj_part_t kk,
    const mj_part_t &num_cuts,
    const double &used_imbalance_tolerance,
    Kokkos::View<mj_scalar_t *, device_t> & current_global_part_weights,
    Kokkos::View<mj_scalar_t *, device_t> & current_local_part_weights,
    Kokkos::View<mj_scalar_t *, device_t> & current_part_target_weights,
    Kokkos::View<bool *, device_t> & current_cut_line_determined,
    Kokkos::View<mj_scalar_t *, device_t> & current_cut_coordinates,
    Kokkos::View<mj_scalar_t *, device_t> & current_cut_upper_bounds,
    Kokkos::View<mj_scalar_t *, device_t> & current_cut_lower_bounds,
    Kokkos::View<mj_scalar_t *, device_t> & current_global_left_closest_points,
    Kokkos::View<mj_scalar_t *, device_t> & current_global_right_closest_points,
    Kokkos::View<mj_scalar_t *, device_t> & current_cut_lower_bound_weights,
    Kokkos::View<mj_scalar_t *, device_t> & current_cut_upper_weights,
    Kokkos::View<mj_scalar_t *, device_t> & new_current_cut_coordinates,
    Kokkos::View<mj_scalar_t *, device_t> &
      current_part_cut_line_weight_to_put_left,
    Kokkos::View<mj_part_t *, device_t> & view_rectilinear_cut_count);

  /*! \brief Function fills up the num_points_in_all_processor_parts, so that
   * it has the number of coordinates in each processor of each part.
   * to access how many points processor i has on part j,
   * num_points_in_all_processor_parts[i * num_parts + j].
   * \param num_procs is the number of processors for migration operation.
   * \param num_parts is the number of parts in the current partitioning.
   * \param num_points_in_all_processor_parts is the output array that holds
   * the number of coordinates in each part in each processor.
   */
  void get_processor_num_points_in_parts(
    mj_part_t num_procs,
    mj_part_t num_parts,
    mj_gno_t *&num_points_in_all_processor_parts);

  /*! \brief Function writes the new permutation arrays after the migration.
   * \param output_num_parts is the number of parts assigned to the processor.
   * \param num_parts is the number of parts right before migration.
   */
  void fill_permutation_array(
    mj_part_t output_num_parts,
    mj_part_t num_parts);

  /*! \brief Function creates consistent chunks for task partitioning. Used only
   * in the case of sequential task partitioning, where consistent handle of the
   * points on the cuts are required.
   * \param num_parts is the number of parts.
   * \param mj_current_dim_coords is 1 dimensional array holding the
   * coordinate values.
   * \param current_concurrent_cut_coordinate is 1 dimensional array holding
   * the cut coordinates.
   * \param coordinate_begin is the start index of the given partition on
   * partitionedPointPermutations.
   * \param coordinate_end is the end index of the given partition on
   * partitionedPointPermutations.
   * \param used_local_cut_line_weight_to_left holds how much weight of the
   * coordinates on the cutline should be put on left side.
   * \param out_part_xadj is the indices of begginning and end of the parts in
   * the output partition.
   * \param coordInd is the index according to which the partitioning is done.
   * \param longest_dim_part  DOCWORK: documentation
   * \param mj_scalar_t DOCWORK: documentation
   * \param p_coord_dimension_range_sorted DOCWORK: documentation
   */
  void create_consistent_chunks(
    mj_part_t num_parts,
    Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
    Kokkos::View<mj_scalar_t *, device_t> & current_concurrent_cut_coordinate,
    mj_lno_t coordinate_begin,
    mj_lno_t coordinate_end,
    Kokkos::View<mj_scalar_t *, device_t> & used_local_cut_line_weight_to_left,
    Kokkos::View<mj_lno_t *, device_t> & out_part_xadj,
    int coordInd,
    bool longest_dim_part,
    uSignedSortItem<int, mj_scalar_t, char> *p_coord_dimension_range_sorted);

  /*! \brief Function checks if should do migration or not.
   * \param current_num_parts is the number of parts in the process.
   * \param output_part_begin_index is the number that will be used as
   * beginning part number
   * \param output_part_boxes is the array that holds the part boxes
   * \param is_data_ever_migrated true if data was migrated
   * if the data is ever migrated during the partitioning.
   */
  void set_final_parts(
    mj_part_t current_num_parts,
    mj_part_t output_part_begin_index,
    RCP<mj_partBoxVector_t> &output_part_boxes,
    bool is_data_ever_migrated);
};

/*! \brief Multi Jagged coordinate partitioning algorithm default constructor.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::AlgMJ():
  mj_env(), mj_problemComm(), comm(), imbalance_tolerance(0),
  recursion_depth(0), coord_dim(0),
  num_weights_per_coord(0), initial_num_loc_coords(0),
  initial_num_glob_coords(0),
  num_local_coords(0), num_global_coords(0),
  sEpsilon(std::numeric_limits<mj_scalar_t>::epsilon() * 100),
  distribute_points_on_cut_lines(true),
  max_concurrent_part_calculation(1),
  mj_run_as_rcb(false), mj_user_recursion_depth(0),
  mj_keep_part_boxes(false),
  check_migrate_avoid_migration_option(0), migration_type(0),
  minimum_migration_imbalance(0.30),
  num_first_level_parts(1),
  total_num_cut(0), total_num_part(0), max_num_part_along_dim(0),
  max_num_cut_along_dim(0),
  max_num_total_part_along_dim(0),
  total_dim_num_reduce_all(0),
  last_dim_num_part(0),
  mj_num_teams(0),
  num_global_parts(1),
  kept_boxes(), global_box(),
  myRank(0), myActualRank(0),
  divide_to_prime_first(false)
{
}

/*! \brief Special function for partitioning for task mapping.
 * Runs sequential, and performs deterministic partitioning for the
 * partitioning the points along a cutline.
 * \param env library configuration and problem parameters
 * \param num_total_coords number of total coordinates
 * \param num_selected_coords : the number of selected coordinates. This is to
 * set, if there are n processors, but only m<n processors are selected for
 * mapping.
 * \param num_target_part: number of target global parts.
 * \param coord_dim_: coordinate dimension for coordinates
 * \param mj_coordinates_: the coordinates
 * \param inital_adjList_output_adjlist: Array allocated by caller, in the size
 * of num_total_coords, first num_selected_coords elements should list the
 * indices of the selected processors. This is output for output permutation
 * array.
 * \param output_xadj: The output part xadj array, pointing beginning and end of
 * each part on output permutation array (inital_adjList_output_adjlist).
 * Returned in CSR format: part i's info in output_xadj[i] : output_xadj[i+1]
 * \param rd: recursion depth
 * \param part_no_array_: possibly null part_no_array, specifying how many parts
 * each should be divided during partitioning.
 * \param partition_along_longest_dim DOCWORK: Documentation
 * \param num_ranks_per_node DOCWORK: Documentation
 * \param divide_to_prime_first_ DOCWORK: Documentation
 *
 *  Nonuniform first level partitioning:
 *  Currently used for Dragonfly task mapping by partitioning Dragonfly RCA
 *  machine coordinates and application coordinates.
 *  An optimization that completely partitions the most important machine
 *  dimension first (i.e. the Dragonfly group coordinate, or RCA's x
 *  coordinate). The standard MJ alg follows after the nonuniform first level
 *  partitioning.
 *
 * \param num_target_first_level_parts: the number of parts requested after
 * the first level of partitioning (resulting parts may be imbalanced)
 * \param first_level_dist: an array requesting the distribution of elements
 * in each part after the first cut (used for nonuniform first cuts)
 *
 *  Ex. (first level partitioning): If we have 120 elements,
 *  num_first_level_parts = 3, first_level_distribution = [4, 10, 6], then
 *  part sizes after first level will be [24, 60, 36]. Standard uniform MJ
 *  continues for all subsequent levels.
*/
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  sequential_task_partitioning(
  const RCP<const Environment> &env,
  mj_lno_t num_total_coords,
  mj_lno_t num_selected_coords,
  size_t num_target_part,
  int coord_dim_,
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t **, Kokkos::LayoutLeft, device_t> &
    mj_coordinates_,
  Kokkos::View<mj_lno_t *, device_t> & initial_adjList_output_adjlist,
  mj_lno_t *output_xadj,
  int recursion_depth_,
  const Kokkos::View<mj_part_t *, Kokkos::HostSpace> & part_no_array_,
  bool partition_along_longest_dim,
  int num_ranks_per_node,
  bool divide_to_prime_first_,
  mj_part_t num_first_level_parts_,
  const Kokkos::View<mj_part_t *, Kokkos::HostSpace> & first_level_distribution_)
{
  this->mj_env = env;
  const RCP<Comm<int> > commN;
  this->mj_problemComm = Teuchos::DefaultComm<int>::getDefaultSerialComm(commN);
  this->comm = Teuchos::rcp_const_cast<Comm<int> >(this->mj_problemComm);
  this->myActualRank = this->myRank = 1;

  this->divide_to_prime_first = divide_to_prime_first_;
  //weights are uniform for task mapping

  //parts are uniform for task mapping
  //as input indices.
  this->imbalance_tolerance = 0;
  this->num_global_parts = num_target_part;
  this->part_no_array = part_no_array_;
  this->recursion_depth = recursion_depth_;

  // If nonuniform first level partitioning, the requested num of parts and the
  // requested distribution of elements for each part
  this->num_first_level_parts = num_first_level_parts_;

  this->first_level_distribution = first_level_distribution_;

  this->coord_dim = coord_dim_;
  this->num_local_coords = num_total_coords;

  this->num_global_coords = num_total_coords;
  this->mj_coordinates = mj_coordinates_;


  this->initial_mj_gnos =
    Kokkos::View<mj_gno_t*, device_t>("gids", this->num_local_coords);

  this->num_weights_per_coord = 0;

  this->mj_uniform_weights = Kokkos::View<bool*, Kokkos::HostSpace>(
    "uniform weights", 1);
  this->mj_uniform_weights(0) = true;

  this->mj_weights = Kokkos::View<mj_scalar_t**, device_t>
    ("weights", 1, 1);

  this->mj_uniform_parts =
    Kokkos::View<bool*, Kokkos::HostSpace>("uniform parts", 1);
  this->mj_uniform_parts(0) = true;

  this->set_part_specifications();

  this->allocate_set_work_memory();

  // Do single init
  auto local_part_xadj = this->part_xadj;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
    KOKKOS_LAMBDA (int dummy) {
    local_part_xadj(0) = static_cast<mj_lno_t>(num_selected_coords);
  });

  Kokkos::deep_copy(coordinate_permutations, initial_adjList_output_adjlist);

  mj_part_t current_num_parts = 1;

  Kokkos::View<mj_scalar_t *, device_t> current_cut_coordinates =
    this->all_cut_coordinates;

  mj_part_t future_num_parts = this->total_num_part;

  std::vector<mj_part_t> *future_num_part_in_parts =
    new std::vector<mj_part_t>();
  std::vector<mj_part_t> *next_future_num_parts_in_parts =
    new std::vector<mj_part_t>();
  next_future_num_parts_in_parts->push_back(this->num_global_parts);
  RCP<mj_partBoxVector_t> t1;
  RCP<mj_partBoxVector_t> t2;

  std::vector <uSignedSortItem<int, mj_scalar_t, char>>
    coord_dimension_range_sorted(this->coord_dim);
  uSignedSortItem<int, mj_scalar_t, char> *p_coord_dimension_range_sorted =
    &(coord_dimension_range_sorted[0]);
  std::vector <mj_scalar_t> coord_dim_mins(this->coord_dim);
  std::vector <mj_scalar_t> coord_dim_maxs(this->coord_dim);

  // Need a device counter - how best to allocate?
  // Putting this allocation in the loops is very costly so moved out here.
  Kokkos::View<mj_part_t*, device_t>
    view_rectilinear_cut_count("view_rectilinear_cut_count", 1);
  Kokkos::View<size_t*, device_t>
    view_total_reduction_size("view_total_reduction_size", 1);

  for(int rd = 0; rd < this->recursion_depth; ++rd) {
    // next_future_num_parts_in_parts will be as the size of outnumParts,
    // and this will hold how many more parts that each output part
    // should be divided. this array will also be used to determine the weight
    // ratios of the parts.
    // swap the arrays to use iteratively..
    std::vector<mj_part_t> *tmpPartVect = future_num_part_in_parts;
    future_num_part_in_parts = next_future_num_parts_in_parts;
    next_future_num_parts_in_parts = tmpPartVect;

    // clear next_future_num_parts_in_parts array as
    // getPartitionArrays expects it to be empty.
    next_future_num_parts_in_parts->clear();

    // returns the total number of output parts for this dimension partitioning.
    mj_part_t output_part_count_in_dimension =
      this->update_part_num_arrays(
        future_num_part_in_parts,
        next_future_num_parts_in_parts,
        future_num_parts,
        current_num_parts,
        rd,
        t1,
        t2, num_ranks_per_node);

    // if the number of obtained parts equal to current number of parts,
    // skip this dimension. For example, this happens when 1 is given in
    // the input part array is given. P=4,5,1,2
    if(output_part_count_in_dimension == current_num_parts) {
      tmpPartVect = future_num_part_in_parts;
      future_num_part_in_parts = next_future_num_parts_in_parts;
      next_future_num_parts_in_parts = tmpPartVect;
      continue;
    }

    //convert i to string to be used for debugging purposes.
    std::string istring = std::to_string(rd);

    // alloc Memory to point the indices
    // of the parts in the permutation array.
    this->new_part_xadj = Kokkos::View<mj_lno_t*, device_t>(
      "new part xadj", output_part_count_in_dimension);

    // the index where in the outtotalCounts will be written.

    mj_part_t output_part_index = 0;

    // whatever is written to outTotalCounts will be added with previousEnd
    // so that the points will be shifted.
    mj_part_t output_coordinate_end_index = 0;

    mj_part_t current_work_part = 0;
    mj_part_t current_concurrent_num_parts = 1;

    mj_part_t obtained_part_index = 0;

    // get the coordinate axis along which the partitioning will be done.
    int coordInd = rd % this->coord_dim;

    Kokkos::View<mj_scalar_t *, device_t> mj_current_dim_coords =
      Kokkos::subview(this->mj_coordinates, Kokkos::ALL, coordInd);

    auto host_process_local_min_max_coord_total_weight =
      Kokkos::create_mirror_view(process_local_min_max_coord_total_weight);
    auto host_global_min_max_coord_total_weight =
      Kokkos::create_mirror_view(global_min_max_coord_total_weight);

    // run for all available parts.
    for(; current_work_part < current_num_parts;
      current_work_part += current_concurrent_num_parts) {

      mj_part_t actual_work_part_count = 0;

      // initialization for 1D partitioning.
      // get the min and max coordinates of each part
      // together with the part weights of each part.
      for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {
        mj_part_t current_work_part_in_concurrent_parts =
          current_work_part + kk;

        // if this part wont be partitioned any further
        // dont do any work for this part.
        mj_part_t partition_count = host_num_partitioning_in_current_dim(
          current_work_part_in_concurrent_parts);
        if(partition_count == 1) {
          continue;
        }
        ++actual_work_part_count;
        if(partition_along_longest_dim) {
          auto local_process_local_min_max_coord_total_weight =
            this->process_local_min_max_coord_total_weight;
          for(int coord_traverse_ind = 0;
            coord_traverse_ind < this->coord_dim; ++coord_traverse_ind) {

            Kokkos::View<mj_scalar_t *, device_t> coords =
              Kokkos::subview(this->mj_coordinates, Kokkos::ALL, coord_traverse_ind);

            this->mj_get_local_min_max_coord_totW(
              current_work_part,
              current_concurrent_num_parts,
              coords);

            coord_dimension_range_sorted[coord_traverse_ind].id =
              coord_traverse_ind;
            coord_dimension_range_sorted[coord_traverse_ind].signbit = 1;

            Kokkos::deep_copy(host_process_local_min_max_coord_total_weight,
              process_local_min_max_coord_total_weight);

            coord_dim_mins[coord_traverse_ind] =
              host_process_local_min_max_coord_total_weight(kk);
            coord_dim_maxs[coord_traverse_ind] =
              host_process_local_min_max_coord_total_weight(
                kk + current_concurrent_num_parts);
            coord_dimension_range_sorted[coord_traverse_ind].val =
              host_process_local_min_max_coord_total_weight(
                kk + current_concurrent_num_parts) -
              host_process_local_min_max_coord_total_weight(kk);
          }

          uqSignsort(this->coord_dim, p_coord_dimension_range_sorted);
          coordInd = p_coord_dimension_range_sorted[this->coord_dim - 1].id;
          auto set_min = coord_dim_mins[coordInd];
          auto set_max = coord_dim_maxs[coordInd];
          Kokkos::parallel_for(
            Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
              (0, 1), KOKKOS_LAMBDA (int dummy) {
            local_process_local_min_max_coord_total_weight(kk) = set_min;
            local_process_local_min_max_coord_total_weight(
              kk + current_concurrent_num_parts) = set_max;
          });

          mj_current_dim_coords =
            Kokkos::subview(this->mj_coordinates, Kokkos::ALL, coordInd);
        }
        else {
            Kokkos::View<mj_scalar_t *, device_t> coords =
              Kokkos::subview(this->mj_coordinates, Kokkos::ALL, coordInd);
            this->mj_get_local_min_max_coord_totW(
              current_work_part,
              current_concurrent_num_parts,
              coords);
        }
      }

      // 1D partitioning
      if(actual_work_part_count > 0) {
        // obtain global Min max of the part.
        this->mj_get_global_min_max_coord_totW(
          current_concurrent_num_parts,
          this->process_local_min_max_coord_total_weight,
          this->global_min_max_coord_total_weight);

        // update host copy
        Kokkos::deep_copy(host_global_min_max_coord_total_weight,
          global_min_max_coord_total_weight);

        // represents the total number of cutlines
        // whose coordinate should be determined.
        mj_part_t total_incomplete_cut_count = 0;

        //Compute weight ratios for parts & cuts:
        //e.g., 0.25  0.25 0.5   0.5  0.75  0.75  1.0
        //      part0 cut0 part1 cut1 part2 cut2  part3
        mj_part_t concurrent_part_cut_shift = 0;
        mj_part_t concurrent_part_part_shift = 0;
        for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {
          mj_scalar_t min_coordinate =
            host_global_min_max_coord_total_weight(kk);
          mj_scalar_t max_coordinate = host_global_min_max_coord_total_weight(
              kk + current_concurrent_num_parts);
          mj_scalar_t global_total_weight = host_global_min_max_coord_total_weight(
              kk + 2*current_concurrent_num_parts);

          mj_part_t concurrent_current_part_index = current_work_part + kk;

          mj_part_t partition_count = host_num_partitioning_in_current_dim(
            concurrent_current_part_index);

          Kokkos::View<mj_scalar_t *, device_t> usedCutCoordinate =
            Kokkos::subview(current_cut_coordinates,
              std::pair<mj_lno_t, mj_lno_t>(
                concurrent_part_cut_shift,
                current_cut_coordinates.size()));
          Kokkos::View<mj_scalar_t *, device_t>
            current_target_part_weights =
            Kokkos::subview(target_part_weights,
              std::pair<mj_lno_t, mj_lno_t>(
                concurrent_part_part_shift,
                target_part_weights.size()));

          // shift the usedCutCoordinate array as noCuts.
          concurrent_part_cut_shift += partition_count - 1;
          // shift the partRatio array as noParts.
          concurrent_part_part_shift += partition_count;
          // calculate only if part is not empty,
          // and part will be further partitioend.
          if(partition_count > 1 && min_coordinate <= max_coordinate) {
            // increase allDone by the number of cuts of the current
            // part's cut line number.
            total_incomplete_cut_count += partition_count - 1;

            this->incomplete_cut_count(kk) = partition_count - 1;

            // When num_first_level_parts != 1 we have
            // nonuniform partitioning on the first level, providing
            // requested number of parts (num_first_level_parts) and
            // requested distribution in parts (first_level_distribution)

            // Get the target part weights given a desired distribution
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
                obtained_part_index,
                rd == 0 ? this->num_first_level_parts : 1,
                this->first_level_distribution);

            mj_lno_t coordinate_end_index =
              host_part_xadj(concurrent_current_part_index);
            mj_lno_t coordinate_begin_index =
              (concurrent_current_part_index==0) ? 0 :
                host_part_xadj[concurrent_current_part_index - 1];

            // get the initial estimated part assignments of the coordinates.
            this->set_initial_coordinate_parts(
              max_coordinate,
              min_coordinate,
              coordinate_begin_index, coordinate_end_index,
              this->coordinate_permutations,
              mj_current_dim_coords,
              this->assigned_part_ids,
              partition_count);
          }
          else {
            // e.g., if have fewer coordinates than parts, don't need to do
            // next dim.
            this->incomplete_cut_count(kk) = 0;
          }
          obtained_part_index += partition_count;
        }

        // used imbalance, it is always 0, as it is difficult
        // to estimate a range.
        double used_imbalance = 0;

        // Determine cut lines for k parts here.
        this->mj_env->timerStart(MACRO_TIMERS,
          mj_timer_base_string + "mj_1D_part()");

        this->mj_1D_part(
          mj_current_dim_coords,
          used_imbalance,
          current_work_part,
          current_concurrent_num_parts,
          current_cut_coordinates,
          total_incomplete_cut_count,
          view_rectilinear_cut_count,
          view_total_reduction_size);

        this->mj_env->timerStop(MACRO_TIMERS,
          mj_timer_base_string + "mj_1D_part()");
      }
      else {
        obtained_part_index += current_concurrent_num_parts;
      }
      // create part chunks
      {
        mj_part_t output_array_shift = 0;
        mj_part_t cut_shift = 0;
        size_t tlr_shift = 0;
        size_t partweight_array_shift = 0;

        for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {
          mj_part_t current_concurrent_work_part = current_work_part + kk;

          mj_part_t num_parts = host_num_partitioning_in_current_dim(
            current_concurrent_work_part);

          // if the part is empty, skip the part.
          int coordinateA_bigger_than_coordinateB =
            host_global_min_max_coord_total_weight(kk) >
            host_global_min_max_coord_total_weight(
              kk + current_concurrent_num_parts);

          if((num_parts != 1) && coordinateA_bigger_than_coordinateB) {
            // we still need to write the begin and end point of the empty part.
            // simply set it zero, the array indices will be shifted later
            auto local_new_part_xadj = this->new_part_xadj;
            Kokkos::parallel_for(
              Kokkos::RangePolicy<typename mj_node_t::execution_space,
                mj_part_t> (0, num_parts), KOKKOS_LAMBDA(mj_part_t jj) {
                local_new_part_xadj(
                  output_part_index + output_array_shift + jj) = 0;
            });

            cut_shift += num_parts - 1;
            tlr_shift += (4 *(num_parts - 1) + 1);
            output_array_shift += num_parts;
            partweight_array_shift += (2 * (num_parts - 1) + 1);
            continue;
          }
          mj_lno_t coordinate_end =
            host_part_xadj(current_concurrent_work_part);
          mj_lno_t coordinate_begin =
            current_concurrent_work_part==0 ? 0 :
              host_part_xadj(current_concurrent_work_part-1);

          Kokkos::View<mj_scalar_t *, device_t>
            current_concurrent_cut_coordinate =
            Kokkos::subview(current_cut_coordinates,
              std::pair<mj_lno_t, mj_lno_t>(
                cut_shift,
                current_cut_coordinates.size()));
          Kokkos::View<mj_scalar_t *, device_t>
            used_local_cut_line_weight_to_left =
            Kokkos::subview(process_cut_line_weight_to_put_left,
              std::pair<mj_lno_t, mj_lno_t>(
                cut_shift,
                process_cut_line_weight_to_put_left.size()));

          this->thread_part_weight_work =
            Kokkos::subview(
              this->thread_part_weights,
              std::pair<mj_lno_t, mj_lno_t>(
                partweight_array_shift,
                this->thread_part_weights.size()));

          if(num_parts > 1) {
            // Rewrite the indices based on the computed cuts.
            Kokkos::View<mj_lno_t *, device_t> subview_new_part_xadj =
              Kokkos::subview(this->new_part_xadj,
                std::pair<mj_lno_t, mj_lno_t>(
                  output_part_index + output_array_shift,
                  this->new_part_xadj.size()));

            this->create_consistent_chunks(
              num_parts,
              mj_current_dim_coords,
              current_concurrent_cut_coordinate,
              coordinate_begin,
              coordinate_end,
              used_local_cut_line_weight_to_left,
              subview_new_part_xadj,
              coordInd,
              partition_along_longest_dim,
              p_coord_dimension_range_sorted);
          }
          else {
            // if this part is partitioned into 1 then just copy
            // the old values.
            mj_lno_t part_size = coordinate_end - coordinate_begin;

            auto local_new_part_xadj = this->new_part_xadj;
            Kokkos::parallel_for(
              Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
                (0, 1), KOKKOS_LAMBDA (int dummy) {
              local_new_part_xadj(output_part_index + output_array_shift)
                = part_size;
            });

            auto subview_new_coordinate_permutations =
              Kokkos::subview(this->new_coordinate_permutations,
                std::pair<mj_lno_t, mj_lno_t>(
                  coordinate_begin,
                  coordinate_begin + part_size));
            auto subview_coordinate_permutations =
              Kokkos::subview(this->coordinate_permutations,
                std::pair<mj_lno_t, mj_lno_t>(
                  coordinate_begin,
                  coordinate_begin + part_size));
            Kokkos::deep_copy(subview_new_coordinate_permutations,
              subview_coordinate_permutations);
          }

          cut_shift += num_parts - 1;
          tlr_shift += (4 *(num_parts - 1) + 1);
          output_array_shift += num_parts;
          partweight_array_shift += (2 * (num_parts - 1) + 1);
        }

        // shift cut coordinates so that all cut coordinates are stored.
        // current_cut_coordinates += cutShift;

        // getChunks from coordinates partitioned the parts and
        // wrote the indices as if there were a single part.
        // now we need to shift the beginning indices.
        for(mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk) {
          mj_part_t num_parts =
            host_num_partitioning_in_current_dim(current_work_part + kk);
          auto local_new_part_xadj = this->new_part_xadj;
          auto local_mj_current_dim_coords = mj_current_dim_coords;
          auto local_new_coordinate_permutations =
            new_coordinate_permutations;
          Kokkos::parallel_for(
            Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t> (
            0, num_parts), KOKKOS_LAMBDA (mj_part_t ii) {
            //shift it by previousCount
            local_new_part_xadj(output_part_index+ii) +=
              output_coordinate_end_index;

            if(ii % 2 == 1) {
              mj_lno_t coordinate_end =
                local_new_part_xadj(output_part_index+ii);
              mj_lno_t coordinate_begin =
                local_new_part_xadj(output_part_index);

              for(mj_lno_t task_traverse = coordinate_begin;
                task_traverse < coordinate_end; ++task_traverse) {
                mj_lno_t l = local_new_coordinate_permutations(task_traverse);
                //MARKER: FLIPPED ZORDER BELOW
                local_mj_current_dim_coords(l) = -local_mj_current_dim_coords(l);
              }
            }
          });

          // increase the previous count by current end.
          mj_part_t get_single;
          Kokkos::parallel_reduce("Read new_part_xadj",
            Kokkos::RangePolicy<typename mj_node_t::execution_space, int>(0, 1),
            KOKKOS_LAMBDA(int dummy, mj_part_t & set_single) {
            set_single = local_new_part_xadj(output_part_index + num_parts - 1);
          }, get_single);;

          output_coordinate_end_index = get_single;
          // increase the current out.
          output_part_index += num_parts;
        }
      }
    }

    // end of this partitioning dimension
    // set the current num parts for next dim partitioning
    current_num_parts = output_part_count_in_dimension;

    //swap the coordinate permutations for the next dimension.
    Kokkos::View<mj_lno_t *, device_t> tmp = this->coordinate_permutations;
    this->coordinate_permutations = this->new_coordinate_permutations;
    this->new_coordinate_permutations = tmp;

    this->part_xadj = this->new_part_xadj;
    this->host_part_xadj = Kokkos::create_mirror_view(part_xadj);
    Kokkos::deep_copy(host_part_xadj, part_xadj); // keep in sync
    this->new_part_xadj = Kokkos::View<mj_lno_t*, device_t>("empty", 0);
  }

  Kokkos::deep_copy(initial_adjList_output_adjlist, coordinate_permutations);

  // Return output_xadj in CSR format
  output_xadj[0] = 0;
  for(size_t i = 0; i < this->num_global_parts ; ++i) {
    output_xadj[i+1] = host_part_xadj(i);
  }

  delete future_num_part_in_parts;
  delete next_future_num_parts_in_parts;
}

/*! \brief Function returns the part boxes stored
 * returns null if boxes are not stored, and prints warning mesage.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
RCP<typename AlgMJ
  <mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t,mj_node_t>::mj_partBox_t>
AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t, mj_node_t>::
  get_global_box() const
{
  return this->global_box;
}

/*! \brief Function call, if the part boxes are intended to be kept.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t,
  mj_node_t>::set_to_keep_part_boxes()
{
  this->mj_keep_part_boxes = true;
}

/* \brief Either the mj array (part_no_array) or num_global_parts should be
 * provided in the input. part_no_array takes
 * precedence if both are provided.
 * Depending on these parameters, total cut/part number,
 * maximum part/cut number along a dimension, estimated number of reduceAlls,
 * and the number of parts before the last dimension is calculated.
 * */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  set_part_specifications()
{
  this->total_num_cut = 0; //how many cuts will be totally
  this->total_num_part = 1;    //how many parts will be totally
  this->max_num_part_along_dim = 0; // maximum part count along a dimension.
  this->total_dim_num_reduce_all = 0; // estimate on #reduceAlls can be done.
  this->last_dim_num_part = 1; //max no of parts that might occur
  //during the partition before the
  //last partitioning dimension.
  this->max_num_cut_along_dim = 0;
  this->max_num_total_part_along_dim = 0;

  if(this->part_no_array.size()) {
    auto local_recursion_depth = this->recursion_depth;

    this->total_dim_num_reduce_all =
      this->total_num_part * this->recursion_depth;

    this->total_num_part = 1;
    for(int i = 0; i < local_recursion_depth; ++i) {
      this->total_num_part *= this->part_no_array(i);
    }

    mj_part_t track_max = 0;
    for(int i = 0; i < local_recursion_depth; ++i) {
      if(part_no_array(i) > track_max) {
        track_max = this->part_no_array(i);
      };
    }

    this->last_dim_num_part = this->total_num_part /
      this->part_no_array(local_recursion_depth-1);

    this->max_num_part_along_dim = track_max;
    this->num_global_parts = this->total_num_part;
  } else {
    mj_part_t future_num_parts = this->num_global_parts;

    // If using nonuniform first level partitioning.
    // initial value max_num_part_along_dim == num_first_level_parts
    if (this->first_level_distribution.size() != 0 &&
        this->num_first_level_parts > 1) {
      this->max_num_part_along_dim = this->num_first_level_parts;
    }

    // we need to calculate the part numbers now, to determine
    // the maximum along the dimensions.
    for(int rd = 0; rd < this->recursion_depth; ++rd) {
      mj_part_t maxNoPartAlongI = 0;
      mj_part_t nfutureNumParts = 0;

      // Nonuniform first level partitioning sets part specificiations for
      // rd == 0 only, given requested num of parts and distribution in parts
      // for the first level.
      if (rd == 0 &&
          this->first_level_distribution.size() != 0 &&
          this->num_first_level_parts > 1) {

        maxNoPartAlongI = this->num_first_level_parts;
        this->max_num_part_along_dim = this->num_first_level_parts;

        mj_part_t sum_first_level_dist = 0;
        mj_part_t max_part = 0;

        // Cumulative sum of distribution of parts and size of largest part
        for (int i = 0; i < this->num_first_level_parts; ++i) {
          sum_first_level_dist += this->first_level_distribution(i);
          if (this->first_level_distribution(i) > max_part)
            max_part = this->first_level_distribution(i);
        }

        // Total parts in largest nonuniform superpart from
        // first level partitioning
        nfutureNumParts =
          this->num_global_parts * max_part / sum_first_level_dist;
      }
      // Standard uniform partitioning this level
      else {
        maxNoPartAlongI = this->get_part_count(future_num_parts,
          1.0f / (this->recursion_depth - rd));
        if (maxNoPartAlongI > this->max_num_part_along_dim)
          this->max_num_part_along_dim = maxNoPartAlongI;
        nfutureNumParts = future_num_parts / maxNoPartAlongI;
        if (future_num_parts % maxNoPartAlongI) {
                ++nfutureNumParts;
        }
      }
      future_num_parts = nfutureNumParts;
    }
    this->total_num_part = this->num_global_parts;

    if(this->divide_to_prime_first) {
      this->total_dim_num_reduce_all = this->num_global_parts * 2;
      this->last_dim_num_part = this->num_global_parts;
    }
    else {
      //this is the lower bound.
      //estimate reduceAll Count here.
      //we find the upperbound instead.
      size_t p = 1;
      for(int i = 0; i < this->recursion_depth; ++i) {
        this->total_dim_num_reduce_all += p;
        p *= this->max_num_part_along_dim;
      }

      if(p / this->max_num_part_along_dim > this->num_global_parts) {
        this->last_dim_num_part = this->num_global_parts;
      }
      else {
        this->last_dim_num_part = p / this->max_num_part_along_dim;
      }
    }
  }

  this->total_num_cut = this->total_num_part - 1;
  this->max_num_cut_along_dim = this->max_num_part_along_dim - 1;
  this->max_num_total_part_along_dim = this->max_num_part_along_dim +
    size_t(this->max_num_cut_along_dim);
  // maxPartNo is P, maxCutNo = P-1, matTotalPartcount = 2P-1

  // refine the concurrent part count, if it is given bigger than the maximum
  // possible part count.
  if(this->max_concurrent_part_calculation > this->last_dim_num_part) {
    if(this->mj_problemComm->getRank() == 0) {
      std::cerr << "Warning: Concurrent part count (" <<
        this->max_concurrent_part_calculation <<
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
  typename mj_part_t, typename mj_node_t>
inline mj_part_t AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  get_part_count(mj_part_t num_total_future, double root)
{
  double fp = pow(num_total_future, root);
  mj_part_t ip = mj_part_t(fp);
  if(fp - ip < std::numeric_limits<float>::epsilon() * 100) {
    return ip;
  }
  else {
    return ip  + 1;
  }
}

/* \brief Function returns how many parts that will be obtained after this
 * dimension partitioning. It sets how many parts each current part will be
 * partitioned into in this dimension to device_num_partitioning_in_current_dim
 * view, sets how many total future parts each obtained part will be
 * partitioned into in next_future_num_parts_in_parts vector. If part boxes are
 * kept, then sets initializes the output_part_boxes as its ancestor.
 * \param future_num_part_in_parts: input, how many future parts each current
 * part will be partitioned into.
 * \param next_future_num_parts_in_parts: output, how many future parts each
 * obtained part will be partitioned into.
 * \param future_num_parts: output, max number of future parts that will be
 * obtained from a single
 * \param current_num_parts: input, how many parts are there currently.
 * \param current_iteration: input, current dimension iteration number.
 * \param input_part_boxes: input, if boxes are kept, current boxes.
 * \param output_part_boxes: output, if boxes are kept, the initial box
 * boundaries for obtained parts.
 * \param atomic_part_count  DOCWORK: Documentation
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
mj_part_t AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  update_part_num_arrays(
  std::vector<mj_part_t> *future_num_part_in_parts,
  std::vector<mj_part_t> *next_future_num_parts_in_parts,
  mj_part_t &future_num_parts,
  mj_part_t current_num_parts,
  int current_iteration,
  RCP<mj_partBoxVector_t> input_part_boxes,
  RCP<mj_partBoxVector_t> output_part_boxes,
  mj_part_t atomic_part_count)
{
  std::vector<mj_part_t> num_partitioning_in_current_dim;

  // how many parts that will be obtained after this dimension.
  mj_part_t output_num_parts = 0;
  if(this->part_no_array.size()) {
    // when the partNo array is provided as input,
    // each current partition will be partition to the same number of parts.
    // we dont need to use the future_num_part_in_parts vector in this case.
    mj_part_t current_part_no_array =
      this->part_no_array(current_iteration);

    if(current_part_no_array < 1) {
      std::cout << "Current recursive iteration: " << current_iteration <<
        " part_no_array[" << current_iteration << "] is given as:" <<
        current_part_no_array << std::endl;
      std::terminate();
    }
    if(current_part_no_array == 1) {
      return current_num_parts;
    }

    // If using part_no_array, ensure compatibility with num_first_level_parts.
    if (this->first_level_distribution.size() != 0 &&
      current_iteration == 0 &&
      current_part_no_array != this->num_first_level_parts) {
      std::cout << "Current recursive iteration: " << current_iteration
        << " part_no_array[" << current_iteration << "] is given as: " <<
        current_part_no_array << " and contradicts num_first_level_parts: " <<
        this->num_first_level_parts << std::endl;
      std::terminate();
    }

    for(mj_part_t ii = 0; ii < current_num_parts; ++ii) {
      num_partitioning_in_current_dim.push_back(current_part_no_array);
    }

/*
    std::cout << "\n\nme: " << this->myRank << " current_iteration: " <<
      current_iteration << " current_num_parts: " <<
      current_num_parts << "\n\n";

    std::cout << "\n\nnum_partitioning_in_current_dim[0]: " <<
      num_partitioning_in_current_dim[0] << "\n\n";

    std::cout << "\n\nfuture_num_parts: " << future_num_parts
     << " num_partitioning_in_current_dim[0]: " <<
     num_partitioning_in_current_dim[0] << " " <<
     future_num_parts / num_partitioning_in_current_dim[0] << "\n\n";
*/

    future_num_parts /= num_partitioning_in_current_dim[0];
    output_num_parts = current_num_parts *
      num_partitioning_in_current_dim[0];
    if(this->mj_keep_part_boxes) {
      for(mj_part_t k = 0; k < current_num_parts; ++k) {
        //initialized the output boxes as its ancestor.
        for(mj_part_t j = 0; j <
          num_partitioning_in_current_dim[0]; ++j) {
          output_part_boxes->push_back((*input_part_boxes)[k]);
        }
      }
    }

    // set the how many more parts each part will be divided.
    // this is obvious when partNo array is provided as input.
    // however, fill this so weights will be calculated according to this array.
    for(mj_part_t ii = 0; ii < output_num_parts; ++ii) {
      next_future_num_parts_in_parts->push_back(future_num_parts);
    }
  }
  else {
    // if partNo array is not provided as input, future_num_part_in_parts
    // holds how many parts each part should be divided. Initially it holds a
    // single number equal to the total number of global parts.

    // calculate the future_num_parts from beginning,
    // since each part might be divided into different number of parts.
    future_num_parts = 1;

    // cout << "i:" << i << std::endl;
    for(mj_part_t ii = 0; ii < current_num_parts; ++ii) {
      // get how many parts a part should be divided.
      mj_part_t future_num_parts_of_part_ii = (*future_num_part_in_parts)[ii];

      // get the ideal number of parts that is close to the
      // (recursion_depth - i) root of the future_num_parts_of_part_ii.
      mj_part_t num_partitions_in_current_dim =
        this->get_part_count(future_num_parts_of_part_ii,
          1.0 / (this->recursion_depth - current_iteration)
                                    );
      if(num_partitions_in_current_dim > this->max_num_part_along_dim) {
        std::cerr << "ERROR: maxPartNo calculation is wrong."
          " num_partitions_in_current_dim: "
          << num_partitions_in_current_dim <<  " this->max_num_part_along_dim: "
          << this->max_num_part_along_dim <<
          " this->recursion_depth: " << this->recursion_depth <<
          " current_iteration:" << current_iteration <<
          " future_num_parts_of_part_ii: " << future_num_parts_of_part_ii <<
          " might need to fix max part no calculation for "
          "largest_prime_first partitioning." <<
          std::endl;
        std::terminate();
      }
      // add this number to vector_num_partitioning_in_current_dim vector.
      // num_partitioning_in_current_dim.push_back(num_partitions_in_current_dim);
      // mj_part_t largest_prime_factor = num_partitions_in_current_dim;

      // Update part num arrays when on current_iteration == 0 and
      // using nonuniform first level partitioning
      // with requested num parts (num_first_level_parts) and
      // a requested distribution in parts (first_level_distribution).
      if (current_iteration == 0 &&
          this->first_level_distribution.size() != 0 &&
          this->num_first_level_parts > 1) {
        // Only 1 current part to begin and partitions into
        // num_first_level_parts many parts
        num_partitioning_in_current_dim.push_back(this->num_first_level_parts);

        // The output number of parts from first level partitioning
        output_num_parts = this->num_first_level_parts;

        // Remaining parts left to partition for all future levels
        future_num_parts /= this->num_first_level_parts;

        mj_part_t max_part = 0;
        mj_part_t sum_first_level_dist = 0;

        // Cumulative sum of distribution of first level parts
        // and size of largest first level part
        for (int i = 0; i < this->num_first_level_parts; ++i) {
          sum_first_level_dist += this->first_level_distribution(i);

          if (this->first_level_distribution(i) > max_part)
            max_part = this->first_level_distribution(i);
        }

        // Maximum # of remaining parts left to partition for all future levels
        future_num_parts = this->num_global_parts * max_part / sum_first_level_dist;

        // Number of parts remaining left to partition for each future_part
        // The sum must exactly equal global_num_parts
        for (int i = 0; i < this->num_first_level_parts; ++i) {
          next_future_num_parts_in_parts->push_back(this->first_level_distribution(i) *
              this->num_global_parts / sum_first_level_dist);
        }
      }
      else if (this->divide_to_prime_first) {
        // Add this number to num_partitioning_in_current_dim vector.
        num_partitioning_in_current_dim.push_back(num_partitions_in_current_dim);

        mj_part_t largest_prime_factor = num_partitions_in_current_dim;

        //increase the output number of parts.
        output_num_parts += num_partitions_in_current_dim;

        if (future_num_parts_of_part_ii == atomic_part_count ||
            future_num_parts_of_part_ii % atomic_part_count != 0) {
          atomic_part_count = 1;
        }

        largest_prime_factor =
          this->find_largest_prime_factor(future_num_parts_of_part_ii / atomic_part_count);

        // We divide to  num_partitions_in_current_dim. But we adjust the weights
        // based on largest prime/ if num_partitions_in_current_dim = 2,
        // largest prime = 5 --> we divide to 2 parts with weights 3x and 2x.
        // if the largest prime is less than part count, we use the part count
        // so that we divide uniformly.
        if (largest_prime_factor < num_partitions_in_current_dim) {
          largest_prime_factor = num_partitions_in_current_dim;
        }
        //ideal number of future partitions for each part.
        mj_part_t ideal_num_future_parts_in_part =
          (future_num_parts_of_part_ii / atomic_part_count) / largest_prime_factor;
        //if num_partitions_in_current_dim = 2, largest prime = 5 then ideal weight is 2x
        mj_part_t ideal_prime_scale = largest_prime_factor / num_partitions_in_current_dim;

/*
        std::cout << "\ncurrent num part: " << ii
          << " largest_prime_factor: " << largest_prime_factor
          << " To Partition: " << future_num_parts_of_part_ii << "\n\n";
*/

        for (mj_part_t iii = 0; iii < num_partitions_in_current_dim; ++iii) {
          //if num_partitions_in_current_dim = 2, largest prime = 5 then ideal weight is 2x
          mj_part_t my_ideal_primescale = ideal_prime_scale;
          //left over weighs. Left side is adjusted to be 3x, right side stays as 2x
          if (iii < (largest_prime_factor) % num_partitions_in_current_dim) {
            ++my_ideal_primescale;
          }
          //scale with 'x';
          mj_part_t num_future_parts_for_part_iii =
            ideal_num_future_parts_in_part * my_ideal_primescale;

           //if there is a remainder in the part increase the part weight.
          if (iii < (future_num_parts_of_part_ii / atomic_part_count) % largest_prime_factor) {
            //if not uniform, add 1 for the extra parts.
            ++num_future_parts_for_part_iii;
          }

          next_future_num_parts_in_parts->push_back(num_future_parts_for_part_iii * atomic_part_count);

          //if part boxes are stored, initialize the box of the parts as the ancestor.
          if (this->mj_keep_part_boxes) {
            output_part_boxes->push_back((*input_part_boxes)[ii]);
          }

          //set num future_num_parts to maximum in this part.
          if (num_future_parts_for_part_iii > future_num_parts)
            future_num_parts = num_future_parts_for_part_iii;

        }
      }
      else {
        // Add this number to num_partitioning_in_current_dim vector.
        num_partitioning_in_current_dim.push_back(num_partitions_in_current_dim);

        //increase the output number of parts.
        output_num_parts += num_partitions_in_current_dim;

        if((future_num_parts_of_part_ii == atomic_part_count) ||
          (future_num_parts_of_part_ii % atomic_part_count != 0)) {
          atomic_part_count = 1;
        }
        //ideal number of future partitions for each part.
        mj_part_t ideal_num_future_parts_in_part =
          (future_num_parts_of_part_ii / atomic_part_count) /
          num_partitions_in_current_dim;
        for(mj_part_t iii = 0; iii < num_partitions_in_current_dim; ++iii) {
          mj_part_t num_future_parts_for_part_iii =
            ideal_num_future_parts_in_part;

          //if there is a remainder in the part increase the part weight.
          if(iii < (future_num_parts_of_part_ii / atomic_part_count) %
            num_partitions_in_current_dim) {
            // if not uniform, add 1 for the extra parts.
            ++num_future_parts_for_part_iii;
          }

          next_future_num_parts_in_parts->push_back(
            num_future_parts_for_part_iii * atomic_part_count);

          // if part boxes are stored, initialize the box of the parts as
          // the ancestor.
          if(this->mj_keep_part_boxes) {
            output_part_boxes->push_back((*input_part_boxes)[ii]);
          }
          //set num future_num_parts to maximum in this part.
          if(num_future_parts_for_part_iii > future_num_parts)
            future_num_parts = num_future_parts_for_part_iii;
        }
      }
    }
  }
  // move temp std::vector to host view
  device_num_partitioning_in_current_dim = Kokkos::View<
    mj_part_t*, device_t>("test", num_partitioning_in_current_dim.size());
  host_num_partitioning_in_current_dim =
    Kokkos::create_mirror_view(device_num_partitioning_in_current_dim);
  for(size_t n = 0; n < num_partitioning_in_current_dim.size(); ++n) {
    host_num_partitioning_in_current_dim(n) =
      num_partitioning_in_current_dim[n];
  }
  // setup device equivalent - this data is used on host and device and it's
  // more efficient to just setup array on both sides now rather than copy
  // values as needed later.
  Kokkos::deep_copy(device_num_partitioning_in_current_dim,
    host_num_partitioning_in_current_dim);
  return output_num_parts;
}

/* \brief Allocates and initializes the work memory that will be used by MJ.
 * */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  allocate_set_work_memory()
{
  // Throughout the partitioning execution,
  // instead of the moving the coordinates, hold a permutation array for parts.
  // coordinate_permutations holds the current permutation.
  this->coordinate_permutations = Kokkos::View<mj_lno_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("coordinate_permutations"),
    this->num_local_coords);
  auto local_coordinate_permutations = coordinate_permutations;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t> (
    0, this->num_local_coords), KOKKOS_LAMBDA (mj_lno_t i) {
      local_coordinate_permutations(i) = i;
  });

  // new_coordinate_permutations holds the current permutation.
  this->new_coordinate_permutations = Kokkos::View<mj_lno_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("num_local_coords"),
    this->num_local_coords);

  this->assigned_part_ids = Kokkos::View<mj_part_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("assigned parts"), 0);
  if(this->num_local_coords > 0) {
    this->assigned_part_ids = Kokkos::View<mj_part_t*, device_t>(
      Kokkos::ViewAllocateWithoutInitializing("assigned part ids"),
      this->num_local_coords);
  }

  // single partition starts at index-0, and ends at numLocalCoords
  // inTotalCounts array holds the end points in coordinate_permutations array
  // for each partition. Initially sized 1, and single element is set to
  // numLocalCoords.
  this->part_xadj = Kokkos::View<mj_lno_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("part xadj"), 1);
  this->host_part_xadj = Kokkos::create_mirror_view(part_xadj);
  host_part_xadj(0) = num_local_coords;
  Kokkos::deep_copy(this->part_xadj, host_part_xadj);

  // the ends points of the output, this is allocated later.
  this->new_part_xadj = Kokkos::View<mj_lno_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("empty"), 0);

  // only store this much if cuts are needed to be stored.
  this->all_cut_coordinates = Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("all cut coordinates"),
    this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

  // how much weight percentage should a MPI put left side of the each cutline
  this->process_cut_line_weight_to_put_left = Kokkos::View<mj_scalar_t*,
    device_t>(Kokkos::ViewAllocateWithoutInitializing("empty"), 0);

  // how much weight percentage should each thread in MPI put left side of
  // each outline
  this->thread_cut_line_weight_to_put_left =
    Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("empty"), 0);

  if(this->distribute_points_on_cut_lines) {
    this->process_cut_line_weight_to_put_left =
      Kokkos::View<mj_scalar_t *, device_t>(
      Kokkos::ViewAllocateWithoutInitializing(
        "process_cut_line_weight_to_put_left"),
      this->max_num_cut_along_dim * this->max_concurrent_part_calculation);
    this->thread_cut_line_weight_to_put_left =
      Kokkos::View<mj_scalar_t *, device_t>(
      Kokkos::ViewAllocateWithoutInitializing(
        "thread_cut_line_weight_to_put_left"),
      this->max_num_cut_along_dim);
    this->process_rectilinear_cut_weight =
      Kokkos::View<mj_scalar_t *, device_t>(
      Kokkos::ViewAllocateWithoutInitializing("process_rectilinear_cut_weight"),
      this->max_num_cut_along_dim);
    this->global_rectilinear_cut_weight =
      Kokkos::View<mj_scalar_t *, device_t>(
      Kokkos::ViewAllocateWithoutInitializing("global_rectilinear_cut_weight"),
      this->max_num_cut_along_dim);
  }

  // work array to manipulate coordinate of cutlines in different iterations.
  // necessary because previous cut line information is used for determining
  // the next cutline information. therefore, cannot update the cut work array
  // until all cutlines are determined.
  this->cut_coordinates_work_array =
    Kokkos::View<mj_scalar_t *, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("cut_coordinates_work_array"),
    this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

  // cumulative part weight array.
  this->target_part_weights = Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("target_part_weights"),
    this->max_num_part_along_dim * this->max_concurrent_part_calculation);

  // upper bound coordinate of a cut line
  this->cut_upper_bound_coordinates =
    Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("cut_upper_bound_coordinates"),
    this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

  // lower bound coordinate of a cut line
  this->cut_lower_bound_coordinates =
    Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("cut_lower_bound_coordinates"),
    this->max_num_cut_along_dim* this->max_concurrent_part_calculation);

  // lower bound weight of a cut line
  this->cut_lower_bound_weights =
    Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("cut_lower_bound_weights"),
    this->max_num_cut_along_dim* this->max_concurrent_part_calculation);

  //upper bound weight of a cut line
  this->cut_upper_bound_weights =
    Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("cut_upper_bound_weights"),
    this->max_num_cut_along_dim* this->max_concurrent_part_calculation);

  // combined array to exchange the min and max coordinate,
  // and total weight of part.
  this->process_local_min_max_coord_total_weight =
    Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing(
      "process_local_min_max_coord_total_weight"),
    3 * this->max_concurrent_part_calculation);

  // global combined array with the results for min, max and total weight.
  this->global_min_max_coord_total_weight =
    Kokkos::View<mj_scalar_t*, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("global_min_max_coord_total_weight"),
    3 * this->max_concurrent_part_calculation);

  // is_cut_line_determined is used to determine if a cutline is
  // determined already. If a cut line is already determined, the next
  // iterations will skip this cut line.
  this->is_cut_line_determined = Kokkos::View<bool *, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("is_cut_line_determined"),
    this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

  // incomplete_cut_count count holds the number of cutlines that have not
  // been finalized for each part when concurrentPartCount>1, using this
  // information, if incomplete_cut_count[x]==0, then no work is done for
  // this part.
  this->device_incomplete_cut_count = Kokkos::View<mj_part_t *, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("device_incomplete_cut_count"),
    this->max_concurrent_part_calculation);
  this->incomplete_cut_count =
    Kokkos::create_mirror_view(device_incomplete_cut_count);

  // local part weights of each thread.
  this->thread_part_weights = Kokkos::View<double *, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("thread_part_weights"),
    this->max_num_total_part_along_dim * this->max_concurrent_part_calculation);

  this->thread_cut_left_closest_point = Kokkos::View<mj_scalar_t *, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("thread_cut_left_closest_point"),
    this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

  // thread_cut_right_closest_point to hold the closest coordinate to a
  // cutline from right (for each thread)
  this->thread_cut_right_closest_point = Kokkos::View<mj_scalar_t *, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("thread_cut_right_closest_point"),
    this->max_num_cut_along_dim * this->max_concurrent_part_calculation);

  // to store how many points in each part a thread has.
  this->thread_point_counts = Kokkos::View<mj_lno_t *, device_t>(
    Kokkos::ViewAllocateWithoutInitializing("thread_point_counts"),
    this->max_num_part_along_dim);

  // for faster communication, concatanation of
  // totalPartWeights sized 2P-1, since there are P parts and P-1 cut lines
  // leftClosest distances sized P-1, since P-1 cut lines
  // rightClosest distances size P-1, since P-1 cut lines.
  this->total_part_weight_left_right_closests =
    Kokkos::View<mj_scalar_t*, device_t>(
      Kokkos::ViewAllocateWithoutInitializing(
        "total_part_weight_left_right_closests"),
      (this->max_num_total_part_along_dim + this->max_num_cut_along_dim * 2) *
      this->max_concurrent_part_calculation);

  this->global_total_part_weight_left_right_closests =
    Kokkos::View<mj_scalar_t*, device_t>(
      Kokkos::ViewAllocateWithoutInitializing(
        "global_total_part_weight_left_right_closests"),
      (this->max_num_total_part_along_dim +
      this->max_num_cut_along_dim * 2) * this->max_concurrent_part_calculation);

  this->current_mj_gnos = Kokkos::View<mj_gno_t*, device_t>(
      Kokkos::ViewAllocateWithoutInitializing("gids"), num_local_coords);

  this->owner_of_coordinate = Kokkos::View<int *, Kokkos::HostSpace>(
    Kokkos::ViewAllocateWithoutInitializing("owner_of_coordinate"),
    num_local_coords);

  // changes owners back to host - so we don't run them on device
  // this improves migration code but means we have to serial init here.
  // Note we might allow this to be OpenMP when available even for CUDA.
  Kokkos::deep_copy(owner_of_coordinate, myActualRank);

  auto local_current_mj_gnos = current_mj_gnos;
  auto local_initial_mj_gnos = initial_mj_gnos;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
    (0, num_local_coords), KOKKOS_LAMBDA (mj_lno_t j) {
    local_current_mj_gnos(j) = local_initial_mj_gnos(j);
  });
}

/* \brief compute the global bounding box
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t,
  mj_node_t>::compute_global_box()
{
  //local min coords
  mj_scalar_t *mins = new mj_scalar_t[this->coord_dim];
  //global min coords
  mj_scalar_t *gmins = new mj_scalar_t[this->coord_dim];
  //local max coords
  mj_scalar_t *maxs = new mj_scalar_t[this->coord_dim];
  //global max coords
  mj_scalar_t *gmaxs = new mj_scalar_t[this->coord_dim];

  auto local_mj_coordinates = this->mj_coordinates;

  // If we are only doing 2 parts then we don't need these values
  // for y and z. Init them all to 0 first
  for(int i = 0; i < this->coord_dim; ++i) {
    mins[i] = 0;
    maxs[i] = 0;
  }

  for(int i = 0; i < std::min(this->recursion_depth, this->coord_dim); ++i) {
    Kokkos::parallel_reduce("MinReduce",
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
        (0, this->num_local_coords),
      KOKKOS_LAMBDA(mj_lno_t j, mj_scalar_t & running_min) {
      if(local_mj_coordinates(j,i) < running_min) {
        running_min = local_mj_coordinates(j,i);
      }
    }, Kokkos::Min<mj_scalar_t>(mins[i]));
    Kokkos::parallel_reduce("MaxReduce",
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
        (0, this->num_local_coords),
      KOKKOS_LAMBDA(mj_lno_t j, mj_scalar_t & running_max) {
      if(local_mj_coordinates(j,i) > running_max) {
        running_max = local_mj_coordinates(j,i);
      }
    }, Kokkos::Max<mj_scalar_t>(maxs[i]));
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
  delete [] mins;
  delete [] gmins;
  delete [] maxs;
  delete [] gmaxs;
}

/* \brief for part communication we keep track of the box boundaries.
 * This is performed when either asked specifically, or when geometric mapping
 * is performed afterwards.
 * This function initializes a single box with all global min, max coordinates.
 * \param initial_partitioning_boxes the input and output vector for boxes.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t,
  mj_node_t>::init_part_boxes(
  RCP<mj_partBoxVector_t> & initial_partitioning_boxes)
{
    mj_partBox_t tmp_box(*global_box);
    initial_partitioning_boxes->push_back(tmp_box);
}

/*! \brief Function to determine the local minimum and maximum coordinate, and
 * local total weight in the given set of local points.
 * DOCWORK: Document parameters
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
          typename mj_part_t,
          typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_get_local_min_max_coord_totW(
  mj_part_t current_work_part,
  mj_part_t current_concurrent_num_parts,
  Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords)
{
  auto local_coordinate_permutations = this->coordinate_permutations;
  auto local_process_local_min_max_coord_total_weight =
    this->process_local_min_max_coord_total_weight;
  auto local_mj_weights = this->mj_weights;

  bool bUniformWeights = mj_uniform_weights(0);

  for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {

    mj_part_t concurrent_current_part = current_work_part + kk;
    mj_lno_t coordinate_begin_index = concurrent_current_part == 0 ? 0 :
      host_part_xadj(concurrent_current_part - 1);
    mj_lno_t coordinate_end_index =
      host_part_xadj(concurrent_current_part);

    mj_scalar_t my_min_coord = 0;
    mj_scalar_t my_max_coord = 0;
    mj_scalar_t my_total_weight;
    //if the part is empty.
    //set the min and max coordinates as reverse.
    if(coordinate_begin_index >= coordinate_end_index)
    {
      my_min_coord = std::numeric_limits<mj_scalar_t>::max();
      my_max_coord = -std::numeric_limits<mj_scalar_t>::max();
      my_total_weight = 0;
    }
    else {
      // get min
      Kokkos::parallel_reduce("get min",
        Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
          (coordinate_begin_index, coordinate_end_index),
        KOKKOS_LAMBDA (mj_lno_t j, mj_scalar_t & running_min) {
        int i = local_coordinate_permutations(j);
        if(mj_current_dim_coords(i) < running_min)
          running_min = mj_current_dim_coords(i);
      }, Kokkos::Min<mj_scalar_t>(my_min_coord));
      // get max
      Kokkos::parallel_reduce("get max",
        Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
          (coordinate_begin_index, coordinate_end_index),
        KOKKOS_LAMBDA (mj_lno_t j, mj_scalar_t & running_max) {
        int i = local_coordinate_permutations(j);
        if(mj_current_dim_coords(i) > running_max)
          running_max = mj_current_dim_coords(i);
      }, Kokkos::Max<mj_scalar_t>(my_max_coord));
      if(bUniformWeights) {
        my_total_weight = coordinate_end_index - coordinate_begin_index;
      }
      else {
        my_total_weight = 0;
        Kokkos::parallel_reduce("get weight",
          Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
            (coordinate_begin_index, coordinate_end_index),
          KOKKOS_LAMBDA (mj_lno_t j, mj_scalar_t & lsum) {
          int i = local_coordinate_permutations(j);
          lsum += local_mj_weights(i,0);
        }, my_total_weight);
      }
    }

    // single write
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
      (0, 1), KOKKOS_LAMBDA (int dummy) {
      local_process_local_min_max_coord_total_weight(kk) =
        my_min_coord;
      local_process_local_min_max_coord_total_weight(
        kk + current_concurrent_num_parts) = my_max_coord;
      local_process_local_min_max_coord_total_weight(
        kk + 2*current_concurrent_num_parts) = my_total_weight;
    });
  }
}

/*! \brief Function that reduces global minimum and maximum coordinates with
 * global total weight from given local arrays.
 * \param current_concurrent_num_parts is the number of parts whose cut lines
 * will be calculated concurrently.
 * \param local_min_max_total is the array holding local min and max coordinate
 * values with local total weight.
 * First concurrentPartCount entries are minimums of the parts, next
 * concurrentPartCount entries are max, and then the total weights.
 * \param global_min_max_total is the output array holding global min and global
 * coordinate values with global total weight.
 * The structure is same as localMinMaxTotal.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t,
  mj_node_t>::mj_get_global_min_max_coord_totW(
  mj_part_t current_concurrent_num_parts,
  Kokkos::View<mj_scalar_t *, device_t> & local_min_max_total,
  Kokkos::View<mj_scalar_t *, device_t> & global_min_max_total) {
  // reduce min for first current_concurrent_num_parts elements, reduce
  // max for next concurrentPartCount elements, reduce sum for the last
  // concurrentPartCount elements.
  if(this->comm->getSize()  > 1) {
    // We're using explicit host here as Spectrum MPI would fail
    // with the prior HostMirror UVMSpace to UVMSpace setup.
    auto host_local_min_max_total =
      Kokkos::create_mirror_view(Kokkos::HostSpace(), local_min_max_total);
    auto host_global_min_max_total =
      Kokkos::create_mirror_view(Kokkos::HostSpace(), global_min_max_total);
    Kokkos::deep_copy(host_local_min_max_total, local_min_max_total);
    Teuchos::MultiJaggedCombinedMinMaxTotalReductionOp<int, mj_scalar_t>
      reductionOp(current_concurrent_num_parts,
        current_concurrent_num_parts, current_concurrent_num_parts);
    try {
      reduceAll<int, mj_scalar_t>(
        *(this->comm),
      reductionOp,
      3 * current_concurrent_num_parts,
      host_local_min_max_total.data(),
      host_global_min_max_total.data());
    }
    Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))
    Kokkos::deep_copy(global_min_max_total, host_global_min_max_total);
  }
  else {
    mj_part_t s = 3 * current_concurrent_num_parts;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
      (0, s), KOKKOS_LAMBDA (mj_part_t i) {
      global_min_max_total(i) = local_min_max_total(i);
    });
  }
}

/*! \brief Function that calculates the new coordinates for the cut lines.
 * Function is called inside the parallel region.
 * \param min_coord minimum coordinate in the range.
 * \param max_coord maximum coordinate in the range.
 * \param num_cuts holds the number of cuts in current partitioning dimension.
 * \param global_weight holds the global total weight in the current part.
 * \param initial_cut_coords is the output array for the initial cut lines.
 * \param target_part_weights is the output array holding the cumulative ratios
 * of parts in current partitioning.
 * For partitioning to 4 uniformly, target_part_weights will be
 * (0.25 * globalTotalWeight, 0.5 *globalTotalWeight,
 *  0.75 * globalTotalWeight, globalTotalWeight).
 * \param future_num_part_in_parts is the vector that holds how many more parts
 * each part will be divided into more for the parts at the beginning of this
 * coordinate partitioning.
 * \param next_future_num_parts_in_parts is the vector that holds how many more
 * parts each part will be divided into more for the parts that will be obtained
 * at the end of this coordinate partitioning.
 * \param concurrent_current_part is the index of the part in the
 * future_num_part_in_parts vector.
 * \param obtained_part_index holds the amount of shift in the
 * next_future_num_parts_in_parts for the output parts.
 *
 * Nonuniform first level partitioning:
 * \param num_target_first_level_parts is the number of parts requested
 * after the first level of partitioning (resulting parts may be imbalanced)
 * \param target_first_level_dist is an array requesting the distribution of
 * elements in each part after the first cut (used for nonuniform first cuts)
 *
 * Ex. If we have num_first_level_parts = 3, first_level_dist = [4, 10, 6], then
 * target_part_weights will be [.20, .70, 1.00] * global_weight
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_get_initial_cut_coords_target_weights(
  mj_scalar_t min_coord,
  mj_scalar_t max_coord,
  mj_part_t num_cuts/*p-1*/ ,
  mj_scalar_t global_weight,
  /*p - 1 sized, coordinate of each cut line*/
  Kokkos::View<mj_scalar_t *, device_t> & initial_cut_coords,
  /*cumulative weights, at left side of each cut line. p-1 sized*/
  Kokkos::View<mj_scalar_t *, device_t> & current_target_part_weights ,
  std::vector <mj_part_t> *future_num_part_in_parts, //the vecto
  std::vector <mj_part_t> *next_future_num_parts_in_parts,
  mj_part_t concurrent_current_part,
  mj_part_t obtained_part_index,
  mj_part_t num_target_first_level_parts,
  const Kokkos::View<mj_part_t *, Kokkos::HostSpace> & target_first_level_dist)
{
  mj_scalar_t coord_range = max_coord - min_coord;

  // We decided we could keep some std::vectors around for now. Eventually
  // it would be nice to have everything just as views with some being device
  // and some host. This particular case needs a bit of work to get setup
  // in a cleaner way so not going to mess with it at the moment.

  bool bUniformPartsCheck =
    num_target_first_level_parts <= 1 && this->mj_uniform_parts(0);

  if(!bUniformPartsCheck) {
    bool bValidNonUniformTargetWeights =
      (num_target_first_level_parts > 1 && target_first_level_dist.size() != 0);
    if(!bValidNonUniformTargetWeights) {
      std::cerr << "MJ does not support non uniform part weights beyond the first partition" << std::endl;
      std::terminate();
    }
  }

  Kokkos::View<mj_scalar_t*, device_t> device_cumulative(
    "device_cumulative", num_cuts);
  auto host_cumulative = Kokkos::create_mirror_view(device_cumulative);

  mj_scalar_t cumulative = 0;

  if(bUniformPartsCheck) {
    // How many total future parts the part will be partitioned into.
    mj_scalar_t total_future_part_count_in_part =
      static_cast<mj_scalar_t>((*future_num_part_in_parts)[concurrent_current_part]);

    // How much each part should weigh in ideal case.
    mj_scalar_t unit_part_weight =
      global_weight / total_future_part_count_in_part;

    for(mj_part_t i = 0; i < num_cuts; ++i) {
      cumulative += unit_part_weight * static_cast<mj_scalar_t>((*next_future_num_parts_in_parts)[i + obtained_part_index]);
      host_cumulative(i) = cumulative;
    }
  }
  else {
    // Sum of entries in the first level partition distribution vector
    mj_scalar_t sum_target_first_level_dist = 0.0;
    for (int i = 0; i < num_target_first_level_parts; ++i) {
      sum_target_first_level_dist += target_first_level_dist(i);
    }

    for(mj_part_t i = 0; i < num_cuts; ++i) {
      cumulative += global_weight * target_first_level_dist(i) /
        sum_target_first_level_dist;
      host_cumulative(i) = cumulative;
    }
  }

  Kokkos::deep_copy(device_cumulative, host_cumulative);

  Kokkos::parallel_for("Write num in parts",
    Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
    (0, num_cuts), KOKKOS_LAMBDA(mj_part_t cut) {
      // set target part weight.
      current_target_part_weights(cut) = device_cumulative(cut);
      initial_cut_coords(cut) = min_coord +
        (coord_range * device_cumulative(cut)) / global_weight;
      // set this multiple times but here for device handling
      current_target_part_weights(num_cuts) = global_weight;
  });

  // round the target part weights.
  // Note need to discuss regarding DragonFly commits and determine if we
  // would not simply check mj_uniform_weights here.
  if (!bUniformPartsCheck || this->mj_uniform_weights[0]) {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
        (0, num_cuts + 1),
      KOKKOS_LAMBDA (mj_part_t i) {
      current_target_part_weights(i) =
        long(current_target_part_weights(i) + 0.5);
    });
  }
}

/*! \brief Function that calculates the new coordinates for the cut lines.
 * Function is called inside the parallel region.
 * \param max_coordinate maximum coordinate in the range.
 * \param min_coordinate minimum coordinate in the range.
 * \param concurrent_current_part_index is the index of the part in the
 * inTotalCounts vector.
 * \param coordinate_begin_index holds the beginning of the coordinates in
 * current part.
 * \param coordinate_end_index holds end of the coordinates in current part.
 * \param mj_current_coordinate_permutations is the permutation array, holds the
 * real indices of coordinates on mj_current_dim_coords array.
 * \param mj_current_dim_coords is the 1D array holding the coordinates.
 * \param mj_part_ids is the array holding the partIds of each coordinate.
 * \param partition_count is the number of parts that the current part will be
 * partitioned into.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  set_initial_coordinate_parts(
  mj_scalar_t &max_coordinate,
  mj_scalar_t &min_coordinate,
  mj_lno_t coordinate_begin_index,
  mj_lno_t coordinate_end_index,
  Kokkos::View<mj_lno_t *, device_t> & mj_current_coordinate_permutations,
  Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
  Kokkos::View<mj_part_t *, device_t> & mj_part_ids,
  mj_part_t &partition_count)
{
  mj_scalar_t coordinate_range = max_coordinate - min_coordinate;

  // if there is single point, or if all points are along a line.
  // set initial part to 0 for all.
  if(std::abs(coordinate_range) < this->sEpsilon ) {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
      (coordinate_begin_index, coordinate_end_index),
      KOKKOS_LAMBDA (mj_lno_t ii) {
      mj_part_ids(mj_current_coordinate_permutations[ii]) = 0;
    });
  }
  else {
    // otherwise estimate an initial part for each coordinate.
    // assuming uniform distribution of points.
    mj_scalar_t slice = coordinate_range / partition_count;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
      (coordinate_begin_index, coordinate_end_index),
      KOKKOS_LAMBDA (mj_lno_t ii) {
      mj_lno_t iii = mj_current_coordinate_permutations[ii];
      mj_part_t pp =
        mj_part_t((mj_current_dim_coords[iii] - min_coordinate) / slice);
      if(pp >= partition_count) {
        pp = partition_count - 1; // don't want last coord in an invalid part
      }
      mj_part_ids[iii] = 2 * pp;
    });
  }
}

/*! \brief Function that is responsible from 1D partitioning of the given range
 * of coordinates.
 * \param mj_current_dim_coords is 1 dimensional array holding
 * coordinate values.
 * \param imbalanceTolerance is the maximum allowed imbalance ratio.
 * \param current_work_part is the beginning index of concurrentPartCount parts.
 * \param current_concurrent_num_parts is the number of parts whose cut lines
 * will be calculated concurrently.
 * \param current_cut_coordinates is array holding the coordinates of the cut.
 * \param total_incomplete_cut_count is the number of cut lines whose positions
 * should be calculated.
 * \param view_rectilinear_cut_count   DOCWORK: Documentation
 * \param view_total_reduction_size  DOCWORK: Documentation
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t,mj_node_t>::mj_1D_part(
  Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
  double used_imbalance_tolerance,
  mj_part_t current_work_part,
  mj_part_t current_concurrent_num_parts,
  Kokkos::View<mj_scalar_t *, device_t> & current_cut_coordinates,
  mj_part_t total_incomplete_cut_count,
  Kokkos::View<mj_part_t *, device_t> & view_rectilinear_cut_count,
  Kokkos::View<size_t*, device_t> & view_total_reduction_size)
{
  this->temp_cut_coords = current_cut_coordinates;

  Teuchos::MultiJaggedCombinedReductionOp<mj_part_t, mj_scalar_t>
               *reductionOp = NULL;

  bool bSingleProcess = (this->comm->getSize() == 1);

  std::vector<mj_part_t> temp(host_num_partitioning_in_current_dim.size());
  if(!bSingleProcess) {
    for(size_t n = 0; n < host_num_partitioning_in_current_dim.size(); ++n) {
      temp[n] = host_num_partitioning_in_current_dim(n);
    }
    reductionOp = new Teuchos::MultiJaggedCombinedReductionOp
      <mj_part_t, mj_scalar_t>(
        &temp,
        current_work_part,
        current_concurrent_num_parts);
  }

  auto local_cut_lower_bound_coordinates =
    cut_lower_bound_coordinates;
  auto local_cut_upper_bound_coordinates =
    cut_upper_bound_coordinates;
  auto local_cut_upper_bound_weights = cut_upper_bound_weights;
  auto local_cut_lower_bound_weights = cut_lower_bound_weights;
  bool local_distribute_points_on_cut_lines = distribute_points_on_cut_lines;
  auto local_process_cut_line_weight_to_put_left =
    process_cut_line_weight_to_put_left;
  auto local_temp_cut_coords = temp_cut_coords;
  auto local_global_total_part_weight_left_right_closests =
    global_total_part_weight_left_right_closests;
  auto local_cut_coordinates_work_array =
    cut_coordinates_work_array;
  auto local_part_xadj = part_xadj;
  auto local_global_min_max_coord_total_weight =
    global_min_max_coord_total_weight;
  auto local_target_part_weights =
    target_part_weights;
  auto local_global_rectilinear_cut_weight =
    global_rectilinear_cut_weight;
  auto local_process_rectilinear_cut_weight =
    process_rectilinear_cut_weight;

  auto local_is_cut_line_determined = this->is_cut_line_determined;
  auto local_device_num_partitioning_in_current_dim =
    device_num_partitioning_in_current_dim;

  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
    KOKKOS_LAMBDA (int dummy) {

    // these need to be initialized
    view_rectilinear_cut_count(0) = 0;
    view_total_reduction_size(0) = 0;

    // initialize the lower and upper bounds of the cuts.
    mj_part_t next = 0;
    for(mj_part_t i = 0; i < current_concurrent_num_parts; ++i) {
      mj_part_t num_part_in_dim =
        local_device_num_partitioning_in_current_dim(current_work_part + i);
      mj_part_t num_cut_in_dim = num_part_in_dim - 1;
      view_total_reduction_size(0) += (4 * num_cut_in_dim + 1);

      for(mj_part_t ii = 0; ii < num_cut_in_dim; ++ii) {
        local_is_cut_line_determined(next) = false;
        // min coordinate
        local_cut_lower_bound_coordinates(next) =
          local_global_min_max_coord_total_weight(i);
        // max coordinate
        local_cut_upper_bound_coordinates(next) =
          local_global_min_max_coord_total_weight(
          i + current_concurrent_num_parts);
        // total weight
        local_cut_upper_bound_weights(next) =
          local_global_min_max_coord_total_weight(
          i + 2 * current_concurrent_num_parts);
        local_cut_lower_bound_weights(next) = 0;
        if(local_distribute_points_on_cut_lines) {
          local_process_cut_line_weight_to_put_left(next) = 0;
        }
        ++next;
      }
    }
  });

  // loop_count allows the kernel to behave differently on the first loop
  // and subsequent loops. First loop we do a binary search and subsequent
  // loops we simply step towards our target.
  int loop_count = 0;
  while (total_incomplete_cut_count != 0) {
    this->mj_1D_part_get_part_weights(
      current_concurrent_num_parts,
      current_work_part,
      mj_current_dim_coords,
      loop_count);
    ++loop_count;

    this->mj_combine_rightleft_and_weights(
      current_work_part,
      current_concurrent_num_parts);

    // now sum up the results of mpi processors.
    if(!bSingleProcess) {
      // We're using explicit host here as Spectrum MPI would fail
      // with the prior HostMirror UVMSpace to UVMSpace setup.
      auto host_total_part_weight_left_right_closests =
        Kokkos::create_mirror_view(Kokkos::HostSpace(),
        total_part_weight_left_right_closests);
      auto host_global_total_part_weight_left_right_closests =
        Kokkos::create_mirror_view(Kokkos::HostSpace(),
        global_total_part_weight_left_right_closests);

      Kokkos::deep_copy(host_total_part_weight_left_right_closests,
        total_part_weight_left_right_closests);

      size_t host_view_total_reduction_size;
      Kokkos::parallel_reduce("Read single",
        Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
        KOKKOS_LAMBDA(int dummy, size_t & set_single) {
        set_single = view_total_reduction_size(0);
      }, host_view_total_reduction_size);

      reduceAll<int, mj_scalar_t>( *(this->comm), *reductionOp,
        host_view_total_reduction_size,
        host_total_part_weight_left_right_closests.data(),
        host_global_total_part_weight_left_right_closests.data());
      Kokkos::deep_copy(global_total_part_weight_left_right_closests,
        host_global_total_part_weight_left_right_closests);
    }
    else {
      local_global_total_part_weight_left_right_closests =
        this->total_part_weight_left_right_closests;
    }

    // how much cut will be shifted for the next part in the concurrent
    // part calculation.
    mj_part_t cut_shift = 0;

    // how much the concantaneted array will be shifted for the next part
    // in concurrent part calculation.
    size_t tlr_shift = 0;

    Kokkos::View<mj_part_t*, Kokkos::HostSpace>
      save_initial_incomplete_cut_count("save_initial_incomplete_cut_count",
      current_concurrent_num_parts);

    for(mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk) {

      mj_part_t num_parts =
        host_num_partitioning_in_current_dim(current_work_part + kk);

      mj_part_t num_cuts = num_parts - 1;
      size_t num_total_part = num_parts + size_t (num_cuts);

      //if the cuts of this cut has already been completed.
      //nothing to do for this part.
      //just update the shift amount and proceed.
      mj_part_t kk_incomplete_cut_count = this->incomplete_cut_count(kk);

      if(kk_incomplete_cut_count == 0) {
        cut_shift += num_cuts;
        tlr_shift += (num_total_part + 2 * num_cuts);
        continue;
      }

      Kokkos::View<mj_scalar_t *, device_t> current_local_part_weights =
        Kokkos::subview(this->total_part_weight_left_right_closests,
          std::pair<mj_lno_t, mj_lno_t>(
            tlr_shift,
            this->total_part_weight_left_right_closests.size()));

      Kokkos::View<mj_scalar_t *, device_t> current_global_tlr =
        Kokkos::subview(
          local_global_total_part_weight_left_right_closests,
          std::pair<mj_lno_t, mj_lno_t>(
            tlr_shift,
            local_global_total_part_weight_left_right_closests.size()));
      Kokkos::View<mj_scalar_t *, device_t>
        current_global_left_closest_points =
        Kokkos::subview(current_global_tlr,
          std::pair<mj_lno_t, mj_lno_t>(
            num_total_part,
            current_global_tlr.size()));
      Kokkos::View<mj_scalar_t *, device_t>
        current_global_right_closest_points =
        Kokkos::subview(current_global_tlr,
          std::pair<mj_lno_t, mj_lno_t>(
            num_total_part + num_cuts,
            current_global_tlr.size()));
      Kokkos::View<mj_scalar_t *, device_t> current_global_part_weights =
        current_global_tlr;

      Kokkos::View<bool *, device_t> current_cut_line_determined =
        Kokkos::subview(this->is_cut_line_determined,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift,
            this->is_cut_line_determined.size()));
      Kokkos::View<mj_scalar_t *, device_t> current_part_target_weights =
        Kokkos::subview(local_target_part_weights,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift + kk,
            local_target_part_weights.size()));
      Kokkos::View<mj_scalar_t *, device_t>
        current_part_cut_line_weight_to_put_left =
        Kokkos::subview(local_process_cut_line_weight_to_put_left,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift,
            local_process_cut_line_weight_to_put_left.size()));

      save_initial_incomplete_cut_count(kk) =
        kk_incomplete_cut_count;

      Kokkos::View<mj_scalar_t *, device_t>
        current_cut_lower_bound_weights =
        Kokkos::subview(local_cut_lower_bound_weights,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift,
            local_cut_lower_bound_weights.size()));
      Kokkos::View<mj_scalar_t *, device_t> current_cut_upper_weights =
        Kokkos::subview(local_cut_upper_bound_weights,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift,
            local_cut_upper_bound_weights.size()));
      Kokkos::View<mj_scalar_t *, device_t> current_cut_upper_bounds =
        Kokkos::subview(local_cut_upper_bound_coordinates,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift,
            local_cut_upper_bound_coordinates.size()));
      Kokkos::View<mj_scalar_t *, device_t> current_cut_lower_bounds =
        Kokkos::subview(local_cut_lower_bound_coordinates,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift,
            local_cut_lower_bound_coordinates.size()));

      // Now compute the new cut coordinates.
      Kokkos::View<mj_scalar_t*, device_t> sub_temp_cut_coords =
        Kokkos::subview(this->temp_cut_coords,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift, this->temp_cut_coords.size()));
      Kokkos::View<mj_scalar_t*, device_t> sub_cut_coordinates_work_array =
        Kokkos::subview(this->cut_coordinates_work_array,
          std::pair<mj_lno_t, mj_lno_t>(
            cut_shift, this->cut_coordinates_work_array.size()));

      this->mj_get_new_cut_coordinates(
        current_concurrent_num_parts,
        kk,
        num_cuts,
        used_imbalance_tolerance,
        current_global_part_weights,
        current_local_part_weights,
        current_part_target_weights,
        current_cut_line_determined,
        sub_temp_cut_coords,
        current_cut_upper_bounds,
        current_cut_lower_bounds,
        current_global_left_closest_points,
        current_global_right_closest_points,
        current_cut_lower_bound_weights,
        current_cut_upper_weights,
        sub_cut_coordinates_work_array,
        current_part_cut_line_weight_to_put_left,
        view_rectilinear_cut_count);

      cut_shift += num_cuts;
      tlr_shift += (num_total_part + 2 * num_cuts);
    } // end of kk loop

    for(mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk) {
      mj_part_t iteration_complete_cut_count =
        save_initial_incomplete_cut_count(kk) - this->incomplete_cut_count(kk);
      total_incomplete_cut_count -= iteration_complete_cut_count;
    }

    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
      (0, local_temp_cut_coords.size()), KOKKOS_LAMBDA(int n) {
      auto t = local_temp_cut_coords(n);
      local_temp_cut_coords(n) = local_cut_coordinates_work_array(n);
      local_cut_coordinates_work_array(n) = t;
    });
  } // end of the while loop

  // Needed only if keep_cuts; otherwise can simply swap array pointers
  // cutCoordinates and cutCoordinatesWork.
  // (at first iteration, cutCoordinates == cutCoorindates_tmp).
  // computed cuts must be in cutCoordinates.
  if(current_cut_coordinates != local_temp_cut_coords) {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
      (0, 1), KOKKOS_LAMBDA(int dummy) {
      mj_part_t next = 0;
      for(mj_part_t i = 0; i < current_concurrent_num_parts; ++i) {
        mj_part_t num_parts = -1;
        num_parts = local_device_num_partitioning_in_current_dim(
          current_work_part + i);
        mj_part_t num_cuts = num_parts - 1;
        for(mj_part_t ii = 0; ii < num_cuts; ++ii) {
          current_cut_coordinates(next + ii) = local_temp_cut_coords(next + ii);
        }
        next += num_cuts;
      }
      for(int n = 0; n <
        static_cast<int>(local_cut_coordinates_work_array.size()); ++n) {
        local_cut_coordinates_work_array(n) = local_temp_cut_coords(n);
      }
    });
  }

  delete reductionOp;
}

template<class scalar_t>
struct Zoltan2_MJArrayType {
  scalar_t * ptr;

  // With new kokkos setup parallel_reduce will call empty constructor and
  // we update the ptr in the init method.
  KOKKOS_INLINE_FUNCTION
  Zoltan2_MJArrayType() : ptr(NULL) {};

  KOKKOS_INLINE_FUNCTION
  Zoltan2_MJArrayType(scalar_t * pSetPtr) : ptr(pSetPtr) {};
};

#ifndef KOKKOS_ENABLE_CUDA

template<class policy_t, class scalar_t, class part_t>
struct ArrayCombinationReducer {

  typedef ArrayCombinationReducer reducer;
  typedef Zoltan2_MJArrayType<scalar_t> value_type;
  scalar_t max_scalar;
  value_type * value;
  int value_count_rightleft;
  int value_count_weights;

  KOKKOS_INLINE_FUNCTION ArrayCombinationReducer(
    scalar_t mj_max_scalar,
    value_type &val,
    int mj_value_count_rightleft,
    int mj_value_count_weights) :
      max_scalar(mj_max_scalar),
      value(&val),
      value_count_rightleft(mj_value_count_rightleft),
      value_count_weights(mj_value_count_weights)
  {}

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src)  const {
    for(int n = 0; n < value_count_weights; ++n) {
      dst.ptr[n] += src.ptr[n];
    }

    for(int n = value_count_weights + 2;
      n < value_count_weights + value_count_rightleft - 2; n += 2) {
      if(src.ptr[n] > dst.ptr[n]) {
        dst.ptr[n] = src.ptr[n];
      }
      if(src.ptr[n+1] < dst.ptr[n+1]) {
        dst.ptr[n+1] = src.ptr[n+1];
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type& dst, const volatile value_type& src) const {
    for(int n = 0; n < value_count_weights; ++n) {
      dst.ptr[n] += src.ptr[n];
    }

    for(int n = value_count_weights + 2;
      n < value_count_weights + value_count_rightleft - 2; n += 2) {
      if(src.ptr[n] > dst.ptr[n]) {
        dst.ptr[n] = src.ptr[n];
      }
      if(src.ptr[n+1] < dst.ptr[n+1]) {
        dst.ptr[n+1] = src.ptr[n+1];
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const {
    dst.ptr = value->ptr; // must update ptr

    for(int n = 0; n < value_count_weights; ++n) {
      dst.ptr[n] = 0;
    }

    for(int n = value_count_weights;
      n < value_count_weights + value_count_rightleft; n += 2) {
      dst.ptr[n]   = -max_scalar;
      dst.ptr[n+1] =  max_scalar;
    }
  }
};
#endif // KOKKOS_ENABLE_CUDA

template<class policy_t, class scalar_t, class part_t, class index_t,
  class device_t, class array_t>
struct ReduceWeightsFunctor {
  typedef typename policy_t::member_type member_type;
  typedef Kokkos::View<scalar_t*> scalar_view_t;

#ifndef KOKKOS_ENABLE_CUDA
  typedef array_t value_type[];
#endif

  int loop_count;
  array_t max_scalar;

  part_t concurrent_current_part;
  part_t num_cuts;
  part_t current_work_part;
  part_t current_concurrent_num_parts;
  int value_count_rightleft;
  int value_count_weights;
  int value_count;
  Kokkos::View<index_t*, device_t> permutations;
  Kokkos::View<scalar_t *, device_t> coordinates;
  Kokkos::View<scalar_t**, device_t> weights;
  Kokkos::View<part_t*, device_t> parts;
  Kokkos::View<scalar_t *, device_t> cut_coordinates;
  Kokkos::View<index_t *, device_t> part_xadj;
  bool uniform_weights0;
  scalar_t sEpsilon;

#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::View<double *, device_t> current_part_weights;
  Kokkos::View<scalar_t *, device_t> current_left_closest;
  Kokkos::View<scalar_t *, device_t> current_right_closest;
#endif // KOKKOS_ENABLE_CUDA

  ReduceWeightsFunctor(
    int mj_loop_count,
    array_t mj_max_scalar,
    part_t mj_concurrent_current_part,
    part_t mj_num_cuts,
    part_t mj_current_work_part,
    part_t mj_current_concurrent_num_parts,
    part_t mj_left_right_array_size,
    part_t mj_weight_array_size,
    Kokkos::View<index_t*, device_t> & mj_permutations,
    Kokkos::View<scalar_t *, device_t> & mj_coordinates,
    Kokkos::View<scalar_t**, device_t> & mj_weights,
    Kokkos::View<part_t*, device_t> & mj_parts,
    Kokkos::View<scalar_t *, device_t> & mj_cut_coordinates,
    Kokkos::View<index_t *, device_t> & mj_part_xadj,
    bool mj_uniform_weights0,
    scalar_t mj_sEpsilon
#ifdef KOKKOS_ENABLE_CUDA
    ,Kokkos::View<double *, device_t> & mj_current_part_weights,
    Kokkos::View<scalar_t *, device_t> & mj_current_left_closest,
    Kokkos::View<scalar_t *, device_t> & mj_current_right_closest
#endif // KOKKOS_ENABLE_CUDA
    ) :
      loop_count(mj_loop_count),
      max_scalar(mj_max_scalar),
      concurrent_current_part(mj_concurrent_current_part),
      num_cuts(mj_num_cuts),
      current_work_part(mj_current_work_part),
      current_concurrent_num_parts(mj_current_concurrent_num_parts),
      value_count_rightleft(mj_left_right_array_size),
      value_count_weights(mj_weight_array_size),
      value_count(mj_weight_array_size+mj_left_right_array_size),
      permutations(mj_permutations),
      coordinates(mj_coordinates),
      weights(mj_weights),
      parts(mj_parts),
      cut_coordinates(mj_cut_coordinates),
      part_xadj(mj_part_xadj),
      uniform_weights0(mj_uniform_weights0),
      sEpsilon(mj_sEpsilon)
#ifdef KOKKOS_ENABLE_CUDA
      ,current_part_weights(mj_current_part_weights),
      current_left_closest(mj_current_left_closest),
      current_right_closest(mj_current_right_closest)
#endif // KOKKOS_ENABLE_CUDA
  {
  }

  size_t team_shmem_size (int team_size) const {
#ifdef KOKKOS_ENABLE_CUDA
    int result = sizeof(array_t) *
      (value_count_weights + value_count_rightleft);
#else
    int result = sizeof(array_t) *
      (value_count_weights + value_count_rightleft) * team_size;
#endif

    // pad this to a multiple of 8 or it will run corrupt
    int remainder = result % 8;
    if(remainder != 0) {
      result += 8 - remainder;
    }
    return result;
  }

  KOKKOS_INLINE_FUNCTION
#ifdef KOKKOS_ENABLE_CUDA
  void operator() (const member_type & teamMember) const {
#else
  void operator() (const member_type & teamMember, value_type teamSum) const {
#endif

    index_t all_begin = (concurrent_current_part == 0) ? 0 :
      part_xadj(concurrent_current_part - 1);
    index_t all_end = part_xadj(concurrent_current_part);

    index_t num_working_points = all_end - all_begin;
    int num_teams = teamMember.league_size();

    index_t stride = num_working_points / num_teams;
    if((num_working_points % num_teams) > 0) {
      stride += 1; // make sure we have coverage for the final points
    }

    // the last team may have less work than the other teams
    // the last team can be empty (begin > end) if num_teams > stride
    // which is true for many teams and small numbers of coords (tests)
    index_t begin = all_begin + stride * teamMember.league_rank();
    index_t end = begin + stride;
    if(end > all_end) {
      end = all_end;
    }

#ifdef KOKKOS_ENABLE_CUDA
    size_t sh_mem_size = sizeof(array_t) * (value_count_weights +
      value_count_rightleft);

    array_t * shared_ptr = (array_t *) teamMember.team_shmem().get_shmem(
      sh_mem_size);

    // init the shared array to 0
    Kokkos::single(Kokkos::PerTeam(teamMember), [=] () {
      for(int n = 0; n < value_count_weights; ++n) {
        shared_ptr[n] = 0;
      }
      for(int n = value_count_weights;
        n < value_count_weights + value_count_rightleft; n += 2) {
        shared_ptr[n]   = -max_scalar;
        shared_ptr[n+1] =  max_scalar;
      }
    });
    teamMember.team_barrier();

    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(teamMember, begin, end),
      [=] (index_t ii) {
#else // KOKKOS_ENABLE_CUDA
    // create the team shared data - each thread gets one of the arrays
    size_t sh_mem_size = sizeof(array_t) * (value_count_weights +
      value_count_rightleft) * teamMember.team_size();

    array_t * shared_ptr = (array_t *) teamMember.team_shmem().get_shmem(
      sh_mem_size);

    // select the array for this thread
    Zoltan2_MJArrayType<array_t> array(&shared_ptr[teamMember.team_rank() *
        (value_count_weights + value_count_rightleft)]);

    // create reducer which handles the Zoltan2_MJArrayType class
    ArrayCombinationReducer<policy_t, array_t, part_t> arraySumReducer(
      max_scalar, array,
      value_count_rightleft,
      value_count_weights);

    Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(teamMember, begin, end),
      [=] (size_t ii, Zoltan2_MJArrayType<array_t>& threadSum) {
#endif // KOKKOS_ENABLE_CUDA

      int i = permutations(ii);
      scalar_t coord = coordinates(i);
      array_t w = uniform_weights0 ? 1 : (array_t) weights(i,0);

      // now check each part and it's right cut
      index_t part = parts(i)/2;

      int upper = num_cuts;
      int lower = 0;

      // binary search - find matching part
      while(true) {
        scalar_t a = (part == 0) ? -max_scalar : cut_coordinates(part-1);
        scalar_t b = (part == num_cuts) ? max_scalar : cut_coordinates(part);

        if(coord >= a + sEpsilon && coord <= b - sEpsilon) {
#ifdef KOKKOS_ENABLE_CUDA
          Kokkos::atomic_add(&shared_ptr[part*2], w);
#else
          threadSum.ptr[part*2] += w;
#endif

          parts(i) = part*2;

          // now handle the left/right closest part
#ifdef KOKKOS_ENABLE_CUDA
          array_t new_value = (array_t) coord;
          array_t prev_value = shared_ptr[value_count_weights + part * 2 + 1];
          while(new_value < prev_value) {
            prev_value = Kokkos::atomic_compare_exchange(
                &shared_ptr[value_count_weights + part * 2 + 1],
                prev_value, new_value);
          }
          prev_value = shared_ptr[value_count_weights + part * 2 + 2];
          while(new_value > prev_value) {
            prev_value = Kokkos::atomic_compare_exchange(
              &shared_ptr[value_count_weights + part * 2 + 2],
              prev_value, new_value);
          }
#else
          // note cut to left needs to set right closest and cut to right needs
          // to set left closest. It's index +1 and +2 instead of -1 and +0
          // because right/left segment is padded with an extra pair at
          // begining and end to avoid branching with if checks.
          if(coord < threadSum.ptr[value_count_weights + part * 2 + 1]) {
            threadSum.ptr[value_count_weights + part * 2 + 1] = coord;
          }
          if(coord > threadSum.ptr[value_count_weights + part * 2 + 2]) {
            threadSum.ptr[value_count_weights + part * 2 + 2] = coord;
          }
#endif

          break;
        }
        else if(part != num_cuts) {
          if(coord < b + sEpsilon && coord > b - sEpsilon) {
            // Note if on cut we set right/left closest to the cut itself
            // but we add +2 because we buffered the area with an extra slot
            // to reduce cuda branching. So it's +2, +3 instead of +0, +1.
#ifdef KOKKOS_ENABLE_CUDA
            Kokkos::atomic_add(&shared_ptr[part*2+1], w);
            shared_ptr[value_count_weights + part * 2 + 2] = b;
            shared_ptr[value_count_weights + part * 2 + 3] = b;
#else
            threadSum.ptr[part*2+1] += w;
            threadSum.ptr[value_count_weights + part * 2 + 2] = b;
            threadSum.ptr[value_count_weights + part * 2 + 3] = b;
#endif

            parts(i) = part*2+1;

            // Need to scan up for any other cuts of same coordinate
            // This is costly but it's only relevant for the fix4785 test
            // which loads a lot of coordinates on the same point, so without
            // this our cuts would all just sit at 0.
            part_t base_b = part;
            scalar_t base_coord = cut_coordinates(base_b);
            part += 1;
            while(part < num_cuts) {
              b = cut_coordinates(part);
              scalar_t delta = b - base_coord;
              if(delta < 0) delta = -delta;
              if(delta < sEpsilon) {
                // Note if on cut we set right/left closest to the cut itself
                // but we add +2 because we buffered the area with an extra slot
                // to reduce cuda branching. So it's +2, +3 instead of +0, +1.
#ifdef KOKKOS_ENABLE_CUDA
                Kokkos::atomic_add(&shared_ptr[part*2+1], w);
                shared_ptr[value_count_weights + part * 2 + 2] = b;
                shared_ptr[value_count_weights + part * 2 + 3] = b;
#else
                threadSum.ptr[part*2+1] += w;
                threadSum.ptr[value_count_weights + part * 2 + 2] = b;
                threadSum.ptr[value_count_weights + part * 2 + 3] = b;
#endif
              }
              else { break; }
              ++part;
            }
            part = base_b - 1;
            while(part >= 0) {
              b = cut_coordinates(part);
              scalar_t delta = b - base_coord;
              if(delta < 0) delta = -delta;
              if(delta < sEpsilon) {
                // Note if on cut we set right/left closest to the cut itself
                // but we add +2 because we buffered the area with an extra slot
                // to reduce cuda branching. So it's +2, +3 instead of +0, +1.
#ifdef KOKKOS_ENABLE_CUDA
                Kokkos::atomic_add(&shared_ptr[part*2+1], w);
                shared_ptr[value_count_weights + part * 2 + 2] = b;
                shared_ptr[value_count_weights + part * 2 + 3] = b;
#else
                threadSum.ptr[part*2+1] += w;
                threadSum.ptr[value_count_weights + part * 2 + 2] = b;
                threadSum.ptr[value_count_weights + part * 2 + 3] = b;
#endif
              }
              else { break; }
              --part;
            }

            break;
          }
        }

        if(loop_count != 0) {
          // subsequent loops can just step towards target
          if(coord < b) {
            part -= 1;
          }
          else {
            part += 1;
          }
        }
        else {
          // initial loop binary search
          if(coord < b) {
            if(part == lower + 1) {
              part = lower;
            }
            else {
              upper = part - 1;
              part -= (part - lower)/2;
            }
          }
          else if(part == upper - 1) {
            part = upper;
          }
          else {
            lower = part + 1;
            part += (upper - part)/2;
          }
        }
      }
#ifdef KOKKOS_ENABLE_CUDA
    });
#else // KOKKOS_ENABLE_CUDA
    }, arraySumReducer);
#endif // KOKKOS_ENABLE_CUDA

    teamMember.team_barrier();

    // collect all the team's results
    Kokkos::single(Kokkos::PerTeam(teamMember), [=] () {
      for(int n = 0; n < value_count_weights; ++n) {
#ifdef KOKKOS_ENABLE_CUDA
        Kokkos::atomic_add(&current_part_weights(n),
          static_cast<double>(shared_ptr[n]));
#else // KOKKOS_ENABLE_CUDA
        teamSum[n] += array.ptr[n];
#endif // KOKKOS_ENABLE_CUDA
      }

#ifdef KOKKOS_ENABLE_CUDA
      int insert_left = 0;
      int insert_right = 0;
#endif

      for(int n = 2 + value_count_weights;
        n < value_count_weights + value_count_rightleft - 2; n += 2) {
#ifdef KOKKOS_ENABLE_CUDA
        scalar_t new_value = shared_ptr[n+1];
        scalar_t prev_value = current_right_closest(insert_right);
        while(new_value < prev_value) {
          prev_value = Kokkos::atomic_compare_exchange(
            &current_right_closest(insert_right), prev_value, new_value);
        }

        new_value = shared_ptr[n];
        prev_value = current_left_closest(insert_left);
        while(new_value > prev_value) {
          prev_value = Kokkos::atomic_compare_exchange(
            &current_left_closest(insert_left), prev_value, new_value);
        }

        ++insert_left;
        ++insert_right;
#else // KOKKOS_ENABLE_CUDA
        if(array.ptr[n] > teamSum[n]) {
          teamSum[n] = array.ptr[n];
        }
        if(array.ptr[n+1] < teamSum[n+1]) {
          teamSum[n+1] = array.ptr[n+1];
        }
#endif // KOKKOS_ENABLE_CUDA
      }
    });

    teamMember.team_barrier();
  }

#ifndef KOKKOS_ENABLE_CUDA
  KOKKOS_INLINE_FUNCTION
  void join(value_type dst, const value_type src)  const {
    for(int n = 0; n < value_count_weights; ++n) {
      dst[n] += src[n];
    }

    for(int n = value_count_weights + 2;
      n < value_count_weights + value_count_rightleft - 2; n += 2) {
      if(src[n] > dst[n]) {
        dst[n] = src[n];
      }
      if(src[n+1] < dst[n+1]) {
        dst[n+1] = src[n+1];
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type dst, const volatile value_type src) const {
    for(int n = 0; n < value_count_weights; ++n) {
      dst[n] += src[n];
    }

    for(int n = value_count_weights + 2;
      n < value_count_weights + value_count_rightleft - 2; n += 2) {
      if(src[n] > dst[n]) {
        dst[n] = src[n];
      }
      if(src[n+1] < dst[n+1]) {
        dst[n+1] = src[n+1];
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type dst) const {
    for(int n = 0; n < value_count_weights; ++n) {
      dst[n] = 0;
    }

    for(int n = value_count_weights;
      n < value_count_weights + value_count_rightleft; n += 2) {
      dst[n]   = -max_scalar;
      dst[n+1] =  max_scalar;
    }
  }
#endif // KOKKOS_ENABLE_CUDA
};

/*! \brief Function that calculates the weights of each part according to given
 * part cut coordinates. Function is called inside the parallel region. Thread
 * specific work arrays are provided as function parameter.
 * \param total_part_count is the sum of number of cutlines and number of parts.
 * Simply it is 2*P - 1.
 * DOCWORK: Documentation
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t,mj_part_t, mj_node_t>::
  mj_1D_part_get_part_weights(
  mj_part_t current_concurrent_num_parts,
  mj_part_t current_work_part,
  Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
  int loop_count)
{
  auto local_is_cut_line_determined = is_cut_line_determined;
  auto local_thread_part_weights = thread_part_weights;
  auto local_thread_cut_left_closest_point = thread_cut_left_closest_point;
  auto local_thread_cut_right_closest_point = thread_cut_right_closest_point;

  // Create some locals so we don't use this inside the kernels
  // which causes problems
  auto local_sEpsilon = this->sEpsilon;
  auto local_assigned_part_ids = this->assigned_part_ids;
  auto local_coordinate_permutations = this->coordinate_permutations;
  auto local_mj_weights = this->mj_weights;
  auto local_part_xadj = this->part_xadj;
  auto local_global_min_max_coord_total_weight =
    this->global_min_max_coord_total_weight;

  typedef Kokkos::TeamPolicy<typename mj_node_t::execution_space> policy_t;

  auto local_device_num_partitioning_in_current_dim =
    device_num_partitioning_in_current_dim;

  Kokkos::deep_copy(device_incomplete_cut_count, this->incomplete_cut_count);
  auto local_device_incomplete_cut_count = device_incomplete_cut_count;

  mj_part_t total_part_shift = 0;

  mj_part_t concurrent_cut_shifts = 0;
  for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {
    Kokkos::View<mj_scalar_t *, device_t> local_temp_cut_coords =
      Kokkos::subview(temp_cut_coords, std::pair<mj_lno_t, mj_lno_t>(
        concurrent_cut_shifts, temp_cut_coords.size()));

    mj_part_t num_parts =
      host_num_partitioning_in_current_dim(current_work_part + kk);
    mj_part_t num_cuts = num_parts - 1;
    mj_part_t total_part_count = num_parts + num_cuts;
    mj_part_t weight_array_length = num_cuts + num_parts;

    // for right/left closest + buffer cut on either side
    mj_part_t right_left_array_length = (num_cuts + 2) * 2;

    if(this->incomplete_cut_count(kk) == 0) {
      total_part_shift += total_part_count;
      concurrent_cut_shifts += num_cuts;
      continue;
    }

    // if not set use 60 - was initial testing amount but somewhat arbitrary
    auto policy_ReduceWeightsFunctor = policy_t(
      mj_num_teams ? mj_num_teams : 60, Kokkos::AUTO);

#ifndef KOKKOS_ENABLE_CUDA
    int total_array_length =
      weight_array_length + right_left_array_length;
#endif

    // Using float here caused some numerical errors for coord on cut calculations.
    // Probably that can be fixed with proper epsilon adjustment but since cuda
    // doesn't reduce right now the shared memory pressure is no longer relevant.
    // Just use scalar_t to match the original algorithm.
    typedef mj_scalar_t array_t;

#ifndef KOKKOS_ENABLE_CUDA
    array_t * reduce_array =
      new array_t[static_cast<size_t>(total_array_length)];
#endif // KOKKOS_ENABLE_CUDA

    int offset_cuts = 0;
    for(int kk2 = 0; kk2 < kk; ++kk2) {
      offset_cuts +=
        host_num_partitioning_in_current_dim(current_work_part + kk2) - 1;
    }
    Kokkos::View<double *, device_t> my_current_part_weights =
      Kokkos::subview(local_thread_part_weights,
        std::pair<mj_lno_t, mj_lno_t>(total_part_shift,
         total_part_shift + total_part_count));
    Kokkos::View<mj_scalar_t *, device_t> my_current_left_closest =
      Kokkos::subview(local_thread_cut_left_closest_point,
      std::pair<mj_lno_t, mj_lno_t>(
        offset_cuts,
        local_thread_cut_left_closest_point.size()));
    Kokkos::View<mj_scalar_t *, device_t> my_current_right_closest =
      Kokkos::subview(local_thread_cut_right_closest_point,
        std::pair<mj_lno_t, mj_lno_t>(
          offset_cuts,
          local_thread_cut_right_closest_point.size()));

    array_t max_scalar = std::numeric_limits<array_t>::max();

#ifdef KOKKOS_ENABLE_CUDA
    // initialize values
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
      KOKKOS_LAMBDA (int dummy) {
      for(int n = 0; n < weight_array_length; ++n) {
        my_current_part_weights(n) = 0;
      }
      for(int n = 0; n < num_cuts; ++n) {
        my_current_left_closest(n) = -max_scalar;
        my_current_right_closest(n) = max_scalar;
      }
    });
#endif

    mj_part_t concurrent_current_part =
      current_work_part + kk;

    ReduceWeightsFunctor<policy_t, mj_scalar_t, mj_part_t, mj_lno_t,
      typename mj_node_t::device_type, array_t>
      teamFunctor(
        loop_count,
        max_scalar,
        concurrent_current_part,
        num_cuts,
        current_work_part,
        current_concurrent_num_parts,
        right_left_array_length,
        weight_array_length,
        coordinate_permutations,
        mj_current_dim_coords,
        mj_weights,
        assigned_part_ids,
        local_temp_cut_coords,
        part_xadj,
        mj_uniform_weights(0), // host and currently only relevant to slot 0
        sEpsilon
#ifdef KOKKOS_ENABLE_CUDA
        ,my_current_part_weights,
        my_current_left_closest,
        my_current_right_closest
#endif
        );

#ifdef KOKKOS_ENABLE_CUDA
    Kokkos::parallel_for(policy_ReduceWeightsFunctor, teamFunctor);
#else
    Kokkos::parallel_reduce(policy_ReduceWeightsFunctor,
      teamFunctor, reduce_array);
#endif

#ifndef KOKKOS_ENABLE_CUDA
    auto hostArray = Kokkos::create_mirror_view(my_current_part_weights);

    for(int i = 0; i < static_cast<int>(total_part_count); ++i) {
      hostArray(i) = reduce_array[i];
    }

    Kokkos::deep_copy(my_current_part_weights, hostArray);

    auto hostLeftArray = Kokkos::create_mirror_view(my_current_left_closest);
    auto hostRightArray = Kokkos::create_mirror_view(my_current_right_closest);
    for(mj_part_t cut = 0; cut < num_cuts; ++cut) {
      hostLeftArray(cut)  = reduce_array[weight_array_length + (cut+1)*2+0];
      hostRightArray(cut) = reduce_array[weight_array_length + (cut+1)*2+1];
    }
    Kokkos::deep_copy(my_current_left_closest, hostLeftArray);
    Kokkos::deep_copy(my_current_right_closest, hostRightArray);

    delete [] reduce_array;
#endif

    total_part_shift += total_part_count;
    concurrent_cut_shifts += num_cuts;
  }

  auto local_temp_cut_coords = temp_cut_coords;

  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
    (0, current_concurrent_num_parts), KOKKOS_LAMBDA(mj_part_t kk) {
    mj_part_t num_parts = local_device_num_partitioning_in_current_dim(
      current_work_part + kk);
    mj_part_t num_cuts = num_parts - 1;
    mj_part_t total_part_count = num_parts + num_cuts;

    if(local_device_incomplete_cut_count(kk) > 0) {
      // get the prefix sum
      // This is an inefficiency but not sure if it matters much
      size_t offset = 0;
      size_t offset_cuts = 0;
      for(mj_part_t kk2 = 0; kk2 < kk; ++kk2) {
        auto num_parts_kk2 = local_device_num_partitioning_in_current_dim(
          current_work_part + kk2);
        offset += num_parts_kk2 * 2 - 1;
        offset_cuts += num_parts_kk2 - 1;
      }

      for(mj_part_t i = 1; i < total_part_count; ++i) {
        // check for cuts sharing the same position; all cuts sharing a position
        // have the same weight == total weight for all cuts sharing the
        // position. Don't want to accumulate that total weight more than once.
        if(i % 2 == 0 && i > 1 && i < total_part_count - 1 &&
          std::abs(local_temp_cut_coords(offset_cuts + i / 2) -
            local_temp_cut_coords(offset_cuts + i /2 - 1))
            < local_sEpsilon) {
          // i % 2 = 0 when part i represents the cut coordinate.
          // if it is a cut, and if next cut also has the same coordinate, then
          // dont addup.
          local_thread_part_weights(offset + i)
            = local_thread_part_weights(offset + i-2);
          continue;
        }

        // otherwise do the prefix sum.
        local_thread_part_weights(offset + i) +=
          local_thread_part_weights(offset + i-1);
      }
    }
  });
}

/*! \brief Function that reduces the result of multiple threads
 * for left and right closest points and part weights in a single mpi process.
 * \param current_work_part holds the index of the first part (important when
 * concurrent parts are used.)
 * \param current_concurrent_num_parts is the number of parts whose cut lines
 * will be calculated concurrently.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_combine_rightleft_and_weights(
  mj_part_t current_work_part,
  mj_part_t current_concurrent_num_parts)
{
  auto local_thread_part_weights = this->thread_part_weights;
  auto local_is_cut_line_determined = this->is_cut_line_determined;
  auto local_thread_cut_left_closest_point =
    this->thread_cut_left_closest_point;
  auto local_thread_cut_right_closest_point =
    this->thread_cut_right_closest_point;
  auto local_total_part_weight_left_right_closests =
    this->total_part_weight_left_right_closests;
  auto local_device_num_partitioning_in_current_dim =
    device_num_partitioning_in_current_dim;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int>(0,1),
    KOKKOS_LAMBDA (int dummy) {

    size_t tlr_array_shift = 0;
    mj_part_t cut_shift = 0;
    size_t total_part_array_shift = 0;

    // iterate for all concurrent parts to find the left and right closest
    // points in the process.
    for(mj_part_t i = 0; i < current_concurrent_num_parts; ++i) {

      mj_part_t num_parts_in_part =
        local_device_num_partitioning_in_current_dim(current_work_part + i);
      mj_part_t num_cuts_in_part = num_parts_in_part - 1;
      size_t num_total_part_in_part =
        num_parts_in_part + size_t (num_cuts_in_part);

      // iterate for cuts in a single part.
      for(int ii = 0; ii < num_cuts_in_part; ++ii) {
        mj_part_t next = tlr_array_shift + ii;
        mj_part_t cut_index = cut_shift + ii;

        if(!local_is_cut_line_determined(cut_index)) {
          mj_scalar_t left_closest_in_process =
            local_thread_cut_left_closest_point(cut_index);
          mj_scalar_t right_closest_in_process =
            local_thread_cut_right_closest_point(cut_index);

          // store the left and right closes points.
          local_total_part_weight_left_right_closests(
            num_total_part_in_part + next) = left_closest_in_process;

          local_total_part_weight_left_right_closests(
            num_total_part_in_part + num_cuts_in_part + next) =
            right_closest_in_process;
        }
      }

      for(size_t j = 0; j < num_total_part_in_part; ++j) {
        mj_part_t cut_ind = j / 2 + cut_shift;

        // need to check j !=  num_total_part_in_part - 1
        // which is same as j/2 != num_cuts_in_part.
        // we cannot check it using cut_ind, because of the concurrent part
        // concantanetion.
        if(j == num_total_part_in_part - 1 ||
          !local_is_cut_line_determined(cut_ind)) {
          double pwj = local_thread_part_weights(total_part_array_shift + j);
          local_total_part_weight_left_right_closests(tlr_array_shift + j) = pwj;
        }
      }

      // set the shift position in the arrays
      cut_shift += num_cuts_in_part;
      tlr_array_shift += num_total_part_in_part + 2 * num_cuts_in_part;
      total_part_array_shift += num_total_part_in_part;
    }
  });
}

/*! \brief
 * Function that calculates the next pivot position,
 * according to given coordinates of upper bound and lower bound, the weights at
 * upper and lower bounds, and the expected weight.
 * \param cut_upper_bound is the upper bound coordinate of the cut.
 * \param cut_lower_bound is the lower bound coordinate of the cut.
 * \param cut_upper_weight is the weights at the upper bound of the cut.
 * \param cut_lower_weight is the weights at the lower bound of the cut.
 * \param expected_weight is the expected weight that should be placed on the
 * left of the cut line.
 * \param new_cut_position  DOCWORK: Documentation
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t,
  mj_node_t>::mj_calculate_new_cut_position(mj_scalar_t cut_upper_bound,
  mj_scalar_t cut_lower_bound,
  mj_scalar_t cut_upper_weight,
  mj_scalar_t cut_lower_weight,
  mj_scalar_t expected_weight,
  mj_scalar_t &new_cut_position) {

  if(std::abs(cut_upper_bound - cut_lower_bound) < this->sEpsilon) {
    new_cut_position = cut_upper_bound; //or lower bound does not matter.
  }

  if(std::abs(cut_upper_weight - cut_lower_weight) < this->sEpsilon) {
    new_cut_position = cut_lower_bound;
  }

  mj_scalar_t coordinate_range = (cut_upper_bound - cut_lower_bound);
  mj_scalar_t weight_range = (cut_upper_weight - cut_lower_weight);
  mj_scalar_t my_weight_diff = (expected_weight - cut_lower_weight);

  mj_scalar_t required_shift = (my_weight_diff / weight_range);
  int scale_constant = 20;
  int shiftint= int (required_shift * scale_constant);
  if(shiftint == 0) shiftint = 1;
  required_shift = mj_scalar_t (shiftint) / scale_constant;
  new_cut_position = coordinate_range * required_shift + cut_lower_bound;
}

#ifndef KOKKOS_ENABLE_CUDA

template<class policy_t, class scalar_t>
struct ArrayReducer {

  typedef ArrayReducer reducer;
  typedef Zoltan2_MJArrayType<scalar_t> value_type;
  value_type * value;
  int value_count;

  KOKKOS_INLINE_FUNCTION ArrayReducer(
    value_type &val,
    int mj_value_count) :
      value(&val),
      value_count(mj_value_count)
  {}

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src)  const {
    for(int n = 0; n < value_count; ++n) {
      dst.ptr[n] += src.ptr[n];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type& dst, const volatile value_type& src) const {
    for(int n = 0; n < value_count; ++n) {
      dst.ptr[n] += src.ptr[n];
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const {
    dst.ptr = value->ptr; // must update ptr
    for(int n = 0; n < value_count; ++n) {
      dst.ptr[n] = 0;
    }
  }
};

#endif

template<class policy_t, class scalar_t, class part_t, class index_t,
  class device_t, class array_t>
struct ReduceArrayFunctor {
  typedef typename policy_t::member_type member_type;
  typedef Kokkos::View<scalar_t*> scalar_view_t;

#ifndef KOKKOS_ENABLE_CUDA
  typedef array_t value_type[];
#endif

  part_t concurrent_current_part;
  int value_count;
  Kokkos::View<index_t*, device_t> permutations;
  Kokkos::View<scalar_t *, device_t> coordinates;
  Kokkos::View<part_t*, device_t> parts;
  Kokkos::View<index_t *, device_t> part_xadj;
  Kokkos::View<index_t *, device_t> track_on_cuts;

#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::View<int *, device_t> local_point_counts;
#endif // KOKKOS_ENABLE_CUDA

  ReduceArrayFunctor(
    part_t mj_concurrent_current_part,
    part_t mj_weight_array_size,
    Kokkos::View<index_t*, device_t> & mj_permutations,
    Kokkos::View<scalar_t *, device_t> & mj_coordinates,
    Kokkos::View<part_t*, device_t> & mj_parts,
    Kokkos::View<index_t *, device_t> & mj_part_xadj,
    Kokkos::View<index_t *, device_t> & mj_track_on_cuts
#ifdef KOKKOS_ENABLE_CUDA
    ,Kokkos::View<int *, device_t> & mj_local_point_counts
#endif // KOKKOS_ENABLE_CUDA
    ) :
      concurrent_current_part(mj_concurrent_current_part),
      value_count(mj_weight_array_size),
      permutations(mj_permutations),
      coordinates(mj_coordinates),
      parts(mj_parts),
      part_xadj(mj_part_xadj),
      track_on_cuts(mj_track_on_cuts)
#ifdef KOKKOS_ENABLE_CUDA
      ,local_point_counts(mj_local_point_counts)
#endif
  {
  }

  size_t team_shmem_size (int team_size) const {
#ifdef KOKKOS_ENABLE_CUDA
    int result = sizeof(array_t) * (value_count);
#else
    int result = sizeof(array_t) * (value_count) * team_size;
#endif

    // pad this to a multiple of 8 or it will run corrupt
    int remainder = result % 8;
    if(remainder != 0) {
      result += 8 - remainder;
    }
    return result;
  }

  KOKKOS_INLINE_FUNCTION
#ifdef KOKKOS_ENABLE_CUDA
  void operator() (const member_type & teamMember) const {
#else
  void operator() (const member_type & teamMember, value_type teamSum) const {
#endif
    index_t all_begin = (concurrent_current_part == 0) ? 0 :
      part_xadj(concurrent_current_part - 1);
    index_t all_end = part_xadj(concurrent_current_part);

    index_t num_working_points = all_end - all_begin;
    int num_teams = teamMember.league_size();

    index_t stride = num_working_points / num_teams;
    if((num_working_points % num_teams) > 0) {
      stride += 1; // make sure we have coverage for the final points
    }

    index_t begin = all_begin + stride * teamMember.league_rank();
    index_t end = begin + stride;
    if(end > all_end) {
      end = all_end; // the last team may have less work than the other teams
    }

    int track_on_cuts_insert_index = track_on_cuts.size() - 1;

    // create the team shared data - each thread gets one of the arrays
#ifdef KOKKOS_ENABLE_CUDA
    size_t sh_mem_size = sizeof(array_t) * (value_count);
#else
    size_t sh_mem_size =
      sizeof(array_t) * (value_count) * teamMember.team_size();
#endif

    array_t * shared_ptr = (array_t *) teamMember.team_shmem().get_shmem(
      sh_mem_size);

#ifdef KOKKOS_ENABLE_CUDA
    // init the shared array to 0
    Kokkos::single(Kokkos::PerTeam(teamMember), [=] () {
      for(int n = 0; n < value_count; ++n) {
        shared_ptr[n] = 0;
      }
    });
    teamMember.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, begin, end),
      [=] (index_t ii) {
#else // KOKKOS_ENABLE_CUDA
    // select the array for this thread
    Zoltan2_MJArrayType<array_t> array(&shared_ptr[teamMember.team_rank() *
        (value_count)]);

    // create reducer which handles the Zoltan2_MJArrayType class
    ArrayReducer<policy_t, array_t> arrayReducer(array, value_count);

    Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(teamMember, begin, end),
      [=] (size_t ii, Zoltan2_MJArrayType<array_t>& threadSum) {
#endif // KOKKOS_ENABLE_CUDA

      index_t coordinate_index = permutations(ii);
      part_t place = parts(coordinate_index);
      part_t part = place / 2;
      if(place % 2 == 0) {
#ifdef KOKKOS_ENABLE_CUDA
        Kokkos::atomic_add(&shared_ptr[part], 1);
#else
        threadSum.ptr[part] += 1;
#endif

        parts(coordinate_index) = part;
      }
      else {
        // fill a tracking array so we can process these slower points
        // in next cycle
        index_t set_index = Kokkos::atomic_fetch_add(
            &track_on_cuts(track_on_cuts_insert_index), 1);
        track_on_cuts(set_index) = ii;
      }
#ifdef KOKKOS_ENABLE_CUDA
    });
#else // KOKKOS_ENABLE_CUDA
    }, arrayReducer);
#endif // KOKKOS_ENABLE_CUDA

    teamMember.team_barrier();

    // collect all the team's results
    Kokkos::single(Kokkos::PerTeam(teamMember), [=] () {
      for(int n = 0; n < value_count; ++n) {
#ifdef KOKKOS_ENABLE_CUDA
        Kokkos::atomic_add(&local_point_counts(n), shared_ptr[n]);
#else // KOKKOS_ENABLE_CUDA
        teamSum[n] += array.ptr[n];
#endif // KOKKOS_ENABLE_CUDA
      }
    });

    teamMember.team_barrier();
  }

#ifndef KOKKOS_ENABLE_CUDA

  KOKKOS_INLINE_FUNCTION
  void join(value_type dst, const value_type src)  const {
    for(int n = 0; n < value_count; ++n) {
      dst[n] += src[n];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type dst, const volatile value_type src) const {
    for(int n = 0; n < value_count; ++n) {
      dst[n] += src[n];
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type dst) const {
    for(int n = 0; n < value_count; ++n) {
      dst[n] = 0;
    }
  }
#endif
};

/*! \brief Function that determines the permutation indices of the coordinates.
 * \param num_parts is the number of parts.
 * \param mj_current_dim_coords is 1 dimensional array holding the
 * coordinate values.
 * \param current_concurrent_cut_coordinate is 1 dimensional array holding the
 * cut coordinates.
 * \param coordinate_begin is the start index of the given partition on
 * partitionedPointPermutations.
 * \param coordinate_end is the end index of the given partition on
 * partitionedPointPermutations.
 * \param used_local_cut_line_weight_to_left holds how much weight of the
 * coordinates on the cutline should be put on left side.
 * \param out_part_xadj is the indices of coordinates calculated for the
 * partition on next dimension.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
mj_create_new_partitions(
  mj_part_t num_parts,
  mj_part_t current_concurrent_work_part,
  Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
  Kokkos::View<mj_scalar_t *, device_t> & current_concurrent_cut_coordinate,
  Kokkos::View<mj_scalar_t *, device_t> & used_local_cut_line_weight_to_left,
  Kokkos::View<mj_lno_t *, device_t> & out_part_xadj)
{
  // Get locals for cuda
  auto local_thread_part_weight_work = this->thread_part_weight_work;
  auto local_point_counts = this->thread_point_counts;
  auto local_distribute_points_on_cut_lines =
    this->distribute_points_on_cut_lines;
  auto local_thread_cut_line_weight_to_put_left =
    this->thread_cut_line_weight_to_put_left;
  auto local_sEpsilon = this->sEpsilon;
  auto local_coordinate_permutations = this->coordinate_permutations;
  auto local_mj_weights = this->mj_weights;
  auto local_assigned_part_ids = this->assigned_part_ids;
  auto local_new_coordinate_permutations = this->new_coordinate_permutations;

  mj_part_t num_cuts = num_parts - 1;

  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
    KOKKOS_LAMBDA(int dummy) {

      if(local_distribute_points_on_cut_lines) {
        for(int i = 0; i < num_cuts; ++i) {
          mj_scalar_t left_weight = used_local_cut_line_weight_to_left(i);
          if(left_weight > local_sEpsilon) {
            // the weight of thread ii on cut.
            mj_scalar_t thread_ii_weight_on_cut =
              local_thread_part_weight_work(i * 2 + 1) -
              local_thread_part_weight_work(i * 2);

            if(thread_ii_weight_on_cut < left_weight) {
              // if left weight is bigger than threads weight on cut.
              local_thread_cut_line_weight_to_put_left(i) =
                thread_ii_weight_on_cut;
            }
            else {
              // if thread's weight is bigger than space, then put only a portion.
              local_thread_cut_line_weight_to_put_left(i) = left_weight;
            }
            left_weight -= thread_ii_weight_on_cut;
          }
          else {
            local_thread_cut_line_weight_to_put_left(i) = 0;
          }
        }

        // this is a special case. If cutlines share the same coordinate,
        // their weights are equal. We need to adjust the ratio for that.
        for(mj_part_t i = num_cuts - 1; i > 0 ; --i) {
          if(std::abs(current_concurrent_cut_coordinate(i) -
            current_concurrent_cut_coordinate(i -1)) < local_sEpsilon) {
              local_thread_cut_line_weight_to_put_left(i) -=
                local_thread_cut_line_weight_to_put_left(i - 1);
          }
          local_thread_cut_line_weight_to_put_left(i) =
            static_cast<long long>((local_thread_cut_line_weight_to_put_left(i) +
            least_signifiance) * significance_mul) /
            static_cast<mj_scalar_t>(significance_mul);
        }
      }

      for(mj_part_t i = 0; i < num_parts; ++i) {
        local_point_counts(i) = 0;
      }
  });

  mj_lno_t coordinate_begin_index =
    current_concurrent_work_part == 0 ? 0 :
    host_part_xadj(current_concurrent_work_part - 1);
  mj_lno_t coordinate_end_index =
    host_part_xadj(current_concurrent_work_part);

  mj_lno_t total_on_cut;
  Kokkos::parallel_reduce("Get total_on_cut",
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (
      coordinate_begin_index, coordinate_end_index),
    KOKKOS_LAMBDA(int ii, mj_lno_t & val) {
    mj_lno_t coordinate_index = local_coordinate_permutations(ii);
    mj_part_t coordinate_assigned_place =
      local_assigned_part_ids(coordinate_index);
    if(coordinate_assigned_place % 2 == 1) {
      val += 1;
    }
  }, total_on_cut);

  Kokkos::View<mj_lno_t *, device_t> track_on_cuts;
  if(total_on_cut > 0) {
    track_on_cuts = Kokkos::View<mj_lno_t *, device_t>(
      "track_on_cuts", // would do WithoutInitialization but need last init to 0
      total_on_cut + 1); // extra index to use for tracking
  }

  // here we need to parallel reduce an array to count coords in each part
  // atomically adding, especially for low part count would kill us
  // in the original setup we kept arrays allocated for each thread but for
  // the cuda version we'd like to avoid allocating N arrays for the number
  // of teams/threads which would be complicated based on running openmp or
  // cuda.
  typedef Kokkos::TeamPolicy<typename mj_node_t::execution_space> policy_t;

  // if not set use 60 - somewhat arbitrary based on initial performance tests
  int use_num_teams = mj_num_teams ? mj_num_teams : 60;

  auto policy_ReduceFunctor = policy_t(use_num_teams, Kokkos::AUTO);
  typedef int array_t;

  // just need parts - on the cuts will be handled in a separate serial
  // call after this.
#ifndef KOKKOS_ENABLE_CUDA
  array_t * reduce_array = new array_t[static_cast<size_t>(num_parts)];
#endif

  ReduceArrayFunctor<policy_t, mj_scalar_t, mj_part_t, mj_lno_t,
    typename mj_node_t::device_type, array_t>teamFunctor(
      current_concurrent_work_part,
      num_parts,
      coordinate_permutations,
      mj_current_dim_coords,
      assigned_part_ids,
      part_xadj,
      track_on_cuts
#ifdef KOKKOS_ENABLE_CUDA
      ,local_point_counts
#endif
      );

#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::parallel_for(policy_ReduceFunctor, teamFunctor);
#else
  Kokkos::parallel_reduce(policy_ReduceFunctor, teamFunctor, reduce_array);
#endif

#ifndef KOKKOS_ENABLE_CUDA
  for(mj_part_t part = 0; part < num_parts; ++part) {
    local_point_counts(part) = reduce_array[part];
  }
  delete [] reduce_array;
#endif

  // the last member is utility used for atomically inserting the values.
  // Sorting here avoids potential indeterminancy in the partitioning results
  if(track_on_cuts.size() > 0) { // size 0 unused, or size is minimum of 2
    auto track_on_cuts_sort = Kokkos::subview(track_on_cuts,
      std::pair<mj_lno_t, mj_lno_t>(0, track_on_cuts.size() - 1)); // do not sort last element
    Kokkos::sort(track_on_cuts_sort);
  }

  bool uniform_weights0 = this->mj_uniform_weights(0);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
    KOKKOS_LAMBDA (int dummy) {

    for(int j = 0; j < total_on_cut; ++j) {
      int ii = track_on_cuts(j);
      mj_lno_t coordinate_index = local_coordinate_permutations(ii);
      mj_scalar_t coordinate_weight = uniform_weights0 ? 1 :
        local_mj_weights(coordinate_index,0);
      mj_part_t coordinate_assigned_place =
        local_assigned_part_ids(coordinate_index);
      mj_part_t coordinate_assigned_part = coordinate_assigned_place / 2;
      // if it is on the cut.
      if(local_distribute_points_on_cut_lines &&
        local_thread_cut_line_weight_to_put_left(
          coordinate_assigned_part) > local_sEpsilon) {
        // if the rectilinear partitioning is allowed,
        // and the thread has still space to put on the left of the cut
        // then thread puts the vertex to left.
        local_thread_cut_line_weight_to_put_left(
          coordinate_assigned_part) -= coordinate_weight;
        // if putting the vertex to left increased the weight more
        // than expected, and if the next cut is on the same coordinate,
        // then we need to adjust how much weight next cut puts to its left as
        // well, in order to take care of the imbalance.
        if(local_thread_cut_line_weight_to_put_left(
            coordinate_assigned_part) < 0 && coordinate_assigned_part <
            num_cuts - 1 &&
            std::abs(current_concurrent_cut_coordinate(
            coordinate_assigned_part+1) -
            current_concurrent_cut_coordinate(
            coordinate_assigned_part)) < local_sEpsilon)
        {
          local_thread_cut_line_weight_to_put_left(
            coordinate_assigned_part + 1) +=
            local_thread_cut_line_weight_to_put_left(
            coordinate_assigned_part);
        }
        ++local_point_counts(coordinate_assigned_part);
        local_assigned_part_ids(coordinate_index) =
          coordinate_assigned_part;
      }
      else {
        // if there is no more space on the left, put the coordinate to the
        // right of the cut.
        ++coordinate_assigned_part;
        // this while loop is necessary when a line is partitioned into more
        // than 2 parts.
        while(local_distribute_points_on_cut_lines &&
          coordinate_assigned_part < num_cuts)
        {
          // traverse all the cut lines having the same partitiong
          if(std::abs(current_concurrent_cut_coordinate(
            coordinate_assigned_part) -
            current_concurrent_cut_coordinate(
              coordinate_assigned_part - 1)) < local_sEpsilon)
          {
            // if line has enough space on left, put it there.
            if(local_thread_cut_line_weight_to_put_left(
              coordinate_assigned_part) > local_sEpsilon &&
              local_thread_cut_line_weight_to_put_left(
                coordinate_assigned_part) >=
                std::abs(local_thread_cut_line_weight_to_put_left(
                  coordinate_assigned_part) - coordinate_weight))
            {
              local_thread_cut_line_weight_to_put_left(
                coordinate_assigned_part) -= coordinate_weight;
              // Again if it put too much on left of the cut,
              // update how much the next cut sharing the same coordinate will
              // put to its left.
              if(local_thread_cut_line_weight_to_put_left(
                coordinate_assigned_part) < 0 &&
                coordinate_assigned_part < num_cuts - 1 &&
                std::abs(current_concurrent_cut_coordinate(
                  coordinate_assigned_part+1) -
                current_concurrent_cut_coordinate(
                  coordinate_assigned_part)) < local_sEpsilon)
              {
                local_thread_cut_line_weight_to_put_left(
                  coordinate_assigned_part + 1) +=
                  local_thread_cut_line_weight_to_put_left(
                    coordinate_assigned_part);
              }
              break;
            }
          }
          else {
            break;
          }
          ++coordinate_assigned_part;
        }
        local_point_counts(coordinate_assigned_part) += 1;
        local_assigned_part_ids(coordinate_index) = coordinate_assigned_part;
      }
    }

    for(int j = 0; j < num_parts; ++j) {
      out_part_xadj(j) = local_point_counts(j);
      local_point_counts(j) = 0;

      if(j != 0) {
        out_part_xadj(j) += out_part_xadj(j - 1);
        local_point_counts(j) += out_part_xadj(j - 1);
      }
    }
  });

  // here we will determine insert indices for N teams
  // then all the teams can fill

#ifdef KOKKOS_ENABLE_CUDA

  // This is the fastest so far - just straight atomic writes for CUDA
  // However this is not a deterministic result since it is atomic.
  // The final result will be deterministic.
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t> (
    coordinate_begin_index, coordinate_end_index),
    KOKKOS_LAMBDA (mj_lno_t ii) {
    mj_lno_t i = local_coordinate_permutations(ii);
    mj_part_t p = local_assigned_part_ids(i);
    mj_lno_t idx = Kokkos::atomic_fetch_add(&local_point_counts(p), 1);
    local_new_coordinate_permutations(coordinate_begin_index + idx) = i;
  });

#else

#ifdef KOKKOS_ENABLE_OPENMP
  // will return and fix this - revert back to 1 for clear auto testing
  const int num_threads = 1; // Kokkos::OpenMP::impl_max_hardware_threads();
#else
  const int num_threads = 1;
#endif

  const int num_teams = 1; // cuda is handled above using a different format

  // allow init - we want all 0's first
  Kokkos::View<mj_lno_t*, device_t>
    point_counter("insert indices", num_teams * num_threads * num_parts);

  // count how many coords per thread
  // then we will fill each independently
  Kokkos::TeamPolicy<typename mj_node_t::execution_space>
    block_policy(num_teams, num_threads);
  typedef typename Kokkos::TeamPolicy<typename mj_node_t::execution_space>::
    member_type member_type;
  mj_lno_t range = coordinate_end_index - coordinate_begin_index;
  mj_lno_t block_size = range / num_teams + 1;
  Kokkos::parallel_for(block_policy, KOKKOS_LAMBDA(member_type team_member) {
    int team = team_member.league_rank();
    int team_offset = team * num_threads * num_parts;
    mj_lno_t begin = coordinate_begin_index + team * block_size;
    mj_lno_t end = begin + block_size;
    if(end > coordinate_end_index) {
      end = coordinate_end_index;
    }

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, begin, end),
      [=] (mj_lno_t ii) {
      int thread = team_member.team_rank();
      mj_lno_t i = local_coordinate_permutations(ii);
      mj_part_t p = local_assigned_part_ids(i);
      int index = team_offset + thread * num_parts + p;
      ++point_counter(index);
    });
  });

  // now prefix sum
  // we currently have the counts in the slots
  // we want the first counter for each part to be 0
  // then the rest should be the sum of all the priors
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
    KOKKOS_LAMBDA (int dummy) {
    int num_sets = point_counter.size() / num_parts;
    for(int set = num_sets - 1; set >= 1; set -=1) {
      int base = set * num_parts;
      for(int part = 0; part < num_parts; ++part) {
        point_counter(base + part) = point_counter(base + part - num_parts);
      }
    }

    for(int part = 0; part < num_parts; ++part) {
      point_counter(part) = 0;
    }

    for(int set = 1; set < num_sets; ++set) {
      int base = set * num_parts;
      for(int part = 0; part < num_parts; ++part) {
        point_counter(base + part) += point_counter(base + part - num_parts);
      }
    }
  });

  // now permute
  Kokkos::parallel_for(block_policy, KOKKOS_LAMBDA(member_type team_member) {
    int team = team_member.league_rank();
    int team_offset = team * num_threads * num_parts;
    mj_lno_t begin = coordinate_begin_index + team * block_size;
    mj_lno_t end = begin + block_size;
    if(end > coordinate_end_index) {
      end = coordinate_end_index;
    }
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, begin, end),
      [=] (mj_lno_t ii) {
      int thread = team_member.team_rank();
      mj_lno_t i = local_coordinate_permutations(ii);
      mj_part_t p = local_assigned_part_ids(i);
      int index = team_offset + thread * num_parts + p;
      int set_counter = (point_counter(index)++) + local_point_counts(p);
      local_new_coordinate_permutations(coordinate_begin_index + set_counter) = i;
    });
  });
#endif
}

/*! \brief Function that calculates the new coordinates for the cut lines.
 * Function is called inside the parallel region.
 * \param num_total_part is the sum of number of cutlines and number of parts.
 * Simply it is 2*P - 1.
 * \param num_cuts is the number of cut lines. P - 1.
 * \param max_coordinate is the maximum coordinate in the current range of
 * coordinates and in the current dimension.
 * \param min_coordinate is the maximum coordinate in the current range of
 * coordinates and in the current dimension.
 * \param global_total_weight is the global total weight in the current range of
 * coordinates.
 * \param used_imbalance_tolerance is the maximum allowed imbalance ratio.
 * \param current_global_part_weights is the array holding the weight of parts.
 * Assumes there are 2*P - 1 parts (cut lines are seperate parts).
 * \param current_local_part_weights is the local totalweight of the processor.
 * \param current_part_target_weights are the desired cumulative part ratios,
 * sized P.
 * \param current_cut_line_determined is the boolean array to determine if the
 * correct position for a cut line is found.
 * \param current_cut_coordinates is the array holding the coordinates of each
 * cut line. Sized P - 1.
 * \param current_cut_upper_bounds is the array holding the upper bound
 * coordinate for each cut line. Sized P - 1.
 * \param current_cut_lower_bounds is the array holding the lower bound
 * coordinate for each cut line. Sized P - 1.
 * \param current_global_left_closest_points is the array holding the closest
 * points to the cut lines from left.
 * \param current_global_right_closest_points is the array holding the closest
 * points to the cut lines from right.
 * \param current_cut_lower_bound_weights is the array holding the weight of the
 * parts at the left of lower bound coordinates.
 * \param current_cut_upper_weights is the array holding the weight of the parts
 * at the left of upper bound coordinates.
 * \param new_current_cut_coordinates is the work array, sized P - 1.
 * \param current_part_cut_line_weight_ratio holds how much weight of the
 * coordinates on the cutline should be put on left side.
 * \param rectilinear_cut_count is the count of cut lines whose balance can be
 * achived via distributing the points in same coordinate to different parts.
 * \param my_num_incomplete_cut is the number of cutlines whose position has not
 * been determined yet. For K > 1 it is the count in a single part (whose cut
 * lines are determined).
 * DOCWORK: Check all documentation in this method
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t,
  mj_node_t>::mj_get_new_cut_coordinates(
  mj_part_t current_concurrent_num_parts,
  mj_part_t kk,
  const mj_part_t &num_cuts,
  const double &used_imbalance_tolerance,
  Kokkos::View<mj_scalar_t *, device_t> & current_global_part_weights,
  Kokkos::View<mj_scalar_t *, device_t> & current_local_part_weights,
  Kokkos::View<mj_scalar_t *, device_t> & current_part_target_weights,
  Kokkos::View<bool *, device_t> & current_cut_line_determined,
  Kokkos::View<mj_scalar_t *, device_t> & current_cut_coordinates,
  Kokkos::View<mj_scalar_t *, device_t> & current_cut_upper_bounds,
  Kokkos::View<mj_scalar_t *, device_t> & current_cut_lower_bounds,
  Kokkos::View<mj_scalar_t *, device_t> & current_global_left_closest_points,
  Kokkos::View<mj_scalar_t *, device_t> & current_global_right_closest_points,
  Kokkos::View<mj_scalar_t *, device_t> & current_cut_lower_bound_weights,
  Kokkos::View<mj_scalar_t *, device_t> & current_cut_upper_weights,
  Kokkos::View<mj_scalar_t *, device_t> & new_current_cut_coordinates,
  Kokkos::View<mj_scalar_t *, device_t> &
    current_part_cut_line_weight_to_put_left,
  Kokkos::View<mj_part_t *, device_t> & view_rectilinear_cut_count)
{
  Kokkos::deep_copy(device_incomplete_cut_count, this->incomplete_cut_count);

  auto local_device_incomplete_cut_count = device_incomplete_cut_count;
  auto local_sEpsilon = sEpsilon;
  auto local_distribute_points_on_cut_lines = distribute_points_on_cut_lines;
  auto local_global_rectilinear_cut_weight = global_rectilinear_cut_weight;
  auto local_process_rectilinear_cut_weight = process_rectilinear_cut_weight;
  auto local_global_min_max_coord_total_weight =
    global_min_max_coord_total_weight;

  // Note for a 22 part system I tried removing the outer loop
  // and doing each sub loop as a simple parallel_for over num_cuts.
  // But that was about twice as slow (10ms) as the current form (5ms)
  // so I think the overhead of launching the new global parallel kernels
  // is costly. This form is just running one team so effectively using
  // a single warp to process the cuts. I expect with a lot of parts this
  // might need changing.
  Kokkos::TeamPolicy<typename mj_node_t::execution_space>
    policy_one_team(1, Kokkos::AUTO());
  typedef typename Kokkos::TeamPolicy<typename mj_node_t::execution_space>::
    member_type member_type;
  Kokkos::parallel_for(policy_one_team, KOKKOS_LAMBDA(member_type team_member) {

    mj_scalar_t min_coordinate =
      local_global_min_max_coord_total_weight(kk);
    mj_scalar_t max_coordinate =
      local_global_min_max_coord_total_weight(
      kk + current_concurrent_num_parts);
    mj_scalar_t global_total_weight =
      local_global_min_max_coord_total_weight(
      kk + current_concurrent_num_parts * 2);

    Kokkos::parallel_for(Kokkos::TeamThreadRange (team_member, num_cuts),
      [=] (mj_part_t i) {
      // if left and right closest points are not set yet,
      // set it to the cut itself.
      if(min_coordinate -
        current_global_left_closest_points(i) > local_sEpsilon) {
        current_global_left_closest_points(i) =
          current_cut_coordinates(i);
      }
      if(current_global_right_closest_points(i) -
        max_coordinate > local_sEpsilon) {
        current_global_right_closest_points(i) =
          current_cut_coordinates(i);
      }
    });
    team_member.team_barrier(); // for end of Kokkos::TeamThreadRange

    Kokkos::parallel_for(Kokkos::TeamThreadRange (team_member, num_cuts),
      [=] (mj_part_t i) {
      // seen weight in the part
      mj_scalar_t seen_weight_in_part = 0;
      // expected weight for part.
      mj_scalar_t expected_weight_in_part = 0;
      // imbalance for the left and right side of the cut.
      double imbalance_on_left = 0, imbalance_on_right = 0;
      if(local_distribute_points_on_cut_lines) {
        // init the weight on the cut.
        local_global_rectilinear_cut_weight(i) = 0;
        local_process_rectilinear_cut_weight(i) = 0;
      }
      bool bContinue = false;
      // if already determined at previous iterations,
      // then just write the coordinate to new array, and proceed.
      if(current_cut_line_determined(i)) {
        new_current_cut_coordinates(i) =
          current_cut_coordinates(i);
        bContinue = true;
      }
      if(!bContinue) {
        //current weight of the part at the left of the cut line.
        seen_weight_in_part = current_global_part_weights(i * 2);

        //expected ratio
        expected_weight_in_part = current_part_target_weights(i);

       //leftImbalance = imbalanceOf(seenW, globalTotalWeight, expected);
        imbalance_on_left = calculate_imbalance(seen_weight_in_part,
          expected_weight_in_part);
        // rightImbalance = imbalanceOf(globalTotalWeight - seenW,
        // globalTotalWeight, 1 - expected);
        imbalance_on_right = calculate_imbalance(global_total_weight -
          seen_weight_in_part, global_total_weight - expected_weight_in_part);
        bool is_left_imbalance_valid = std::abs(imbalance_on_left) -
          used_imbalance_tolerance < local_sEpsilon ;
        bool is_right_imbalance_valid = std::abs(imbalance_on_right) -
          used_imbalance_tolerance < local_sEpsilon;
        //if the cut line reaches to desired imbalance.
        if(is_left_imbalance_valid && is_right_imbalance_valid) {
          current_cut_line_determined(i) = true;
          Kokkos::atomic_add(&local_device_incomplete_cut_count(kk), -1);
          new_current_cut_coordinates(i) = current_cut_coordinates(i);
        }
        else if(imbalance_on_left < 0) {
          //if left imbalance < 0 then we need to move the cut to right.
          if(local_distribute_points_on_cut_lines) {
            // if it is okay to distribute the coordinate on
            // the same coordinate to left and right.
            // then check if we can reach to the target weight by including the
            // coordinates in the part.
            if(current_global_part_weights(i * 2 + 1) ==
              expected_weight_in_part) {
              // if it is we are done.
              current_cut_line_determined(i) = true;
              Kokkos::atomic_add(&local_device_incomplete_cut_count(kk), -1);

              //then assign everything on the cut to the left of the cut.
              new_current_cut_coordinates(i) =
                current_cut_coordinates(i);
              //for this cut all the weight on cut will be put to left.
              current_part_cut_line_weight_to_put_left(i) =
                current_local_part_weights(i * 2 + 1) -
                current_local_part_weights(i * 2);
              bContinue = true;
            }
            else if(current_global_part_weights(i * 2 + 1) >
              expected_weight_in_part) {
              // if the weight is larger than the expected weight,
              // then we need to distribute some points to left, some to right.
              current_cut_line_determined(i) = true;
              Kokkos::atomic_add(&view_rectilinear_cut_count(0), 1);

              // increase the num cuts to be determined with rectilinear
              // partitioning.
              Kokkos::atomic_add(&local_device_incomplete_cut_count(kk), -1);
              new_current_cut_coordinates(i) =
                current_cut_coordinates(i);
              local_process_rectilinear_cut_weight[i] =
                current_local_part_weights(i * 2 + 1) -
                current_local_part_weights(i * 2);
              bContinue = true;
            }
          }

          if(!bContinue) {

            // we need to move further right,so set lower bound to current line,
            // and shift it to the closes point from right.
            current_cut_lower_bounds(i) =
              current_global_right_closest_points(i);

            //set the lower bound weight to the weight we have seen.
            current_cut_lower_bound_weights(i) = seen_weight_in_part;

            // compare the upper bound with what has been found in the
            // last iteration.
            // we try to make more strict bounds for the cut here.
            for(mj_part_t ii = i + 1; ii < num_cuts ; ++ii) {
              mj_scalar_t p_weight = current_global_part_weights(ii * 2);
              mj_scalar_t line_weight =
                current_global_part_weights(ii * 2 + 1);
              if(p_weight >= expected_weight_in_part) {
                // if a cut on the right has the expected weight, then we found
                // our cut position. Set up and low coordiantes to this
                // new cut coordinate, but we need one more iteration to
                // finalize the cut position, as wee need to update the part ids.
                if(p_weight == expected_weight_in_part) {
                  current_cut_upper_bounds(i) =
                    current_cut_coordinates(ii);
                  current_cut_upper_weights(i) = p_weight;
                  current_cut_lower_bounds(i) =
                    current_cut_coordinates(ii);
                  current_cut_lower_bound_weights(i) = p_weight;
                } else if(p_weight < current_cut_upper_weights(i)) {
                  // if a part weight is larger then my expected weight,
                  // but lower than my upper bound weight, update upper bound.
                  current_cut_upper_bounds(i) =
                    current_global_left_closest_points(ii);
                  current_cut_upper_weights(i) = p_weight;
                }
                break;
              }
              // if comes here then pw < ew
              // then compare the weight against line weight.
              if(line_weight >= expected_weight_in_part) {
                // if the line is larger than the expected weight, then we need
                // to reach to the balance by distributing coordinates on
                // this line.
                current_cut_upper_bounds(i) =
                  current_cut_coordinates(ii);
                current_cut_upper_weights(i) = line_weight;
                current_cut_lower_bounds(i) =
                  current_cut_coordinates(ii);
                current_cut_lower_bound_weights(i) = p_weight;
                break;
              }
              // if a stricter lower bound is found,
              // update the lower bound.
              if(p_weight <= expected_weight_in_part && p_weight >=
                current_cut_lower_bound_weights(i)) {
                current_cut_lower_bounds(i) =
                  current_global_right_closest_points(ii);
                current_cut_lower_bound_weights(i) = p_weight;
              }
            }

            mj_scalar_t new_cut_position = 0;
            this->mj_calculate_new_cut_position(
              current_cut_upper_bounds(i),
              current_cut_lower_bounds(i),
              current_cut_upper_weights(i),
              current_cut_lower_bound_weights(i),
              expected_weight_in_part, new_cut_position);

            // if cut line does not move significantly.
            // then finalize the search.
            if(std::abs(current_cut_coordinates(i) -
              new_cut_position) < local_sEpsilon) {
              current_cut_line_determined(i) = true;
              Kokkos::atomic_add(&local_device_incomplete_cut_count(kk), -1);

              //set the cut coordinate and proceed.
              new_current_cut_coordinates(i) =
                current_cut_coordinates(i);
            } else {
              new_current_cut_coordinates(i) = new_cut_position;
            }
          } // bContinue
        } else {
          // need to move the cut line to left.
          // set upper bound to current line.
          current_cut_upper_bounds(i) =
            current_global_left_closest_points(i);
          current_cut_upper_weights(i) =
            seen_weight_in_part;
          // compare the current cut line weights with
          // previous upper and lower bounds.
          for(int ii = i - 1; ii >= 0; --ii) {
            mj_scalar_t p_weight =
              current_global_part_weights(ii * 2);
            mj_scalar_t line_weight =
              current_global_part_weights(ii * 2 + 1);
            if(p_weight <= expected_weight_in_part) {
              if(p_weight == expected_weight_in_part) {
                // if the weight of the part is my expected weight
                // then we find the solution.
                current_cut_upper_bounds(i) =
                  current_cut_coordinates(ii);
                current_cut_upper_weights(i) = p_weight;
                current_cut_lower_bounds(i) =
                  current_cut_coordinates(ii);
                current_cut_lower_bound_weights(i) = p_weight;
              }
              else if(p_weight > current_cut_lower_bound_weights(i)) {
                // if found weight is bigger than the lower bound
                // then update the lower bound.
                current_cut_lower_bounds(i) =
                  current_global_right_closest_points(ii);
                current_cut_lower_bound_weights(i) = p_weight;

                // at the same time, if weight of line is bigger than the
                // expected weight, then update the upper bound as well.
                // in this case the balance will be obtained by distributing
                // weights on this cut position.
                if(line_weight > expected_weight_in_part) {
                  current_cut_upper_bounds(i) =
                    current_global_right_closest_points(ii);
                  current_cut_upper_weights(i) = line_weight;
                }
              }
              break;
            }
            // if the weight of the cut on the left is still bigger than
            // my weight, and also if the weight is smaller than the current
            // upper weight, or if the weight is equal to current upper
            // weight, but on the left of the upper weight, then update
            // upper bound.
            if(p_weight >= expected_weight_in_part &&
              (p_weight < current_cut_upper_weights(i) ||
              (p_weight == current_cut_upper_weights(i) &&
                current_cut_upper_bounds(i) >
                  current_global_left_closest_points(ii)))) {
              current_cut_upper_bounds(i) =
                current_global_left_closest_points(ii);
              current_cut_upper_weights(i) = p_weight;
            }
          }
          mj_scalar_t new_cut_position = 0;
          this->mj_calculate_new_cut_position(
            current_cut_upper_bounds(i),
            current_cut_lower_bounds(i),
            current_cut_upper_weights(i),
            current_cut_lower_bound_weights(i),
            expected_weight_in_part,
            new_cut_position);

            // if cut line does not move significantly.
            if(std::abs(current_cut_coordinates(i) -
              new_cut_position) < local_sEpsilon) {
              current_cut_line_determined(i) = true;
              Kokkos::atomic_add(&local_device_incomplete_cut_count(kk), -1);
              //set the cut coordinate and proceed.
              new_current_cut_coordinates(i) =
                current_cut_coordinates(i);
            } else {
              new_current_cut_coordinates(i) =
                new_cut_position;
            }
          }
        }; // bContinue
      });

      team_member.team_barrier(); // for end of Kokkos::TeamThreadRange
  });

  // view_rectilinear_cut_count
  mj_part_t rectilinear_cut_count;
  Kokkos::parallel_reduce("Read bDoingWork",
    Kokkos::RangePolicy<typename mj_node_t::execution_space, int>(0, 1),
    KOKKOS_LAMBDA(int dummy, int & set_single) {
    set_single = view_rectilinear_cut_count(0);
  }, rectilinear_cut_count);

  if(rectilinear_cut_count > 0) {
    auto host_local_process_rectilinear_cut_weight =
      Kokkos::create_mirror_view(Kokkos::HostSpace(),
      local_process_rectilinear_cut_weight);
    auto host_local_global_rectilinear_cut_weight =
      Kokkos::create_mirror_view(Kokkos::HostSpace(),
      local_global_rectilinear_cut_weight);
    Kokkos::deep_copy(host_local_process_rectilinear_cut_weight,
      local_process_rectilinear_cut_weight);
    Kokkos::deep_copy(host_local_global_rectilinear_cut_weight,
      local_global_rectilinear_cut_weight);
    Teuchos::scan<int,mj_scalar_t>(
      *comm, Teuchos::REDUCE_SUM,
      num_cuts,
      host_local_process_rectilinear_cut_weight.data(),
      host_local_global_rectilinear_cut_weight.data());
    Kokkos::deep_copy(local_process_rectilinear_cut_weight,
      host_local_process_rectilinear_cut_weight);
    Kokkos::deep_copy(local_global_rectilinear_cut_weight,
      host_local_global_rectilinear_cut_weight);

    Kokkos::parallel_for("finish up mj_get_new_cut_coordinates",
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
      KOKKOS_LAMBDA(int dummy) {
      for(mj_part_t i = 0; i < num_cuts; ++i) {
        // if cut line weight to be distributed.
        if(local_global_rectilinear_cut_weight(i) > 0) {
          // expected weight to go to left of the cut.
          mj_scalar_t expected_part_weight = current_part_target_weights(i);
          // the weight that should be put to left of the cut.
          mj_scalar_t necessary_weight_on_line_for_left =
            expected_part_weight - current_global_part_weights(i * 2);

          // the weight of the cut in the process
          mj_scalar_t my_weight_on_line =
            local_process_rectilinear_cut_weight(i);

          // the sum of the cut weights upto this process,
          // including the weight of this process.
          mj_scalar_t weight_on_line_upto_process_inclusive =
            local_global_rectilinear_cut_weight(i);
          // the space on the left side of the cut after all processes
          // before this process (including this process)
          // puts their weights on cut to left.
          mj_scalar_t space_to_put_left =
            necessary_weight_on_line_for_left -
            weight_on_line_upto_process_inclusive;
          // add my weight to this space to find out how much space
          // is left to me.
          mj_scalar_t space_left_to_me =
            space_to_put_left + my_weight_on_line;

          /*
          cout << "expected_part_weight:" << expected_part_weight
            << " necessary_weight_on_line_for_left:"
            << necessary_weight_on_line_for_left
            << " my_weight_on_line" << my_weight_on_line
            << " weight_on_line_upto_process_inclusive:"
            << weight_on_line_upto_process_inclusive
            << " space_to_put_left:" << space_to_put_left
            << " space_left_to_me" << space_left_to_me << endl;
           */

          if(space_left_to_me < 0) {
            // space_left_to_me is negative and i dont need to put
            // anything to left.
            current_part_cut_line_weight_to_put_left(i) = 0;
          }
          else if(space_left_to_me >= my_weight_on_line) {
            // space left to me is bigger than the weight of the
            // processor on cut.
            // so put everything to left.
            current_part_cut_line_weight_to_put_left(i) =
              my_weight_on_line;
            // cout << "setting current_part_cut_line_weight_to_put_left
            // to my_weight_on_line:" << my_weight_on_line << endl;
          }
          else {
            // put only the weight as much as the space.
            current_part_cut_line_weight_to_put_left(i) =
              space_left_to_me;
            // cout << "setting current_part_cut_line_weight_to_put_left
            // to space_left_to_me:" << space_left_to_me << endl;
          }
        }
      }
      view_rectilinear_cut_count(0) = 0;
    });
  }

  Kokkos::deep_copy(this->incomplete_cut_count, device_incomplete_cut_count);
}

/*! \brief Function fills up the num_points_in_all_processor_parts, so that
 * it has the number of coordinates in each processor of each part.
 * to access how many points processor i has on part j,
 * num_points_in_all_processor_parts[i * num_parts + j].
 * \param num_procs is the number of processor attending to migration operation.
 * \param num_parts is the number of parts that exist in current partitioning.
 * \param num_points_in_all_processor_parts is the output array that holds
 * the number of coordinates in each part in each processor.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  get_processor_num_points_in_parts(
  mj_part_t num_procs,
  mj_part_t num_parts,
  mj_gno_t *&num_points_in_all_processor_parts)
{
  // initially allocation_size is num_parts
  size_t allocation_size = num_parts * (num_procs + 1);

  // this will be output
  // holds how many each processor has in each part.
  // last portion is the sum of all processor points in each part.

  // allocate memory for the local num coordinates in each part.
  mj_gno_t *num_local_points_in_each_part_to_reduce_sum =
    new mj_gno_t[allocation_size];

  // this is the portion of the memory which will be used
  // at the summation to obtain total number of processors' points in each part.
  mj_gno_t *my_local_points_to_reduce_sum =
    num_local_points_in_each_part_to_reduce_sum + num_procs * num_parts;

  // this is the portion of the memory where each stores its local number.
  // this information is needed by other processors.
  mj_gno_t *my_local_point_counts_in_each_part =
    num_local_points_in_each_part_to_reduce_sum + this->myRank * num_parts;

  // initialize the array with 0's.
  memset(num_local_points_in_each_part_to_reduce_sum, 0,
    sizeof(mj_gno_t)*allocation_size);

  auto local_new_part_xadj = this->new_part_xadj;
  Kokkos::View<mj_gno_t *, typename mj_node_t::device_type> points_per_part(
    Kokkos::ViewAllocateWithoutInitializing("points per part"), num_parts);
  Kokkos::parallel_for("get vals on device",
    Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_gno_t>
    (0, num_parts), KOKKOS_LAMBDA(mj_gno_t i) {
    points_per_part(i) =
      local_new_part_xadj(i) - ((i == 0) ? 0 : local_new_part_xadj(i-1));
  });
  auto host_points_per_part = Kokkos::create_mirror_view(points_per_part);
  Kokkos::deep_copy(host_points_per_part, points_per_part);
  for(int i = 0; i < num_parts; ++i) {
    my_local_points_to_reduce_sum[i] = host_points_per_part(i);
  }

  // copy the local num parts to the last portion of array, so that this portion
  // will represent the global num points in each part after the reduction.
  memcpy (my_local_point_counts_in_each_part, my_local_points_to_reduce_sum,
    sizeof(mj_gno_t) * (num_parts) );

  // reduceAll operation.
  // the portion that belongs to a processor with index p
  // will start from myRank * num_parts.
  // the global number of points will be held at the index
  try{
    reduceAll<int, mj_gno_t>(
      *(this->comm),
      Teuchos::REDUCE_SUM,
      allocation_size,
      num_local_points_in_each_part_to_reduce_sum,
      num_points_in_all_processor_parts);
  }
  Z2_THROW_OUTSIDE_ERROR(*(this->mj_env))

  delete [] num_local_points_in_each_part_to_reduce_sum;
}

/*! \brief Function checks if should do migration or not.
 * It returns true to point that migration should be done when
 * -migration_reduce_all_population are higher than a predetermined value
 * -num_coords_for_last_dim_part that left for the last dimension partitioning
 * is less than a predetermined value - the imbalance of the processors on the
 * parts are higher than given threshold.
 * \param migration_reduce_all_population is the multiplication of the number of
 * reduceall operations estimated and the number of processors.
 * \param num_coords_for_last_dim_part is the estimated number of coordinates in
 * a part per processor in the last dimension partitioning.
 * \param num_procs is number of processors attending to migration operation.
 * \param num_parts is number of parts that exist in the current partitioning.
 * \param num_points_in_all_processor_parts is the input array that holds
 * the number of coordinates in each part in each processor.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
bool AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_check_to_migrate(
  size_t migration_reduce_all_population,
  mj_lno_t num_coords_for_last_dim_part,
  mj_part_t num_procs,
  mj_part_t num_parts,
  mj_gno_t *num_points_in_all_processor_parts)
{
  // if reduce all count and population in the last dim is too high
  if(migration_reduce_all_population > future_reduceall_cutoff) {
    return true;
  }

  // if the work in a part per processor in the last dim is too low.
  if(num_coords_for_last_dim_part < min_work_last_dim) {
    return true;
  }

  // if migration is to be checked and the imbalance is too high
  if(this->check_migrate_avoid_migration_option == 0) {
    double global_imbalance = 0;
    // global shift to reach the sum of coordiante count in each part.
    size_t global_shift = num_procs * num_parts;

    for(mj_part_t ii = 0; ii < num_procs; ++ii) {
      for(mj_part_t i = 0; i < num_parts; ++i) {
       double ideal_num = num_points_in_all_processor_parts[global_shift + i]
         / double(num_procs);

       global_imbalance += std::abs(ideal_num -
         num_points_in_all_processor_parts[ii * num_parts + i]) /  (ideal_num);
      }
    }
    global_imbalance /= num_parts;
    global_imbalance /= num_procs;

    if(global_imbalance <= this->minimum_migration_imbalance) {
      return false;
    }
    else {
      return true;
    }
  }
  else {
    // if migration is forced
    return true;
  }
}

/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent.
 * \param num_parts is the number of parts that exist in the current
 * partitioning.
 * \param part_assignment_proc_begin_indices ([i]) points to the first processor
 * index that part i will be sent to.
 * \param processor_chains_in_parts the array that holds the linked list
 * structure, started from part_assignment_proc_begin_indices ([i]).
 * \param send_count_to_each_proc array array storing the number of points to
 * be sent to each part.
 * \param coordinate_destinations is the output array that holds which part
 * each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  assign_send_destinations(
  mj_part_t num_parts,
  mj_part_t *part_assignment_proc_begin_indices,
  mj_part_t *processor_chains_in_parts,
  mj_lno_t *send_count_to_each_proc,
  int *coordinate_destinations) {

  auto host_new_part_xadj = Kokkos::create_mirror_view(this->new_part_xadj);
  deep_copy(host_new_part_xadj, this->new_part_xadj);

  auto host_new_coordinate_permutations =
    Kokkos::create_mirror_view(this->new_coordinate_permutations);
  deep_copy(host_new_coordinate_permutations, this->new_coordinate_permutations);

  for(mj_part_t p = 0; p < num_parts; ++p) {
    mj_lno_t part_begin = 0;
    if(p > 0) part_begin = host_new_part_xadj(p - 1);
    mj_lno_t part_end = host_new_part_xadj(p);
    // get the first part that current processor will send its part-p.
    mj_part_t proc_to_sent = part_assignment_proc_begin_indices[p];
    // initialize how many point I sent to this processor.
    mj_lno_t num_total_send = 0;
    for(mj_lno_t j=part_begin; j < part_end; j++) {
      mj_lno_t local_ind = host_new_coordinate_permutations(j);
      while (num_total_send >= send_count_to_each_proc[proc_to_sent]) {
        // then get the next processor to send the points in part p.
        num_total_send = 0;
        // assign new processor to part_assign_begin[p]
        part_assignment_proc_begin_indices[p] =
          processor_chains_in_parts[proc_to_sent];
        // remove the previous processor
        processor_chains_in_parts[proc_to_sent] = -1;
        // choose the next processor as the next one to send.
        proc_to_sent = part_assignment_proc_begin_indices[p];
      }
      // write the gno index to corresponding position in sendBuf.
      coordinate_destinations[local_ind] = proc_to_sent;
      ++num_total_send;
    }
  }
}

/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent.
 * \param num_points_in_all_processor_parts is the array holding the num points
 * in each part in each proc.
 * \param num_parts is the number of parts that exist in the
 * current partitioning.
 * \param num_procs is the number of processor attending to migration operation.
 * \param send_count_to_each_proc array array storing the number of points to
 * be sent to each part.
 * \param processor_ranks_for_subcomm is the ranks of the processors that will
 * be in the subcommunicator with me.
 * \param next_future_num_parts_in_parts is the vector, how many more parts
 * each part will be divided into in the future.
 * \param out_part_index is the index of the part to which the processor
 * is assigned.
 * \param output_part_numbering_begin_index is how much the numbers should
 * be shifted when numbering the result parts.
 * \param coordinate_destinations is the output array that holds which part
 * each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_assign_proc_to_parts(
  mj_gno_t * num_points_in_all_processor_parts,
  mj_part_t num_parts,
  mj_part_t num_procs,
  mj_lno_t *send_count_to_each_proc,
  std::vector<mj_part_t> &processor_ranks_for_subcomm,
  std::vector<mj_part_t> *next_future_num_parts_in_parts,
  mj_part_t &out_part_index,
  mj_part_t &output_part_numbering_begin_index,
  int * coordinate_destinations) {
  mj_gno_t *global_num_points_in_parts =
    num_points_in_all_processor_parts + num_procs * num_parts;
  mj_part_t *num_procs_assigned_to_each_part = new mj_part_t[num_parts];

  // boolean variable if the process finds its part to be assigned.
  bool did_i_find_my_group = false;

  mj_part_t num_free_procs = num_procs;
  mj_part_t minimum_num_procs_required_for_rest_of_parts = num_parts - 1;

  double max_imbalance_difference = 0;
  mj_part_t max_differing_part = 0;

  // find how many processor each part requires.
  for(mj_part_t i = 0; i < num_parts; i++) {

    // scalar portion of the required processors
    double scalar_required_proc = num_procs *
      (double (global_num_points_in_parts[i]) /
      double (this->num_global_coords));

    // round it to closest integer; make sure have at least one proc.
    mj_part_t required_proc =
      static_cast<mj_part_t> (0.5 + scalar_required_proc);
    if(required_proc == 0) required_proc = 1;

    // if assigning the required num procs, creates problems for the rest
    // of the parts, then only assign {num_free_procs -
    // (minimum_num_procs_required_for_rest_of_parts)} procs to this part.
    if(num_free_procs -
      required_proc < minimum_num_procs_required_for_rest_of_parts)  {
        required_proc = num_free_procs -
          (minimum_num_procs_required_for_rest_of_parts);
      }

      // reduce the free processor count
      num_free_procs -= required_proc;

      // reduce the free minimum processor count required for the rest of the
      // part by 1.
      --minimum_num_procs_required_for_rest_of_parts;

      // part (i) is assigned to (required_proc) processors.
      num_procs_assigned_to_each_part[i] = required_proc;

      // because of the roundings some processors might be left as unassigned.
      // we want to assign those processors to the part with most imbalance.
      // find the part with the maximum imbalance here.
      double imbalance_wrt_ideal =
        (scalar_required_proc - required_proc) /  required_proc;
      if(imbalance_wrt_ideal > max_imbalance_difference) {
        max_imbalance_difference = imbalance_wrt_ideal;
        max_differing_part = i;
      }
    }

    // assign extra processors to the part with maximum imbalance
    // than the ideal.
    if(num_free_procs > 0) {
      num_procs_assigned_to_each_part[max_differing_part] +=  num_free_procs;
    }

    // now find what are the best processors with least migration for each part.

    // part_assignment_proc_begin_indices ([i]) is the array that holds the
    // beginning index of a processor that processor sends its data for part - i
    mj_part_t *part_assignment_proc_begin_indices = new mj_part_t[num_parts];

    // the next processor send is found in processor_chains_in_parts,
    // in linked list manner.
    mj_part_t *processor_chains_in_parts = new mj_part_t [num_procs];
    mj_part_t *processor_part_assignments = new mj_part_t[num_procs];

    // initialize the assignment of each processor.
    // this has a linked list implementation.
    // the beginning of processors assigned
    // to each part is hold at  part_assignment_proc_begin_indices[part].
    // then the next processor assigned to that part is located at
    // proc_part_assignments[part_assign_begins[part]], this is a chain
    // until the value of -1 is reached.
    for(int i = 0; i < num_procs; ++i ) {
      processor_part_assignments[i] = -1;
      processor_chains_in_parts[i] = -1;
    }
    for(int i = 0; i < num_parts; ++i ) {
      part_assignment_proc_begin_indices[i] = -1;
    }

    // std::cout << "Before migration: mig type:" <<
    //   this->migration_type << std::endl;
    // Allocate memory for sorting data structure.
    uSignedSortItem<mj_part_t, mj_gno_t, char> *
      sort_item_num_part_points_in_procs =
       new uSignedSortItem<mj_part_t, mj_gno_t, char>[num_procs];

    for(mj_part_t i = 0; i < num_parts; ++i) {
      // the algorithm tries to minimize the cost of migration, by assigning the
      // processors with highest number of coordinates on that part.
      // here we might want to implement a maximum weighted bipartite matching
      // algorithm.
      for(mj_part_t ii = 0; ii < num_procs; ++ii) {
        sort_item_num_part_points_in_procs[ii].id = ii;
        // if processor is not assigned yet.
        // add its num points to the sort data structure.
        if(processor_part_assignments[ii] == -1) {
          sort_item_num_part_points_in_procs[ii].val =
            num_points_in_all_processor_parts[ii * num_parts + i];
          // indicate that the processor has positive weight.
          sort_item_num_part_points_in_procs[ii].signbit = 1;
        }
        else {
          // if processor is already assigned, insert -nLocal - 1 so that it
          // won't be selected again.
          // would be same if we simply set it to -1, but more information with
          // no extra cost (which is used later) is provided.
          // sort_item_num_part_points_in_procs[ii].val =
          // -num_points_in_all_processor_parts[ii * num_parts + i] - 1;

          // UPDATE: Since above gets warning when unsigned is used to
          // represent, we added extra bit to as sign bit to the sort item.
          // It is 1 for positives, 0 for negatives.
          sort_item_num_part_points_in_procs[ii].val =
            num_points_in_all_processor_parts[ii * num_parts + i];
          sort_item_num_part_points_in_procs[ii].signbit = 0;
        }
      }

      // sort the processors in the part.
      uqSignsort<mj_part_t, mj_gno_t,char>
        (num_procs, sort_item_num_part_points_in_procs);

      /*
      for(mj_part_t ii = 0; ii < num_procs; ++ii) {
        std::cout << "ii:" << ii << " " <<
          sort_item_num_part_points_in_procs[ii].id <<
          " " << sort_item_num_part_points_in_procs[ii].val <<
          " " << int(sort_item_num_part_points_in_procs[ii].signbit) <<
          std::endl;
      }
      */

      mj_part_t required_proc_count = num_procs_assigned_to_each_part[i];
      mj_gno_t total_num_points_in_part = global_num_points_in_parts[i];
      mj_gno_t ideal_num_points_in_a_proc = Teuchos::as<mj_gno_t>(
        ceil(total_num_points_in_part / double (required_proc_count)));

      // starts sending to least heaviest part.
      mj_part_t next_proc_to_send_index = num_procs - required_proc_count;
      mj_part_t next_proc_to_send_id =
        sort_item_num_part_points_in_procs[next_proc_to_send_index].id;
      mj_lno_t space_left_in_sent_proc = ideal_num_points_in_a_proc -
        sort_item_num_part_points_in_procs[next_proc_to_send_index].val;

      // find the processors that will be assigned to this part, which are the
      // heaviest non assigned processors.
      for(mj_part_t ii = num_procs - 1;
        ii >= num_procs - required_proc_count; --ii) {
        mj_part_t proc_id = sort_item_num_part_points_in_procs[ii].id;
        // assign processor to part - i.
        processor_part_assignments[proc_id] = i;
      }

      bool did_change_sign = false;
      // if processor has a minus count, reverse it.
      for(mj_part_t ii = 0; ii < num_procs; ++ii) {
        // TODO:  THE LINE BELOW PRODUCES A WARNING IF gno_t IS UNSIGNED
        // TODO:  SEE BUG 6194
        if(sort_item_num_part_points_in_procs[ii].signbit == 0) {
          did_change_sign = true;
          sort_item_num_part_points_in_procs[ii].signbit = 1;
        }
        else {
          break;
        }
      }

      if(did_change_sign) {
        // resort the processors in the part for the rest of the processors that
        // is not assigned.
        uqSignsort<mj_part_t, mj_gno_t>(num_procs - required_proc_count,
          sort_item_num_part_points_in_procs);
      }

      /*
      for(mj_part_t ii = 0; ii < num_procs; ++ii) {
        std::cout << "after resort ii:" << ii << " " <<
          sort_item_num_part_points_in_procs[ii].id <<
          " " << sort_item_num_part_points_in_procs[ii].val <<
          " " << int(sort_item_num_part_points_in_procs[ii].signbit ) <<
          std::endl;
      }
      */

      // check if this processors is one of the procs assigned to this part.
      // if it is, then get the group.
      if(!did_i_find_my_group) {
        for(mj_part_t ii = num_procs - 1; ii >=
          num_procs - required_proc_count; --ii) {

        mj_part_t proc_id_to_assign = sort_item_num_part_points_in_procs[ii].id;

        // add the proc to the group.
        processor_ranks_for_subcomm.push_back(proc_id_to_assign);

        if(proc_id_to_assign == this->myRank) {
          // if the assigned process is me, then I find my group.
          did_i_find_my_group = true;

          // set the beginning of part i to my rank.
          part_assignment_proc_begin_indices[i] = this->myRank;
          processor_chains_in_parts[this->myRank] = -1;

          // set send count to myself to the number of points that I have
          // in part i.
          send_count_to_each_proc[this->myRank] =
            sort_item_num_part_points_in_procs[ii].val;

          // calculate the shift required for the
          // output_part_numbering_begin_index
          for(mj_part_t in = 0; in < i; ++in) {
            output_part_numbering_begin_index +=
              (*next_future_num_parts_in_parts)[in];
          }
          out_part_index = i;
        }
      }

      // if these was not my group,
      // clear the subcomminicator processor array.
      if(!did_i_find_my_group) {
        processor_ranks_for_subcomm.clear();
      }
    }

    // send points of the nonassigned coordinates to the assigned coordinates.
    // starts from the heaviest nonassigned processor.
    // TODO we might want to play with this part, that allows more
    // computational imbalance but having better communication balance.
    for(mj_part_t ii = num_procs - required_proc_count - 1; ii >= 0; --ii) {
        mj_part_t nonassigned_proc_id =
          sort_item_num_part_points_in_procs[ii].id;
        mj_lno_t num_points_to_sent =
          sort_item_num_part_points_in_procs[ii].val;

      // we set number of points to -to_sent - 1 for the assigned processors.
      // we reverse it here. This should not happen, as we have already
      // reversed them above.
#ifdef MJ_DEBUG
      if(num_points_to_sent < 0) {
        cout << "Migration - processor assignments - for part:" << i
          << "from proc:" << nonassigned_proc_id << " num_points_to_sent:"
          << num_points_to_sent << std::endl;
        std::terminate();
      }
#endif

	    switch (migration_type) {
	      case 0:
	      {
          // now sends the points to the assigned processors.
          while (num_points_to_sent > 0) {
            // if the processor has enough space.
            if(num_points_to_sent <= space_left_in_sent_proc) {
                // reduce the space left in the processor.
                space_left_in_sent_proc -= num_points_to_sent;
                // if my rank is the one that is sending the coordinates.
                if(this->myRank == nonassigned_proc_id) {
                  // set my sent count to the sent processor.
                  send_count_to_each_proc[next_proc_to_send_id] =
                    num_points_to_sent;
                  // save the processor in the list (processor_chains_in_parts
                  // and part_assignment_proc_begin_indices)
                  // that the processor will send its point in part-i.
                  mj_part_t prev_begin = part_assignment_proc_begin_indices[i];
                  part_assignment_proc_begin_indices[i] = next_proc_to_send_id;
                  processor_chains_in_parts[next_proc_to_send_id] = prev_begin;
                }
                num_points_to_sent = 0;
            }
            else {
              // there might be no space left in the processor.
              if(space_left_in_sent_proc > 0) {
                num_points_to_sent -= space_left_in_sent_proc;

                //send as the space left in the processor.
                if(this->myRank == nonassigned_proc_id) {
                  // send as much as the space in this case.
                  send_count_to_each_proc[next_proc_to_send_id] =
                    space_left_in_sent_proc;
                  mj_part_t prev_begin = part_assignment_proc_begin_indices[i];
                  part_assignment_proc_begin_indices[i] = next_proc_to_send_id;
                  processor_chains_in_parts[next_proc_to_send_id] = prev_begin;
                }
              }
              // change the sent part
              ++next_proc_to_send_index;

#ifdef MJ_DEBUG
              if(next_part_to_send_index <  nprocs - required_proc_count ) {
                  cout << "Migration - processor assignments - for part:"
                    << i
                    <<  " next_part_to_send :" << next_part_to_send_index
                    << " nprocs:" << nprocs
                    << " required_proc_count:" << required_proc_count
                    << " Error: next_part_to_send_index <" <<
                    << " nprocs - required_proc_count" << std::endl;
                std::terminate();
              }
#endif
              // send the new id.
              next_proc_to_send_id =
                sort_item_num_part_points_in_procs[next_proc_to_send_index].id;
              // set the new space in the processor.
              space_left_in_sent_proc = ideal_num_points_in_a_proc -
                sort_item_num_part_points_in_procs[next_proc_to_send_index].val;
            }
          }
	      }
	      break;
	      default:
	      {
          // to minimize messages, we want each processor to send its
          // coordinates to only a single point.
          // we do not respect imbalances here, we send all points to the
          // next processor.
		      if(this->myRank == nonassigned_proc_id) {
            // set my sent count to the sent processor.
            send_count_to_each_proc[next_proc_to_send_id] = num_points_to_sent;
            // save the processor in the list (processor_chains_in_parts and
            // part_assignment_proc_begin_indices)
            // that the processor will send its point in part-i.
            mj_part_t prev_begin = part_assignment_proc_begin_indices[i];
            part_assignment_proc_begin_indices[i] = next_proc_to_send_id;
            processor_chains_in_parts[next_proc_to_send_id] = prev_begin;
          }
          num_points_to_sent = 0;
          ++next_proc_to_send_index;

          // if we made it to the heaviest processor we round robin and
          // go to beginning
          if(next_proc_to_send_index == num_procs) {
            next_proc_to_send_index = num_procs - required_proc_count;
          }
          // send the new id.
          next_proc_to_send_id =
            sort_item_num_part_points_in_procs[next_proc_to_send_index].id;
          // set the new space in the processor.
          space_left_in_sent_proc = ideal_num_points_in_a_proc -
            sort_item_num_part_points_in_procs[next_proc_to_send_index].val;
	      }
      }
    }
  }

  /*
  for(int i = 0; i < num_procs;++i) {
    std::cout << "me:" << this->myRank << " to part:" << i << " sends:" <<
      send_count_to_each_proc[i] << std::endl;
  }
  */

  this->assign_send_destinations(
    num_parts,
    part_assignment_proc_begin_indices,
    processor_chains_in_parts,
    send_count_to_each_proc,
    coordinate_destinations);
  delete [] part_assignment_proc_begin_indices;
  delete [] processor_chains_in_parts;
  delete [] processor_part_assignments;
  delete [] sort_item_num_part_points_in_procs;
  delete [] num_procs_assigned_to_each_part;
}

/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent. In addition
 * it calculates the shift amount (output_part_numbering_begin_index) to be
 * done when final numberings of the parts are performed.
 * \param num_parts is the number of parts that exist in the
 * current partitioning.
 * \param sort_item_part_to_proc_assignment is the sorted parts with respect
 * to the assigned processors.
 * \param coordinate_destinations is the output array that holds which part
 * each coordinate should be sent.
 * \param output_part_numbering_begin_index is how much the numbers should be
 * shifted when numbering the result parts.
 * \param next_future_num_parts_in_parts is the vector, how many more parts
 * each part will be divided into in the future.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  assign_send_destinations2(
  mj_part_t num_parts,
  uSortItem<mj_part_t, mj_part_t> * sort_item_part_to_proc_assignment,
  int *coordinate_destinations,
  mj_part_t &output_part_numbering_begin_index,
  std::vector<mj_part_t> *next_future_num_parts_in_parts)
{
  mj_part_t part_shift_amount = output_part_numbering_begin_index;
  mj_part_t previous_processor = -1;

  auto local_new_part_xadj = Kokkos::create_mirror_view(this->new_part_xadj);
  Kokkos::deep_copy(local_new_part_xadj, this->new_part_xadj);

  auto local_new_coordinate_permutations =
    Kokkos::create_mirror_view(this->new_coordinate_permutations);
  Kokkos::deep_copy(local_new_coordinate_permutations,
    this->new_coordinate_permutations);

  for(mj_part_t i = 0; i < num_parts; ++i) {
    mj_part_t p = sort_item_part_to_proc_assignment[i].id;

    // assigned processors are sorted.
    mj_lno_t part_begin_index = 0;

    if(p > 0) {
      part_begin_index = local_new_part_xadj(p - 1);
    }

    mj_lno_t part_end_index = local_new_part_xadj(p);

    mj_part_t assigned_proc = sort_item_part_to_proc_assignment[i].val;
    if(this->myRank == assigned_proc && previous_processor != assigned_proc) {
        output_part_numbering_begin_index =  part_shift_amount;
    }
    previous_processor = assigned_proc;
    part_shift_amount += (*next_future_num_parts_in_parts)[p];

    for(mj_lno_t j= part_begin_index; j < part_end_index; j++) {
      mj_lno_t localInd = local_new_coordinate_permutations(j);
      coordinate_destinations[localInd] = assigned_proc;
    }
  }
}

/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent. In addition it
 * calculates the shift amount (output_part_numbering_begin_index) to be done
 * when final numberings of the parts are performed.
 * \param num_points_in_all_processor_parts is the array holding the num points
 * in each part in each proc.
 * \param num_parts is the number of parts that exist in the
 * current partitioning.
 * \param num_procs is the number of processor attending to migration operation.
 * \param send_count_to_each_proc array array storing the number of points to
 * be sent to each part.
 * \param next_future_num_parts_in_parts is the vector, how many more parts
 * each part will be divided into in the future.
 * \param out_num_part is the number of parts assigned to the process.
 * \param out_part_indices is the indices of the part to which the processor
 * is assigned.
 * \param output_part_numbering_begin_index is how much the numbers should be
 * shifted when numbering the result parts.
 * \param coordinate_destinations is the output array that holds which part
 * each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_assign_parts_to_procs(
  mj_gno_t * num_points_in_all_processor_parts,
  mj_part_t num_parts,
  mj_part_t num_procs,
  mj_lno_t *send_count_to_each_proc,
  std::vector<mj_part_t> *next_future_num_parts_in_parts,
  mj_part_t &out_num_part,
  std::vector<mj_part_t> &out_part_indices,
  mj_part_t &output_part_numbering_begin_index,
  int *coordinate_destinations) {

  out_num_part = 0;
  mj_gno_t *global_num_points_in_parts =
    num_points_in_all_processor_parts + num_procs * num_parts;
  out_part_indices.clear();

  // to sort the parts that is assigned to the processors.
  // id is the part number, sort value is the assigned processor id.
  uSortItem<mj_part_t, mj_part_t> * sort_item_part_to_proc_assignment =
    new uSortItem<mj_part_t, mj_part_t>[num_parts];
  uSortItem<mj_part_t, mj_gno_t> * sort_item_num_points_of_proc_in_part_i =
    new uSortItem<mj_part_t, mj_gno_t>[num_procs];

  // calculate the optimal number of coordinates that should be assigned
  // to each processor.
  mj_lno_t work_each =
    mj_lno_t (this->num_global_coords / (double (num_procs)) + 0.5f);

  // to hold the left space as the number of coordinates to the optimal
  // number in each proc.
  mj_lno_t *space_in_each_processor = new mj_lno_t[num_procs];

  // initialize left space in each.
  for(mj_part_t i = 0; i < num_procs; ++i) {
    space_in_each_processor[i] = work_each;
  }

  // we keep track of how many parts each processor is assigned to.
  // because in some weird inputs, it might be possible that some
  // processors is not assigned to any part. Using these variables,
  // we force each processor to have at least one part.
  mj_part_t *num_parts_proc_assigned = new mj_part_t[num_procs];
  memset(num_parts_proc_assigned, 0, sizeof(mj_part_t) * num_procs);
  int empty_proc_count = num_procs;

  // to sort the parts with decreasing order of their coordiantes.
  // id are the part numbers, sort value is the number of points in each.
  uSortItem<mj_part_t, mj_gno_t> * sort_item_point_counts_in_parts =
    new uSortItem<mj_part_t, mj_gno_t>[num_parts];

  // initially we will sort the parts according to the number of coordinates
  // they have, so that we will start assigning with the part that has the most
  // number of coordinates.
  for(mj_part_t i = 0; i < num_parts; ++i) {
    sort_item_point_counts_in_parts[i].id = i;
    sort_item_point_counts_in_parts[i].val = global_num_points_in_parts[i];
  }

  // sort parts with increasing order of loads.
  uqsort<mj_part_t, mj_gno_t>(num_parts, sort_item_point_counts_in_parts);

  // assigning parts to the processors
  // traverse the part with decreasing order of load.
  // first assign the heaviest part.
  for(mj_part_t j = 0; j < num_parts; ++j) {
    // sorted with increasing order, traverse inverse.
    mj_part_t i = sort_item_point_counts_in_parts[num_parts - 1 - j].id;

    // load of the part
    mj_gno_t load = global_num_points_in_parts[i];

    // assigned processors
    mj_part_t assigned_proc = -1;

    // sort processors with increasing number of points in this part.
    for(mj_part_t ii = 0; ii < num_procs; ++ii) {
      sort_item_num_points_of_proc_in_part_i[ii].id = ii;

      // if there are still enough parts to fill empty processors, than proceed
      // normally, but if empty processor count is equal to the number of part,
      // then we force to part assignments only to empty processors.
      if(empty_proc_count < num_parts - j ||
        num_parts_proc_assigned[ii] == 0) {
        // how many points processor ii has in part i?
        sort_item_num_points_of_proc_in_part_i[ii].val =
          num_points_in_all_processor_parts[ii * num_parts + i];
      }
      else {
        sort_item_num_points_of_proc_in_part_i[ii].val = -1;
      }
    }

    uqsort<mj_part_t, mj_gno_t>(num_procs,
      sort_item_num_points_of_proc_in_part_i);

    // traverse all processors with decreasing load.
    for(mj_part_t iii = num_procs - 1; iii >= 0; --iii) {
      mj_part_t ii = sort_item_num_points_of_proc_in_part_i[iii].id;
      if(assigned_proc == -1 ||
        (space_in_each_processor[ii] > space_in_each_processor[assigned_proc])) {
        assigned_proc = ii;
      }
      else if(space_in_each_processor[ii] == space_in_each_processor[assigned_proc]) {
        if(ii < assigned_proc) {
          // ties go to lower proc
          // not necessary for a valid result but allows testing to compare
          // MPI results and have parts numbers assigned to the same boxes.
          // We don't break here because we may have more ties still to check.
          // The indeterminate state before this is due to Cuda using
          // atomics to refill the permutation array. So non-cuda runs don't
          // actualy need this since they will always have the same pattern.
          assigned_proc = ii;
        }
      }
      else {
        break; // now we can break - we have our part and no more ties.
      }
    }

    if(num_parts_proc_assigned[assigned_proc]++ == 0) {
      --empty_proc_count;
    }

    space_in_each_processor[assigned_proc] -= load;
    //to sort later, part-i is assigned to the proccessor - assignment.
    sort_item_part_to_proc_assignment[j].id = i; //part i

    // assigned to processor - assignment.
    sort_item_part_to_proc_assignment[j].val = assigned_proc;

    // if assigned processor is me, increase the number.
    if(assigned_proc == this->myRank) {
        out_num_part++;//assigned_part_count;
        out_part_indices.push_back(i);
    }

    // increase the send to that processor by the number of points in that
    // part, as everyone send their coordiantes in this part to the
    // processor assigned to this part.
    send_count_to_each_proc[assigned_proc] +=
      num_points_in_all_processor_parts[this->myRank * num_parts + i];
  }

  delete [] num_parts_proc_assigned;
  delete [] sort_item_num_points_of_proc_in_part_i;
  delete [] sort_item_point_counts_in_parts;
  delete [] space_in_each_processor;

  // sort assignments with respect to the assigned processors.
  uqsort<mj_part_t, mj_part_t>(num_parts, sort_item_part_to_proc_assignment);

  // fill sendBuf.
  this->assign_send_destinations2(
    num_parts,
    sort_item_part_to_proc_assignment,
    coordinate_destinations,
    output_part_numbering_begin_index,
    next_future_num_parts_in_parts);

  delete [] sort_item_part_to_proc_assignment;
}


/*! \brief Function fills up coordinate_destinations is the output array
 * that holds which part each coordinate should be sent. In addition it
 * calculates the shift amount (output_part_numbering_begin_index) to be done
 * when final numberings of the parts are performed.
 * \param num_points_in_all_processor_parts is the array holding the num points
 * in each part in each proc.
 * \param num_parts is the number of parts that exist in the current
 * partitioning.
 * \param num_procs is the number of processor attending to migration operation.
 * \param send_count_to_each_proc array array storing the number of points to
 * be sent to each part.
 * \param processor_ranks_for_subcomm is the ranks of the processors that will
 * be in the subcommunicator with me.
 * \param next_future_num_parts_in_parts is the vector, how many more parts
 * each part will be divided into in the future.
 * \param out_num_part is the number of parts assigned to the process.
 * \param out_part_indices is the indices of the part to which the processor
 * is assigned.
 * \param output_part_numbering_begin_index is how much the numbers should be
 * shifted when numbering the result parts.
 * \param coordinate_destinations is the output array that holds which part
 * each coordinate should be sent.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_migration_part_proc_assignment(
  mj_gno_t * num_points_in_all_processor_parts,
  mj_part_t num_parts,
  mj_part_t num_procs,
  mj_lno_t *send_count_to_each_proc,
  std::vector<mj_part_t> &processor_ranks_for_subcomm,
  std::vector<mj_part_t> *next_future_num_parts_in_parts,
  mj_part_t &out_num_part,
  std::vector<mj_part_t> &out_part_indices,
  mj_part_t &output_part_numbering_begin_index,
  int *coordinate_destinations)
{
  processor_ranks_for_subcomm.clear();
  // if(this->num_local_coords > 0)
  if(num_procs > num_parts) {
    // if there are more processors than the number of current part
    // then processors share the existing parts.
    // at the end each processor will have a single part,
    // but a part will be shared by a group of processors.
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

    // there are more parts than the processors.
    // therefore a processor will be assigned multiple parts,
    // the subcommunicators will only have a single processor.
    processor_ranks_for_subcomm.push_back(this->myRank);

    // since there are more parts then procs,
    // assign multiple parts to processors.

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
 * that holds which part each coordinate should be sent. In addition it
 * calculates the shift amount (output_part_numbering_begin_index) to be done
 * when final numberings of the parts are performed.
 * \param num_procs is the number of processor attending to migration operation.
 * \param num_new_local_points is the output to represent the new number of
 * local points.
 * \param iteration is the string for the current iteration.
 * \param coordinate_destinations is the output array that holds which part
 * each coordinate should be sent.
 * \param num_parts is the number of parts that exist in the current
 * partitioning.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_migrate_coords(
  mj_part_t num_procs,
  mj_lno_t &num_new_local_points,
  std::string iteration,
  int *coordinate_destinations,
  mj_part_t num_parts)
{

#ifdef ZOLTAN2_MJ_ENABLE_ZOLTAN_MIGRATION
  if(sizeof(mj_lno_t) <= sizeof(int)) {
    // Cannot use Zoltan_Comm with local ordinals larger than ints.
    // In Zoltan_Comm_Create, the cast int(this->num_local_coords)
    // may overflow.
    ZOLTAN_COMM_OBJ *plan = NULL;
    MPI_Comm mpi_comm = Teuchos::getRawMpiComm(*(this->comm));
    int num_incoming_gnos = 0;
    int message_tag = 7859;

    this->mj_env->timerStart(MACRO_TIMERS,
      mj_timer_base_string + "Migration Z1PlanCreating-" + iteration);
    int ierr = Zoltan_Comm_Create(
      &plan,
      int(this->num_local_coords),
      coordinate_destinations,
      mpi_comm,
      message_tag,
      &num_incoming_gnos);

    Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
    this->mj_env->timerStop(MACRO_TIMERS,
      mj_timer_base_string + "Migration Z1PlanCreating-" + iteration);

    this->mj_env->timerStart(MACRO_TIMERS,
      mj_timer_base_string + "Migration Z1Migration-" + iteration);

    // MPI Buffers should be on Kokkos::HostSpace not Kokkos::CudaUVMSpace

    // migrate gnos.
    {
      auto host_current_mj_gnos = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), this->current_mj_gnos);
      Kokkos::deep_copy(host_current_mj_gnos, this->current_mj_gnos);
      Kokkos::View<mj_gno_t*, device_t> dst_gnos(
        Kokkos::ViewAllocateWithoutInitializing("dst_gnos"), num_incoming_gnos);
      auto host_dst_gnos = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), dst_gnos);
      message_tag++;
      ierr = Zoltan_Comm_Do(
        plan,
        message_tag,
        (char *) host_current_mj_gnos.data(),
        sizeof(mj_gno_t),
        (char *) host_dst_gnos.data());
      Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
      Kokkos::deep_copy(dst_gnos, host_dst_gnos);
      this->current_mj_gnos = dst_gnos;
    }

    //migrate coordinates
    {
      // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
      auto host_src_coordinates = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), this->mj_coordinates);
      Kokkos::deep_copy(host_src_coordinates, this->mj_coordinates);
      Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t>
        dst_coordinates(Kokkos::ViewAllocateWithoutInitializing("mj_coordinates"),
        num_incoming_gnos, this->coord_dim);
      auto host_dst_coordinates = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), dst_coordinates);
      for(int i = 0; i < this->coord_dim; ++i) {
        Kokkos::View<mj_scalar_t*, Kokkos::HostSpace> sub_host_src_coordinates;
        // view could be size 0 if graph was not distributed
        if(host_src_coordinates.extent(0) != 0) {
          sub_host_src_coordinates =
            Kokkos::subview(host_src_coordinates, Kokkos::ALL, i);
        }
        Kokkos::View<mj_scalar_t *, Kokkos::HostSpace> sub_host_dst_coordinates
          = Kokkos::subview(host_dst_coordinates, Kokkos::ALL, i);
        // Note Layout Left means we can do these in contiguous blocks
        message_tag++;
        ierr = Zoltan_Comm_Do(
          plan,
          message_tag,
          (char *) sub_host_src_coordinates.data(),
          sizeof(mj_scalar_t),
          (char *) sub_host_dst_coordinates.data());
        Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
      }
      deep_copy(dst_coordinates, host_dst_coordinates);
      this->mj_coordinates = dst_coordinates;
    }

    // migrate weights.
    {
      auto host_src_weights = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), this->mj_weights);
      Kokkos::deep_copy(host_src_weights, this->mj_weights);
      Kokkos::View<mj_scalar_t**, device_t> dst_weights(
        Kokkos::ViewAllocateWithoutInitializing("mj_weights"),
        num_incoming_gnos, this->num_weights_per_coord);
      auto host_dst_weights = Kokkos::create_mirror_view(dst_weights);
      for(int i = 0; i < this->num_weights_per_coord; ++i) {
        auto sub_host_src_weights
          = Kokkos::subview(host_src_weights, Kokkos::ALL, i);
        auto sub_host_dst_weights
          = Kokkos::subview(host_dst_weights, Kokkos::ALL, i);
        ArrayRCP<mj_scalar_t> sent_weight(this->num_local_coords);
        // Copy because of layout
        for(mj_lno_t n = 0; n < this->num_local_coords; ++n) {
          sent_weight[n] = sub_host_src_weights(n);
        }
        ArrayRCP<mj_scalar_t> received_weight(num_incoming_gnos);
        message_tag++;
        ierr = Zoltan_Comm_Do(
          plan,
          message_tag,
          (char *) sent_weight.getRawPtr(),
          sizeof(mj_scalar_t),
          (char *) received_weight.getRawPtr());
        Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
        // Again we copy by index due to layout
        for(mj_lno_t n = 0; n < num_incoming_gnos; ++n) {
          sub_host_dst_weights(n) = received_weight[n];
        }
      }
      deep_copy(dst_weights, host_dst_weights);
      this->mj_weights = dst_weights;
    }

    // migrate owners.
    {
      // Note that owners we kept on Serial
      Kokkos::View<int *, Kokkos::HostSpace> dst_owners_of_coordinate(
        Kokkos::ViewAllocateWithoutInitializing("owner_of_coordinate"),
        num_incoming_gnos);
      message_tag++;
      ierr = Zoltan_Comm_Do(
        plan,
        message_tag,
        (char *) owner_of_coordinate.data(),
        sizeof(int),
        (char *) dst_owners_of_coordinate.data());
      Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
      this->owner_of_coordinate = dst_owners_of_coordinate;
    }

    // if num procs is less than num parts,
    // we need the part assigment arrays as well, since
    // there will be multiple parts in processor.
    {
      auto host_src_assigned_part_ids = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), this->assigned_part_ids);
      Kokkos::deep_copy(host_src_assigned_part_ids, this->assigned_part_ids);
      Kokkos::View<int *, device_t> dst_assigned_part_ids(
        Kokkos::ViewAllocateWithoutInitializing("assigned_part_ids"),
        num_incoming_gnos);
      auto host_dst_assigned_part_ids = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), dst_assigned_part_ids);
      mj_part_t *new_parts = new mj_part_t[num_incoming_gnos];
      if(num_procs < num_parts) {
        message_tag++;
        ierr = Zoltan_Comm_Do(
          plan,
          message_tag,
          (char *) host_src_assigned_part_ids.data(),
          sizeof(mj_part_t),
          (char *) host_dst_assigned_part_ids.data());
        Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
        Kokkos::deep_copy(dst_assigned_part_ids, host_dst_assigned_part_ids);
      }
      // In original code this would just assign to an uninitialized array
      // if num_procs < num_parts. We're doing the same here.
      this->assigned_part_ids = dst_assigned_part_ids;
    }

    ierr = Zoltan_Comm_Destroy(&plan);
    Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
    num_new_local_points = num_incoming_gnos;
    this->mj_env->timerStop(MACRO_TIMERS,
      mj_timer_base_string + "Migration Z1Migration-" + iteration);
  }
  else
#endif  // ZOLTAN2_MJ_ENABLE_ZOLTAN_MIGRATION
  {
    this->mj_env->timerStart(MACRO_TIMERS, mj_timer_base_string +
      "Migration DistributorPlanCreating-" + iteration);

    Tpetra::Distributor distributor(this->comm);
    ArrayView<const mj_part_t> destinations( coordinate_destinations,
      this->num_local_coords);
    mj_lno_t num_incoming_gnos = distributor.createFromSends(destinations);
    this->mj_env->timerStop(MACRO_TIMERS, mj_timer_base_string +
      "Migration DistributorPlanCreating-" + iteration);
    this->mj_env->timerStart(MACRO_TIMERS, mj_timer_base_string +
      "Migration DistributorMigration-" + iteration);

    // note MPI buffers should all be on Kokkos::HostSpace and not
    // Kokkos::CudaUVMSpace.

    // migrate gnos.
    {
      ArrayRCP<mj_gno_t> received_gnos(num_incoming_gnos);
      auto src_host_current_mj_gnos =
        Kokkos::create_mirror_view(Kokkos::HostSpace(), this->current_mj_gnos);
      Kokkos::deep_copy(src_host_current_mj_gnos, this->current_mj_gnos);
      ArrayView<mj_gno_t> sent_gnos(
        src_host_current_mj_gnos.data(), this->num_local_coords);
      distributor.doPostsAndWaits<mj_gno_t>(sent_gnos, 1, received_gnos());
      this->current_mj_gnos = Kokkos::View<mj_gno_t*, device_t>(
        Kokkos::ViewAllocateWithoutInitializing("gids"), num_incoming_gnos);
      auto host_current_mj_gnos = Kokkos::create_mirror_view(
        this->current_mj_gnos);
      memcpy(host_current_mj_gnos.data(),
        received_gnos.getRawPtr(), num_incoming_gnos * sizeof(mj_gno_t));
      Kokkos::deep_copy(this->current_mj_gnos, host_current_mj_gnos);
    }

    // migrate coordinates
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t>
      dst_coordinates("mj_coordinates", num_incoming_gnos, this->coord_dim);
    auto host_dst_coordinates = Kokkos::create_mirror_view(dst_coordinates);
    auto host_src_coordinates = Kokkos::create_mirror_view(
      Kokkos::HostSpace(), this->mj_coordinates);
    Kokkos::deep_copy(host_src_coordinates, this->mj_coordinates);
    for(int i = 0; i < this->coord_dim; ++i) {
      Kokkos::View<mj_scalar_t*, Kokkos::HostSpace> sub_host_src_coordinates;
      // view could be size 0 if graph was not distributed
      if(host_src_coordinates.extent(0) != 0) {
        sub_host_src_coordinates =
          Kokkos::subview(host_src_coordinates, Kokkos::ALL, i);
      }

      auto sub_host_dst_coordinates
        = Kokkos::subview(host_dst_coordinates, Kokkos::ALL, i);
      // Note Layout Left means we can do these in contiguous blocks
      // This form was causing problems on cuda 10 pascal nodes, issue #6422
      // Doing a manual copy clears the error though it seems this is probably
      // just shifting some kind of race condition or UVM issue around. The
      // bug can be sensitive to simple changes like adding a printf log.

      // Using this form will segfault on cuda 10 pascal node
      //ArrayView<mj_scalar_t> sent_coord(
      //  sub_host_src_coordinates.data(), this->num_local_coords);

      // Manual copy will clear the error but this is probably just due to
      // shifting some kind of race condition.
      ArrayRCP<mj_scalar_t> sent_coord(this->num_local_coords);
      for(int n = 0; n < this->num_local_coords; ++n) {
        sent_coord[n] = sub_host_src_coordinates[n];
      }

      ArrayRCP<mj_scalar_t> received_coord(num_incoming_gnos);
      distributor.doPostsAndWaits<mj_scalar_t>(
        sent_coord(), 1, received_coord());
      memcpy(sub_host_dst_coordinates.data(),
        received_coord.getRawPtr(), num_incoming_gnos * sizeof(mj_scalar_t));
    }
    deep_copy(dst_coordinates, host_dst_coordinates);
    this->mj_coordinates = dst_coordinates;

    // migrate weights.
    Kokkos::View<mj_scalar_t**, device_t> dst_weights(
     "mj_weights", num_incoming_gnos, this->num_weights_per_coord);
    auto host_dst_weights = Kokkos::create_mirror_view(dst_weights);
    auto host_src_weights = Kokkos::create_mirror_view(
      Kokkos::HostSpace(), this->mj_weights);
    Kokkos::deep_copy(host_src_weights, this->mj_weights);
    for(int i = 0; i < this->num_weights_per_coord; ++i) {
      auto sub_host_src_weights
        = Kokkos::subview(host_src_weights, Kokkos::ALL, i);
      auto sub_host_dst_weights
        = Kokkos::subview(host_dst_weights, Kokkos::ALL, i);
      ArrayRCP<mj_scalar_t> sent_weight(this->num_local_coords);

      // TODO: Layout Right means these are not contiguous
      // However we don't have any systems setup with more than 1 weight so
      // really I have not tested any of this code with num weights > 1.
      // I think this is the right thing to do.
      for(mj_lno_t n = 0; n < this->num_local_coords; ++n) {
        sent_weight[n] = sub_host_src_weights(n);
      }
      ArrayRCP<mj_scalar_t> received_weight(num_incoming_gnos);
      distributor.doPostsAndWaits<mj_scalar_t>(
        sent_weight(), 1, received_weight());

      // Again we copy by index due to layout
      for(mj_lno_t n = 0; n < num_incoming_gnos; ++n) {
        sub_host_dst_weights(n) = received_weight[n];
      }
    }
    Kokkos::deep_copy(dst_weights, host_dst_weights);
    this->mj_weights = dst_weights;

    // migrate owners
    {
      // Note owners we kept on Serial
      ArrayView<int> sent_owners(
        owner_of_coordinate.data(), this->num_local_coords);
      ArrayRCP<int> received_owners(num_incoming_gnos);
      distributor.doPostsAndWaits<int>(sent_owners, 1, received_owners());
      this->owner_of_coordinate = Kokkos::View<int *, Kokkos::HostSpace>
        ("owner_of_coordinate", num_incoming_gnos);
      memcpy(this->owner_of_coordinate.data(),
        received_owners.getRawPtr(), num_incoming_gnos * sizeof(int));
    }

    // if num procs is less than num parts,
    // we need the part assigment arrays as well, since
    // there will be multiple parts in processor.
    if(num_procs < num_parts) {
      auto src_host_assigned_part_ids =
        Kokkos::create_mirror_view(Kokkos::HostSpace(), this->assigned_part_ids);
      Kokkos::deep_copy(src_host_assigned_part_ids, assigned_part_ids);
      ArrayView<mj_part_t> sent_partids(
        src_host_assigned_part_ids.data(), this->num_local_coords);
      ArrayRCP<mj_part_t> received_partids(num_incoming_gnos);
      distributor.doPostsAndWaits<mj_part_t>(
        sent_partids, 1, received_partids());
      this->assigned_part_ids = Kokkos::View<mj_part_t *, device_t>
        ("assigned_part_ids", num_incoming_gnos);
      auto host_assigned_part_ids = Kokkos::create_mirror_view(
        this->assigned_part_ids);
      memcpy(
        host_assigned_part_ids.data(),
        received_partids.getRawPtr(),
        num_incoming_gnos * sizeof(mj_part_t));
      Kokkos::deep_copy(this->assigned_part_ids, host_assigned_part_ids);
    }
    else {
      this->assigned_part_ids = Kokkos::View<mj_part_t *, device_t>
        ("assigned_part_ids", num_incoming_gnos);
    }
    this->mj_env->timerStop(MACRO_TIMERS, "" + mj_timer_base_string +
      "Migration DistributorMigration-" + iteration);

    num_new_local_points = num_incoming_gnos;
  }
}

/*! \brief Function creates the new subcomminicator for the processors
 * given in processor_ranks_for_subcomm.
 * \param processor_ranks_for_subcomm is the vector that has the ranks of
 * the processors that will be in the same group.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  create_sub_communicator(std::vector<mj_part_t> &processor_ranks_for_subcomm)
{
  mj_part_t group_size = processor_ranks_for_subcomm.size();
  mj_part_t *ids = new mj_part_t[group_size];
  for(mj_part_t i = 0; i < group_size; ++i) {
    ids[i] = processor_ranks_for_subcomm[i];
  }
  ArrayView<const mj_part_t> idView(ids, group_size);
  this->comm = this->comm->createSubcommunicator(idView);
  delete [] ids;
}

/*! \brief Function writes the new permutation arrays after the migration.
 * \param output_num_parts is the number of parts that is assigned to
 * the processor.
 * \param num_parts is the number of parts right before migration.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  fill_permutation_array(
  mj_part_t output_num_parts,
  mj_part_t num_parts)
{
  // if there is single output part, then simply fill the permutation array.
  if(output_num_parts == 1) {
    auto local_new_coordinate_permutations = this->new_coordinate_permutations;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_lno_t>
        (0, this->num_local_coords),
      KOKKOS_LAMBDA(mj_lno_t i) {
      local_new_coordinate_permutations(i) = i;
    });
    auto local_new_part_xadj = this->new_part_xadj;
    auto local_num_local_coords = this->num_local_coords;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0,1),
      KOKKOS_LAMBDA(int dummy) {
      local_new_part_xadj(0) = local_num_local_coords;
    });
  }
  else {
    auto local_num_local_coords = this->num_local_coords;
    auto local_assigned_part_ids = this->assigned_part_ids;
    auto local_new_part_xadj = this->new_part_xadj;
    auto local_new_coordinate_permutations = this->new_coordinate_permutations;

    // part shift holds the which part number an old part number corresponds to.
    Kokkos::View<mj_part_t*, device_t> part_shifts("part_shifts", num_parts);

    // otherwise we need to count how many points are there in each part.
    // we allocate here as num_parts, because the sent partids are up to
    // num_parts, although there are outout_num_parts different part.
    Kokkos::View<mj_lno_t*, device_t> num_points_in_parts(
      "num_points_in_parts", num_parts);

    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0,1),
      KOKKOS_LAMBDA(int dummy) {

      for(mj_lno_t i = 0; i < local_num_local_coords; ++i) {
        mj_part_t ii = local_assigned_part_ids(i);
        ++num_points_in_parts(ii);
      }

      // write the end points of the parts.
      mj_part_t p = 0;
      mj_lno_t prev_index = 0;
      for(mj_part_t i = 0; i < num_parts; ++i) {
        if(num_points_in_parts(i) > 0) {
          local_new_part_xadj(p) = prev_index + num_points_in_parts(i);
          prev_index += num_points_in_parts(i);
          part_shifts(i) = p++;
        }
      }

      // for the rest of the parts write the end index as end point.
      mj_part_t assigned_num_parts = p - 1;
      for(;p < num_parts; ++p) {
        local_new_part_xadj(p) =
          local_new_part_xadj(assigned_num_parts);
      }
      for(mj_part_t i = 0; i < output_num_parts; ++i) {
        num_points_in_parts(i) = local_new_part_xadj(i);
      }

      // write the permutation array here.
      // get the part of the coordinate i, shift it to obtain the new part number.
      // assign it to the end of the new part numbers pointer.
      for(mj_lno_t i = local_num_local_coords - 1; i >= 0; --i) {
        mj_part_t part =
          part_shifts[mj_part_t(local_assigned_part_ids(i))];
        local_new_coordinate_permutations(--num_points_in_parts[part]) = i;
      }
    });
  }
}

/*! \brief Function checks if should do migration or not.
 * It returns true to point that migration should be done when
 * -migration_reduce_all_population are higher than a predetermined value
 * -num_coords_for_last_dim_part that left for the last dimension partitioning
 * is less than a predetermined value - the imbalance of the processors on the
 * parts are higher than given threshold.
 * \param input_num_parts is the number of parts when migration is called.
 * \param output_num_parts is the output number of parts after migration.
 * \param next_future_num_parts_in_parts is the number of total future parts
 * each part is partitioned into. This will be updated for migration.
 * \param output_part_begin_index is the number that will be used as beginning
 * part number when final solution part numbers are assigned.
 * \param migration_reduce_all_population is the estimated total number of
 * reduceall operations multiplied with number of processors to be used for
 * determining migration.
 * \param num_coords_for_last_dim_part is the estimated number of points in each
 * part, when last dimension partitioning is performed.
 * \param iteration is the string that gives information about the dimension
 * for printing purposes.
 * \param input_part_boxes is the array that holds the part boxes after
 * the migration. (swapped)
 * \param output_part_boxes is the array that holds the part boxes before
 * the migration. (swapped)
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
bool AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  mj_perform_migration(
  mj_part_t input_num_parts,
  mj_part_t &output_num_parts,
  std::vector<mj_part_t> *next_future_num_parts_in_parts,
  mj_part_t &output_part_begin_index,
  size_t migration_reduce_all_population,
  mj_lno_t num_coords_for_last_dim_part,
  std::string iteration,
  RCP<mj_partBoxVector_t> &input_part_boxes,
  RCP<mj_partBoxVector_t> &output_part_boxes)
{
  mj_part_t num_procs = this->comm->getSize();
  this->myRank = this->comm->getRank();

  // this array holds how many points each processor has in each part.
  // to access how many points processor i has on part j,
  // num_points_in_all_processor_parts[i * num_parts + j]
  mj_gno_t *num_points_in_all_processor_parts =
    new mj_gno_t[input_num_parts * (num_procs + 1)];

  // get the number of coordinates in each part in each processor.
  this->get_processor_num_points_in_parts(
    num_procs,
    input_num_parts,
    num_points_in_all_processor_parts);

  // check if migration will be performed or not.
  if(!this->mj_check_to_migrate(
    migration_reduce_all_population,
    num_coords_for_last_dim_part,
    num_procs,
    input_num_parts,
    num_points_in_all_processor_parts)) {
    delete [] num_points_in_all_processor_parts;
    return false;
  }

  mj_lno_t *send_count_to_each_proc = NULL;
  int *coordinate_destinations = new int[this->num_local_coords];
  send_count_to_each_proc = new mj_lno_t[num_procs];

  for(int i = 0; i < num_procs; ++i) {
    send_count_to_each_proc[i] = 0;
  }

  std::vector<mj_part_t> processor_ranks_for_subcomm;
  std::vector<mj_part_t> out_part_indices;

  // determine which processors are assigned to which parts
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

  delete [] send_count_to_each_proc;
  std::vector <mj_part_t> tmpv;

  std::sort (out_part_indices.begin(), out_part_indices.end());
  mj_part_t outP = out_part_indices.size();
  mj_gno_t new_global_num_points = 0;
  mj_gno_t *global_num_points_in_parts =
    num_points_in_all_processor_parts + num_procs * input_num_parts;

  if(this->mj_keep_part_boxes) {
    input_part_boxes->clear();
  }

  // now we calculate the new values for next_future_num_parts_in_parts.
  // same for the part boxes.
  for(mj_part_t i = 0; i < outP; ++i) {
    mj_part_t ind = out_part_indices[i];
    new_global_num_points += global_num_points_in_parts[ind];
    tmpv.push_back((*next_future_num_parts_in_parts)[ind]);
    if(this->mj_keep_part_boxes) {
      input_part_boxes->push_back((*output_part_boxes)[ind]);
    }
  }

  // swap the input and output part boxes.
  if(this->mj_keep_part_boxes) {
    RCP<mj_partBoxVector_t> tmpPartBoxes = input_part_boxes;
    input_part_boxes = output_part_boxes;
    output_part_boxes = tmpPartBoxes;
  }
  next_future_num_parts_in_parts->clear();
  for(mj_part_t i = 0; i < outP; ++i) {
    mj_part_t p = tmpv[i];
    next_future_num_parts_in_parts->push_back(p);
  }

  delete [] num_points_in_all_processor_parts;

  mj_lno_t num_new_local_points = 0;
  //perform the actual migration operation here.
  this->mj_migrate_coords(
    num_procs,
    num_new_local_points,
    iteration,
    coordinate_destinations,
    input_num_parts);

  delete [] coordinate_destinations;
  if(this->num_local_coords != num_new_local_points) {
    this->new_coordinate_permutations = Kokkos::View<mj_lno_t*, device_t>
      (Kokkos::ViewAllocateWithoutInitializing("new_coordinate_permutations"),
      num_new_local_points);
    this->coordinate_permutations = Kokkos::View<mj_lno_t*, device_t>
      (Kokkos::ViewAllocateWithoutInitializing("coordinate_permutations"),
      num_new_local_points);
  }
  this->num_local_coords = num_new_local_points;
  this->num_global_coords = new_global_num_points;

  // create subcommunicator.
  this->create_sub_communicator(processor_ranks_for_subcomm);

  processor_ranks_for_subcomm.clear();

  // fill the new permutation arrays.
  this->fill_permutation_array(output_num_parts, input_num_parts);

  return true;
}

/*! \brief Function creates consistent chunks for task partitioning. Used only
 * in the case of sequential task partitioning, where consistent handle of the
 * points on the cuts are required.
 * \param num_parts is the number of parts.
 * \param mj_current_dim_coords is 1 dimensional array holding the
 * coordinate values.
 * \param current_concurrent_cut_coordinate is 1 dimensional array holding
 * the cut coordinates.
 * \param coordinate_begin is the start index of the given partition on
 * partitionedPointPermutations.
 * \param coordinate_end is the end index of the given partition on
 * partitionedPointPermutations.
 * \param used_local_cut_line_weight_to_left holds how much weight of the
 * coordinates on the cutline should be put on left side.
 * \param out_part_xadj is the indices of begginning and end of the parts in
 * the output partition.
 * \param coordInd is the index according to which the partitioning is done.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  create_consistent_chunks(
  mj_part_t num_parts,
  Kokkos::View<mj_scalar_t *, device_t> & mj_current_dim_coords,
  Kokkos::View<mj_scalar_t *, device_t> & current_concurrent_cut_coordinate,
  mj_lno_t coordinate_begin,
  mj_lno_t coordinate_end,
  Kokkos::View<mj_scalar_t *, device_t> & used_local_cut_line_weight_to_left,
  Kokkos::View<mj_lno_t *, device_t> & out_part_xadj,
  int coordInd,
  bool longest_dim_part,
  uSignedSortItem<int, mj_scalar_t, char> * p_coord_dimension_range_sorted)
{
  // Note that this method is only used by task mapper
  // All code in this file has been verified to run with UVM off by running
  // mj tests and task mapper tests with UVM off. However for this particular
  // method I did not do much for UVM off. I heavily use device to host copies
  // and more or less preserve the original logic. Due to the handling of
  // arrays it will be a bit of work to convert this to as better form.
  // Since it's only relevant to task mapper and I wasn't sure how much priority
  // to give it, I put that on hold until further discussion.
  mj_part_t no_cuts = num_parts - 1;

  // now if the rectilinear partitioning is allowed we decide how
  // much weight each thread should put to left and right.
  if(this->distribute_points_on_cut_lines) {
    auto local_thread_cut_line_weight_to_put_left =
      this->thread_cut_line_weight_to_put_left;
    auto local_thread_part_weight_work =
      this->thread_part_weight_work;
    auto local_sEpsilon = this->sEpsilon;

    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space,
      mj_part_t> (0, no_cuts), KOKKOS_LAMBDA (mj_part_t i) {
      // the left to be put on the left of the cut.
      mj_scalar_t left_weight = used_local_cut_line_weight_to_left(i);
      if(left_weight > local_sEpsilon) {
        // the weight of thread ii on cut.
        mj_scalar_t thread_ii_weight_on_cut =
          local_thread_part_weight_work(i * 2 + 1) -
          local_thread_part_weight_work(i * 2);
        if(thread_ii_weight_on_cut < left_weight) {
          local_thread_cut_line_weight_to_put_left(i) =
            thread_ii_weight_on_cut;
        }
        else {
          local_thread_cut_line_weight_to_put_left(i) = left_weight;
        }
      }
      else {
        local_thread_cut_line_weight_to_put_left(i) = 0;
      }
    });

    if(no_cuts > 0) {
      auto local_least_signifiance = least_signifiance;
      auto local_significance_mul = significance_mul;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
        (0, 1), KOKKOS_LAMBDA (int dummy) {
        // this is a special case. If cutlines share the same coordinate,
        // their weights are equal.
        // we need to adjust the ratio for that.
        for(mj_part_t i = no_cuts - 1; i > 0 ; --i) {
          mj_scalar_t cut1 = current_concurrent_cut_coordinate(i-1);
          mj_scalar_t cut2 = current_concurrent_cut_coordinate(i);
          mj_scalar_t delta = cut2 - cut1;
          mj_scalar_t abs_delta = (delta > 0) ? delta : -delta;
          if(abs_delta < local_sEpsilon) {
            local_thread_cut_line_weight_to_put_left(i) -=
              local_thread_cut_line_weight_to_put_left(i - 1);
          }
          local_thread_cut_line_weight_to_put_left(i) =
            static_cast<long long>((local_thread_cut_line_weight_to_put_left(i) +
            local_least_signifiance) * local_significance_mul) /
            static_cast<mj_scalar_t>(local_significance_mul);
        }
      });
    }
  }

  auto local_thread_point_counts = this->thread_point_counts;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
    (0, num_parts), KOKKOS_LAMBDA (mj_part_t i) {
    local_thread_point_counts(i) = 0;
  });

  // for this specific case we dont want to distribute the points along the
  // cut position randomly, as we need a specific ordering of them. Instead,
  // we put the coordinates into a sort item, where we sort those
  // using the coordinates of points on other dimensions and the index.

  // some of the cuts might share the same position.
  // in this case, if cut i and cut j share the same position
  // cut_map[i] = cut_map[j] = sort item index.
  mj_part_t *cut_map = new mj_part_t[no_cuts];

  typedef uMultiSortItem<mj_lno_t, int, mj_scalar_t> multiSItem;
  typedef std::vector< multiSItem > multiSVector;
  typedef std::vector<multiSVector> multiS2Vector;

  // to keep track of the memory allocated.
  std::vector<mj_scalar_t *>allocated_memory;

  // vector for which the coordinates will be sorted.
  multiS2Vector sort_vector_points_on_cut;

  // the number of cuts that have different coordinates.
  mj_part_t different_cut_count = 1;
  cut_map[0] = 0;

  // now we insert 1 sort vector for all cuts on the different
  // positins.if multiple cuts are on the same position,
  // they share sort vectors.
  multiSVector tmpMultiSVector;
  sort_vector_points_on_cut.push_back(tmpMultiSVector);

  auto local_current_concurrent_cut_coordinate =
    current_concurrent_cut_coordinate;
  auto host_current_concurrent_cut_coordinate =
    Kokkos::create_mirror_view(local_current_concurrent_cut_coordinate);
  Kokkos::deep_copy(host_current_concurrent_cut_coordinate,
    local_current_concurrent_cut_coordinate);

  for(mj_part_t i = 1; i < no_cuts ; ++i) {
    // if cuts share the same cut coordinates
    // set the cutmap accordingly.
    if(std::abs(host_current_concurrent_cut_coordinate(i) -
      host_current_concurrent_cut_coordinate(i-1)) < this->sEpsilon) {
      cut_map[i] = cut_map[i-1];
    }
    else {
      cut_map[i] = different_cut_count++;
      multiSVector tmp2MultiSVector;
      sort_vector_points_on_cut.push_back(tmp2MultiSVector);
    }
  }
  Kokkos::deep_copy(current_concurrent_cut_coordinate,
    host_current_concurrent_cut_coordinate);

  // now the actual part assigment.
  auto host_coordinate_permutations =
    Kokkos::create_mirror_view(coordinate_permutations);
  Kokkos::deep_copy(host_coordinate_permutations, coordinate_permutations);

  auto host_assigned_part_ids = Kokkos::create_mirror_view(assigned_part_ids);
  Kokkos::deep_copy(host_assigned_part_ids, assigned_part_ids);

  auto host_mj_coordinates = Kokkos::create_mirror_view(mj_coordinates);
  Kokkos::deep_copy(host_mj_coordinates, mj_coordinates);

  auto host_thread_point_counts = Kokkos::create_mirror_view(thread_point_counts);
  Kokkos::deep_copy(host_thread_point_counts, thread_point_counts);

  auto local_coord_dim = this->coord_dim;

  for(mj_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii) {
    mj_lno_t i = host_coordinate_permutations(ii);
    mj_part_t pp = host_assigned_part_ids(i);
    mj_part_t p = pp / 2;
    // if the coordinate is on a cut.
    if(pp % 2 == 1 ) {
      mj_scalar_t *vals = new mj_scalar_t[local_coord_dim -1];
      allocated_memory.push_back(vals);

      // we insert the coordinates to the sort item here.
      int val_ind = 0;

      if(longest_dim_part) {
        // std::cout << std::endl << std::endl;
        for(int dim = local_coord_dim - 2; dim >= 0; --dim) {
          // uSignedSortItem<int, mj_scalar_t, char>
          //   *p_coord_dimension_range_sorted
          int next_largest_coord_dim = p_coord_dimension_range_sorted[dim].id;
          // std::cout << "next_largest_coord_dim: " <<
          //   next_largest_coord_dim << " ";
          // Note refactor in progress
          vals[val_ind++] =
            host_mj_coordinates(i,next_largest_coord_dim);
        }
      }
      else {
        for(int dim = coordInd + 1; dim < local_coord_dim; ++dim) {
          vals[val_ind++] = host_mj_coordinates(i,dim);
        }
        for(int dim = 0; dim < coordInd; ++dim) {
          vals[val_ind++] = host_mj_coordinates(i,dim);
        }
      }

      multiSItem tempSortItem(i, local_coord_dim -1, vals);
      //insert the point to the sort vector pointed by the cut_map[p].
      mj_part_t cmap = cut_map[p];
      sort_vector_points_on_cut[cmap].push_back(tempSortItem);
    }
    else {
      //if it is not on the cut, simple sorting.
      ++host_thread_point_counts(p);
      host_assigned_part_ids(i) = p;
    }
  }

  // sort all the sort vectors.
  for(mj_part_t i = 0; i < different_cut_count; ++i) {
    std::sort (sort_vector_points_on_cut[i].begin(),
      sort_vector_points_on_cut[i].end());
  }

  mj_part_t previous_cut_map = cut_map[0];

  auto host_thread_cut_line_weight_to_put_left =
    Kokkos::create_mirror_view(thread_cut_line_weight_to_put_left);
  Kokkos::deep_copy(host_thread_cut_line_weight_to_put_left,
    thread_cut_line_weight_to_put_left);

  auto host_mj_weights = Kokkos::create_mirror_view(mj_weights);
  Kokkos::deep_copy(host_mj_weights, mj_weights);

  // this is how much previous part owns the weight of the current part.
  // when target part weight is 1.6, and the part on the left is given 2,
  // the left has an extra 0.4, while the right has missing 0.4 from the
  // previous cut.
  // This parameter is used to balance this issues.
  // in the above example weight_stolen_from_previous_part will be 0.4.
  // if the left part target is 2.2 but it is given 2,
  // then weight_stolen_from_previous_part will be -0.2.
  mj_scalar_t weight_stolen_from_previous_part = 0;
  for(mj_part_t p = 0; p < no_cuts; ++p) {
    mj_part_t mapped_cut = cut_map[p];

    // if previous cut map is done, and it does not have the same index,
    // then assign all points left on that cut to its right.
    if(previous_cut_map != mapped_cut) {
      mj_lno_t sort_vector_end = (mj_lno_t)
        sort_vector_points_on_cut[previous_cut_map].size() - 1;
      for(; sort_vector_end >= 0; --sort_vector_end) {
        multiSItem t =
          sort_vector_points_on_cut[previous_cut_map][sort_vector_end];
        mj_lno_t i = t.index;
        ++host_thread_point_counts(p);
        host_assigned_part_ids(i) = p;
      }
      sort_vector_points_on_cut[previous_cut_map].clear();
    }

    // TODO: MD: I dont remember why I have it reverse order here.
    mj_lno_t sort_vector_end = (mj_lno_t)
      sort_vector_points_on_cut[mapped_cut].size() - 1;
    // mj_lno_t sort_vector_begin= 0;
    // mj_lno_t sort_vector_size =
    //   (mj_lno_t)sort_vector_points_on_cut[mapped_cut].size();

    // TODO commented for reverse order
    for(; sort_vector_end >= 0; --sort_vector_end) {
      // for(; sort_vector_begin < sort_vector_size; ++sort_vector_begin) {
      // TODO COMMENTED FOR REVERSE ORDER
      multiSItem t = sort_vector_points_on_cut[mapped_cut][sort_vector_end];
      //multiSItem t = sort_vector_points_on_cut[mapped_cut][sort_vector_begin];
      mj_lno_t i = t.index;
      mj_scalar_t w = this->mj_uniform_weights(0) ? 1 :
        this->mj_weights(i,0);
      // part p has enough space for point i, then put it to point i.
      if(host_thread_cut_line_weight_to_put_left(p) +
        weight_stolen_from_previous_part> this->sEpsilon &&
        host_thread_cut_line_weight_to_put_left(p) +
        weight_stolen_from_previous_part -
        std::abs(host_thread_cut_line_weight_to_put_left(p) +
        weight_stolen_from_previous_part - w)> this->sEpsilon)
      {
        host_thread_cut_line_weight_to_put_left(p) -= w;

        sort_vector_points_on_cut[mapped_cut].pop_back();

        ++host_thread_point_counts(p);
        host_assigned_part_ids(i) = p;
        // if putting this weight to left overweights the left cut, then
        // increase the space for the next cut using
        // weight_stolen_from_previous_part.
        if(p < no_cuts - 1 &&
          host_thread_cut_line_weight_to_put_left(p) < this->sEpsilon) {
            if(mapped_cut == cut_map[p + 1] ) {
              // if the cut before the cut indexed at p was also at the same
              // position special case, as we handle the weight differently here.
              if(previous_cut_map != mapped_cut) {
                weight_stolen_from_previous_part =
                  host_thread_cut_line_weight_to_put_left(p);
              }
              else {
                // if the cut before the cut indexed at p was also at the same
                // position we assign extra weights cumulatively in this case.
                weight_stolen_from_previous_part +=
                  host_thread_cut_line_weight_to_put_left(p);
              }
            }
            else{
              weight_stolen_from_previous_part =
                -host_thread_cut_line_weight_to_put_left(p);
            }
            // end assignment for part p
            break;
        }
      } else {
        // if part p does not have enough space for this point
        // and if there is another cut sharing the same positon,
        // again increase the space for the next
        if(p < no_cuts - 1 && mapped_cut == cut_map[p + 1]) {
          if(previous_cut_map != mapped_cut) {
            weight_stolen_from_previous_part =
              host_thread_cut_line_weight_to_put_left(p);
          }
          else {
            weight_stolen_from_previous_part +=
              host_thread_cut_line_weight_to_put_left(p);
          }
        }
        else{
          weight_stolen_from_previous_part =
            -host_thread_cut_line_weight_to_put_left(p);
        }
        // end assignment for part p
        break;
      }
    }
    previous_cut_map = mapped_cut;
  }

  // TODO commented for reverse order
  // put everything left on the last cut to the last part.
  mj_lno_t sort_vector_end = (mj_lno_t)sort_vector_points_on_cut[
    previous_cut_map].size() - 1;

  // mj_lno_t sort_vector_begin= 0;
  // mj_lno_t sort_vector_size = (mj_lno_t)
  //   sort_vector_points_on_cut[previous_cut_map].size();
  // TODO commented for reverse order
  for(; sort_vector_end >= 0; --sort_vector_end) {
    // TODO commented for reverse order
    multiSItem t = sort_vector_points_on_cut[previous_cut_map][sort_vector_end];
    // multiSItem t =
    //   sort_vector_points_on_cut[previous_cut_map][sort_vector_begin];
    mj_lno_t i = t.index;
    ++host_thread_point_counts(no_cuts);
    host_assigned_part_ids(i) = no_cuts;
  }

  sort_vector_points_on_cut[previous_cut_map].clear();
  delete [] cut_map;

  //free the memory allocated for vertex sort items .
  mj_lno_t vSize = (mj_lno_t) allocated_memory.size();
  for(mj_lno_t i = 0; i < vSize; ++i) {
    delete [] allocated_memory[i];
  }

  auto local_out_part_xadj = out_part_xadj;
  auto host_out_part_xadj = Kokkos::create_mirror_view(local_out_part_xadj);
  Kokkos::deep_copy(host_out_part_xadj, out_part_xadj);

  // creation of part_xadj as in usual case.
  for(mj_part_t j = 0; j < num_parts; ++j) {
    host_out_part_xadj(j) = host_thread_point_counts(j);
    host_thread_point_counts(j) = 0;
  }

  // perform prefix sum for num_points in parts.
  for(mj_part_t j = 1; j < num_parts; ++j) {
    host_out_part_xadj(j) += host_out_part_xadj(j - 1);
  }

  // shift the num points in threads thread to obtain the
  // beginning index of each thread's private space.
  for(mj_part_t j = 1; j < num_parts; ++j) {
    host_thread_point_counts(j) += host_out_part_xadj(j - 1);
  }

  auto host_new_coordinate_permutations =
    Kokkos::create_mirror_view(new_coordinate_permutations);
  Kokkos::deep_copy(host_new_coordinate_permutations,
    new_coordinate_permutations);

  // now thread gets the coordinate and writes the index of coordinate to
  // the permutation array using the part index we calculated.
  for(mj_lno_t ii = coordinate_begin; ii < coordinate_end; ++ii) {
    mj_lno_t i = host_coordinate_permutations(ii);
    mj_part_t p = host_assigned_part_ids(i);
    host_new_coordinate_permutations(coordinate_begin +
      host_thread_point_counts(p)++) = i;
  }

  Kokkos::deep_copy(thread_point_counts, host_thread_point_counts);
  Kokkos::deep_copy(new_coordinate_permutations,
    host_new_coordinate_permutations);
  Kokkos::deep_copy(local_out_part_xadj, host_out_part_xadj);
}

/*! \brief Function sends the found partids to the owner of the coordinates, if
 * the data is ever migrated. otherwise, it seets the part numbers and returns.
 * \param current_num_parts is the number of parts in the process.
 * \param output_part_begin_index is the number that will be used as beginning
 * part number
 * \param output_part_boxes is the array that holds the part boxes
 * \param is_data_ever_migrated true if data was migrated
 * if the data is ever migrated during the partitioning.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  set_final_parts(
  mj_part_t current_num_parts,
  mj_part_t output_part_begin_index,
  RCP<mj_partBoxVector_t> &output_part_boxes,
  bool is_data_ever_migrated)
{
    this->mj_env->timerStart(MACRO_TIMERS,
      mj_timer_base_string + "Part_Assignment");

    auto local_part_xadj = part_xadj;
    auto local_mj_keep_part_boxes = mj_keep_part_boxes;
    auto local_coordinate_permutations = coordinate_permutations;
    auto local_assigned_part_ids = assigned_part_ids;

    if(local_mj_keep_part_boxes) {
      for(int i = 0; i < current_num_parts; ++i) {
        (*output_part_boxes)[i].setpId(i + output_part_begin_index);
      }
    }

    Kokkos::TeamPolicy<typename mj_node_t::execution_space> policy(
      current_num_parts, Kokkos::AUTO());
    typedef typename Kokkos::TeamPolicy<typename mj_node_t::execution_space>::
      member_type member_type;
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(member_type team_member) {
      int i = team_member.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange (team_member, (i != 0) ?
        local_part_xadj(i-1) : 0, local_part_xadj(i)),
        [=] (mj_lno_t ii) {
        mj_lno_t k = local_coordinate_permutations(ii);
        local_assigned_part_ids(k) = i + output_part_begin_index;
      });
    });

    if(is_data_ever_migrated) {
#ifdef ZOLTAN2_MJ_ENABLE_ZOLTAN_MIGRATION
    if(sizeof(mj_lno_t) <=  sizeof(int)) {

      // Cannot use Zoltan_Comm with local ordinals larger than ints.
      // In Zoltan_Comm_Create, the cast int(this->num_local_coords)
      // may overflow.

      // if data is migrated, then send part numbers to the original owners.
      ZOLTAN_COMM_OBJ *plan = NULL;
      MPI_Comm mpi_comm = Teuchos::getRawMpiComm(*(this->mj_problemComm));

      int incoming = 0;
      int message_tag = 7856;

      this->mj_env->timerStart(MACRO_TIMERS,
        mj_timer_base_string + "Final Z1PlanCreating");

      // setup incoming count
      int ierr = Zoltan_Comm_Create( &plan, int(this->num_local_coords),
        this->owner_of_coordinate.data(), mpi_comm, message_tag, &incoming);

      Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
      this->mj_env->timerStop(MACRO_TIMERS,
        mj_timer_base_string + "Final Z1PlanCreating" );

      this->mj_env->timerStart(MACRO_TIMERS,
        mj_timer_base_string + "Final Z1PlanComm");

      // MPI Buffers should be on Kokkos::HostSpace not Kokkos::CudaUVMSpace

      // migrate gnos to actual owners.
      auto host_current_mj_gnos = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), this->current_mj_gnos);
      deep_copy(host_current_mj_gnos, this->current_mj_gnos);
      Kokkos::View<mj_gno_t*, device_t> dst_gnos(
        Kokkos::ViewAllocateWithoutInitializing("dst_gnos"), incoming);
      auto host_dst_gnos = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), dst_gnos);
      message_tag++;
      ierr = Zoltan_Comm_Do( plan, message_tag,
        (char *) host_current_mj_gnos.data(),
        sizeof(mj_gno_t), (char *) host_dst_gnos.data());
      Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
      Kokkos::deep_copy(dst_gnos, host_dst_gnos);
      this->current_mj_gnos = dst_gnos;

      // migrate part ids to actual owners.
      auto host_src_part_ids = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), this->assigned_part_ids);
      deep_copy(host_src_part_ids, this->assigned_part_ids);
      Kokkos::View<mj_part_t*, device_t> dst_part_ids(
        Kokkos::ViewAllocateWithoutInitializing("dst_part_ids"), incoming);
      auto host_dst_part_ids = Kokkos::create_mirror_view(
        Kokkos::HostSpace(), dst_part_ids);
      message_tag++;
      ierr = Zoltan_Comm_Do( plan, message_tag,
        (char *) host_src_part_ids.data(),
        sizeof(mj_part_t), (char *) host_dst_part_ids.data());
      Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);
      Kokkos::deep_copy(dst_part_ids, host_dst_part_ids);
      this->assigned_part_ids = dst_part_ids;

      ierr = Zoltan_Comm_Destroy(&plan);
      Z2_ASSERT_VALUE(ierr, ZOLTAN_OK);

      this->num_local_coords = incoming;

      this->mj_env->timerStop(MACRO_TIMERS,
        mj_timer_base_string + "Final Z1PlanComm");
    }
    else
#endif  // ZOLTAN2_MJ_ENABLE_ZOLTAN_MIGRATION
    {
      // setup incoming count
      this->mj_env->timerStart(MACRO_TIMERS,
        mj_timer_base_string + "Final DistributorPlanCreating");
      Tpetra::Distributor distributor(this->mj_problemComm);
      ArrayView<const mj_part_t> owners_of_coords(
        this->owner_of_coordinate.data(), this->num_local_coords);
      mj_lno_t incoming = distributor.createFromSends(owners_of_coords);
      this->mj_env->timerStop(MACRO_TIMERS,
        mj_timer_base_string + "Final DistributorPlanCreating" );

      this->mj_env->timerStart(MACRO_TIMERS,
        mj_timer_base_string + "Final DistributorPlanComm");

      // MPI buffers should be Kokkos::HostSpace, not Kokkos::CudaUVMSpace

      // migrate gnos to actual owners.
      auto src_host_current_mj_gnos =
        Kokkos::create_mirror_view(Kokkos::HostSpace(), this->current_mj_gnos);
      Kokkos::deep_copy(src_host_current_mj_gnos, this->current_mj_gnos);
      ArrayRCP<mj_gno_t> received_gnos(incoming);
      ArrayView<mj_gno_t> sent_gnos(src_host_current_mj_gnos.data(),
        this->num_local_coords);
      distributor.doPostsAndWaits<mj_gno_t>(sent_gnos, 1, received_gnos());
      this->current_mj_gnos = Kokkos::View<mj_gno_t*, device_t>(
        Kokkos::ViewAllocateWithoutInitializing("current_mj_gnos"), incoming);
      auto host_current_mj_gnos = Kokkos::create_mirror_view(
        this->current_mj_gnos);
      memcpy(host_current_mj_gnos.data(),
        received_gnos.getRawPtr(), incoming * sizeof(mj_gno_t));
      Kokkos::deep_copy(this->current_mj_gnos, host_current_mj_gnos);

      // migrate part ids to actual owners.
      auto src_host_assigned_part_ids =
        Kokkos::create_mirror_view(Kokkos::HostSpace(), this->assigned_part_ids);
      Kokkos::deep_copy(src_host_assigned_part_ids, this->assigned_part_ids);
      ArrayView<mj_part_t> sent_partids(src_host_assigned_part_ids.data(),
        this->num_local_coords);
      ArrayRCP<mj_part_t> received_partids(incoming);
      distributor.doPostsAndWaits<mj_part_t>(
        sent_partids, 1, received_partids());
      this->assigned_part_ids =
        Kokkos::View<mj_part_t*, device_t>(
          Kokkos::ViewAllocateWithoutInitializing("assigned_part_ids"),
          incoming);
      auto host_assigned_part_ids = Kokkos::create_mirror_view(
        this->assigned_part_ids);
      memcpy( host_assigned_part_ids.data(),
        received_partids.getRawPtr(), incoming * sizeof(mj_part_t));
      deep_copy(this->assigned_part_ids, host_assigned_part_ids);
      this->num_local_coords = incoming;

      this->mj_env->timerStop(MACRO_TIMERS,
        mj_timer_base_string + "Final DistributorPlanComm");
    }
  }

  this->mj_env->timerStop(MACRO_TIMERS,
    mj_timer_base_string + "Part_Assignment");

  this->mj_env->timerStart(MACRO_TIMERS,
    mj_timer_base_string + "Solution_Part_Assignment");

  // ArrayRCP<mj_part_t> partId;
  // partId = arcp(this->assigned_part_ids, 0, this->num_local_coords, true);

  if(this->mj_keep_part_boxes) {
    this->kept_boxes = compute_global_box_boundaries(output_part_boxes);
  }

  this->mj_env->timerStop(MACRO_TIMERS,
    mj_timer_base_string + "Solution_Part_Assignment");
}

/*!\brief Multi Jagged  coordinate partitioning algorithm.
 * \param distribute_points_on_cut_lines_ :  if partitioning can distribute
 * points on same coordinate to different parts.
 * \param max_concurrent_part_calculation_ : how many parts we can calculate
 * concurrently.
 * \param check_migrate_avoid_migration_option_ : whether to migrate=1,
 * avoid migrate=2, or leave decision to MJ=0
 * \param minimum_migration_imbalance_  : when MJ decides whether to migrate,
 * the minimum imbalance for migration.
 * \param migration_type : whether to migrate for perfect load imbalance (0) or
 * less messages.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  set_partitioning_parameters(
  bool distribute_points_on_cut_lines_,
  int max_concurrent_part_calculation_,
  int check_migrate_avoid_migration_option_,
  double minimum_migration_imbalance_,
  int migration_type_)
{
  this->distribute_points_on_cut_lines = distribute_points_on_cut_lines_;
  this->max_concurrent_part_calculation = max_concurrent_part_calculation_;
  this->check_migrate_avoid_migration_option =
    check_migrate_avoid_migration_option_;
  this->minimum_migration_imbalance = minimum_migration_imbalance_;
  this->migration_type = migration_type_;
}

/*! \brief Multi Jagged  coordinate partitioning algorithm.
 * \param env   library configuration and problem parameters
 * \param problemComm the communicator for the problem
 * \param imbalance_tolerance : the input provided imbalance tolerance.
 * \param num_teams : number of teams for CUDA kernels
 * \param num_global_parts: number of target global parts.
 * \param part_no_array: part no array, if provided this will be used for
 * partitioning.
 * \param recursion_depth: if part no array is provided, it is the length of
 * part no array, if part no is not provided than it is the number of steps that
 * algorithm will divide into num_global_parts parts.
 * \param coord_dim: coordinate dimension
 * \param num_local_coords: number of local coordinates
 * \param num_global_coords: number of global coordinates
 * \param initial_mj_gnos: the list of initial global id's
 * \param mj_coordinates: the two dimensional coordinate array.
 * \param num_weights_per_coord: number of weights per coordinate
 * \param mj_uniform_weights: if weight index [i] has uniform weight or not.
 * \param mj_weights: the two dimensional array for weights
 * \param mj_uniform_parts: if the target partitioning aims uniform parts
 * \param result_assigned_part_ids: Output - 1D pointer, should be provided as
 * null. Memory is given in the function. the result partids corresponding to
 * the coordinates given in result_mj_gnos.
 * \param result_mj_gnos: Output - 1D pointer, should be provided as null.
 * Memory is given in the function. the result coordinate global id's
 * corresponding to the part_ids array.
 */
template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
void AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t>::
  multi_jagged_part(
  const RCP<const Environment> &env,
  RCP<const Comm<int> > &problemComm,
  double imbalance_tolerance_,
  int num_teams_,
  size_t num_global_parts_,
  Kokkos::View<mj_part_t*, Kokkos::HostSpace> & part_no_array_,
  int recursion_depth_,
  int coord_dim_,
  mj_lno_t num_local_coords_,
  mj_gno_t num_global_coords_,
  Kokkos::View<const mj_gno_t*, device_t> & initial_mj_gnos_,
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t> & mj_coordinates_,
  int num_weights_per_coord_,
  Kokkos::View<bool*, Kokkos::HostSpace> & mj_uniform_weights_,
  Kokkos::View<mj_scalar_t**, device_t> & mj_weights_,
  Kokkos::View<bool*, Kokkos::HostSpace> & mj_uniform_parts_,
  Kokkos::View<mj_part_t *, device_t> & result_assigned_part_ids_,
  Kokkos::View<mj_gno_t*, device_t> & result_mj_gnos_)
{

  // see comment above for Zoltan2_AlgMJ_TrackCallsCounter
  int execute_counter = Zoltan2_AlgMJ_TrackCallsCounter::get_counter_AlgMJ();
  this->mj_timer_base_string = "MJ(" + std::to_string(execute_counter) + ") - ";

  this->mj_env = env;
  this->mj_problemComm = problemComm;
  this->myActualRank = this->myRank = this->mj_problemComm->getRank();
  this->mj_env->timerStart(MACRO_TIMERS,
    mj_timer_base_string + "Total");
  this->mj_env->debug(3, "In MultiJagged Jagged");
  this->imbalance_tolerance = imbalance_tolerance_;
  this->mj_num_teams = num_teams_;
  this->num_global_parts = num_global_parts_;
  this->part_no_array = part_no_array_;
  this->recursion_depth = recursion_depth_;
  this->coord_dim = coord_dim_;
  this->num_local_coords = num_local_coords_;
  this->num_global_coords = num_global_coords_;
  this->mj_coordinates = mj_coordinates_;
  this->initial_mj_gnos = initial_mj_gnos_;
  this->num_weights_per_coord = num_weights_per_coord_;
  this->mj_uniform_weights = mj_uniform_weights_;
  this->mj_weights = mj_weights_;
  this->mj_uniform_parts = mj_uniform_parts_;

  // this->set_input_data();

  this->set_part_specifications();

  this->mj_env->timerStart(MACRO_TIMERS,
    mj_timer_base_string + "Allocate Views");
  this->allocate_set_work_memory();
  this->mj_env->timerStop(MACRO_TIMERS,
    mj_timer_base_string + "Allocate Views");

  // We duplicate the comm as we create subcommunicators during migration.
  // We keep the problemComm as it is, while comm changes after each migration.
  this->comm = this->mj_problemComm->duplicate();

#ifdef print_debug
  if(comm->getRank() == 0) {
    std::cout << "size of gno:" << sizeof(mj_gno_t) << std::endl;
    std::cout << "size of lno:" << sizeof(mj_lno_t) << std::endl;
    std::cout << "size of mj_scalar_t:" << sizeof(mj_scalar_t) << std::endl;
  }
#endif

  // initially there is a single partition
  mj_part_t current_num_parts = 1;
  Kokkos::View<mj_scalar_t *, device_t> current_cut_coordinates =
    this->all_cut_coordinates;
  this->mj_env->timerStart(MACRO_TIMERS,
    mj_timer_base_string + "Problem_Partitioning");
  mj_part_t output_part_begin_index = 0;
  mj_part_t future_num_parts = this->total_num_part;
  bool is_data_ever_migrated = false;

  std::vector<mj_part_t> *future_num_part_in_parts =
    new std::vector<mj_part_t> ();
  std::vector<mj_part_t> *next_future_num_parts_in_parts =
    new std::vector<mj_part_t> ();

  next_future_num_parts_in_parts->push_back(this->num_global_parts);

  RCP<mj_partBoxVector_t> input_part_boxes;
  RCP<mj_partBoxVector_t> output_part_boxes;

  if(this->mj_keep_part_boxes) {
    input_part_boxes = RCP<mj_partBoxVector_t>(new mj_partBoxVector_t(), true);
    output_part_boxes = RCP<mj_partBoxVector_t>(new mj_partBoxVector_t(), true);
    compute_global_box();
    this->init_part_boxes(output_part_boxes);
  }

  auto local_part_xadj = this->part_xadj;

  // Need a device counter - how best to allocate?
  // Putting this allocation in the loops is very costly so moved out here.
  Kokkos::View<mj_part_t*, device_t>
    view_rectilinear_cut_count("view_rectilinear_cut_count", 1);
  Kokkos::View<size_t*, device_t>
    view_total_reduction_size("view_total_reduction_size", 1);

  for(int i = 0; i < this->recursion_depth; ++i) {

    // convert i to string to be used for debugging purposes.
    std::string istring = std::to_string(i);

    // next_future_num_parts_in_parts will be as the size of outnumParts,
    // and this will hold how many more parts that each output part
    // should be divided. this array will also be used to determine the weight
    // ratios of the parts. swap the arrays to use iteratively.
    std::vector<mj_part_t> *tmpPartVect= future_num_part_in_parts;
    future_num_part_in_parts = next_future_num_parts_in_parts;
    next_future_num_parts_in_parts = tmpPartVect;

    // clear next_future_num_parts_in_parts array as
    // getPartitionArrays expects it to be empty.
    next_future_num_parts_in_parts->clear();
    if(this->mj_keep_part_boxes) {
      RCP<mj_partBoxVector_t> tmpPartBoxes = input_part_boxes;
      input_part_boxes = output_part_boxes;
      output_part_boxes = tmpPartBoxes;
      output_part_boxes->clear();
    }

    // returns the total no. of output parts for this dimension partitioning.
    mj_part_t output_part_count_in_dimension =
      this->update_part_num_arrays(
        future_num_part_in_parts,
        next_future_num_parts_in_parts,
        future_num_parts,
        current_num_parts,
        i,
        input_part_boxes,
        output_part_boxes, 1);

    // if the number of obtained parts equal to current number of parts,
    // skip this dimension. For example, this happens when 1 is given in the
    // input part array is given. P=4,5,1,2
    if(output_part_count_in_dimension == current_num_parts) {
      //still need to swap the input output arrays.
      tmpPartVect= future_num_part_in_parts;
      future_num_part_in_parts = next_future_num_parts_in_parts;
      next_future_num_parts_in_parts = tmpPartVect;

      if(this->mj_keep_part_boxes) {
        RCP<mj_partBoxVector_t> tmpPartBoxes = input_part_boxes;
        input_part_boxes = output_part_boxes;
        output_part_boxes = tmpPartBoxes;
      }
      continue;
    }

    // get the coordinate axis along which the partitioning will be done.
    int coordInd = i % this->coord_dim;

    Kokkos::View<mj_scalar_t *, device_t> mj_current_dim_coords;
    // view could be size 0 if graph was not distributed
    if(this->mj_coordinates.extent(0) != 0) {
      mj_current_dim_coords =
        Kokkos::subview(this->mj_coordinates, Kokkos::ALL, coordInd);
    }

    this->mj_env->timerStart(MACRO_TIMERS,
      mj_timer_base_string + "Problem_Partitioning_" + istring);

    // alloc Memory to point the indices
    // of the parts in the permutation array.
    this->new_part_xadj = Kokkos::View<mj_lno_t*, device_t>(
      "new part xadj", output_part_count_in_dimension);

    // the index where in the new_part_xadj will be written.
    mj_part_t output_part_index = 0;

    // whatever is written to output_part_index will be added with
    // output_coordinate_end_index so that the points will be shifted.
    mj_part_t output_coordinate_end_index = 0;

    mj_part_t current_work_part = 0;
    mj_part_t current_concurrent_num_parts =
      std::min(current_num_parts - current_work_part,
      this->max_concurrent_part_calculation);

    mj_part_t obtained_part_index = 0;

    auto host_process_local_min_max_coord_total_weight =
      Kokkos::create_mirror_view(process_local_min_max_coord_total_weight);
    auto host_global_min_max_coord_total_weight =
      Kokkos::create_mirror_view(global_min_max_coord_total_weight);

    // run for all available parts.
    for(; current_work_part < current_num_parts;
      current_work_part += current_concurrent_num_parts) {

      current_concurrent_num_parts =
        std::min(current_num_parts - current_work_part,
        this->max_concurrent_part_calculation);

      int bDoingWork_int; // Can't reduce on bool so use int
      auto local_device_num_partitioning_in_current_dim =
        device_num_partitioning_in_current_dim;
      Kokkos::parallel_reduce("Read bDoingWork",
        Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
        KOKKOS_LAMBDA(int dummy, int & set_single) {
        set_single = 0;
        for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {
          if(local_device_num_partitioning_in_current_dim(
            current_work_part + kk) != 1) {
            set_single = 1;
            break;
          }
        }
      }, bDoingWork_int);
      bool bDoingWork = (bDoingWork_int != 0) ? true : false;

      this->mj_get_local_min_max_coord_totW(
        current_work_part,
        current_concurrent_num_parts,
        mj_current_dim_coords);

      // 1D partitioning
      if(bDoingWork) {
        // obtain global Min max of the part.
        this->mj_get_global_min_max_coord_totW(
          current_concurrent_num_parts,
          this->process_local_min_max_coord_total_weight,
          this->global_min_max_coord_total_weight);

        // represents the total number of cutlines
        // whose coordinate should be determined.
        mj_part_t total_incomplete_cut_count = 0;

        // Compute weight ratios for parts & cuts:
        // e.g., 0.25  0.25  0.5    0.5  0.75 0.75  1
        // part0  cut0  part1 cut1 part2 cut2 part3
        mj_part_t concurrent_part_cut_shift = 0;
        mj_part_t concurrent_part_part_shift = 0;

        for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {

            Kokkos::deep_copy(host_global_min_max_coord_total_weight,
              global_min_max_coord_total_weight);

            mj_scalar_t min_coordinate =
              host_global_min_max_coord_total_weight(kk);
            mj_scalar_t max_coordinate =
              host_global_min_max_coord_total_weight(
                kk + current_concurrent_num_parts);

            mj_scalar_t global_total_weight =
              host_global_min_max_coord_total_weight(
                kk + 2 * current_concurrent_num_parts);

          mj_part_t concurrent_current_part_index = current_work_part + kk;

          mj_part_t partition_count = host_num_partitioning_in_current_dim(
              concurrent_current_part_index);

          Kokkos::View<mj_scalar_t *, device_t> usedCutCoordinate =
            Kokkos::subview(current_cut_coordinates,
              std::pair<mj_lno_t, mj_lno_t>(
                concurrent_part_cut_shift, current_cut_coordinates.size()));
          Kokkos::View<mj_scalar_t *, device_t>
            current_target_part_weights =
            Kokkos::subview(target_part_weights,
              std::pair<mj_lno_t, mj_lno_t>(
                concurrent_part_part_shift, target_part_weights.size()));

          // shift the usedCutCoordinate array as noCuts.
          concurrent_part_cut_shift += partition_count - 1;
          // shift the partRatio array as noParts.
          concurrent_part_part_shift += partition_count;

          // calculate only if part is not empty,
          // and part will be further partitioned.
          if(partition_count > 1 && min_coordinate <= max_coordinate) {

            // increase num_cuts_do_be_determined by the number of cuts of the
            // current part's cut line number.
            total_incomplete_cut_count += partition_count - 1;

            this->incomplete_cut_count(kk) = partition_count - 1;

            // get the target weights of the parts
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

            mj_lno_t coordinate_end_index =
              host_part_xadj(concurrent_current_part_index);
            mj_lno_t coordinate_begin_index =
              concurrent_current_part_index==0 ? 0 :
                host_part_xadj(concurrent_current_part_index - 1);

            this->set_initial_coordinate_parts(
              max_coordinate,
              min_coordinate,
              coordinate_begin_index, coordinate_end_index,
              this->coordinate_permutations,
              mj_current_dim_coords,
              this->assigned_part_ids,
              partition_count);
          }
          else {
            // e.g., if have fewer coordinates than parts, don't need to do
            // next dim.
            this->incomplete_cut_count(kk) = 0;
          }

          obtained_part_index += partition_count;
        }

        // used imbalance, it is always 0, as it is difficult to
        // estimate a range.
        double used_imbalance = 0;
        // Determine cut lines for all concurrent parts parts here.
        this->mj_env->timerStart(MACRO_TIMERS,
          mj_timer_base_string + "Problem_Partitioning Get Part Weights");

        this->mj_1D_part(
          mj_current_dim_coords,
          used_imbalance,
          current_work_part,
          current_concurrent_num_parts,
          current_cut_coordinates,
          total_incomplete_cut_count,
          view_rectilinear_cut_count,
          view_total_reduction_size);

        this->mj_env->timerStop(MACRO_TIMERS,
          mj_timer_base_string + "Problem_Partitioning Get Part Weights");
      }

      // create new part chunks
      {
        mj_part_t output_array_shift = 0;
        mj_part_t cut_shift = 0;
        size_t tlr_shift = 0;
        size_t partweight_array_shift = 0;
        for(int kk = 0; kk < current_concurrent_num_parts; ++kk) {

          mj_part_t current_concurrent_work_part = current_work_part + kk;

          mj_part_t num_parts = host_num_partitioning_in_current_dim(
            current_concurrent_work_part);

          // if the part is empty, skip the part.
          int coordinateA_bigger_than_coordinateB =
            host_global_min_max_coord_total_weight(kk) >
            host_global_min_max_coord_total_weight(
              kk + current_concurrent_num_parts);

          if((num_parts != 1) && coordinateA_bigger_than_coordinateB) {
            // we still need to write the begin and end point of the empty part.
            // simply set it zero, the array indices will be shifted later
            auto local_new_part_xadj = this->new_part_xadj;
            Kokkos::parallel_for(
              Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
                (0, num_parts), KOKKOS_LAMBDA (mj_part_t jj) {
                local_new_part_xadj(
                  output_part_index + output_array_shift + jj) = 0;
            });

            cut_shift += num_parts - 1;
            tlr_shift += (4 *(num_parts - 1) + 1);
            output_array_shift += num_parts;
            partweight_array_shift += (2 * (num_parts - 1) + 1);
            continue;
          }

          Kokkos::View<mj_scalar_t *, device_t>
            current_concurrent_cut_coordinate =
            Kokkos::subview(current_cut_coordinates,
              std::pair<mj_lno_t, mj_lno_t>(
                cut_shift,
                current_cut_coordinates.size()));
          Kokkos::View<mj_scalar_t *, device_t>
            used_local_cut_line_weight_to_left =
            Kokkos::subview(process_cut_line_weight_to_put_left,
              std::pair<mj_lno_t, mj_lno_t>(
                cut_shift,
                process_cut_line_weight_to_put_left.size()));

          this->thread_part_weight_work =
            Kokkos::subview(
              this->thread_part_weights,
              std::pair<mj_lno_t, mj_lno_t>(
                partweight_array_shift,
                this->thread_part_weights.extent(0)));

          if(num_parts > 1) {
            if(this->mj_keep_part_boxes) {
              // if part boxes are to be stored update the boundaries.
              for(mj_part_t j = 0; j < num_parts - 1; ++j) {
                mj_scalar_t temp_get_val;
                Kokkos::parallel_reduce("Read single",
                  Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
                  KOKKOS_LAMBDA(int dummy, mj_scalar_t & set_single) {
                  set_single = current_concurrent_cut_coordinate(j);
                }, temp_get_val);
                (*output_part_boxes)
                  [output_array_shift + output_part_index + j].
                  updateMinMax(temp_get_val, 1 /*update max*/, coordInd);
                (*output_part_boxes)
                  [output_array_shift + output_part_index + j + 1].
                  updateMinMax(temp_get_val, 0 /*update max*/, coordInd);
              }
            }

            // Rewrite the indices based on the computed cuts.
            Kokkos::View<mj_lno_t*, device_t> sub_new_part_xadj =
              Kokkos::subview(this->new_part_xadj,
                std::pair<mj_lno_t, mj_lno_t>(
                  output_part_index + output_array_shift,
                  this->new_part_xadj.size()));

            this->mj_create_new_partitions(
              num_parts,
              current_concurrent_work_part,
              mj_current_dim_coords,
              current_concurrent_cut_coordinate,
              used_local_cut_line_weight_to_left,
              sub_new_part_xadj);
          }
          else {

            mj_lno_t coordinate_end = host_part_xadj(
              current_concurrent_work_part);
            mj_lno_t coordinate_begin =
              current_concurrent_work_part==0 ? 0 : host_part_xadj(
              current_concurrent_work_part - 1);

            // if this part is partitioned into 1 then just copy
            // the old values.
            mj_lno_t part_size = coordinate_end - coordinate_begin;

            // Awkward here to set one value - need some broader
            // refactoring to improve this one.
            auto local_new_part_xadj = this->new_part_xadj;
            Kokkos::parallel_for(
              Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
              (0, 1), KOKKOS_LAMBDA (int dummy) {
              local_new_part_xadj(
                output_part_index + output_array_shift) = part_size;
            });

            auto subview_new_coordinate_permutations =
              Kokkos::subview(this->new_coordinate_permutations,
                std::pair<mj_lno_t, mj_lno_t>(
                  coordinate_begin,
                  coordinate_begin + part_size));
            auto subview_coordinate_permutations =
              Kokkos::subview(this->coordinate_permutations,
                std::pair<mj_lno_t, mj_lno_t>(
                  coordinate_begin,
                  coordinate_begin + part_size));
            Kokkos::deep_copy(subview_new_coordinate_permutations,
              subview_coordinate_permutations);
          }
          cut_shift += num_parts - 1;
          output_array_shift += num_parts;
          partweight_array_shift += (2 * (num_parts - 1) + 1);
        }

        // shift cut coordinates so that all cut coordinates are stored.
        // no shift now because we dont keep the cuts.
        // current_cut_coordinates += cut_shift;
        // mj_create_new_partitions from coordinates partitioned the parts
        // and write the indices as if there were a single part.
        // now we need to shift the beginning indices.
        for(mj_part_t kk = 0; kk < current_concurrent_num_parts; ++kk) {
          mj_part_t num_parts =
            host_num_partitioning_in_current_dim(current_work_part + kk);

          // These two kernels are a bit awkward but need broader redesign to
          // avoid this situation.
          auto local_new_part_xadj = this->new_part_xadj;
          Kokkos::parallel_for(
            Kokkos::RangePolicy<typename mj_node_t::execution_space, mj_part_t>
            (0, num_parts), KOKKOS_LAMBDA (mj_part_t ii) {
            local_new_part_xadj(output_part_index+ii) +=
              output_coordinate_end_index;
          });

          // increase the previous count by current end.
          mj_part_t temp_get;
          Kokkos::parallel_reduce("Read single",
            Kokkos::RangePolicy<typename mj_node_t::execution_space, int> (0, 1),
            KOKKOS_LAMBDA(int dummy, mj_part_t & set_single) {
            set_single =
              local_new_part_xadj(output_part_index + num_parts - 1);
          }, temp_get);
          output_coordinate_end_index = temp_get;
          //increase the current out.
          output_part_index += num_parts;
        }
      }
    }

    // end of this partitioning dimension
    int current_world_size = this->comm->getSize();
    long migration_reduce_all_population =
      this->total_dim_num_reduce_all * current_world_size;
    bool is_migrated_in_current_dimension = false;

    // we migrate if there are more partitionings to be done after this step
    // and if the migration is not forced to be avoided.
    // and the operation is not sequential.
    if(future_num_parts > 1 &&
      this->check_migrate_avoid_migration_option >= 0 &&
      current_world_size > 1) {
      this->mj_env->timerStart(MACRO_TIMERS,
        mj_timer_base_string + "Problem_Migration-" + istring);
      mj_part_t num_parts = output_part_count_in_dimension;

      if(this->mj_perform_migration(
        num_parts,
        current_num_parts, //output
        next_future_num_parts_in_parts, //output
        output_part_begin_index,
        migration_reduce_all_population,
        this->num_global_coords / (future_num_parts * current_num_parts),
        istring,
        input_part_boxes, output_part_boxes) )
      {
        is_migrated_in_current_dimension = true;
        is_data_ever_migrated = true;
        this->mj_env->timerStop(MACRO_TIMERS,
          mj_timer_base_string + "Problem_Migration-" + istring);
        // since data is migrated, we reduce the number of reduceAll
        // operations for the last part.
        this->total_dim_num_reduce_all /= num_parts;
      }
      else {
        is_migrated_in_current_dimension = false;
        this->mj_env->timerStop(MACRO_TIMERS,
          mj_timer_base_string + "Problem_Migration-" + istring);
      }
    }

    // swap the coordinate permutations for the next dimension.
    Kokkos::View<mj_lno_t*, device_t> tmp =
      this->coordinate_permutations;
    this->coordinate_permutations =
      this->new_coordinate_permutations;

    this->new_coordinate_permutations = tmp;
    if(!is_migrated_in_current_dimension) {
      this->total_dim_num_reduce_all -= current_num_parts;
      current_num_parts = output_part_count_in_dimension;
    }

    {
      this->part_xadj = this->new_part_xadj;
      local_part_xadj = this->new_part_xadj;
      this->host_part_xadj = Kokkos::create_mirror_view(part_xadj);
      Kokkos::deep_copy(host_part_xadj, part_xadj); // keep in sync

      this->new_part_xadj = Kokkos::View<mj_lno_t*, device_t>("empty", 0);
      this->mj_env->timerStop(MACRO_TIMERS,
        mj_timer_base_string + "Problem_Partitioning_" + istring);
    }
  }

  // Partitioning is done
  delete future_num_part_in_parts;
  delete next_future_num_parts_in_parts;
  this->mj_env->timerStop(MACRO_TIMERS,
    mj_timer_base_string + "Problem_Partitioning");
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
  this->mj_env->timerStop(MACRO_TIMERS,
    mj_timer_base_string + "Total");
  this->mj_env->debug(3, "Out of MultiJagged");
}

template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
RCP<typename AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t, mj_node_t>::
  mj_partBoxVector_t>
AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t, mj_node_t>::
  get_kept_boxes() const
{
  if(this->mj_keep_part_boxes) {
    return this->kept_boxes;
  }
  else {
    throw std::logic_error("Error: part boxes are not stored.");
  }
}

template <typename mj_scalar_t, typename mj_lno_t, typename mj_gno_t,
  typename mj_part_t, typename mj_node_t>
RCP<typename AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t, mj_node_t>::
  mj_partBoxVector_t>
AlgMJ<mj_scalar_t,mj_lno_t,mj_gno_t,mj_part_t, mj_node_t>::
  compute_global_box_boundaries(RCP<mj_partBoxVector_t> &localPartBoxes) const
{
  typedef typename Zoltan2::coordinateModelPartBox::coord_t coord_t;
  mj_part_t ntasks = this->num_global_parts;
  int dim = (*localPartBoxes)[0].getDim();
  coord_t *localPartBoundaries = new coord_t[ntasks * 2 *dim];

  memset(localPartBoundaries, 0, sizeof(coord_t) * ntasks * 2 *dim);

  coord_t *globalPartBoundaries = new coord_t[ntasks * 2 *dim];
  memset(globalPartBoundaries, 0, sizeof(coord_t) * ntasks * 2 *dim);

  coord_t *localPartMins = localPartBoundaries;
  coord_t *localPartMaxs = localPartBoundaries + ntasks * dim;

  coord_t *globalPartMins = globalPartBoundaries;
  coord_t *globalPartMaxs = globalPartBoundaries + ntasks * dim;

  mj_part_t boxCount = localPartBoxes->size();
  for(mj_part_t i = 0; i < boxCount; ++i) {
    mj_part_t pId = (*localPartBoxes)[i].getpId();

    // cout << "me:" << comm->getRank() << " has:" << pId << endl;

    coord_t *lmins = (*localPartBoxes)[i].getlmins();
    coord_t *lmaxs = (*localPartBoxes)[i].getlmaxs();

    for(int j = 0; j < dim; ++j) {
      localPartMins[dim * pId + j] = lmins[j];
      localPartMaxs[dim * pId + j] = lmaxs[j];

      /*
      std::cout << "me:" << comm->getRank()  <<
              " dim * pId + j:"<< dim * pId + j <<
              " localMin:" << localPartMins[dim * pId + j] <<
              " localMax:" << localPartMaxs[dim * pId + j] << std::endl;
      */
    }
  }

  Teuchos::Zoltan2_BoxBoundaries<int, coord_t> reductionOp(ntasks * 2 *dim);

  reduceAll<int, coord_t>(*mj_problemComm, reductionOp,
    ntasks * 2 *dim, localPartBoundaries, globalPartBoundaries);

  RCP<mj_partBoxVector_t> pB(new mj_partBoxVector_t(),true);
  for(mj_part_t i = 0; i < ntasks; ++i) {
    Zoltan2::coordinateModelPartBox tpb(i, dim,
      globalPartMins + dim * i,
      globalPartMaxs + dim * i);

    /*
    for(int j = 0; j < dim; ++j) {
        std::cout << "me:" << comm->getRank()  <<
                " dim * pId + j:"<< dim * i + j <<
                " globalMin:" << globalPartMins[dim * i + j] <<
                " globalMax:" << globalPartMaxs[dim * i + j] << std::endl;
    }
    */

    pB->push_back(tpb);
  }
  delete []localPartBoundaries;
  delete []globalPartBoundaries;
  //RCP <mj_partBoxVector_t> tmpRCPBox(pB, true);
  return pB;
}

/*! \brief Multi Jagged coordinate partitioning algorithm.
 */
template <typename Adapter>
class Zoltan2_AlgMJ : public Algorithm<Adapter>
{

private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef CoordinateModel<typename Adapter::base_adapter_t> coordinateModel_t;

  // For coordinates and weights, MJ needs floats or doubles
  // But Adapter can provide other scalars, e.g., ints.
  // So have separate scalar_t for MJ and adapter.
  typedef typename Adapter::scalar_t adapter_scalar_t;

  // Provide a default type for mj_scalar_t;
  typedef float default_mj_scalar_t;

  // If Adapter provided float or double scalar_t, use it (prevents copies).
  // Otherwise, use the default type of mj_scalar_t;
  typedef typename
          std::conditional<
               (std::is_same<adapter_scalar_t, float>::value ||
                std::is_same<adapter_scalar_t, double>::value),
               adapter_scalar_t, default_mj_scalar_t>::type   mj_scalar_t;

  typedef typename Adapter::gno_t mj_gno_t;
  typedef typename Adapter::lno_t mj_lno_t;
  typedef typename Adapter::part_t mj_part_t;
  typedef typename Adapter::node_t mj_node_t;
  typedef coordinateModelPartBox mj_partBox_t;
  typedef std::vector<mj_partBox_t> mj_partBoxVector_t;
  typedef typename mj_node_t::device_type device_t;
#endif

  AlgMJ<mj_scalar_t, mj_lno_t, mj_gno_t, mj_part_t, mj_node_t> mj_partitioner;

  RCP<const Environment> mj_env; // the environment object
  RCP<const Comm<int> > mj_problemComm; // initial comm object
  RCP<const coordinateModel_t> mj_coords; // coordinate adapter

  // PARAMETERS
  double imbalance_tolerance; // input imbalance tolerance.

  int num_teams; // how many teams to run main loop with

  size_t num_global_parts; // the targeted number of parts

  // input part array specifying num part to divide along each dim.
  Kokkos::View<mj_part_t*, Kokkos::HostSpace> part_no_array;

  // the number of steps that partitioning will be solved in.
  int recursion_depth;

  int coord_dim; // coordinate dimension.
  mj_lno_t num_local_coords; //number of local coords.
  mj_gno_t num_global_coords; //number of global coords.

  // initial global ids of the coordinates.
  Kokkos::View<const mj_gno_t*, device_t> initial_mj_gnos;

  // two dimension coordinate array.
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t>
    mj_coordinates;

  int num_weights_per_coord; // number of weights per coordinate

  // if the target parts are uniform.
  Kokkos::View<bool*, Kokkos::HostSpace> mj_uniform_weights;

  // two dimensional weight array.
  Kokkos::View<mj_scalar_t**, device_t> mj_weights;

  // if the target parts are uniform
  Kokkos::View<bool*, Kokkos::HostSpace> mj_uniform_parts;

  // Nonuniform first level partitioning
  // Currently used for Dragonfly task mapping by partitioning Dragonfly RCA
  // machine coordinates and application coordinates.
  // An optimization that completely partitions the most important machine
  // dimension first (i.e. the Dragonfly group coordinate, or RCA's x
  // coordinate). The standard MJ alg follows after the nonuniform first level
  // partitioning.
  // If used, number of parts for the first level partitioning
  mj_part_t num_first_level_parts;

  // If used, the distribution of parts for the nonuniform
  // first level partitioning
  Kokkos::View<mj_part_t*, Kokkos::HostSpace> first_level_distribution;

  // if partitioning can distribute points on same coordiante to
  // different parts.
  bool distribute_points_on_cut_lines;

  // how many parts we can calculate concurrently.
  mj_part_t max_concurrent_part_calculation;

  // whether to migrate=1, avoid migrate=2, or leave decision to MJ=0
  int check_migrate_avoid_migration_option;

  // when doing the migration, 0 will aim for perfect load-imbalance,
  int migration_type;

  // 1 for minimized messages

  // when MJ decides whether to migrate, the minimum imbalance for migration.
  double minimum_migration_imbalance;
  bool mj_keep_part_boxes; //if the boxes need to be kept.

  // if this is set, then recursion depth is adjusted to its maximum value.
  bool mj_run_as_rcb;
  int mj_premigration_option;
  int min_coord_per_rank_for_premigration;

  // communication graph xadj
  ArrayRCP<mj_part_t> comXAdj_;

  // communication graph adj.
  ArrayRCP<mj_part_t> comAdj_;

  void copy(
    const RCP<PartitioningSolution<Adapter> >&solution);

  void set_input_parameters(const Teuchos::ParameterList &p);

  RCP<mj_partBoxVector_t> getGlobalBoxBoundaries() const;

  bool mj_premigrate_to_subset(
    int used_num_ranks,
    int migration_selection_option,
    RCP<const Environment> mj_env_,
    RCP<const Comm<int> > mj_problemComm_,
    int coord_dim_,
    mj_lno_t num_local_coords_,
    mj_gno_t num_global_coords_,  size_t num_global_parts_,
    Kokkos::View<const mj_gno_t*, device_t> & initial_mj_gnos_,
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t> &
      mj_coordinates_,
    int num_weights_per_coord_,
    Kokkos::View<mj_scalar_t**, device_t> & mj_weights_,
    //results
    RCP<const Comm<int> > &result_problemComm_,
    mj_lno_t & result_num_local_coords_,
    Kokkos::View<mj_gno_t*, device_t> & result_initial_mj_gnos_,
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t> &
      result_mj_coordinates_,
    Kokkos::View<mj_scalar_t**, device_t> & result_mj_weights_,
    int * &result_actual_owner_rank_);

public:

  Zoltan2_AlgMJ(const RCP<const Environment> &env,
    RCP<const Comm<int> > &problemComm,
    const RCP<const coordinateModel_t> &coords) :
      mj_partitioner(),
      mj_env(env),
      mj_problemComm(problemComm),
      mj_coords(coords),
      imbalance_tolerance(0),
      num_teams(0),
      num_global_parts(1),
      recursion_depth(0),
      coord_dim(0),
      num_local_coords(0),
      num_global_coords(0),
      num_weights_per_coord(0),
      num_first_level_parts(1),
      distribute_points_on_cut_lines(true),
      max_concurrent_part_calculation(1),
      check_migrate_avoid_migration_option(0),
      migration_type(0),
      minimum_migration_imbalance(0.30),
      mj_keep_part_boxes(false),
      mj_run_as_rcb(false),
      mj_premigration_option(0),
      min_coord_per_rank_for_premigration(32000),
      comXAdj_(),
      comAdj_()
  {
  }

  ~Zoltan2_AlgMJ()
  {
  }

  /*! \brief Set up validators specific to this algorithm
   */
  static void getValidParameters(ParameterList & pl)
  {
    const bool bUnsorted = true; // this clarifies the flag is for unsrorted
    RCP<Zoltan2::IntegerRangeListValidator<int>> mj_parts_Validator =
    Teuchos::rcp( new Zoltan2::IntegerRangeListValidator<int>(bUnsorted) );
    pl.set("mj_parts", "0", "list of parts for multiJagged partitioning "
      "algorithm. As many as the dimension count.", mj_parts_Validator);

    pl.set("mj_concurrent_part_count", 1, "The number of parts whose cut "
      "coordinates will be calculated concurently.",
      Environment::getAnyIntValidator());

    pl.set("mj_minimum_migration_imbalance", 1.1,
      "mj_minimum_migration_imbalance, the minimum imbalance of the "
      "processors to avoid migration",
      Environment::getAnyDoubleValidator());

    RCP<Teuchos::EnhancedNumberValidator<int>> mj_migration_option_validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(0, 2) );
    pl.set("mj_migration_option", 1, "Migration option, 0 for decision "
      "depending on the imbalance, 1 for forcing migration, 2 for "
      "avoiding migration", mj_migration_option_validator);

    RCP<Teuchos::EnhancedNumberValidator<int>> mj_migration_type_validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(0, 1) );
      pl.set("mj_migration_type", 0,
      "Migration type, 0 for migration to minimize the imbalance "
      "1 for migration to minimize messages exchanged the migration.",
      mj_migration_option_validator);

    // bool parameter
    pl.set("mj_keep_part_boxes", false, "Keep the part boundaries of the "
      "geometric partitioning.", Environment::getBoolValidator());

    // bool parameter
    pl.set("mj_enable_rcb", false, "Use MJ as RCB.",
      Environment::getBoolValidator());

    pl.set("mj_recursion_depth", -1, "Recursion depth for MJ: Must be "
      "greater than 0.", Environment::getAnyIntValidator());

    RCP<Teuchos::EnhancedNumberValidator<int>>
      mj_num_teams_validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(
      0, Teuchos::EnhancedNumberTraits<int>::max()) );
    pl.set("mj_num_teams", 0,
      "How many teams for the main kernel loop"
      , mj_num_teams_validator);

    RCP<Teuchos::EnhancedNumberValidator<int>>
      mj_premigration_option_validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(0, 1024) );

    pl.set("mj_premigration_option", 0,
      "Whether to do premigration or not. 0 for no migration "
      "x > 0 for migration to consecutive processors, "
      "the subset will be 0,x,2x,3x,...subset ranks."
      , mj_premigration_option_validator);

    pl.set("mj_premigration_coordinate_count", 32000, "How many coordinate to "
      "assign each rank in multijagged after premigration"
      , Environment::getAnyIntValidator());
  }

  /*! \brief Multi Jagged  coordinate partitioning algorithm.
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

  mj_part_t pointAssign(int dim, adapter_scalar_t *point) const;

  void boxAssign(int dim, adapter_scalar_t *lower, adapter_scalar_t *upper,
    size_t &nPartsFound, mj_part_t **partsFound) const;

  /*! \brief returns communication graph resulting from MJ partitioning.
   */
  void getCommunicationGraph(
    const PartitioningSolution<Adapter> *solution,
    ArrayRCP<mj_part_t> &comXAdj,
    ArrayRCP<mj_part_t> &comAdj);

  void set_up_partitioning_data( // public for CUDA
    const RCP<PartitioningSolution<Adapter> >&solution);

  private:
    std::string timer_base_string; // used for making timers

    // After loading views from coordinate adapter we may need to copy them
    // if mj type is different, but otherwise we just want to assign the view.
    // So purpose of this code is to make that assign only happen when the types
    // match. The empty case would otherwise not compile.
    // If they don't match the internal code handles allocating the new view
    // and copying the elements. See the test Zoltan2_mj_int_coordinates.
    template<class dst_t, class src_t> // version for same types
    typename std::enable_if<std::is_same<typename dst_t::value_type,
      typename src_t::value_type>::value>::type
    assign_if_same(dst_t & dst, const src_t & src) {
      dst = src;
    }
    template<class dst_t, class src_t> // version for different types
    typename std::enable_if<!std::is_same<typename dst_t::value_type,
      typename src_t::value_type>::value>::type
    assign_if_same(dst_t & dst, const src_t & src) {
      // do nothing - handled manually
    }
};

template <typename Adapter>
bool Zoltan2_AlgMJ<Adapter>::mj_premigrate_to_subset(
  int used_num_ranks,
  int migration_selection_option,
  RCP<const Environment> mj_env_,
  RCP<const Comm<int> > mj_problemComm_,
  int coord_dim_,
  mj_lno_t num_local_coords_,
  mj_gno_t num_global_coords_, size_t num_global_parts_,
  Kokkos::View<const mj_gno_t*, device_t> & initial_mj_gnos_,
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t> & mj_coordinates_,
  int num_weights_per_coord_,
  Kokkos::View<mj_scalar_t**, device_t> & mj_weights_,
  //results
  RCP<const Comm<int> > & result_problemComm_,
  mj_lno_t &result_num_local_coords_,
  Kokkos::View<mj_gno_t*, device_t> & result_initial_mj_gnos_,
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t> &
    result_mj_coordinates_,
  Kokkos::View<mj_scalar_t**, device_t> & result_mj_weights_,
  int * &result_actual_owner_rank_)
{
  mj_env_->timerStart(MACRO_TIMERS,
    timer_base_string + "PreMigration DistributorPlanCreating");

  int myRank = mj_problemComm_->getRank();
  int worldSize = mj_problemComm_->getSize();

  mj_part_t groupsize = worldSize / used_num_ranks;

  std::vector<mj_part_t> group_begins(used_num_ranks + 1, 0);

  mj_part_t i_am_sending_to = 0;
  bool am_i_a_receiver = false;

  for(int i = 0; i < used_num_ranks; ++i) {
    group_begins[i+ 1]  = group_begins[i] + groupsize;
    if(worldSize % used_num_ranks > i) group_begins[i+ 1] += 1;
    if(i == used_num_ranks) group_begins[i+ 1] = worldSize;
    if(myRank >= group_begins[i] && myRank < group_begins[i + 1]) {
      i_am_sending_to = group_begins[i];
    }
    if(myRank == group_begins[i]) {
      am_i_a_receiver = true;
    }
  }

  ArrayView<const mj_part_t> idView(&(group_begins[0]), used_num_ranks );
  result_problemComm_ = mj_problemComm_->createSubcommunicator(idView);

  Tpetra::Distributor distributor(mj_problemComm_);

  std::vector<mj_part_t>
    coordinate_destinations(num_local_coords_, i_am_sending_to);

  ArrayView<const mj_part_t>
    destinations(&(coordinate_destinations[0]), num_local_coords_);
  mj_lno_t num_incoming_gnos = distributor.createFromSends(destinations);
  result_num_local_coords_ = num_incoming_gnos;
  mj_env_->timerStop(MACRO_TIMERS,
    timer_base_string + "PreMigration DistributorPlanCreating");

  mj_env_->timerStart(MACRO_TIMERS,
    timer_base_string + "PreMigration DistributorMigration");

  // MPI Buffers should be on Kokkos::HostSpace not Kokkos::CudaUVMSpace

  // migrate gnos.
  {
    ArrayRCP<mj_gno_t> received_gnos(num_incoming_gnos);
    Kokkos::View<mj_gno_t*, Kokkos::HostSpace> host_initial_mj_gnos(
      Kokkos::ViewAllocateWithoutInitializing("host_initial_mj_gnos"),
      initial_mj_gnos_.size()); // initial_mj_gnos_ is const mj_gno_t *
    Kokkos::deep_copy(host_initial_mj_gnos, initial_mj_gnos_);
    ArrayView<const mj_gno_t> sent_gnos(host_initial_mj_gnos.data(),
      num_local_coords_);
    distributor.doPostsAndWaits<mj_gno_t>(sent_gnos, 1, received_gnos());
    result_initial_mj_gnos_ = Kokkos::View<mj_gno_t*, device_t>(
      Kokkos::ViewAllocateWithoutInitializing("result_initial_mj_gnos_"),
      num_incoming_gnos);
    auto host_result_initial_mj_gnos_ = Kokkos::create_mirror_view(
      result_initial_mj_gnos_);
    memcpy(host_result_initial_mj_gnos_.data(),
      received_gnos.getRawPtr(), num_incoming_gnos * sizeof(mj_gno_t));
    Kokkos::deep_copy(result_initial_mj_gnos_, host_result_initial_mj_gnos_);
  }

  // migrate coordinates
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t> dst_coordinates(
    Kokkos::ViewAllocateWithoutInitializing("mj_coordinates"),
    num_incoming_gnos, this->coord_dim);
  auto host_dst_coordinates = Kokkos::create_mirror_view(
    dst_coordinates);
  auto host_src_coordinates =
    Kokkos::create_mirror_view(Kokkos::HostSpace(), this->mj_coordinates);
  Kokkos::deep_copy(host_src_coordinates, this->mj_coordinates);
  for(int i = 0; i < this->coord_dim; ++i) {
    auto sub_host_src_coordinates
      = Kokkos::subview(host_src_coordinates, Kokkos::ALL, i);
    auto sub_host_dst_coordinates
      = Kokkos::subview(host_dst_coordinates, Kokkos::ALL, i);
    // Note Layout Left means we can do these in contiguous blocks
    ArrayView<mj_scalar_t> sent_coord(
      sub_host_src_coordinates.data(), this->num_local_coords);
    ArrayRCP<mj_scalar_t> received_coord(num_incoming_gnos);
    distributor.doPostsAndWaits<mj_scalar_t>(
      sent_coord, 1, received_coord());
    memcpy(sub_host_dst_coordinates.data(),
      received_coord.getRawPtr(), num_incoming_gnos * sizeof(mj_scalar_t));
  }
  deep_copy(dst_coordinates, host_dst_coordinates);
  result_mj_coordinates_ = dst_coordinates;

  // migrate weights.
  Kokkos::View<mj_scalar_t**, device_t> dst_weights(
    Kokkos::ViewAllocateWithoutInitializing("mj_weights"),
    num_incoming_gnos, this->num_weights_per_coord);
  auto host_dst_weights = Kokkos::create_mirror_view(dst_weights);
  auto host_src_weights = Kokkos::create_mirror_view(this->mj_weights);
  Kokkos::deep_copy(host_src_weights, this->mj_weights);
  for(int i = 0; i < this->num_weights_per_coord; ++i) {
    auto sub_host_src_weights
      = Kokkos::subview(host_src_weights, Kokkos::ALL, i);
    auto sub_host_dst_weights
      = Kokkos::subview(host_dst_weights, Kokkos::ALL, i);
    ArrayRCP<mj_scalar_t> sent_weight(this->num_local_coords);

    // Layout Right means these are not contiguous
    // However we don't have any systems setup with more than 1 weight so
    // really I have not tested any of this code with num weights > 1.
    // I think this is the right thing to do. Note that there are other
    // places in the code which don't handle the possibility of more weights.
    // So evaluating all that and adding tests would be another project.
    for(mj_lno_t n = 0; n < this->num_local_coords; ++n) {
      sent_weight[n] = sub_host_src_weights(n);
    }
    ArrayRCP<mj_scalar_t> received_weight(num_incoming_gnos);
    distributor.doPostsAndWaits<mj_scalar_t>(
      sent_weight(), 1, received_weight());

    // Again we copy by index due to layout
    for(mj_lno_t n = 0; n < num_incoming_gnos; ++n) {
      sub_host_dst_weights(n) = received_weight[n];
    }
  }
  Kokkos::deep_copy(dst_weights, host_dst_weights);
  result_mj_weights_ = dst_weights;

  // migrate the owners of the coordinates
  {
    std::vector<int> owner_of_coordinate(num_local_coords_, myRank);
    ArrayView<int> sent_owners(&(owner_of_coordinate[0]), num_local_coords_);
    ArrayRCP<int> received_owners(num_incoming_gnos);
    distributor.doPostsAndWaits<int>(sent_owners, 1, received_owners());
    result_actual_owner_rank_ = new int[num_incoming_gnos];
    memcpy(
	    result_actual_owner_rank_,
	    received_owners.getRawPtr(),
	    num_incoming_gnos * sizeof(int));
  }

  mj_env_->timerStop(MACRO_TIMERS,
    timer_base_string + "PreMigration DistributorMigration");
  return am_i_a_receiver;
}

/*! \brief Multi Jagged  coordinate partitioning algorithm.
 * \param env   library configuration and problem parameters
 * \param problemComm the communicator for the problem
 * \param coords    a CoordinateModel with user data
 * \param solution  a PartitioningSolution, on input it contains part
 * information, on return it also contains the solution and quality metrics.
 */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution)
{
  // purpose of this code is to validate node and UVM status for the tests
  // std::cout << "Memory Space: " << mj_node_t::memory_space::name() << "  "
  //           << "Execution Space: " << mj_node_t::execution_space::name()
  //           << std::endl;

  int execute_counter =
    Zoltan2_AlgMJ_TrackCallsCounter::get_counter_Zoltan2_AlgMJ();
  timer_base_string = "partition(" + std::to_string(execute_counter) + ") - ";

  this->mj_env->timerStart(MACRO_TIMERS, timer_base_string + "all");
  {
    this->mj_env->timerStart(MACRO_TIMERS, timer_base_string + "setup");

    this->set_up_partitioning_data(solution);

    this->set_input_parameters(this->mj_env->getParameters());
    if(this->mj_keep_part_boxes) {
        this->mj_partitioner.set_to_keep_part_boxes();
    }

    this->mj_partitioner.set_partitioning_parameters(
      this->distribute_points_on_cut_lines,
      this->max_concurrent_part_calculation,
      this->check_migrate_avoid_migration_option,
      this->minimum_migration_imbalance, this->migration_type);

    RCP<const Comm<int> > result_problemComm = this->mj_problemComm;
    mj_lno_t result_num_local_coords = this->num_local_coords;
    Kokkos::View<mj_gno_t*, device_t> result_initial_mj_gnos;
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<mj_scalar_t**, Kokkos::LayoutLeft, device_t>
      result_mj_coordinates = this->mj_coordinates;
    Kokkos::View<mj_scalar_t**, device_t> result_mj_weights =
      this->mj_weights;
    int *result_actual_owner_rank = NULL;

    Kokkos::View<const mj_gno_t*, device_t> result_initial_mj_gnos_ =
      this->initial_mj_gnos;

    // TODO: MD 08/2017: Further discussion is required.
    // MueLu calls MJ when it has very few coordinates per processors,
    // such as 10. For example, it begins with 1K processor with 1K coordinate
    // in each. Then with coarsening this reduces to 10 coordinate per procesor.
    // It calls MJ to repartition these to 10 coordinates.
    // MJ runs with 1K processor, 10 coordinate in each, and partitions to
    // 10 parts.  As expected strong scaling is problem here, because
    // computation is almost 0, and communication cost of MJ linearly increases.
    // Premigration option gathers the coordinates to 10 parts before MJ starts
    // therefore MJ will run with a smalller subset of the problem.
    // Below, I am migrating the coordinates if mj_premigration_option is set,
    // and the result parts are less than the current part count, and the
    // average number of local coordinates is less than some threshold.
    // For example, premigration may not help if 1000 processors are
    // partitioning data to 10, but each of them already have 1M coordinate.
    // In that case, we premigration would not help.
    int current_world_size = this->mj_problemComm->getSize();
    mj_lno_t threshold_num_local_coords =
      this->min_coord_per_rank_for_premigration;
    bool is_pre_migrated = false;
    bool am_i_in_subset = true;

    // Note that we need to add testing for migration and should also cover the
    // zoltan case when ZOLTAN2_MJ_ENABLE_ZOLTAN_MIGRATION is defined.
    // Currently did a minimal test of this code by running mjTest with
    // PM=1, TB=0 then run again with C=3 instead of C=4 (numProcs is 4).
    if(mj_premigration_option > 0 &&
        size_t (current_world_size) > this->num_global_parts &&
        this->num_global_coords < mj_gno_t (
        current_world_size * threshold_num_local_coords))
    {
      if(this->mj_keep_part_boxes) {
        throw std::logic_error("Multijagged: mj_keep_part_boxes and "
          "mj_premigration_option are not supported together yet.");
      }

      is_pre_migrated =true;
      int migration_selection_option = mj_premigration_option;
      if(migration_selection_option * this->num_global_parts >
        (size_t) (current_world_size)) {
        migration_selection_option =
          current_world_size / this->num_global_parts;
      }

      int used_num_ranks = int (this->num_global_coords /
        float (threshold_num_local_coords) + 0.5);

      if(used_num_ranks == 0) {
        used_num_ranks = 1;
      }

      am_i_in_subset = this->mj_premigrate_to_subset(
      used_num_ranks,
        migration_selection_option,
        this->mj_env,
        this->mj_problemComm,
        this->coord_dim,
        this->num_local_coords,
        this->num_global_coords,
        this->num_global_parts,
        this->initial_mj_gnos,
        this->mj_coordinates,
        this->num_weights_per_coord,
        this->mj_weights,
        //results
        result_problemComm,
        result_num_local_coords,
        result_initial_mj_gnos,
        result_mj_coordinates,
        result_mj_weights,
        result_actual_owner_rank);

       result_initial_mj_gnos_ = result_initial_mj_gnos;
     }

    Kokkos::View<mj_part_t *, device_t> result_assigned_part_ids;
    Kokkos::View<mj_gno_t*, device_t> result_mj_gnos;

    this->mj_env->timerStop(MACRO_TIMERS, timer_base_string + "setup");

    if(am_i_in_subset) {
      this->mj_partitioner.multi_jagged_part(
        this->mj_env,
        result_problemComm, //this->mj_problemComm,
        this->imbalance_tolerance,
        this->num_teams,
        this->num_global_parts,
        this->part_no_array,
        this->recursion_depth,
        this->coord_dim,
        result_num_local_coords, //this->num_local_coords,
        this->num_global_coords,
        result_initial_mj_gnos_,
        result_mj_coordinates,
        this->num_weights_per_coord,
        this->mj_uniform_weights,
        result_mj_weights,
        this->mj_uniform_parts,
        result_assigned_part_ids,
        result_mj_gnos
      );
    }

    this->mj_env->timerStart(MACRO_TIMERS, timer_base_string + "cleanup");

    // Reorder results so that they match the order of the input
    std::unordered_map<mj_gno_t, mj_lno_t> localGidToLid;
    localGidToLid.reserve(result_num_local_coords);
    Kokkos::View<mj_gno_t*, Kokkos::HostSpace> host_result_initial_mj_gnos(
      Kokkos::ViewAllocateWithoutInitializing("host_result_initial_mj_gnos"),
      result_initial_mj_gnos_.size());
    Kokkos::deep_copy(host_result_initial_mj_gnos, result_initial_mj_gnos_);
    for(mj_lno_t i = 0; i < result_num_local_coords; i++) {
      localGidToLid[host_result_initial_mj_gnos(i)] = i;
    }

    ArrayRCP<mj_part_t> partId = arcp(new mj_part_t[result_num_local_coords],
        0, result_num_local_coords, true);
    auto host_result_assigned_part_ids =
      Kokkos::create_mirror_view(result_assigned_part_ids);
    Kokkos::deep_copy(host_result_assigned_part_ids, result_assigned_part_ids);
    auto host_result_mj_gnos = Kokkos::create_mirror_view(result_mj_gnos);
    Kokkos::deep_copy(host_result_mj_gnos, result_mj_gnos);
    for(mj_lno_t i = 0; i < result_num_local_coords; i++) {
      mj_lno_t origLID = localGidToLid[host_result_mj_gnos(i)];
      partId[origLID] = host_result_assigned_part_ids(i);
    }

    //now the results are reordered. but if premigration occured,
    //then we need to send these ids to actual owners again.
    if(is_pre_migrated) {
      this->mj_env->timerStart(MACRO_TIMERS, timer_base_string +
        "PostMigration DistributorPlanCreating");
      Tpetra::Distributor distributor(this->mj_problemComm);
      ArrayView<const mj_part_t> actual_owner_destinations(
        result_actual_owner_rank , result_num_local_coords);
      mj_lno_t num_incoming_gnos = distributor.createFromSends(
        actual_owner_destinations);
      if(num_incoming_gnos != this->num_local_coords) {
        throw std::logic_error("Zoltan2 - Multijagged Post Migration - "
          "num incoming is not equal to num local coords");
      }

      mj_env->timerStop(MACRO_TIMERS, timer_base_string +
        "PostMigration DistributorPlanCreating");
      mj_env->timerStart(MACRO_TIMERS, timer_base_string +
        "PostMigration DistributorMigration");
      ArrayRCP<mj_gno_t> received_gnos(num_incoming_gnos);
      ArrayRCP<mj_part_t> received_partids(num_incoming_gnos);
      {
        ArrayView<const mj_gno_t> sent_gnos(host_result_initial_mj_gnos.data(),
         result_num_local_coords);
        distributor.doPostsAndWaits<mj_gno_t>(sent_gnos, 1, received_gnos());
      }

      {
        ArrayView<mj_part_t> sent_partnos(partId());
        distributor.doPostsAndWaits<mj_part_t>(sent_partnos, 1,
         received_partids());
      }

      partId = arcp(new mj_part_t[this->num_local_coords],
        0, this->num_local_coords, true);

      {
        std::unordered_map<mj_gno_t, mj_lno_t> localGidToLid2;
        localGidToLid2.reserve(this->num_local_coords);
        auto host_initial_mj_gnos =
          Kokkos::create_mirror_view(this->initial_mj_gnos);
        Kokkos::deep_copy(host_initial_mj_gnos,
          this->initial_mj_gnos);
        for(mj_lno_t i = 0; i < this->num_local_coords; i++) {
          localGidToLid2[host_initial_mj_gnos(i)] = i;
        }

        for(mj_lno_t i = 0; i < this->num_local_coords; i++) {
          mj_lno_t origLID = localGidToLid2[received_gnos[i]];
          partId[origLID] = received_partids[i];
        }
      }

      {
        delete [] result_actual_owner_rank;
      }
      mj_env->timerStop(MACRO_TIMERS,
        timer_base_string + "PostMigration DistributorMigration");
    }
    solution->setParts(partId);
    this->mj_env->timerStop(MACRO_TIMERS, timer_base_string + "cleanup");
  }

  this->mj_env->timerStop(MACRO_TIMERS, timer_base_string + "all");
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
  int criteria_dim = (this->num_weights_per_coord ?
    this->num_weights_per_coord : 1);
  // From the Solution we get part information.
  // If the part sizes for a given criteria are not uniform,
  // then they are values that sum to 1.0.
  this->num_global_parts = solution->getTargetGlobalNumberOfParts();
  // allocate only two dimensional pointer.
  // raw pointer addresess will be obtained from multivector.
  this->mj_uniform_parts = Kokkos::View<bool *, Kokkos::HostSpace>(
    "uniform parts", criteria_dim);
  this->mj_uniform_weights = Kokkos::View<bool *, Kokkos::HostSpace>(
    "uniform weights", criteria_dim);

  Kokkos::View<const mj_gno_t *, device_t> gnos;
  Kokkos::View<adapter_scalar_t **, Kokkos::LayoutLeft, device_t> xyz_adapter;
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<adapter_scalar_t **, device_t> wgts_adapter;
  this->mj_coords->getCoordinatesKokkos(gnos, xyz_adapter, wgts_adapter);
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<mj_scalar_t **, Kokkos::LayoutLeft, device_t> xyz;
  Kokkos::View<mj_scalar_t **, device_t> wgts;

  // Now we must get the data from the adapter.
  // If the types match we point to the view but if not, we must copy.
  if(std::is_same<mj_scalar_t, adapter_scalar_t>()) {
    // we can just point the views but we must specialize because this code
    // only compiles in this case - for is_same false assign does nothing.
    assign_if_same(xyz, xyz_adapter);
    assign_if_same(wgts, wgts_adapter);
  }
  else {
    // we only allocate a new view if we are going to copy
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    xyz = Kokkos::View<mj_scalar_t **, Kokkos::LayoutLeft, device_t>
      (Kokkos::ViewAllocateWithoutInitializing(
      "xyz"), xyz_adapter.extent(0), xyz_adapter.extent(1));
    wgts = Kokkos::View<mj_scalar_t **, device_t>(
      Kokkos::ViewAllocateWithoutInitializing("wgts"),
      wgts_adapter.extent(0), wgts_adapter.extent(1));

    typedef typename Kokkos::View<mj_scalar_t **, device_t>::size_type view_size_t;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
      (0, xyz_adapter.extent(0)), KOKKOS_LAMBDA (int i) {
      for(view_size_t n = 0; n < xyz_adapter.extent(1); ++n) {
        xyz(i, n) = static_cast<mj_scalar_t>(xyz_adapter(i, n));
      }
    });
    Kokkos::parallel_for(
      Kokkos::RangePolicy<typename mj_node_t::execution_space, int>
      (0, wgts.extent(0)), KOKKOS_LAMBDA (int i) {
      for(view_size_t n = 0; n < wgts.extent(1); ++n) {
        wgts(i, n) = static_cast<mj_scalar_t>(wgts_adapter(i, n));
      }
    });
  }

  // obtain global ids.
  this->initial_mj_gnos = gnos;
  // extract coordinates from multivector.
  this->mj_coordinates = xyz;
  // if no weights are provided set uniform weight.

  if(this->num_weights_per_coord == 0) {
    this->mj_uniform_weights(0) = true;
    Kokkos::resize(this->mj_weights, 0, 0);
  }
  else{
    this->mj_weights = wgts;
    for(int wdim = 0; wdim < this->num_weights_per_coord; ++wdim) {
      this->mj_uniform_weights(wdim) = false;
    }
  }

  for(int wdim = 0; wdim < criteria_dim; ++wdim) {
    if(solution->criteriaHasUniformPartSizes(wdim)) {
      this->mj_uniform_parts(wdim) = true;
    }
    else {
      printf("Error: MJ does not support non uniform target part weights\n");
      std::terminate();
    }
  }
}

/* \brief Sets the partitioning parameters for multijagged algorithm.
 * \param pl: is the parameter list provided to zoltan2 call
 * */
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::set_input_parameters(
  const Teuchos::ParameterList &pl)
{
  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("imbalance_tolerance");
  if(pe) {
    double tol;
    tol = pe->getValue(&tol);
    this->imbalance_tolerance = tol - 1.0;
  }

  // TODO: May be a more relaxed tolerance is needed. RCB uses 10%
  if(this->imbalance_tolerance <= 0) {
    this->imbalance_tolerance= 10e-4;
  }

  // if an input partitioning array is provided.
  Kokkos::resize(this->part_no_array, 0);

  // the length of the input partitioning array.
  this->recursion_depth = 0;

  if(pl.getPtr<int>("mj_num_teams")) {
    this->num_teams = pl.get<int>("mj_num_teams");
  }

  if(pl.getPtr<Array <mj_part_t> >("mj_parts")) {
    auto mj_parts = pl.get<Array <mj_part_t> >("mj_parts");
    int mj_parts_size = static_cast<int>(mj_parts.size());

    // build the view we'll have data on and copy values from host
    this->part_no_array = Kokkos::View<mj_part_t*, Kokkos::HostSpace>(
      "part_no_array", mj_parts_size);
    for(int i = 0; i < mj_parts_size; ++i) {
      this->part_no_array(i) = mj_parts.getRawPtr()[i];
    }

    this->recursion_depth = mj_parts_size - 1;
    this->mj_env->debug(2, "mj_parts provided by user");
  }

  // get mj specific parameters.
  this->distribute_points_on_cut_lines = true;
  this->max_concurrent_part_calculation = 1;

  this->mj_run_as_rcb = false;
  this->mj_premigration_option = 0;
	this->min_coord_per_rank_for_premigration = 32000;

  int mj_user_recursion_depth = -1;
  this->mj_keep_part_boxes = false;
  this->check_migrate_avoid_migration_option = 0;
  this->migration_type = 0;
	this->minimum_migration_imbalance = 0.35;

  pe = pl.getEntryPtr("mj_minimum_migration_imbalance");
  if(pe) {
    double imb;
    imb = pe->getValue(&imb);
    this->minimum_migration_imbalance = imb - 1.0;
  }

  pe = pl.getEntryPtr("mj_migration_option");
  if(pe) {
    this->check_migrate_avoid_migration_option =
      pe->getValue(&this->check_migrate_avoid_migration_option);
  } else {
    this->check_migrate_avoid_migration_option = 0;
  }
  if(this->check_migrate_avoid_migration_option > 1) {
    this->check_migrate_avoid_migration_option = -1;
  }

	///
  pe = pl.getEntryPtr("mj_migration_type");
  if(pe) {
    this->migration_type = pe->getValue(&this->migration_type);
  } else {
    this->migration_type = 0;
  }

	//std::cout << " this->migration_type:" <<  this->migration_type << std::endl;
	///

  pe = pl.getEntryPtr("mj_concurrent_part_count");
  if(pe) {
    this->max_concurrent_part_calculation =
      pe->getValue(&this->max_concurrent_part_calculation);
  } else {
    this->max_concurrent_part_calculation = 1; // Set to 1 if not provided.
  }

  pe = pl.getEntryPtr("mj_keep_part_boxes");
  if(pe) {
    this->mj_keep_part_boxes = pe->getValue(&this->mj_keep_part_boxes);
  } else {
    this->mj_keep_part_boxes = false; // Set to invalid value
  }

  // For now, need keep_part_boxes to do pointAssign and boxAssign.
  // pe = pl.getEntryPtr("keep_cuts");
  // if(pe) {
  //      int tmp = pe->getValue(&tmp);
  //      if(tmp) this->mj_keep_part_boxes = true;
  // }

  //need to keep part boxes if mapping type is geometric.
  if(this->mj_keep_part_boxes == false) {
    pe = pl.getEntryPtr("mapping_type");
    if(pe) {
      int mapping_type = -1;
      mapping_type = pe->getValue(&mapping_type);
      if(mapping_type == 0) {
        mj_keep_part_boxes  = true;
      }
    }
  }

  // need to keep part boxes if mapping type is geometric.
  pe = pl.getEntryPtr("mj_enable_rcb");
  if(pe) {
    this->mj_run_as_rcb = pe->getValue(&this->mj_run_as_rcb);
  } else {
    this->mj_run_as_rcb = false; // Set to invalid value
  }

  pe = pl.getEntryPtr("mj_premigration_option");
  if(pe) {
    mj_premigration_option = pe->getValue(&mj_premigration_option);
  } else {
     mj_premigration_option = 0;
  }

  pe = pl.getEntryPtr("mj_premigration_coordinate_count");
  if(pe) {
    min_coord_per_rank_for_premigration = pe->getValue(&mj_premigration_option);
  } else {
    min_coord_per_rank_for_premigration = 32000;
  }

  pe = pl.getEntryPtr("mj_recursion_depth");
  if(pe) {
    mj_user_recursion_depth = pe->getValue(&mj_user_recursion_depth);
  } else {
    mj_user_recursion_depth = -1; // Set to invalid value
  }

  bool val = false;
  pe = pl.getEntryPtr("rectilinear");
  if(pe) {
    val = pe->getValue(&val);
  }
  if(val) {
    this->distribute_points_on_cut_lines = false;
  } else {
    this->distribute_points_on_cut_lines = true;
  }

  if(this->mj_run_as_rcb) {
    mj_user_recursion_depth =
      (int)(ceil(log ((this->num_global_parts)) / log (2.0)));
  }
  if(this->recursion_depth < 1) {
    if(mj_user_recursion_depth > 0) {
      this->recursion_depth = mj_user_recursion_depth;
    }
    else {
      this->recursion_depth = this->coord_dim;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void Zoltan2_AlgMJ<Adapter>::boxAssign(
  int dim,
  adapter_scalar_t *lower,
  adapter_scalar_t *upper,
  size_t &nPartsFound,
  typename Adapter::part_t **partsFound) const
{
  // TODO:  Implement with cuts rather than boxes to reduce algorithmic
  // TODO:  complexity.  Or at least do a search through the boxes, using
  // TODO:  p x q x r x ... if possible.

  nPartsFound = 0;
  *partsFound = NULL;

  if(this->mj_keep_part_boxes) {

    // Get vector of part boxes
    RCP<mj_partBoxVector_t> partBoxes = this->getGlobalBoxBoundaries();

    size_t nBoxes = (*partBoxes).size();
    if(nBoxes == 0) {
      throw std::logic_error("no part boxes exist");
    }

    // Determine whether the box overlaps the globalBox at all
    RCP<mj_partBox_t> globalBox = this->mj_partitioner.get_global_box();

    if(globalBox->boxesOverlap(dim, lower, upper)) {

      std::vector<typename Adapter::part_t> partlist;

      // box overlaps the global box; find specific overlapping boxes
      for(size_t i = 0; i < nBoxes; i++) {
        try {
          if((*partBoxes)[i].boxesOverlap(dim, lower, upper)) {
            nPartsFound++;
            partlist.push_back((*partBoxes)[i].getpId());
            /*
            std::cout << "Given box (";
            for(int j = 0; j < dim; j++)
              std::cout << lower[j] << " ";
            std::cout << ") x (";
            for(int j = 0; j < dim; j++)
              std::cout << upper[j] << " ";
            std::cout << ") overlaps PartBox "
                      << (*partBoxes)[i].getpId() << " (";
            for(int j = 0; j < dim; j++)
              std::cout << (*partBoxes)[i].getlmins()[j] << " ";
            std::cout << ") x (";
            for(int j = 0; j < dim; j++)
              std::cout << (*partBoxes)[i].getlmaxs()[j] << " ";
            std::cout << ")" << std::endl;
            */
          }
        }
        Z2_FORWARD_EXCEPTIONS;
      }
      if(nPartsFound) {
        *partsFound = new mj_part_t[nPartsFound];
        for(size_t i = 0; i < nPartsFound; i++)
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
  adapter_scalar_t *point) const
{
  // TODO:  Implement with cuts rather than boxes to reduce algorithmic
  // TODO:  complexity.  Or at least do a search through the boxes, using
  // TODO:  p x q x r x ... if possible.

  if(this->mj_keep_part_boxes) {
    typename Adapter::part_t foundPart = -1;

    // Get vector of part boxes
    RCP<mj_partBoxVector_t> partBoxes = this->getGlobalBoxBoundaries();

    size_t nBoxes = (*partBoxes).size();
    if(nBoxes == 0) {
      throw std::logic_error("no part boxes exist");
    }

    // Determine whether the point is within the global domain
    RCP<mj_partBox_t> globalBox = this->mj_partitioner.get_global_box();

    if(globalBox->pointInBox(dim, point)) {

      // point is in the global domain; determine in which part it is.
      size_t i;
      for(i = 0; i < nBoxes; i++) {
        try {
          if((*partBoxes)[i].pointInBox(dim, point)) {
            foundPart = (*partBoxes)[i].getpId();
            // std::cout << "Point (";
            // for(int j = 0; j < dim; j++) std::cout << point[j] << " ";
            //   std::cout << ") found in box " << i << " part " << foundPart
            //     << std::endl;
            // (*partBoxes)[i].print();
            break;
          }
        }
        Z2_FORWARD_EXCEPTIONS;
      }

      if(i == nBoxes) {
        // This error should never occur
        std::ostringstream oss;
        oss << "Point (";
        for(int j = 0; j < dim; j++) oss << point[j] << " ";
        oss << ") not found in domain";
        throw std::logic_error(oss.str());
      }
    }

    else {
      // Point is outside the global domain.
      // Determine to which part it is closest.
      // TODO:  with cuts, would not need this special case

      typedef typename Zoltan2::coordinateModelPartBox::coord_t coord_t;
      size_t closestBox = 0;
      coord_t minDistance = std::numeric_limits<coord_t>::max();
      coord_t *centroid = new coord_t[dim];
      for(size_t i = 0; i < nBoxes; i++) {
        (*partBoxes)[i].computeCentroid(centroid);
        coord_t sum = 0.;
        coord_t diff;
        for(int j = 0; j < dim; j++) {
          diff = centroid[j] - point[j];
          sum += diff * diff;
        }
        if(sum < minDistance) {
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
  if(comXAdj_.getRawPtr() == NULL && comAdj_.getRawPtr() == NULL) {
    RCP<mj_partBoxVector_t> pBoxes = this->getGlobalBoxBoundaries();
    mj_part_t ntasks =  (*pBoxes).size();
    int dim = (*pBoxes)[0].getDim();
    GridHash grid(pBoxes, ntasks, dim);
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
} // namespace Zoltan2

#endif
