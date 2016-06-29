
#ifndef _KOKKOSSPGEMMIMPL_HPP
#define _KOKKOSSPGEMMIMPL_HPP
//#define HASHTRACK

//#define TRACK_INSERTS
#define GPU_EXPERIMENTAL
//#define NUMERIC_USE_STATICMEM
//#define twostep
#include <KokkosKernels_Utils.hpp>
#include <KokkosKernels_HashmapAccumulator.hpp>
#include "KokkosKernels_Uniform_Initialized_MemoryPool.hpp"
#include "KokkosKernels_GraphColor.hpp"
#include <fstream>
#include <sstream>
#include <string>
//#include <bitset>

#define KOKKOSKERNELS_SPGEMM_BLOCKSIZE 256
#define KOKKOSKERNELS_SPGEMM_HASHSIZE 16
#define KOKKOSKERNELS_SPGEMM_ANDHASHFUNCTION 15
#define KOKKOSKERNELS_SPGEMM_SHMEMSIZE 16128//12032//16128//16384

namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{

template <typename HandleType,
  typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
  typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_
  >
class KokkosSPGEMM{
public:

  typedef a_row_view_t_ a_row_view_t;
  typedef a_lno_nnz_view_t_ a_in_lno_nnz_view_t;
  typedef a_scalar_nnz_view_t_ a_in_scalar_nnz_view_t;

  typedef b_lno_row_view_t_ b_in_lno_row_view_t;
  typedef b_lno_nnz_view_t_ b_in_lno_nnz_view_t;
  typedef b_scalar_nnz_view_t_ b_in_scalar_nnz_view_t;



  typedef typename a_row_view_t::non_const_value_type size_type;
  typedef typename a_row_view_t::const_value_type const_size_type;


  typedef typename a_in_lno_nnz_view_t::non_const_value_type nnz_lno_t;
  typedef typename a_in_lno_nnz_view_t::const_value_type const_nnz_lno_t;

  typedef typename a_in_scalar_nnz_view_t::non_const_value_type scalar_t;
  typedef typename a_in_scalar_nnz_view_t::const_value_type const_scalar_t;


  typedef typename a_row_view_t::const_type const_a_lno_row_view_t;
  typedef typename a_row_view_t::non_const_type non_const_a_lno_row_view_t;

  typedef typename a_in_lno_nnz_view_t::const_type const_a_lno_nnz_view_t;
  typedef typename a_in_lno_nnz_view_t::non_const_type non_const_a_lno_nnz_view_t;

  typedef typename a_in_scalar_nnz_view_t::const_type const_a_scalar_nnz_view_t;
  typedef typename a_in_scalar_nnz_view_t::non_const_type non_const_a_scalar_nnz_view_t;


  typedef typename b_in_lno_row_view_t::const_type const_b_lno_row_view_t;
  typedef typename b_in_lno_row_view_t::non_const_type non_const_b_lno_row_view_t;

  typedef typename b_in_lno_nnz_view_t::const_type const_b_lno_nnz_view_t;
  typedef typename b_in_lno_nnz_view_t::non_const_type non_const_b_lno_nnz_view_t;

  typedef typename b_in_scalar_nnz_view_t::const_type const_b_scalar_nnz_view_t;
  typedef typename b_in_scalar_nnz_view_t::non_const_type non_const_b_scalar_nnz_view_t;

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t; //Host view type


  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_host_view_t nnz_lno_persistent_work_host_view_t; //Host view type


  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;


  typedef typename HandleType::bool_persistent_view_t bool_persistent_view_t;
  typedef typename HandleType::bool_temp_view_t bool_temp_view_t;


  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t ;
  typedef typename team_policy_t::member_type team_member_t ;

  struct CountTag{};
  struct CountandFillTag{};

  struct SortedCountTag{};
  struct SortedFillTag{};
  struct SortedFillExtraTag{};

  struct FillTag{};
  struct FillExtraTag{};

  struct MultiCoreTag{};
  struct MultiCoreTag2{};
  struct MultiCoreTag3{};

  struct GPUTag{};
  struct GPUTag2{};
  struct GPUTag3{};

  struct MultiCoreTagCount{};
  struct MultiCoreMultiHashTag{};
  struct MultiCoreDenseAccumulatorTag{};


  struct Numeric1Tag{};
  struct Numeric2Tag{};
  struct Numeric3Tag{};

  //typedef Kokkos::RangePolicy<MultiCoreTag, MyExecSpace> my_count_multicore_space;
  typedef Kokkos::TeamPolicy<MultiCoreTagCount, MyExecSpace> multicore_team_count_policy_t ;
  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag, MyExecSpace> multicore_dense_team_count_policy_t ;

  typedef Kokkos::TeamPolicy<MultiCoreTag, MyExecSpace> multicore_team_policy_t ;
  typedef Kokkos::TeamPolicy<MultiCoreTag2, MyExecSpace> multicore2_team_policy_t ;
  typedef Kokkos::TeamPolicy<MultiCoreTag3, MyExecSpace> multicore3_team_policy_t ;
  typedef Kokkos::TeamPolicy<GPUTag, MyExecSpace> gpu_team_policy_t ;
  typedef Kokkos::TeamPolicy<GPUTag2, MyExecSpace> gpu2_team_policy_t ;
  typedef Kokkos::TeamPolicy<GPUTag3, MyExecSpace> gpu3_team_policy_t ;


  typedef Kokkos::RangePolicy<CountTag, MyExecSpace> my_count_exec_space;
  typedef Kokkos::RangePolicy<FillTag, MyExecSpace> my_fill_exec_space;
  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace> team_count_policy_t ;
  typedef Kokkos::TeamPolicy<CountandFillTag, MyExecSpace> team_countandfill_policy_t ;

  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace> team_fill_policy_t ;
  typedef Kokkos::TeamPolicy<FillExtraTag, MyExecSpace> team_fillextra_policy_t ;
  typedef Kokkos::RangePolicy<FillExtraTag, MyExecSpace> my_fillextra_exec_space;


  typedef Kokkos::TeamPolicy<Numeric1Tag, MyExecSpace> team_numeric1_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric2Tag, MyExecSpace> team_numeric2_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric3Tag, MyExecSpace> team_numeric3_policy_t ;

private:
  HandleType *handle;
  nnz_lno_t m;
  nnz_lno_t n;
  nnz_lno_t k;


  const_a_lno_row_view_t row_mapA;
  const_a_lno_nnz_view_t entriesA;
  const_a_scalar_nnz_view_t valsA;
  bool transposeA;

  const_b_lno_row_view_t row_mapB;
  const_b_lno_nnz_view_t entriesB;
  const_b_scalar_nnz_view_t valsB;
  bool transposeB;

  const int verbose = 1;
public:
  KokkosSPGEMM(
      HandleType *handle_,
      nnz_lno_t m_,
      nnz_lno_t n_,
      nnz_lno_t k_,
      const_a_lno_row_view_t row_mapA_,
      const_a_lno_nnz_view_t entriesA_,
      bool transposeA_,
      const_b_lno_row_view_t row_mapB_,
      const_b_lno_nnz_view_t entriesB_,
      bool transposeB_):handle (handle_), m(m_), n(n_), k(k_),
          row_mapA(row_mapA_), entriesA(entriesA_), valsA(), transposeA(transposeA_),
          row_mapB(row_mapB_), entriesB(entriesB_), valsB(), transposeB(transposeB_)
          //,row_mapC(), entriesC(), valsC()
          {}

  KokkosSPGEMM(
      HandleType *handle_,
      nnz_lno_t m_,
      nnz_lno_t n_,
      nnz_lno_t k_,
        const_a_lno_row_view_t row_mapA_,
        const_a_lno_nnz_view_t entriesA_,
        const_a_scalar_nnz_view_t valsA_,
        bool transposeA_,
        const_b_lno_row_view_t row_mapB_,
        const_b_lno_nnz_view_t entriesB_,
        const_b_scalar_nnz_view_t valsB_,
        bool transposeB_):handle (handle_), m(m_), n(n_), k(k_),
            row_mapA(row_mapA_), entriesA(entriesA_), valsA(valsA_), transposeA(transposeA_),
            row_mapB(row_mapB_), entriesB(entriesB_), valsB(valsB_), transposeB(transposeB_)
            //,row_mapB(), entriesC(), valsC()
            {}

  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename c_row_view_t>
  struct PredicMaxRowNNZ{

    //typedef typename a_row_view_t::non_const_value_type size_type;
    //typedef typename a_nnz_view_t::non_const_value_type nnz_lno_t;


    nnz_lno_t m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;
    b_row_view_t row_mapB;
    c_row_view_t rough_row_mapC;

    const size_type min_val;

    PredicMaxRowNNZ(
        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        c_row_view_t rough_row_mapC_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_mapB(row_mapB_), rough_row_mapC(rough_row_mapC_),
          min_val(KOKKOSKERNELS_MACRO_MIN(-std::numeric_limits<size_type>::max(), 0)) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember, size_type &overal_max) const {
      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //check ii is out of range. if it is, just return.
      if (row_index >= m) return;
      const size_type col_begin = row_mapA[row_index];
      const size_type col_end = row_mapA[row_index + 1];
      const nnz_lno_t left_work = col_end - col_begin;

      size_type max_num_results_in_row = 0;

      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, left_work),
          [&] (nnz_lno_t i, size_type & valueToUpdate) {
        const size_type adjind = i + col_begin;
        const nnz_lno_t colIndex = entriesA[adjind];
        valueToUpdate += row_mapB [colIndex + 1] - row_mapB[colIndex];
      },
      max_num_results_in_row);
      if (overal_max < max_num_results_in_row) { overal_max = max_num_results_in_row;}
    }

    KOKKOS_INLINE_FUNCTION
    void join (volatile size_type& dst,const volatile size_type& src) const {
      if (dst < src) { dst = src;}
    }


    KOKKOS_INLINE_FUNCTION
    void init (size_type& dst) const
    {
      dst = min_val;
    }
  };

  template <typename row_view_t, typename new_row_view_t, typename new_nnz_view_t>
  struct copyMatrix{

    const nnz_lno_t numrows;

    row_view_t row_map;
    new_nnz_view_t set_index_entries;
    new_nnz_view_t set_entries;


    new_row_view_t new_row_map;
    new_nnz_view_t out_set_index_entries;
    new_nnz_view_t out_set_entries;

    copyMatrix(
        row_view_t row_map_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_,

        new_row_view_t new_row_map_,
        new_row_view_t out_set_index_entries_,
        new_row_view_t out_set_entries_
        ):
      row_map(row_map_),

      new_row_map(new_row_map_),
      set_index_entries(set_index_entries_),
      set_entries(set_entries_),

      out_set_index_entries(out_set_index_entries_),
      out_set_entries(out_set_entries_),
      numrows(row_map_.dimension_0() - 1)
      {}


    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {

      nnz_lno_t i = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (i >= numrows) return;


      size_type rowBegin = new_row_map(i);
      nnz_lno_t left_work = new_row_map(i + 1) - rowBegin;
      size_type oldRowBegin= row_map(i);

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, left_work),
          [&] (nnz_lno_t i) {
        const size_type adjind = i + rowBegin;
        const size_type oldadjind = i + oldRowBegin;

        out_set_index_entries(adjind) = set_index_entries(oldadjind);
        out_set_entries (adjind) =  set_entries(oldadjind);
      });
    }
  };

  template <typename row_view_t, typename nnz_view_t, typename new_row_view_t, typename new_nnz_view_t>
  struct SingleStepZipMatrixSingleVector{

    row_view_t row_map;
    nnz_view_t entries;
    const nnz_lno_t compression_bit_mask;
    const int compression_bit_divide_shift;
    int vector_size;

    new_row_view_t new_row_map;
    //new_row_view_t out_set_index_entries;
    //new_nnz_view_t out_set_entries;


    new_nnz_view_t set_index_begins;
    new_nnz_view_t set_index_nexts;

    new_nnz_view_t set_index_entries;
    new_nnz_view_t set_entries;

    nnz_lno_t *pset_index_begins;
    nnz_lno_t *pset_index_nexts;

    nnz_lno_t * pset_index_entries;
    nnz_lno_t * pset_entries;
    const nnz_lno_t numrows;

    SingleStepZipMatrixSingleVector(
        row_view_t row_map_,
        nnz_view_t entries_,
        nnz_lno_t compression_bit_mask_,
        int compression_bit_divide_shift_,
        int vector_size_,
        new_row_view_t new_row_map_,
        /*Global memory hash space. in the size of nnz*/
        new_nnz_view_t set_index_begins_,
        new_nnz_view_t set_index_nexts_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_
        ):
      row_map(row_map_),
      entries(entries_),
      compression_bit_mask(compression_bit_mask_),
      compression_bit_divide_shift(compression_bit_divide_shift_),
      vector_size(vector_size_),
      new_row_map(new_row_map_),
      //out_set_index_entries(), out_set_entries(),

      /*Global memory hashes*/
      set_index_begins(set_index_begins_),
      set_index_nexts(set_index_nexts_),
      set_index_entries(set_index_entries_),
      set_entries(set_entries_),
      pset_index_begins(set_index_begins_.ptr_on_device()),
      pset_index_nexts(set_index_nexts_.ptr_on_device()),
      pset_index_entries(set_index_entries_.ptr_on_device()),
      pset_entries(set_entries_.ptr_on_device()),
      numrows(row_map_.dimension_0() - 1)
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {
      nnz_lno_t row_ind = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //std::cout << "i:" << i << std::endl;
      if (row_ind >= numrows) return;
      size_type rowBegin = row_map(row_ind);

      nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

      //same as second level hash map.
      //second level hashmap.
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm2(left_work, left_work, pset_index_begins + rowBegin, pset_index_nexts+ rowBegin, pset_index_entries+ rowBegin, pset_entries+ rowBegin);

      nnz_lno_t used_size = 0;
      for (nnz_lno_t i = 0; i < left_work; ++i){
        nnz_lno_t n_set = 1;
        const size_type adjind = i + rowBegin;
        const nnz_lno_t n = entries[adjind];
        nnz_lno_t n_set_index = n >> compression_bit_divide_shift;
        n_set = n_set << (n & compression_bit_mask);
        nnz_lno_t hash = n_set_index % left_work;

        int insertion = hm2.sequential_insert_into_hash_mergeOr(
            hash, n_set_index,n_set, &used_size, hm2.max_value_size);
      }
      new_row_map(row_ind) = used_size;
    }
  };

  template <typename row_view_t, typename nnz_view_t, typename new_row_view_t, typename new_nnz_view_t>
  struct SingleStepZipMatrix{


    const nnz_lno_t numrows;
    row_view_t row_map;
    nnz_view_t entries;
    const nnz_lno_t compression_bit_mask;
    const int compression_bit_divide_shift;
    int vector_size;

    new_row_view_t new_row_map;
    //new_row_view_t out_set_index_entries;
    //new_nnz_view_t out_set_entries;


    new_nnz_view_t set_index_begins;
    new_nnz_view_t set_index_nexts;

    new_nnz_view_t set_index_entries;
    new_nnz_view_t set_entries;

    nnz_lno_t *pset_index_begins;
    nnz_lno_t *pset_index_nexts;

    nnz_lno_t * pset_index_entries;
    nnz_lno_t * pset_entries;


    //pool_memory_space m_space;
    const int shared_memory_size;

    SingleStepZipMatrix(
        row_view_t row_map_,
        nnz_view_t entries_,
        nnz_lno_t compression_bit_mask_,
        int compression_bit_divide_shift_,
        int vector_size_,
        new_row_view_t new_row_map_,
        /*Global memory hash space. in the size of nnz*/
        new_nnz_view_t set_index_begins_,
        new_nnz_view_t set_index_nexts_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_
        ):
      row_map(row_map_),
      entries(entries_),
      compression_bit_mask(compression_bit_mask_),
      compression_bit_divide_shift(compression_bit_divide_shift_),
      vector_size(vector_size_),
      new_row_map(new_row_map_),
      //out_set_index_entries(), out_set_entries(),

      /*Global memory hashes*/
      set_index_begins(set_index_begins_),
      set_index_nexts(set_index_nexts_),
      set_index_entries(set_index_entries_),
      set_entries(set_entries_),
      pset_index_begins(set_index_begins_.ptr_on_device()),
      pset_index_nexts(set_index_nexts_.ptr_on_device()),
      pset_index_entries(set_index_entries_.ptr_on_device()),
      pset_entries(set_entries_.ptr_on_device()),
      numrows(row_map_.dimension_0() - 1),
      shared_memory_size(KOKKOSKERNELS_SPGEMM_SHMEMSIZE)
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {

      nnz_lno_t row_ind = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //std::cout << "i:" << i << std::endl;
      if (row_ind >= numrows) return;

      //std::cout << "i:" << i << std::endl;
      //how much shared memory a thread will have in team
      int thread_memory = shared_memory_size / teamMember.team_size();
      //allocate all shared memory
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      //allocate memory in the size of vectors
      nnz_lno_t *result_keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * vector_size;
      nnz_lno_t *result_vals = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * vector_size;

      thread_memory -= vector_size * sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;

      //calculate the memory needed for 1 single hash entry.
      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2; //begins, nexts, and keys. No need for vals yet.
      //how many hash entries can be held in shared memory.
      nnz_lno_t shared_memory_hash_size = thread_memory / unit_memory;


      //points to the beginning of hashes
      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      nnz_lno_t * vals = (nnz_lno_t *) (all_shared_memory);;


      //this is a hash for individual row elements. therefore, the size can be at most
      //the number of elements in a row, so the nnz_lno_t can be used instead of size_type here.

      //first level hashmap
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      size_type rowBegin = row_map(row_ind);
      size_type rowBeginP = rowBegin;


      nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

      //same as second level hash map.
      //second level hashmap.
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm2(left_work, left_work, pset_index_begins + rowBegin, pset_index_nexts+ rowBegin, pset_index_entries+ rowBegin, pset_entries+ rowBegin);



      //initialize begins.
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
          [&] (int i) {
        begins[i] = -1;
      });
      //initialize begins.
      //these are initialized before the loop.
      /*
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, left_work),
          [&] (nnz_lno_t i) {
        pset_index_begins[rowBegin + i] = -1;
      });
      */

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
      });

      //lno_t neighbor_set_count = 0;

      while (left_work){

        //std::cout << "left_work:" << left_work << std::endl;
        nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

        //first get the portion of the work for the vector lane.
        nnz_lno_t n_set_index = -1;
        nnz_lno_t n_set = 1;

        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, work_to_handle),
            [&] (nnz_lno_t i) {
          const size_type adjind = i + rowBegin;
          const nnz_lno_t n = entries[adjind];
          n_set_index = n >> compression_bit_divide_shift;
          n_set = n_set << (n & compression_bit_mask);
        });



        //it is possible that multiple threads have same values.
        //first merge them, as a result of this operation we will have the n_sets merged,
        //if a thread's value merged to some other threads we have n_set = -1.
        hm.vector_mergeOr_MEM(teamMember, vector_size, n_set_index,n_set, result_keys, result_vals);



        nnz_lno_t hash = n_set_index % shared_memory_hash_size;
        if (n_set_index == -1) hash = -1;
        int overall_num_unsuccess = 0;

        int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(
                                          teamMember, vector_size, hash,n_set_index, n_set, used_hash_sizes, shared_memory_hash_size);

        Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &overall_num_unsuccess) {
          overall_num_unsuccess += num_unsuccess;
        }, overall_num_unsuccess);


        //if one of the inserts was successfull, which means we run out shared memory
        if (overall_num_unsuccess){
          nnz_lno_t hash = -1;
          if (num_unsuccess) hash = n_set_index % hm2.hash_key_size;

          int insertion = hm2.vector_atomic_insert_into_hash_mergeOr(
              teamMember, vector_size, hash,n_set_index,n_set, used_hash_sizes + 1, hm2.max_value_size);
        }

        left_work -= work_to_handle;
        rowBegin += work_to_handle;
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
        new_row_map(row_ind) = used_hash_sizes[0] + used_hash_sizes[1];
      });


      size_type written_index = used_hash_sizes[1];
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, used_hash_sizes[0]),
          [&] (nnz_lno_t i) {
        pset_index_entries[rowBeginP + written_index + i] = keys[i];
        pset_entries[rowBeginP + written_index + i] = vals[i];
      });
    }

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }

  };

  template <typename row_view_t, typename nnz_view_t, typename new_row_view_t, typename new_nnz_view_t>
  struct zipMatrix{
    typedef typename row_view_t::non_const_value_type index_t;
    typedef typename nnz_view_t::non_const_value_type lno_t;
    row_view_t row_map;
    nnz_view_t entries;
    const lno_t compression_bit_mask;
    const int compression_bit_divide_shift;

    new_row_view_t new_row_map;
    new_nnz_view_t set_index_entries;
    new_nnz_view_t set_entries;
    new_nnz_view_t set_entries_begin;
    new_nnz_view_t set_entries_next;
    const lno_t numrows;

    zipMatrix(
        row_view_t row_map_,nnz_view_t entries_,
        lno_t compression_bit_mask_, int compression_bit_divide_shift_,
        new_row_view_t new_row_map_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_ ):
      row_map(row_map_),
      entries(entries_),
      compression_bit_mask(compression_bit_mask_),
      compression_bit_divide_shift(compression_bit_divide_shift_),
      new_row_map(new_row_map_),
      set_index_entries(set_index_entries_),
      set_entries(set_entries_),
      numrows(row_map_.dimension_0() - 1){}

    zipMatrix(
        row_view_t row_map_,nnz_view_t entries_,
        lno_t compression_bit_mask_, int compression_bit_divide_shift_,
        new_row_view_t new_row_map_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_ ,
        new_nnz_view_t set_entries_begin_,
        new_nnz_view_t set_entries_next_):
      row_map(row_map_),
      entries(entries_),
      compression_bit_mask(compression_bit_mask_),
      compression_bit_divide_shift(compression_bit_divide_shift_),
      new_row_map(new_row_map_),
      set_index_entries(set_index_entries_),
      set_entries(set_entries_),
      set_entries_begin(set_entries_begin_),
      set_entries_next(set_entries_next_),
      numrows(row_map_.dimension_0() - 1){}


    KOKKOS_INLINE_FUNCTION
    void operator()(const FillExtraTag&, const size_t & i) const {
      index_t newrowBegin = new_row_map(i);
      const index_t rowBegin = row_map(i);
      const index_t rowEnd = row_map(i + 1);
      index_t ii = rowBegin;
      //lno_t col_set_index = 0;
      lno_t n = entries[ii];

      lno_t prev_n_set_index = n >> compression_bit_divide_shift;
      lno_t prev_n_set = 1; lno_t n_mask = n & compression_bit_mask;  prev_n_set = prev_n_set << n_mask;

      //lno_t prev_set_begin_index = ii;
      //set_entries_begin(prev_n_set_index) = ii++;

      set_entries_begin(newrowBegin) = ii++;

      for (; ii < rowEnd; ++ii){
        lno_t n = entries[ii];
        lno_t n_set_index = n >> compression_bit_divide_shift;

        //set_entries_next(ii - 1) = ii;
        if (n_set_index != prev_n_set_index){
          set_entries(newrowBegin) = prev_n_set;
          set_index_entries(newrowBegin++) = prev_n_set_index;
          prev_n_set_index = n_set_index;
          prev_n_set = 0;
          set_entries_begin(newrowBegin) = ii;
          //set_entries_next(ii - 1) = -1;
        }

        lno_t set_mark = 1; lno_t n_mask = n & compression_bit_mask;  set_mark = set_mark << n_mask;
        prev_n_set = prev_n_set | set_mark;
      }
      //set_entries_next(ii - 1) = -1;
      set_entries(newrowBegin) = prev_n_set;
      set_index_entries(newrowBegin) = prev_n_set_index;
      set_entries_begin(newrowBegin + 1) = ii;
    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const FillTag&, const size_t & i) const {
      index_t newrowBegin = new_row_map(i);
      const index_t rowBegin = row_map(i);
      const index_t rowEnd = row_map(i + 1);
      index_t ii = rowBegin;
      //lno_t col_set_index = 0;
      lno_t n = entries[ii++];
      lno_t prev_n_set_index = n >> compression_bit_divide_shift;
      lno_t prev_n_set = 1; lno_t n_mask = n & compression_bit_mask;  prev_n_set = prev_n_set << n_mask;
      for (; ii < rowEnd; ++ii){
        lno_t n = entries[ii];
        lno_t n_set_index = n >> compression_bit_divide_shift;
        if (n_set_index != prev_n_set_index){
          set_entries(newrowBegin) = prev_n_set;
          set_index_entries(newrowBegin++) = prev_n_set_index;
          prev_n_set_index = n_set_index;
          prev_n_set = 0;
        }
        lno_t set_mark = 1; lno_t n_mask = n & compression_bit_mask;  set_mark = set_mark << n_mask;
        prev_n_set = prev_n_set | set_mark;
      }
      set_entries(newrowBegin) = prev_n_set;
      set_index_entries(newrowBegin) = prev_n_set_index;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const lno_t & i, lno_t & overall_size) const {
      const index_t rowBegin = row_map(i);
      const index_t rowEnd = row_map(i + 1);
      lno_t prev_n_set = -1;
      lno_t neighbor_set_count = 0;
      for (index_t ii = rowBegin; ii < rowEnd; ++ii){
        lno_t n = entries(ii);
        lno_t n_set = n >> compression_bit_divide_shift;
        neighbor_set_count += (prev_n_set != n_set);
        prev_n_set = n_set;
      }
      new_row_map(i) = neighbor_set_count;
      overall_size += neighbor_set_count;
    }

    size_t team_shmem_size (int team_size) const {
      return 16384;
    }
  };



  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t,typename c_row_view_t>
  size_t getMaxRoughRowNNZ(
      nnz_lno_t m,
      a_row_view_t row_mapA,
      a_nnz_view_t entriesA,
      b_row_view_t row_mapB,
      c_row_view_t row_mapC){
    int teamSizeMax = 0;
    int vector_size = 0;

    PredicMaxRowNNZ<a_row_view_t, a_nnz_view_t, b_row_view_t, c_row_view_t> pcnnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        row_mapC);

    int max_allowed_team_size = team_policy_t::team_size_max(pcnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<size_type, MyExecSpace>(
        max_allowed_team_size, vector_size, teamSizeMax, m, entriesA.dimension_0());

    typename c_row_view_t::non_const_value_type rough_size = 0;
    Kokkos::parallel_reduce( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), pcnnnz, rough_size);
    MyExecSpace::fence();
    return rough_size;
  }

#if 0
  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
            typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
            typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t>
  struct calculateC{
    size_type m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;
    a_scalar_view_t valuesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesB;
    b_scalar_view_t valuesB;


    row_lno_temp_work_view_t compressed_b_rowmap;
    row_lno_temp_work_view_t compressed_b_set_indices;
    row_lno_temp_work_view_t compressed_b_sets;
    row_lno_temp_work_view_t compressed_b_set_begins;
    row_lno_temp_work_view_t compressed_b_set_nexts;

    c_row_view_t row_mapC;
    c_nnz_view_t entriesC;
    c_scalar_view_t valuesC;

    typedef typename a_row_view_t::non_const_value_type c_row_lno_t;
    typedef typename b_nnz_view_t::non_const_value_type c_nnz_lno_t;
    typedef typename c_scalar_view_t::non_const_value_type c_scalar_t;

    //const int KOKKOSKERNELS_SPGEMM_HASHSIZE;
    //const int KOKKOSKERNELS_SPGEMM_BLOCKSIZE;

    const c_nnz_lno_t hash_size;
    const c_nnz_lno_t pow2_hash_func;
    const c_nnz_lno_t maxUncompressedNonZero;
    const c_nnz_lno_t maxCompressedNonZero;
    //const c_nnz_lno_t MaxFlops;

    const size_t requiredMem;
    const size_t multicore_chunk_size;

#ifdef TRACK_INSERTS
    //c_row_view_t overal_inserts, overall_hash_op, wasted_ops;
    Kokkos::View<long, Kokkos::OpenMP> overal_inserts, overall_hash_op, wasted_ops, used_hashes;
    //c_scalar_view_t hash_op_per_insert;
#endif

    calculateC(
        size_type m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,
        a_scalar_view_t valuesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesB_,
        b_scalar_view_t valuesB_,


        c_row_view_t row_mapC_,
        c_nnz_view_t entriesC_,
        c_scalar_view_t valuesC_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_), valuesA(valuesA_),
          row_mapB(row_mapB_), entriesB(entriesB_), valuesB(valuesB_),
          row_mapC(row_mapC_), entriesC(entriesC_), valuesC(valuesC_),
          hash_size(0), pow2_hash_func(0),
          maxUncompressedNonZero(16384), maxCompressedNonZero(0), requiredMem(KOKKOSKERNELS_SPGEMM_SHMEMSIZE),
          multicore_chunk_size(0)
          //,KOKKOSKERNELS_SPGEMM_HASHSIZE(16)//, KOKKOSKERNELS_SPGEMM_BLOCKSIZE(block_size)
#ifdef TRACK_INSERTS
          //,hash_op_per_insert("hash", m)
    ,overal_inserts("inserts", 1)
          ,overall_hash_op ("hash", 1)
    ,wasted_ops("wasted",1)
    ,used_hashes("hashes", 1)
#endif
        {
        }

    calculateC(
        size_type m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,
        a_scalar_view_t valuesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesB_,
        b_scalar_view_t valuesB_,


        c_row_view_t row_mapC_,
        c_nnz_view_t entriesC_,
        c_scalar_view_t valuesC_,
        c_nnz_lno_t min_pow2_hash_size,
        c_nnz_lno_t MaxRoughNonZero_,
        size_t multicore_chunk_size_
        //,c_nnz_lno_t MaxFlops_
        ):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_), valuesA(valuesA_),
          row_mapB(row_mapB_), entriesB(entriesB_), valuesB(valuesB_),
          row_mapC(row_mapC_), entriesC(entriesC_), valuesC(valuesC_),
          hash_size(min_pow2_hash_size), pow2_hash_func(min_pow2_hash_size -1),
          maxUncompressedNonZero(MaxRoughNonZero_),maxCompressedNonZero(0),
          requiredMem(MaxRoughNonZero_ * 2 * sizeof(c_nnz_lno_t) +
              MaxRoughNonZero_ * sizeof(c_scalar_t) +
              40 + hash_size * 2 * sizeof(c_nnz_lno_t) ),
          multicore_chunk_size(multicore_chunk_size_)
          //,MaxFlops(MaxFlops_)

#ifdef TRACK_INSERTS
          //,hash_op_per_insert("hash", m)
    ,overal_inserts("inserts", 1)
          ,overall_hash_op ("hash", 1)
    ,wasted_ops("wasted",1)
    ,used_hashes("hashes", 1)
#endif
        {
        }


    calculateC(
        size_type m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,
        a_scalar_view_t valuesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesB_,
        b_scalar_view_t valuesB_,

        row_lno_temp_work_view_t compressed_b_rowmap_,
        row_lno_temp_work_view_t compressed_b_set_indices_,
        row_lno_temp_work_view_t compressed_b_sets_,
        row_lno_temp_work_view_t compressed_b_set_begins_,
        row_lno_temp_work_view_t compressed_b_set_nexts_,

        c_row_view_t row_mapC_,
        c_nnz_view_t entriesC_,
        c_scalar_view_t valuesC_,
        c_nnz_lno_t min_pow2_hash_size,
        c_nnz_lno_t maxUncompressedNNZ,
        c_nnz_lno_t maxCompressedNNZ,
        size_t multicore_chunk_size_
        //,c_nnz_lno_t MaxFlops_
        ):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_), valuesA(valuesA_),
          row_mapB(row_mapB_), entriesB(entriesB_), valuesB(valuesB_),

          compressed_b_rowmap(compressed_b_rowmap_),
          compressed_b_set_indices(compressed_b_set_indices_),
          compressed_b_sets(compressed_b_sets_),
          compressed_b_set_begins(compressed_b_set_begins_),
          compressed_b_set_nexts(compressed_b_set_nexts_),

          row_mapC(row_mapC_), entriesC(entriesC_), valuesC(valuesC_),
          hash_size(min_pow2_hash_size), pow2_hash_func(min_pow2_hash_size -1),
          maxUncompressedNonZero(maxUncompressedNNZ),maxCompressedNonZero(maxCompressedNNZ),
          requiredMem(
              maxCompressedNonZero * 3 * sizeof(c_nnz_lno_t) +
              maxCompressedNonZero * 8 * sizeof(c_scalar_t) *  sizeof(c_nnz_lno_t) +
              100 + hash_size * 2 * sizeof(c_nnz_lno_t)
              ),
          multicore_chunk_size(multicore_chunk_size_)
          //,MaxFlops(MaxFlops_)

#ifdef TRACK_INSERTS
          //,hash_op_per_insert("hash", m)
    ,overal_inserts("inserts", 1)
          ,overall_hash_op ("hash", 1)
    ,wasted_ops("wasted",1)
    ,used_hashes("hashes", 1)
#endif
        {
        }
    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreDenseAccumulatorTag&, const team_member_t & teamMember) const {
      size_type row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);

      c_nnz_lno_t *used_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(maxUncompressedNonZero * sizeof (c_nnz_lno_t));
      c_scalar_t *dense_accumulator = (c_scalar_t *) teamMember.team_shmem().get_shmem(maxUncompressedNonZero * sizeof (c_scalar_t));


      for (int i = 0; i < maxUncompressedNonZero; ++i){
        dense_accumulator[i] = 0;
      }

      for (; row_index < row_max; ++row_index){
        c_nnz_lno_t used_index_count = 0;
        //c_nnz_lno_t next_available_entry = 0;


        //check ii is out of range. if it is, just return.
        size_type col_begin = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];

        for (; col_begin < col_end; ++col_begin){
          //printf("col_begin:%d\n", col_begin);
          const nnz_lno_t rowb = entriesA[col_begin];
          const c_scalar_t aval = valuesA[col_begin];

          //std::cout << "  for column:" << rowb << std::endl;
          size_type colb_ind = row_mapB[rowb];
          const size_type col_b_end = row_mapB[rowb + 1];
          for (; colb_ind < col_b_end; ++colb_ind){
              const c_nnz_lno_t b_col_ind = entriesB[colb_ind];
              if (dense_accumulator[b_col_ind] == 0){
                used_indices[used_index_count++] = b_col_ind;
              }
              dense_accumulator[b_col_ind] +=valuesB[colb_ind] * aval;
          }
        }
        c_nnz_lno_t  my_nonzeros = row_mapC(row_index);
        for (c_nnz_lno_t i = 0; i < used_index_count; ++i){
          c_nnz_lno_t index = used_indices[i];
          entriesC(my_nonzeros) = index;
          valuesC(my_nonzeros++) = dense_accumulator[index];
          dense_accumulator[index] = 0;
        }

      }

    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const FillExtraTag&, const team_member_t & teamMember) const {
      size_type row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);


      int lnot_size = sizeof (c_nnz_lno_t) * 8;

      int multiplier_shift = 0;
      //const c_nnz_lno_t remainder_and = lnot_size - 1;

      while (lnot_size > 1) {
        ++multiplier_shift;
        lnot_size = lnot_size >> 1;
      }
      lnot_size = sizeof (c_nnz_lno_t) * 8;


      //std::cout << "multiplier shift:" << multiplier_shift << " remainder_and:" << remainder_and  << std::endl;

      const size_t hash_begin_array_size = hash_size * sizeof(c_nnz_lno_t);
      c_nnz_lno_t *used_hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);

      const size_t hash_index_array_size = maxCompressedNonZero * sizeof(c_nnz_lno_t);
      c_nnz_lno_t *hash_nexts = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_set_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);

      const size_t multiplication_val_size = maxCompressedNonZero * sizeof(c_nnz_lno_t) * 8;
      c_scalar_t *set_multiplication_vals = (c_scalar_t *) teamMember.team_shmem().get_shmem(
                  sizeof(c_scalar_t) * multiplication_val_size);

      /*
      for (size_t i = 0; i < maxCompressedNonZero; ++i){
        hash_nexts[i] = -1;
      }
      */

      for (size_t i = 0; i < hash_size; ++i){
        hash_begins[i] = -1;
      }

      for (size_t i = 0; i < multiplication_val_size; ++i){
        set_multiplication_vals[i] = 0;
      }
#ifdef HASHTRACK
      size_t numhashes = 0;
      size_t numpoints = 0;
      size_t num_comp = 0;
      size_t num_mul = 0;
#endif
      for (; row_index < row_max; ++row_index){
        c_nnz_lno_t used_hash_counts = 0;
        c_nnz_lno_t next_available_entry = 0;

        //check ii is out of range. if it is, just return.
        size_type col_begin = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];
        for (; col_begin < col_end; ++col_begin){

          //printf("col_begin:%d\n", col_begin);
          const nnz_lno_t rowb = entriesA[col_begin];
          const c_scalar_t aval = valuesA[col_begin];

          size_type colb_ind = compressed_b_rowmap[rowb];

          const size_type col_b_end = compressed_b_rowmap[rowb + 1];

          for (; colb_ind < col_b_end; ){

            const c_nnz_lno_t b_col_set_ind = compressed_b_set_indices[colb_ind];
            const c_nnz_lno_t b_col_set = compressed_b_sets[colb_ind];
            //const c_scalar_t bmul = valuesB[colb_ind] * aval;
            {
              //printf("b_col_set_index:%d", b_col_set_index)
              const c_nnz_lno_t hash = (b_col_set_ind ) & pow2_hash_func;
              c_nnz_lno_t hash_ind = hash_begins[hash];
#ifdef HASHTRACK
              ++num_mul;
#endif
              if (hash_ind == -1){
                hash_begins[hash] = next_available_entry;
                hash_set_values[next_available_entry] = b_col_set;

                const c_nnz_lno_t accumulator_index = next_available_entry * sizeof (c_nnz_lno_t) * 8;
                {
                  c_nnz_lno_t b_adj_index = compressed_b_set_begins[colb_ind];
                  const c_nnz_lno_t b_adj_index_end = compressed_b_set_begins[colb_ind + 1];

                  for (; b_adj_index < b_adj_index_end; ++b_adj_index){
                    const c_scalar_t b_scalar = valuesB [b_adj_index] * aval;
                    const c_nnz_lno_t b_col_index = entriesB[b_adj_index];
                    set_multiplication_vals[accumulator_index +  b_col_index % (sizeof (c_nnz_lno_t) * 8) ] += b_scalar;
                  }

                  /*
                  do{

                    const c_nnz_lno_t b_col_index = entriesB[b_adj_index];
                    const c_scalar_t b_scalar = valuesB [b_adj_index] * aval;
                    set_multiplication_vals[accumulator_index +  b_adj_index % (sizeof (c_nnz_lno_t) * 8) ] += b_scalar;
                    b_adj_index = compressed_b_set_nexts[b_adj_index];
                  }while (b_adj_index != -1);
                  */
                }
                hash_indices[next_available_entry++] = b_col_set_ind;
                used_hash_indices[used_hash_counts++] = hash;
                //std::cout << "\t\tDONE INSERTEDED for the first time" << std::endl;
              }
              else {
                do {
                  if(b_col_set_ind ^ hash_indices[hash_ind]){ //if they are not equal
                    hash_ind = hash_nexts[hash_ind];
#ifdef HASHTRACK
                    ++num_comp;
#endif
                  }
                  else {
                    hash_set_values[hash_ind] = hash_set_values[hash_ind] | b_col_set;
                    const c_nnz_lno_t accumulator_index = hash_ind * sizeof (c_nnz_lno_t) * 8;
                    {

                      c_nnz_lno_t b_adj_index = compressed_b_set_begins[colb_ind];
                      const c_nnz_lno_t b_adj_index_end = compressed_b_set_begins[colb_ind + 1];

                      for (; b_adj_index < b_adj_index_end; ++b_adj_index){
                        const c_scalar_t b_scalar = valuesB [b_adj_index] * aval;
                        const c_nnz_lno_t b_col_index = entriesB[b_adj_index];
                        set_multiplication_vals[accumulator_index +  b_col_index % (sizeof (c_nnz_lno_t) * 8) ] += b_scalar;
                      }
                      /*
                      do{
                        const c_nnz_lno_t b_col_index = entriesB[b_adj_index];
                        const c_scalar_t b_scalar = valuesB [b_adj_index] * aval;
                        set_multiplication_vals[accumulator_index + b_adj_index % (sizeof (c_nnz_lno_t) * 8) ] += b_scalar;
                        b_adj_index = compressed_b_set_nexts[b_adj_index];
                      }while (b_adj_index != -1);
                       */
                    }
                    //hash_set_values[hash_ind] = hash_set_values[hash_ind] +  valuesB[colb_ind] * aval;
                    goto endloop; //break;
                  }
                }
                while(hash_ind != -1);

                hash_nexts[next_available_entry] = hash_begins[hash];
                hash_begins[hash] = next_available_entry;
                hash_set_values[next_available_entry] = b_col_set;
                const c_nnz_lno_t accumulator_index = next_available_entry * sizeof (c_nnz_lno_t) * 8;
                {

                  c_nnz_lno_t b_adj_index = compressed_b_set_begins[colb_ind];
                  const c_nnz_lno_t b_adj_index_end = compressed_b_set_begins[colb_ind + 1];

                  for (; b_adj_index < b_adj_index_end; ++b_adj_index){
                    const c_scalar_t b_scalar = valuesB [b_adj_index] * aval;
                    const c_nnz_lno_t b_col_index = entriesB[b_adj_index];
                    set_multiplication_vals[accumulator_index +  b_col_index % (sizeof (c_nnz_lno_t) * 8) ] += b_scalar;
                  }
                  /*
                  do{

                    const c_nnz_lno_t b_col_index = entriesB[b_adj_index];

                    const c_scalar_t b_scalar = valuesB [b_adj_index] * aval;

                    set_multiplication_vals[accumulator_index +  b_adj_index % (sizeof (c_nnz_lno_t) * 8) ] += b_scalar;
                    b_adj_index = compressed_b_set_nexts[b_adj_index];
                  }while (b_adj_index != -1);
                  */
                }
                hash_indices[next_available_entry++] = b_col_set_ind;

              }
            }
            endloop:
            ++colb_ind;
          }
        }

        c_nnz_lno_t  my_nonzeros = row_mapC(row_index);
        for (c_nnz_lno_t i = 0; i < used_hash_counts; ++i){

          const int hash = used_hash_indices[i];
          hash_begins[hash] = -1;
        }

        for (c_nnz_lno_t i = 0; i < next_available_entry; ++i){
            c_nnz_lno_t b_col_set_index = hash_indices[i] * sizeof(c_nnz_lno_t) * 8;
            c_nnz_lno_t b_col_set= hash_set_values[i];

            c_nnz_lno_t in_set_index = 0;
            c_nnz_lno_t shift_val = 1;

            while (b_col_set){
              if (b_col_set & shift_val){
                const c_nnz_lno_t accumulator_index =
                    i * sizeof (c_nnz_lno_t) * 8 +  in_set_index;
                entriesC(my_nonzeros) = b_col_set_index + in_set_index;
                valuesC(my_nonzeros++) = set_multiplication_vals[accumulator_index];
                set_multiplication_vals[accumulator_index] =0;
                b_col_set = b_col_set ^ shift_val;
              }
              shift_val = shift_val << 1;
              ++in_set_index;
            }
        }
      }

#ifdef HASHTRACK
      printf("numhash:%ld np:%ld numComparison:%ld numOP:%ld numrows:%ld\n", numhashes, numpoints, num_comp, num_mul,row_max - teamMember.league_rank() * multicore_chunk_size);
#endif
    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTagCount&, const team_member_t & teamMember) const {
      size_type row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);


      const size_t hash_index_array_size = maxUncompressedNonZero * sizeof(c_nnz_lno_t);
      const size_t hash_begin_array_size = hash_size * sizeof(c_nnz_lno_t);
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);
      c_nnz_lno_t *hash_nexts = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *used_hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);

      /*
      for (size_t i = 0; i < maxUncompressedNonZero; ++i){
        hash_nexts[i] = -1;
      }
      */

      for (size_t i = 0; i < hash_size; ++i){
        hash_begins[i] = -1;
      }
#ifdef HASHTRACK
      size_t numhashes = 0;
      size_t numpoints = 0;
      size_t num_comp = 0;
      size_t num_mul = 0;
#endif
      for (; row_index < row_max; ++row_index){
        c_nnz_lno_t used_hash_counts = 0;
        c_nnz_lno_t next_available_entry = 0;


        c_nnz_view_t hash_indices = Kokkos::subview(entriesC, Kokkos::make_pair(row_mapC(row_index), row_mapC(row_index + 1)));
        c_scalar_view_t hash_values = Kokkos::subview(valuesC, Kokkos::make_pair(row_mapC(row_index), row_mapC(row_index + 1)));

        size_type col_begin = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];


        for (; col_begin < col_end; ++col_begin){
          //printf("col_begin:%d\n", col_begin);
          const nnz_lno_t rowb = entriesA[col_begin];
          const c_scalar_t aval = valuesA[col_begin];

          //std::cout << "  for column:" << rowb << std::endl;
          size_type colb_ind = row_mapB[rowb];
          const size_type col_b_end = row_mapB[rowb + 1];
          for (; colb_ind < col_b_end; ){
            //printf("colb_ind:%d\n", colb_ind);

            const c_nnz_lno_t b_col_ind = entriesB[colb_ind];
            //const c_scalar_t bmul = valuesB[colb_ind] * aval;
            {
              //printf("b_col_set_index:%d", b_col_set_index)
              const c_nnz_lno_t hash = (b_col_ind ) & pow2_hash_func;
              //std::cout << "\tinserting:" << b_col_set_index << " hash:" << hash << std::endl;
              c_nnz_lno_t hash_ind = hash_begins[hash];
              //std::cout << "\t\thash begin:" << hash_ind << std::endl;
#ifdef HASHTRACK
              ++num_mul;
#endif
              if (hash_ind == -1){
                //std::cout << "\t\tinserting for the first time" << std::endl;
                hash_begins[hash] = next_available_entry;
                hash_indices[next_available_entry] = b_col_ind;
                hash_values[next_available_entry++] = valuesB[colb_ind] * aval;
                used_hash_indices[used_hash_counts++] = hash;
              }
              else {
                do {
                  if(b_col_ind ^ hash_indices[hash_ind]){ //if they are not equal
                    hash_ind = hash_nexts[hash_ind];
#ifdef HASHTRACK
                    ++num_comp;
#endif
                  }
                  else {
                    hash_values[hash_ind] = hash_values[hash_ind] +  valuesB[colb_ind] * aval;
                    goto endloop; //break;
                  }
                }
                while (hash_ind != -1);

                //std::cout << "\t\tinserting for the first time" << std::endl;

                hash_nexts[next_available_entry] = hash_begins[hash];
                hash_begins[hash] = next_available_entry;
                hash_indices[next_available_entry] = b_col_ind;
                hash_values[next_available_entry++] = valuesB[colb_ind] * aval;
              }
            }
            endloop:
            ++colb_ind;
          }
        }

        for (c_nnz_lno_t i = 0; i < used_hash_counts; ++i){
#ifdef HASHTRACK
          ++numhashes;
#endif
          const int hash = used_hash_indices[i];
          hash_begins[hash] = -1;

#ifdef HASHTRACK
            ++numpoints;
#endif

        }

      }
#ifdef HASHTRACK
      printf("numhash:%ld np:%ld numComparison:%ld numOP:%ld numrows:%ld\n", numhashes, numpoints, num_comp, num_mul,row_max - teamMember.league_rank() * multicore_chunk_size);
#endif
    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const FillTag&, const team_member_t & teamMember) const {

      const int team_size = teamMember.team_size();
      const int team_rank = teamMember.team_rank();

      size_type row_index = (teamMember.league_rank()  * team_size + team_rank) * KOKKOSKERNELS_SPGEMM_BLOCKSIZE;
      //if (row_index >= m) return;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + KOKKOSKERNELS_SPGEMM_BLOCKSIZE, m);
      //int row_max = row_index + KOKKOSKERNELS_SPGEMM_BLOCKSIZE;
      //if(row_max > m) row_max = m;



      const int hash_array_bytes = KOKKOSKERNELS_SPGEMM_HASHSIZE * sizeof(c_nnz_lno_t) * team_size;
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_bytes);
      c_nnz_lno_t *hash_ends = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_bytes);

      const int overall_set_size = ((KOKKOSKERNELS_SPGEMM_SHMEMSIZE - hash_array_bytes * 2) - 16) / (sizeof(c_nnz_lno_t) + sizeof(c_scalar_t));


      c_nnz_lno_t * team_shared_keys = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(
          overall_set_size * sizeof(c_nnz_lno_t));
      c_scalar_t * team_shared_values = (c_scalar_t *) teamMember.team_shmem().get_shmem(
          overall_set_size * sizeof(c_scalar_t));

      const int team_set_size = overall_set_size / team_size ;
      hash_ends += team_rank * KOKKOSKERNELS_SPGEMM_HASHSIZE;
      hash_begins += team_rank * KOKKOSKERNELS_SPGEMM_HASHSIZE;
      team_shared_keys += team_rank * team_set_size;
      team_shared_values += team_rank * team_set_size;

      const int team_hash_set_size = team_set_size / KOKKOSKERNELS_SPGEMM_HASHSIZE;
      for (int i = 0; i < KOKKOSKERNELS_SPGEMM_HASHSIZE; ++i){
        hash_ends[i] = hash_begins[i] = team_hash_set_size * i;
      }

      for (; row_index < row_max; ++row_index)
      {
        //unsigned int used_hashes = 0;
        //check ii is out of range. if it is, just return.
        size_type col_a_index = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];

#ifdef TRACK_INSERTS
        long hash_ops = 0, inserts = 0;
        long wasted_op_count = 0;
#endif
        //c_nnz_lno_t  current_heap_size = 0;

        if (col_a_index < col_end){
          const nnz_lno_t rowb = entriesA[col_a_index];
          const c_scalar_t aval = valuesA[col_a_index];
          size_type col_b_begin = row_mapB[rowb];
          const size_type col_b_end = row_mapB[rowb + 1];

          for (; col_b_begin < col_b_end; ++col_b_begin){
            const c_nnz_lno_t b_col_set_index = entriesB[col_b_begin];
            const int hash = b_col_set_index & KOKKOSKERNELS_SPGEMM_ANDHASHFUNCTION;//b_col_set_index % HASHSIZE;
            //used_hashes = used_hashes | 1u << hash;
            const c_nnz_lno_t hash_ind = hash_ends[hash]++;
            team_shared_keys[hash_ind] = b_col_set_index;
            team_shared_values[hash_ind] = valuesB[col_b_begin] * aval;
          }
        }
        ++col_a_index;


        for (; col_a_index < col_end; ++col_a_index){


          const nnz_lno_t rowb = entriesA[col_a_index];
          const c_scalar_t aval = valuesA[col_a_index];

          size_type col_b_begin = row_mapB[rowb];
          const size_type col_b_end = row_mapB[rowb + 1];

#ifdef TRACK_INSERTS
          inserts += col_b_end - col_b_begin;
          int hashes_before[HASHSIZE];
          for (int h = 0; h < HASHSIZE; ++h){
            hashes_before[h] = hash_ends[h];
          }
#endif



          for (; col_b_begin < col_b_end; ){
            const c_nnz_lno_t b_col_set_index = entriesB[col_b_begin];
            const int hash = b_col_set_index & KOKKOSKERNELS_SPGEMM_ANDHASHFUNCTION;//b_col_set_index % HASHSIZE;
            //used_hashes = used_hashes | 1u << hash;
            c_nnz_lno_t hash_ind = hash_begins[hash];//team_hash_set_size * hash;//hash_begins[hash];
            const c_nnz_lno_t hashend = hash_ends[hash];

            for (; hash_ind < hashend; ){
#ifdef TRACK_INSERTS
              ++hash_ops;
              if (hash_ind >= hashes_before[hash]){
                wasted_op_count++;
              }
#endif
              if(b_col_set_index ^ team_shared_keys[hash_ind]){ //if they are not equal
                ++hash_ind;
              }
              else {
                team_shared_values[hash_ind] = team_shared_values[hash_ind] + valuesB[col_b_begin] * aval;
                goto endloop; //break;
              }
            }

#ifdef TRACK_INSERTS
              ++hash_ops;
#endif
            team_shared_keys[hash_ind] = b_col_set_index;
            team_shared_values[hash_ind] = valuesB[col_b_begin] * aval;
            ++hash_ends[hash];// = i + 1;

            endloop:
            ++col_b_begin;
          }
        }

#ifdef TRACK_INSERTS
        overal_inserts(0) += inserts;
        overall_hash_op(0) += hash_ops;
        wasted_ops(0) += wasted_op_count;
#endif

        c_nnz_lno_t adjInd = row_mapC(row_index);


        for (int i = 0; i < KOKKOSKERNELS_SPGEMM_HASHSIZE; ++i){
          //if (used_hashes & (1u << i))
          {


            const c_nnz_lno_t hashend = hash_ends[i];

            /*
            if (i+1 < KOKKOSKERNELS_SPGEMM_HASHSIZE && qhashend > hash_begins[i + 1]){
              std::cout << "row:" << row_index << " i:" <<i <<" hashend:" << hashend << " hash_begins[i + 1]:" << hash_begins[i + 1] << std::endl;
            }
            */
#ifdef TRACK_INSERTS
            if (hash_begins[i] < hashend){
              used_hashes(0) += 1;
            }
#endif
            c_nnz_lno_t current_heap_ind = hash_ends[i] = hash_begins[i];
            for ( /*team_hash_set_size * i*/; current_heap_ind < hashend; ++current_heap_ind){
              entriesC(adjInd) = team_shared_keys[current_heap_ind];
              valuesC(adjInd++) = team_shared_values[current_heap_ind];
            }
          }
        }
      }
    }

    // Provide the shared memory capacity.
    // This function takes the team_size as an argument ,
    // which allows team_size dependent allocations.
    size_t team_shmem_size (int team_size) const {
      return requiredMem;
    }


  };

  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t>
  struct SymbolicC{
    size_type m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
    c_row_view_t compressedrowmapC;
    c_nnz_view_t entriesSetIndicesC;
    c_nnz_view_t entriesSetsC;
    typedef typename a_row_view_t::non_const_value_type c_row_lno_t;
    typedef typename b_nnz_view_t::non_const_value_type c_nnz_lno_t;

    const c_nnz_lno_t hash_size;
    const c_nnz_lno_t pow2_hash_func;
    const c_nnz_lno_t MaxRoughNonZero;

    const size_t requiredMem;
    const size_t multicore_chunk_size;

    //const int KOKKOSKERNELS_SPGEMM_HASHSIZE;
    //const int KOKKOSKERNELS_SPGEMM_BLOCKSIZE;


    SymbolicC(
        size_type m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,
        c_row_view_t rowmapC_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_mapB(row_mapB_), entriesSetIndicesB(entriesSetIndicesB_), entriesSetsB(entriesSetsB_),
          rowmapC(rowmapC_),
          hash_size(0), pow2_hash_func(0),
          MaxRoughNonZero(16384), requiredMem(KOKKOSKERNELS_SPGEMM_SHMEMSIZE),
          multicore_chunk_size(0),compressedrowmapC()
          //,KOKKOSKERNELS_SPGEMM_HASHSIZE(16),BLOCKSIZE (block_size)
          {}

    SymbolicC(
        size_type m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,
        c_row_view_t rowmapC_,
        c_nnz_lno_t min_pow2_hash_size,
        c_nnz_lno_t MaxRoughNonZero_,
        size_t multicore_chunk_size_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_mapB(row_mapB_), entriesSetIndicesB(entriesSetIndicesB_), entriesSetsB(entriesSetsB_),
          rowmapC(rowmapC_),
          hash_size(min_pow2_hash_size), pow2_hash_func(min_pow2_hash_size -1),
          MaxRoughNonZero(MaxRoughNonZero_),
          requiredMem(MaxRoughNonZero_ * 3 * sizeof(c_nnz_lno_t) + 40 + hash_size * 2 * sizeof(c_nnz_lno_t) ),
          multicore_chunk_size(multicore_chunk_size_),compressedrowmapC()
          {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const FillExtraTag&, const team_member_t & teamMember) const {
      size_type row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);
      //printf("myRowMin:%d rowMax:%d\n", row_index, row_max);


      const size_t hash_index_array_size = MaxRoughNonZero * sizeof(c_nnz_lno_t);
      const size_t hash_begin_array_size = hash_size * sizeof(c_nnz_lno_t);
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);
      c_nnz_lno_t *hash_nexts = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *used_hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);


      /*
      for (size_t i = 0; i < MaxRoughNonZero; ++i){
        hash_nexts[i] = -1;
      }
      */

      for (size_t i = 0; i < hash_size; ++i){
        hash_begins[i] = -1;
      }
#ifdef HASHTRACK
      size_t numhashes = 0;
      size_t numpoints = 0;
      size_t num_comp = 0;
      size_t num_mul = 0;
#endif
      for (; row_index < row_max; ++row_index){
        c_nnz_lno_t used_hash_counts = 0;
        c_nnz_lno_t next_available_entry = 0;

        //check ii is out of range. if it is, just return.
        size_type col_begin = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];
        for (; col_begin < col_end; ++col_begin){
          //printf("col_begin:%d\n", col_begin);
          const nnz_lno_t rowb = entriesA[col_begin];

          //std::cout << "  for column:" << rowb << std::endl;
          size_type colb_ind = row_mapB[rowb];
          const size_type col_b_end = row_mapB[rowb + 1];
          for (; colb_ind < col_b_end; ){
            //printf("colb_ind:%d\n", colb_ind);
            const c_nnz_lno_t b_col_set_index = entriesSetIndicesB[colb_ind];
            const c_nnz_lno_t b_col_set = entriesSetsB[colb_ind];
            {
              //printf("b_col_set_index:%d", b_col_set_index)
              const c_nnz_lno_t hash = b_col_set_index & pow2_hash_func;
              //std::cout << "\tinserting:" << b_col_set_index << " hash:" << hash << std::endl;
              c_nnz_lno_t hash_ind = hash_begins[hash];
              //std::cout << "\t\thash begin:" << hash_ind << std::endl;
#ifdef HASHTRACK
              ++num_mul;
#endif
              if (hash_ind == -1){
                //std::cout << "\t\tinserting for the first time" << std::endl;
                hash_nexts[next_available_entry] = hash_begins[hash];
                hash_begins[hash] = next_available_entry;
                hash_indices[next_available_entry] = b_col_set_index;
                hash_values[next_available_entry++] = b_col_set;
                used_hash_indices[used_hash_counts++] = hash;
              }
              else {
                while (hash_nexts[hash_ind] != -1){
                  if(b_col_set_index ^ hash_indices[hash_ind]){ //if they are not equal
                    hash_ind = hash_nexts[hash_ind];
#ifdef HASHTRACK
                    ++num_comp;
#endif
                  }
                  else {
                    hash_values[hash_ind] = hash_values[hash_ind] | b_col_set;
                    goto endloop; //break;
                  }
                }
                if(b_col_set_index ^ hash_indices[hash_ind]){ //if they are not equal

                  hash_nexts[next_available_entry] = hash_begins[hash];
                  hash_begins[hash] = next_available_entry;
                  hash_indices[next_available_entry] = b_col_set_index;
                  hash_values[next_available_entry++] = b_col_set;
                }
                else {
                  hash_values[hash_ind] = hash_values[hash_ind] | b_col_set;
                }
              }
            }
            endloop:
            ++colb_ind;
          }
        }
        c_nnz_lno_t  my_nonzeros = 0;
        for (c_nnz_lno_t i = 0; i < used_hash_counts; ++i){
          const int hash = used_hash_indices[i];
          hash_begins[hash] = -1;
        }

        for (c_nnz_lno_t i = 0; i < next_available_entry; ++i){
            c_nnz_lno_t c_rows = hash_values[i];
            for (; c_rows; my_nonzeros++)
            {
              c_rows &= c_rows - 1; // clear the least significant bit set
            }
        }
        rowmapC(row_index) = my_nonzeros;
        compressedrowmapC(row_index) = next_available_entry;
      }

    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTagCount&, const team_member_t & teamMember) const {
      size_type row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);
      //printf("myRowMin:%d rowMax:%d\n", row_index, row_max);


      const size_t hash_index_array_size = MaxRoughNonZero * sizeof(c_nnz_lno_t);
      const size_t hash_begin_array_size = hash_size * sizeof(c_nnz_lno_t);
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);
      c_nnz_lno_t *hash_nexts = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *used_hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);

      /*
      for (size_t i = 0; i < MaxRoughNonZero; ++i){
        hash_nexts[i] = -1;
      }
      */

      for (size_t i = 0; i < hash_size; ++i){
        hash_begins[i] = -1;
      }
#ifdef HASHTRACK
      size_t numhashes = 0;
      size_t numpoints = 0;
      size_t num_comp = 0;
      size_t num_mul = 0;
#endif
      for (; row_index < row_max; ++row_index){
        c_nnz_lno_t used_hash_counts = 0;
        c_nnz_lno_t next_available_entry = 0;

        //check ii is out of range. if it is, just return.
        size_type col_begin = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];
        for (; col_begin < col_end; ++col_begin){
          //printf("col_begin:%d\n", col_begin);
          const nnz_lno_t rowb = entriesA[col_begin];

          //std::cout << "  for column:" << rowb << std::endl;
          size_type colb_ind = row_mapB[rowb];
          const size_type col_b_end = row_mapB[rowb + 1];
          for (; colb_ind < col_b_end; ){
            //printf("colb_ind:%d\n", colb_ind);
            const c_nnz_lno_t b_col_set_index = entriesSetIndicesB[colb_ind];
            const c_nnz_lno_t b_col_set = entriesSetsB[colb_ind];
            {
              //printf("b_col_set_index:%d", b_col_set_index)
              const c_nnz_lno_t hash = b_col_set_index & pow2_hash_func;
              //std::cout << "\tinserting:" << b_col_set_index << " hash:" << hash << std::endl;
              c_nnz_lno_t hash_ind = hash_begins[hash];
              //std::cout << "\t\thash begin:" << hash_ind << std::endl;
#ifdef HASHTRACK
              ++num_mul;
#endif
              if (hash_ind == -1){
                //std::cout << "\t\tinserting for the first time" << std::endl;
                hash_nexts[next_available_entry] = hash_begins[hash];
                hash_begins[hash] = next_available_entry;
                hash_indices[next_available_entry] = b_col_set_index;
                hash_values[next_available_entry++] = b_col_set;
                used_hash_indices[used_hash_counts++] = hash;
              }
              else {
                while (hash_nexts[hash_ind] != -1){
                  if(b_col_set_index ^ hash_indices[hash_ind]){ //if they are not equal
                    hash_ind = hash_nexts[hash_ind];
#ifdef HASHTRACK
                    ++num_comp;
#endif
                  }
                  else {
                    hash_values[hash_ind] = hash_values[hash_ind] | b_col_set;
                    goto endloop; //break;
                  }
                }
                if(b_col_set_index ^ hash_indices[hash_ind]){ //if they are not equal
                  hash_nexts[next_available_entry] = hash_begins[hash];
                  hash_begins[hash] = next_available_entry;
                  hash_indices[next_available_entry] = b_col_set_index;
                  hash_values[next_available_entry++] = b_col_set;
                }
                else {
                  hash_values[hash_ind] = hash_values[hash_ind] | b_col_set;
                }
              }
            }
            endloop:
            ++colb_ind;
          }
        }

        c_nnz_lno_t  my_nonzeros = 0;
        for (c_nnz_lno_t i = 0; i < used_hash_counts; ++i){
          const int hash = used_hash_indices[i];
          hash_begins[hash] = -1;
        }

        for (c_nnz_lno_t i = 0; i < next_available_entry; ++i){
            c_nnz_lno_t c_rows = hash_values[i];
            for (; c_rows; my_nonzeros++)
            {
              c_rows &= c_rows - 1; // clear the least significant bit set
            }
        }
        rowmapC(row_index) = my_nonzeros;

      }
#ifdef HASHTRACK
      printf("numhash:%ld np:%ld numComparison:%ld numOP:%ld\n", numhashes, numpoints, num_comp, num_mul);
#endif

    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {
      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      size_type row_index = (teamMember.league_rank()  * team_size + team_rank) * KOKKOSKERNELS_SPGEMM_BLOCKSIZE;

      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + KOKKOSKERNELS_SPGEMM_BLOCKSIZE, m);

      const int hash_array_size = KOKKOSKERNELS_SPGEMM_HASHSIZE * sizeof(c_nnz_lno_t) * team_size;

      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_size);
      c_nnz_lno_t *hash_ends = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_size);

      const int set_size = (8184 - hash_array_size );

      c_nnz_lno_t * team_shared_keys = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(set_size);
      c_nnz_lno_t * team_shared_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(set_size);

      const int team_set_size = (set_size / team_size) / (sizeof (c_nnz_lno_t)) ;

      //set_size = (set_size / team_size) / sizeof(c_nnz_lno_t);
      //hash_array_size = (hash_array_size / team_size) / sizeof(c_nnz_lno_t);

      hash_ends += team_rank * KOKKOSKERNELS_SPGEMM_HASHSIZE;
      hash_begins += team_rank * KOKKOSKERNELS_SPGEMM_HASHSIZE;
      team_shared_keys += team_rank * team_set_size;
      team_shared_values += team_rank * team_set_size;
      const int team_hash_set_size = team_set_size / KOKKOSKERNELS_SPGEMM_HASHSIZE;

      for (int i = 0; i < KOKKOSKERNELS_SPGEMM_HASHSIZE; ++i){
        hash_ends[i] = hash_begins[i] = team_hash_set_size * i;
      }


      for (; row_index < row_max; ++row_index){
        //check ii is out of range. if it is, just return.
        const size_type col_begin = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];

        for (size_type col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
          const nnz_lno_t rowb = entriesA[col_a_index];
          size_type colb_ind = row_mapB[rowb];
          const size_type col_b_end = row_mapB[rowb + 1];
          for (; colb_ind < col_b_end; ){
            const c_nnz_lno_t b_col_set_index = entriesSetIndicesB[colb_ind];
            const c_nnz_lno_t b_col_set = entriesSetsB[colb_ind];
            //insert_to_heap(b_col_set_index, b_col_set, team_shared_keys, team_shared_values, hash_begins, hash_ends);

            {
              const int hash = b_col_set_index & KOKKOSKERNELS_SPGEMM_ANDHASHFUNCTION;
              c_nnz_lno_t hash_ind = hash_begins[hash];
              const c_nnz_lno_t hashend = hash_ends[hash];

              for (; hash_ind < hashend; ){
                if(b_col_set_index ^ team_shared_keys[hash_ind]){ //if they are not equal
                  ++hash_ind;
                }
                else {
                  team_shared_values[hash_ind] = team_shared_values[hash_ind] | b_col_set;
                  goto endloop; //break;
                }
              }
              team_shared_keys[hash_ind] = b_col_set_index;
              team_shared_values[hash_ind] = b_col_set;
              ++hash_ends[hash];// = i + 1;
            }
            endloop:
            ++colb_ind;
          }
        }
        c_nnz_lno_t  my_nonzeros = 0;
        for (int i = 0; i < KOKKOSKERNELS_SPGEMM_HASHSIZE; ++i){
          const c_nnz_lno_t hashend = hash_ends[i];
          c_nnz_lno_t current_heap_ind = hash_ends[i] = hash_begins[i];
          for (; current_heap_ind < hashend; ++current_heap_ind){
            c_nnz_lno_t c_rows = team_shared_values[current_heap_ind];
            for (; c_rows; my_nonzeros++)
            {
              c_rows &= c_rows - 1; // clear the least significant bit set
            }
          }
        }
        rowmapC(row_index) = my_nonzeros;
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const FillTag&, const team_member_t & teamMember) const {

      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      size_type row_index = teamMember.league_rank()  * team_size + team_rank;
      if (row_index >= m) return;


      int hash_array_size = KOKKOSKERNELS_SPGEMM_HASHSIZE * sizeof(c_nnz_lno_t) * team_size;
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_size);
      c_nnz_lno_t *hash_ends = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_size);

      int set_size = (KOKKOSKERNELS_SPGEMM_SHMEMSIZE - hash_array_size * 2) / 2 - 8;
      c_nnz_lno_t * team_shared_keys = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(set_size);
      c_nnz_lno_t * team_shared_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(set_size);
      set_size = set_size / team_size / sizeof(c_nnz_lno_t);
      hash_array_size = hash_array_size / team_size / sizeof(c_nnz_lno_t);

      hash_ends += team_rank * hash_array_size;
      hash_begins += team_rank * hash_array_size;
      team_shared_keys += team_rank * set_size;
      team_shared_values += team_rank * set_size;

      set_size = set_size / KOKKOSKERNELS_SPGEMM_HASHSIZE;
      for (int i = 0; i < KOKKOSKERNELS_SPGEMM_HASHSIZE; ++i){
        hash_begins[i] = hash_ends[i] = set_size * i;


      }



      //check ii is out of range. if it is, just return.
      const size_type col_begin = row_mapA[row_index];
      const size_type col_end = row_mapA[row_index + 1];

      c_nnz_lno_t  current_heap_size = 0;
      for (size_type col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t rowb = entriesA[col_a_index];
        const size_type col_b_begin = row_mapB[rowb];
        const size_type col_b_end = row_mapB[rowb + 1];
        for (size_type colb_ind = col_b_begin; colb_ind < col_b_end; ++colb_ind){
          c_nnz_lno_t b_col_set_index = entriesSetIndicesB[colb_ind];
          c_nnz_lno_t b_col_set = entriesSetsB[colb_ind];
          //insert_to_heap(b_col_set_index, b_col_set, team_shared_keys, team_shared_values, current_heap_size);
          insert_to_heap(b_col_set_index, b_col_set, team_shared_keys, team_shared_values, hash_begins, hash_ends);
        }
      }

      c_nnz_lno_t  current_heap_ind = 0;
      c_nnz_lno_t adjInd = rowmapC(row_index);
      c_nnz_lno_t my_nonzeros = 0;

      c_nnz_lno_t adj = 0;
      for (int i = 0; i < KOKKOSKERNELS_SPGEMM_HASHSIZE; ++i){
        for (c_nnz_lno_t current_heap_ind = hash_begins[i]; current_heap_ind < hash_ends[i]; ++current_heap_ind){
          entriesSetIndicesC(adj + adjInd) = team_shared_keys[current_heap_ind];
          c_nnz_lno_t c_rows = team_shared_values[current_heap_ind];
          entriesSetsC(adj + adjInd) = c_rows;

          /*
          if (row_index == 0) {
            std::cout << "row_index:" << row_index << " key:" << team_shared_keys[current_heap_ind] << " c_rows:" << c_rows<< std::endl;
          }
          */
          ++adj;
          nnz_lno_t unit = 1;
          while (c_rows){
            if (c_rows & unit){
              my_nonzeros++;
            }
            c_rows = c_rows & ~unit;
            unit = unit << 1;
          }
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void insert_to_heap(
        const c_nnz_lno_t &key,
        const c_nnz_lno_t &val,
        c_nnz_lno_t *keys,
        c_nnz_lno_t *vals,
        c_nnz_lno_t *hash_begins,
        c_nnz_lno_t *hash_ends
#ifdef TRACK_INSERTS_SYM
        ,int &hash_ops
#endif
        ) const {

      //int hash = key % HASHSIZE;
      const int hash = key & KOKKOSKERNELS_SPGEMM_ANDHASHFUNCTION;
      //int hash = 0;

      int i = hash_begins[hash];
      for (; i < hash_ends[hash]; ++i){

        if (key == keys[i]){
          vals[i] = vals[i] | val;
#ifdef TRACK_INSERTS_SYM
          hash_ops += i -  hash_begins[hash] + 1;
#endif
          return;
        }
      }
#ifdef TRACK_INSERTS_SYM
      hash_ops += i -  hash_begins[hash] + 1;
#endif
      keys[i] = key;
      vals[i] = val;
      hash_ends[hash] = i + 1;
    }

    // Provide the shared memory capacity.
    // This function takes the team_size as an argument ,
    // which allows team_size dependent allocations.
    size_t team_shmem_size (int team_size) const {
      return requiredMem;
    }


  };
#endif


  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t>
  struct symbolic_sorted{
    size_type m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
    c_row_view_t compressedrowmapC;
    c_nnz_view_t entriesSetIndicesC;
    c_nnz_view_t entriesSetsC;
    typedef typename a_row_view_t::non_const_value_type c_row_lno_t;
    typedef typename b_nnz_view_t::non_const_value_type c_nnz_lno_t;

    const c_nnz_lno_t max_a_nnz;
    const size_t requiredMem;
    const size_t multicore_chunk_size;

    symbolic_sorted(
        size_type m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,
        c_row_view_t rowmapC_,
        c_nnz_lno_t a_max_nnz,
        size_t multicore_chunk_size_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_mapB(row_mapB_), entriesSetIndicesB(entriesSetIndicesB_), entriesSetsB(entriesSetsB_),
          rowmapC(rowmapC_),
          max_a_nnz(a_max_nnz),
          requiredMem(a_max_nnz * 5 * sizeof(c_nnz_lno_t) + 40 ),
          multicore_chunk_size(multicore_chunk_size_),
          compressedrowmapC()
          {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const FillExtraTag&, const team_member_t & teamMember) const {
      size_type row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);
      //printf("myRowMin:%d rowMax:%d\n", row_index, row_max);


      const size_t heap_size = max_a_nnz * sizeof(c_nnz_lno_t);

      c_nnz_lno_t *heap_col_set_index_vals = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(heap_size);
      c_nnz_lno_t *heap_col_sets_vals = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(heap_size);
      c_nnz_lno_t *heap_brow_index_vals = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(heap_size);
      c_nnz_lno_t *b_row_map_begin = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(heap_size);
      c_nnz_lno_t *b_row_map_end = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(heap_size);

#ifdef HASHTRACK
      size_t numhashes = 0;
      size_t num_del = 0;
      size_t num_comp = 0;
      size_t num_add = 0;
#endif

      for (; row_index < row_max; ++row_index){

        c_nnz_lno_t col_begin = row_mapA[row_index];
        const c_nnz_lno_t col_end = row_mapA[row_index + 1];
        //const c_nnz_lno_t num_cols = col_end - col_begin;

        //std::cout << "row:" << row_index << " num_cols :" << num_cols  << std::endl;
        c_nnz_lno_t heap_size = 0;

        for (; col_begin < col_end; ++col_begin){
          const c_nnz_lno_t rowb = entriesA[col_begin];

          c_nnz_lno_t b_adj = row_mapB[rowb];
          c_nnz_lno_t b_adj_end = row_mapB[rowb + 1];

          c_nnz_lno_t col_adj_ind = col_begin - row_mapA[row_index];
          b_row_map_begin[col_adj_ind] = b_adj + 1;
          b_row_map_end[col_adj_ind] = b_adj_end;
          c_nnz_lno_t inserted_index = heap_size++;
          if (b_adj < b_adj_end){

            //std::cout << "\tinserting " << entriesSetIndicesB[b_adj] << " to heap pos:" << inserted_index << std::endl;
            heap_col_set_index_vals[inserted_index] = entriesSetIndicesB[b_adj];
            //heap_col_sets_vals[inserted_index] = entriesSetsB[b_adj];
            heap_brow_index_vals[inserted_index] = col_adj_ind;

            c_nnz_lno_t parent = (inserted_index - 1) / 2;
            //compoare with parent
            while (heap_col_set_index_vals[parent] > heap_col_set_index_vals[inserted_index]){
#ifdef HASHTRACK
              ++num_add;
              ++num_comp;
#endif

              c_nnz_lno_t tmp = heap_col_set_index_vals[inserted_index];
              heap_col_set_index_vals[inserted_index] = heap_col_set_index_vals[parent];
              heap_col_set_index_vals[parent] = tmp;

              //tmp = heap_col_sets_vals[inserted_index];
              //heap_col_sets_vals[inserted_index] = heap_col_sets_vals[parent];
              //heap_col_sets_vals[parent] = tmp;

              tmp = heap_brow_index_vals[inserted_index];
              heap_brow_index_vals[inserted_index] = heap_brow_index_vals[parent];
              heap_brow_index_vals[parent] = tmp;

              inserted_index = parent;
              parent = (inserted_index - 1) / 2;
            }
          }
        }

        c_nnz_lno_t my_compressed_nonzeros = 0;
        c_nnz_lno_t my_uncompressed_nonzeros = 0;
        c_nnz_lno_t my_latest_set_ind = heap_col_set_index_vals[0];
        c_nnz_lno_t my_latest_set_val = 0;

        while (heap_size > 0){

          c_nnz_lno_t top_val = heap_col_set_index_vals[0];
          c_nnz_lno_t top_row = heap_brow_index_vals[0];
          c_nnz_lno_t top_set = entriesSetsB[b_row_map_begin[top_row] - 1];

          if (top_val == my_latest_set_ind) {
            my_latest_set_val = my_latest_set_val | top_set;
          }
          else {
            ++my_compressed_nonzeros;
            for (; my_latest_set_val; my_uncompressed_nonzeros++)
            {
              my_latest_set_val &= my_latest_set_val - 1; // clear the least significant bit set
            }

            my_latest_set_ind = top_val;
            my_latest_set_val = top_set;
          }

          c_nnz_lno_t row_col_ind = b_row_map_begin[top_row]++;
          if (row_col_ind < b_row_map_end[top_row]){
            //c_nnz_lno_t inserted_index = heap_size++;
            heap_col_set_index_vals[0] = entriesSetIndicesB[row_col_ind];
            //heap_col_sets_vals[inserted_index] = entriesSetsB[row_col_ind];
            heap_brow_index_vals[0] = top_row;
          }
          else {
            heap_col_set_index_vals[0] = heap_col_set_index_vals[--heap_size];
            heap_brow_index_vals[0] = heap_brow_index_vals[heap_size];
          }
          {
            c_nnz_lno_t heap_index = 0;
            while(1){
#ifdef HASHTRACK
              ++num_del;
              ++num_comp;
#endif

              c_nnz_lno_t child1 = heap_index * 2 + 1;
              if (child1 >= heap_size)break;

              c_nnz_lno_t child2 = child1 + 1;
              c_nnz_lno_t min_child = child1;


              if (child2 < heap_size && heap_col_set_index_vals[child2] < heap_col_set_index_vals[child1]){
                min_child = child2;
              }
              if (heap_col_set_index_vals[heap_index] < heap_col_set_index_vals[min_child]) break;

              c_nnz_lno_t tmp = heap_col_set_index_vals[min_child];
              heap_col_set_index_vals[min_child] = heap_col_set_index_vals[heap_index];
              heap_col_set_index_vals[heap_index] = tmp;

              //tmp = heap_col_sets_vals[min_child];
              //heap_col_sets_vals[min_child] = heap_col_sets_vals[heap_index];
              //heap_col_sets_vals[heap_index] = tmp;

              tmp = heap_brow_index_vals[min_child];
              heap_brow_index_vals[min_child] = heap_brow_index_vals[heap_index];
              heap_brow_index_vals[heap_index] = tmp;

              heap_index = min_child;
            }
          }


#if notdefined
          heap_col_set_index_vals[0] = heap_col_set_index_vals[--heap_size];
          //heap_col_sets_vals[0] = heap_col_sets_vals[heap_size];
          heap_brow_index_vals[0] = heap_brow_index_vals[heap_size];

          c_nnz_lno_t heap_index = 0;

          //delete from heap
          while(1){
            ++num_del;
            ++num_comp;

            c_nnz_lno_t child1 = heap_index * 2 + 1;
            if (child1 >= heap_size)break;

            c_nnz_lno_t child2 = child1 + 1;
            c_nnz_lno_t min_child = child1;


            if (child2 < heap_size && heap_col_set_index_vals[child2] < heap_col_set_index_vals[child1]){
              min_child = child2;
            }
            if (heap_col_set_index_vals[heap_index] < heap_col_set_index_vals[min_child]) break;

            c_nnz_lno_t tmp = heap_col_set_index_vals[min_child];
            heap_col_set_index_vals[min_child] = heap_col_set_index_vals[heap_index];
            heap_col_set_index_vals[heap_index] = tmp;

            //tmp = heap_col_sets_vals[min_child];
            //heap_col_sets_vals[min_child] = heap_col_sets_vals[heap_index];
            //heap_col_sets_vals[heap_index] = tmp;

            tmp = heap_brow_index_vals[min_child];
            heap_brow_index_vals[min_child] = heap_brow_index_vals[heap_index];
            heap_brow_index_vals[heap_index] = tmp;

            heap_index = min_child;
          }
          //now insert heap from the row with min val
          c_nnz_lno_t row_col_ind = b_row_map_begin[top_row]++;
          if (row_col_ind < b_row_map_end[top_row]){
            c_nnz_lno_t inserted_index = heap_size++;
            heap_col_set_index_vals[inserted_index] = entriesSetIndicesB[row_col_ind];
            //heap_col_sets_vals[inserted_index] = entriesSetsB[row_col_ind];
            heap_brow_index_vals[inserted_index] = top_row;

            c_nnz_lno_t parent = (inserted_index - 1) / 2;
            //compoare with parent
            while (heap_col_set_index_vals[parent] > heap_col_set_index_vals[inserted_index]){
              ++num_add;
              ++num_comp;

              c_nnz_lno_t tmp = heap_col_set_index_vals[inserted_index];
              heap_col_set_index_vals[inserted_index] = heap_col_set_index_vals[parent];
              heap_col_set_index_vals[parent] = tmp;

              //tmp = heap_col_sets_vals[inserted_index];
              //heap_col_sets_vals[inserted_index] = heap_col_sets_vals[parent];
              //heap_col_sets_vals[parent] = tmp;

              tmp = heap_brow_index_vals[inserted_index];
              heap_brow_index_vals[inserted_index] = heap_brow_index_vals[parent];
              heap_brow_index_vals[parent] = tmp;

              inserted_index = parent;
              parent = (inserted_index - 1) / 2;
            }
          }
#endif
        }

        ++my_compressed_nonzeros;
        for (; my_latest_set_val; my_uncompressed_nonzeros++)
        {
          my_latest_set_val &= my_latest_set_val - 1; // clear the least significant bit set
        }


        /*
        std::cout << "row:" << row_index <<
                      " my_uncompressed_nonzeros:" << my_uncompressed_nonzeros <<
                      " my_compressed_nonzeros:" << my_compressed_nonzeros << std::endl;
        */
        rowmapC(row_index) = my_uncompressed_nonzeros;
        compressedrowmapC(row_index) = my_compressed_nonzeros;
      }
#ifdef HASHTRACK
      std::cout << "num_comp:"<< num_comp << " num_del:" << num_del << " num_add:" << num_add <<  std::endl;
#endif


    }

    // Provide the shared memory capacity.
    // This function takes the team_size as an argument ,
    // which allows team_size dependent allocations.
    size_t team_shmem_size (int team_size) const {
      return requiredMem;
    }


  };


  template <
      typename out_row_view_t,
      typename out_nnz_view_t,
      typename in_row_map_t,
      typename in_nnz_t>
  struct unzipMatrix{
    typedef typename out_row_view_t::non_const_type non_const_c_lno_row_view_t;
    size_type m;
    in_row_map_t in_row_map;
    in_nnz_t in_set_index_entries;
    in_nnz_t in_set_entries;
    out_row_view_t out_rowmap;
    out_nnz_view_t out_entries;
    int set_size;
    size_t shared_memory_size;

    unzipMatrix(size_type m_,
                in_row_map_t row_map_c_copy_,
                in_nnz_t c_set_index_entries_,
                in_nnz_t c_set_entries_,
                out_row_view_t rowmapC_,
                out_row_view_t entriesC_, size_t shared_memory_size_): m(m_), in_row_map(row_map_c_copy_),
                    in_set_index_entries(c_set_index_entries_), in_set_entries(c_set_entries_),
                    out_rowmap(rowmapC_), out_entries(entriesC_), set_size(sizeof(typename in_nnz_t::value_type) * 8),
                    shared_memory_size(shared_memory_size_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {
      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      size_type row_index = teamMember.league_rank()  * team_size + team_rank;

      if (row_index >= m) return;
      //printf("row_index:%d m:%d\n", row_index, m);
      //check ii is out of range. if it is, just return.
      const size_type col_begin = in_row_map[row_index];
      const size_type col_end = in_row_map[row_index + 1];



      size_type nnz = 0;

      Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
                    [&] (size_type i, size_type &nnz) {
        const size_type col_a_index = i + col_begin;
        nnz_lno_t c_rows = in_set_entries[col_a_index];
                size_type c = 0;
        while (c_rows){
          c_rows &= (c_rows-1) ;
          c++;
        }
        nnz += c;
      }, nnz);

      //out_rowmap(row_index) = nnz;
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        out_rowmap(row_index) = nnz;
      });

    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const FillTag&, const team_member_t & teamMember) const {

      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      size_type row_index = teamMember.league_rank()  * team_size + team_rank;


      size_type *adj_index = (size_type *) teamMember.team_shmem().get_shmem(team_size);



      if (row_index >= m) return;

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        adj_index[team_rank] = out_rowmap(row_index);
      });
      //check ii is out of range. if it is, just return.
      const size_type col_begin = in_row_map[row_index];
      const size_type col_end = in_row_map[row_index + 1];
      size_type c = 0;
      //size_type out_row_map_index = out_rowmap(row_index);

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
          [&] (size_type i) {
        const size_type col_a_index = i + col_begin;

        nnz_lno_t c_rows = in_set_entries[col_a_index];
        nnz_lno_t c_rows_set_index = in_set_index_entries[col_a_index];
        int current_row = 0;
        nnz_lno_t unit = 1;


        while (c_rows){
          if (c_rows & unit){

            size_type wind = Kokkos::atomic_fetch_add(adj_index + team_rank , 1);
            out_entries(wind) = set_size * c_rows_set_index + current_row;
          }
          current_row++;
          c_rows = c_rows & ~unit;
          unit = unit << 1;
        }

      });
    }


    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }
  };


  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
            typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
            typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t,
            typename mpool_type>
  struct NumericCMEM_CPU{
    nnz_lno_t numrows;
    nnz_lno_t numcols;

    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;
    a_scalar_view_t valuesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesB;
    b_scalar_view_t valuesB;

    c_row_view_t rowmapC;
    c_nnz_view_t entriesC;
    c_scalar_view_t valuesC;
    mpool_type memory_space;

    nnz_lno_t *pEntriesC;
    scalar_t *pVals;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
    nnz_lno_t team_work_size;


    NumericCMEM_CPU(
        nnz_lno_t m_,
        nnz_lno_t k_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,
        a_scalar_view_t valuesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesB_,
        b_scalar_view_t valuesB_,

        c_row_view_t rowmapC_,
        c_nnz_view_t entriesC_,
        c_scalar_view_t valuesC_,

        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):
          numrows(m_),
          numcols(k_),
          row_mapA (row_mapA_),
          entriesA(entriesA_),
          valuesA(valuesA_),

          row_mapB(row_mapB_),
          entriesB(entriesB_),
          valuesB(valuesB_),

          rowmapC(rowmapC_),
          entriesC(entriesC_),
          valuesC(valuesC_), memory_space(), pEntriesC(entriesC_.ptr_on_device()), pVals(valuesC.ptr_on_device()),
          my_exec_space(my_exec_space_),
          team_work_size(){
          }


    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {

      nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

      scalar_t * dense_accum= NULL;
      size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
      while (dense_accum == NULL){
        dense_accum = (scalar_t * )( memory_space.allocate_chunk(tid));
      }
      char *marker = (char *) (dense_accum + numcols);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {

        const size_type c_row_begin = rowmapC[row_index];
        nnz_lno_t *myentries = pEntriesC + c_row_begin;
        scalar_t *myvals = pVals + c_row_begin;

        nnz_lno_t current_col_index = 0;
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t nnza = nnz_lno_t(row_mapA[row_index + 1] - col_begin);

        for (nnz_lno_t colind = 0; colind < nnza; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;
          for (int i = 0; i < left_work; ++i){
            const size_type adjind = i + rowBegin;
            nnz_lno_t b_col_ind = entriesB[adjind];
            scalar_t b_val = valuesB[adjind] * valA;
            if (marker[b_col_ind] == 0){
              marker[b_col_ind] = 1;
              myentries[current_col_index++] = b_col_ind;
            }
            dense_accum[b_col_ind] += b_val;
          }
        }
        for (nnz_lno_t i = 0; i < current_col_index; ++i){
          nnz_lno_t ind = myentries[i];
          myvals[i] = dense_accum[ind];
          dense_accum[ind] = 0;
          marker [ind] = 0;
        }
      });

    }

  };

  template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__,
            typename b_row_view_t__, typename b_nnz_view_t__, typename b_scalar_view_t__,
            typename c_row_view_t__, typename c_nnz_view_t__, typename c_scalar_view_t__>
  struct NumericCMEM{
    nnz_lno_t numrows;

    a_row_view_t__ row_mapA;
    a_nnz_view_t__ entriesA;
    a_scalar_view_t__ valuesA;

    b_row_view_t__ row_mapB;
    b_nnz_view_t__ entriesB;
    b_scalar_view_t__ valuesB;

    c_row_view_t__ rowmapC;
    c_nnz_view_t__ entriesC;
    c_scalar_view_t__ valuesC;

    c_nnz_view_t__ beginsC;
    c_nnz_view_t__ nextsC;

    nnz_lno_t *pbeginsC, *pnextsC, *pEntriesC;
    scalar_t *pvaluesC;

    const size_t shared_memory_size;
    int vector_size;
    nnz_lno_t team_work_size;

    NumericCMEM(
        nnz_lno_t m_,
        a_row_view_t__ row_mapA_,
        a_nnz_view_t__ entriesA_,
        a_scalar_view_t__ valuesA_,

        b_row_view_t__ row_mapB_,
        b_nnz_view_t__ entriesB_,
        b_scalar_view_t__ valuesB_,

        c_row_view_t__ rowmapC_,
        c_nnz_view_t__ entriesC_,
        c_scalar_view_t__ valuesC_,

        c_nnz_view_t__ beginsC_,
        c_nnz_view_t__ nextsC_,

        const size_type sharedMemorySize_):
          numrows(m_),
          row_mapA (row_mapA_),
          entriesA(entriesA_),
          valuesA(valuesA_),

          row_mapB(row_mapB_),
          entriesB(entriesB_),
          valuesB(valuesB_),

          rowmapC(rowmapC_),
          entriesC(entriesC_),
          valuesC(valuesC_),
          beginsC(beginsC_),
          nextsC(nextsC_),
          pbeginsC(beginsC_.ptr_on_device()), pnextsC(nextsC_.ptr_on_device()),
          pEntriesC(entriesC_.ptr_on_device()), pvaluesC(valuesC_.ptr_on_device()),
          shared_memory_size(sharedMemorySize_),
          vector_size (), team_work_size()
          {
          }

    KOKKOS_INLINE_FUNCTION
    void operator()(const GPUTag&, const team_member_t & teamMember) const {


      nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

      //nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      int thread_memory = (shared_memory_size /8 / teamMember.team_size()) * 8;
      /*
      if (team_row_begin == 0){
        Kokkos::single(Kokkos::PerTeam(teamMember),[&] () {
        printf("teamMember.team_size():%d thread_memory:%d\n", teamMember.team_size(), thread_memory);
      });
      }
      */

      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));


      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;
      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t) ; //begins, nexts, and keys. No need for vals yet.
      nnz_lno_t shared_memory_hash_size = (thread_memory - sizeof(nnz_lno_t) * 2) / unit_memory ;
      shared_memory_hash_size = (shared_memory_hash_size >> 1) << 1;

      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      scalar_t * vals = (scalar_t *) (all_shared_memory);


      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm2(0, 0,
          NULL, NULL, NULL, NULL);
      /*
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm2(global_memory_hash_size, global_memory_hash_size,
          pbeginsC + c_row_begin, pnextsC + c_row_begin, pEntriesC + c_row_begin, pvaluesC + c_row_begin);
          */


      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
        const size_type c_row_begin = rowmapC[row_index];
        const nnz_lno_t global_memory_hash_size = nnz_lno_t(rowmapC[row_index + 1] - c_row_begin);

        hm2.hash_key_size = global_memory_hash_size;
        hm2.max_value_size = global_memory_hash_size;
        hm2.keys = pEntriesC + c_row_begin;
        hm2.values = pvaluesC + c_row_begin;
        hm2.hash_begins = pbeginsC + c_row_begin;
        hm2.hash_nexts = pnextsC + c_row_begin;

        //initialize begins.
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
            [&] (int i) {
          begins[i] = -1;
        });

        //initialize hash usage sizes
        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          used_hash_sizes[0] = 0;
          used_hash_sizes[1] = 0;
        });

        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);

        for (nnz_lno_t colind = 0; colind < left_work; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          while (left_work){
            nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);
            nnz_lno_t b_col_ind = -1;
            scalar_t b_val = -1;
            nnz_lno_t hash = -1;
            Kokkos::parallel_for(
                Kokkos::ThreadVectorRange(teamMember, work_to_handle),
                [&] (nnz_lno_t i) {
              const size_type adjind = i + rowBegin;
              b_col_ind = entriesB[adjind];
              b_val = valuesB[adjind] * valA;
              hash = b_col_ind % shared_memory_hash_size;
            });

            int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeAdd(
                teamMember, vector_size,
                hash, b_col_ind, b_val,
                used_hash_sizes,
                shared_memory_hash_size);

            int overall_num_unsuccess = 0;

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
                [&] (const int threadid, int &overall_num_unsuccess) {
              overall_num_unsuccess += num_unsuccess;
            }, overall_num_unsuccess);

            if (overall_num_unsuccess){
              nnz_lno_t hash = -1;
              if (num_unsuccess) {
                hash = b_col_ind % global_memory_hash_size;
              }



              int insertion = hm2.vector_atomic_insert_into_hash_mergeAdd(
                  teamMember, vector_size,
                  hash,b_col_ind,b_val,
                  used_hash_sizes + 1, hm2.max_value_size
              );

            }
            left_work -= work_to_handle;
            rowBegin += work_to_handle;
          }
        }

        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
        });


        size_type num_elements = used_hash_sizes[0];


        size_type written_index = used_hash_sizes[1];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, num_elements),
            [&] (size_type i) {
          pEntriesC[c_row_begin + written_index + i] = keys[i];
          pvaluesC[c_row_begin + written_index + i] = vals[i];
        });
      });
    }

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }
  };


  template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__,
            typename b_row_view_t__, typename b_nnz_view_t__, typename b_scalar_view_t__,
            typename c_row_view_t__, typename c_nnz_view_t__, typename c_scalar_view_t__>
  struct NumericCCOLOR{
    nnz_lno_t numrows;
    nnz_lno_t numcols;

    a_row_view_t__ row_mapA;
    a_nnz_view_t__ entriesA;
    a_scalar_view_t__ valuesA;

    b_row_view_t__ row_mapB;
    b_nnz_view_t__ entriesB;
    b_scalar_view_t__ valuesB;

    c_row_view_t__ rowmapC;
    c_nnz_view_t__ entriesC;
    c_scalar_view_t__ valuesC;
    nnz_lno_t *pEntriesC ;
    scalar_t *pVals;

    scalar_temp_work_view_t denseAccumulator; //initially all zeroes
    scalar_t *pdenseAccumulator; //initially all zeroes


    //bool_temp_view_t denseAccumulatorFlags; //initially all false.
    //bool * pdenseAccumulatorFlags; //initially all false.

    nnz_lno_t team_work_size;

    nnz_lno_t color_begin;
    nnz_lno_t color_end;
    nnz_lno_persistent_work_view_t color_adj;
    nnz_lno_persistent_work_view_t vertex_colors;

    nnz_lno_t consecutive_chunk_size;
    nnz_lno_t consecutive_all_color_chunk_size;

    nnz_lno_t chunk_divison;
    nnz_lno_t chunk_and;
    NumericCCOLOR(

        nnz_lno_t m_,
        nnz_lno_t k_,


        a_row_view_t__ row_mapA_,
        a_nnz_view_t__ entriesA_,
        a_scalar_view_t__ valuesA_,

        b_row_view_t__ row_mapB_,
        b_nnz_view_t__ entriesB_,
        b_scalar_view_t__ valuesB_,

        c_row_view_t__ rowmapC_,
        c_nnz_view_t__ entriesC_,
        c_scalar_view_t__ valuesC_,
        scalar_temp_work_view_t denseAccumulator_, //initially all zeroes
        //bool_temp_view_t denseAccumulatorFlags_, //initially all false.



        nnz_lno_t team_row_work_size_):
          numrows(m_), numcols(k_),
          row_mapA (row_mapA_),
          entriesA(entriesA_),
          valuesA(valuesA_),

          row_mapB(row_mapB_),
          entriesB(entriesB_),
          valuesB(valuesB_),

          rowmapC(rowmapC_),
          entriesC(entriesC_),
          valuesC(valuesC_), pEntriesC(entriesC_.ptr_on_device()), pVals(valuesC.ptr_on_device()),
          denseAccumulator (denseAccumulator_), pdenseAccumulator(denseAccumulator_.ptr_on_device()),
          //denseAccumulatorFlags (denseAccumulatorFlags_), pdenseAccumulatorFlags(denseAccumulatorFlags_.ptr_on_device()),
          team_work_size(team_row_work_size_),
          color_begin(0),
          color_end(0),
          color_adj(),
          vertex_colors(), consecutive_chunk_size(numcols), consecutive_all_color_chunk_size(numcols), chunk_divison (0), chunk_and(0)
          {
          }

    //one color at a time.
    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric1Tag&, const team_member_t & teamMember) const {
      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& color_index_index) {
        nnz_lno_t row_index = color_adj(color_index_index);

        //nnz_lno_t color = vertex_colors(row_index);
        scalar_t *mydenseAccumulator = pdenseAccumulator;// + numcols * color;
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
        const size_type c_row_begin = rowmapC[row_index];
        nnz_lno_t *my_entries = pEntriesC + c_row_begin;
        for (nnz_lno_t colind = 0; colind < left_work; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, left_work),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            nnz_lno_t acc_index = entriesB[adjind];
            //nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);

            scalar_t b_val = valuesB[adjind] * valA;
            mydenseAccumulator[acc_index] += b_val;
          });
        }
        scalar_t *my_vals = pVals + c_row_begin;

        nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, row_size),
            [&] (nnz_lno_t i) {
          nnz_lno_t acc_index = my_entries[i];
          //nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);
          my_vals[i] = mydenseAccumulator[acc_index];
          mydenseAccumulator[acc_index] = 0;
        });
      });
    }


    //multi-color minimized writes
    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric2Tag&, const team_member_t & teamMember) const {


      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& color_index_index) {
        nnz_lno_t row_index = color_adj(color_index_index);

        nnz_lno_t color = vertex_colors(row_index);
        scalar_t *mydenseAccumulator = pdenseAccumulator + numcols * color;

        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
        const size_type c_row_begin = rowmapC[row_index];
        nnz_lno_t *my_entries = pEntriesC + c_row_begin;
        for (nnz_lno_t colind = 0; colind < left_work; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, left_work),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            nnz_lno_t acc_index = entriesB[adjind];
            //nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);
            scalar_t b_val = valuesB[adjind] * valA;
            mydenseAccumulator[acc_index] += b_val;
          });
        }
        scalar_t *my_vals = pVals + c_row_begin;

        nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, row_size),
            [&] (nnz_lno_t i) {
          nnz_lno_t acc_index = my_entries[i];
          //nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);
          my_vals[i] = mydenseAccumulator[acc_index];
          mydenseAccumulator[acc_index] = 0;
        });
      });
    }

    //multi-color minimized reads
    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric3Tag&, const team_member_t & teamMember) const {


      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& color_index_index) {
        nnz_lno_t row_index = color_adj(color_index_index);

        nnz_lno_t color = vertex_colors(row_index);
        scalar_t *mydenseAccumulator = pdenseAccumulator;// + numcols * color;
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
        const size_type c_row_begin = rowmapC[row_index];
        nnz_lno_t *my_entries = pEntriesC + c_row_begin;
        for (nnz_lno_t colind = 0; colind < left_work; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, left_work),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            nnz_lno_t col_ind = entriesB[adjind];
            nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);
            scalar_t b_val = valuesB[adjind] * valA;
            mydenseAccumulator[acc_index] += b_val;
          });
        }
        scalar_t *my_vals = pVals + c_row_begin;

        nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, row_size),
            [&] (nnz_lno_t i) {
          nnz_lno_t col_ind = my_entries[i];
          nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);
          my_vals[i] = mydenseAccumulator[acc_index];
          mydenseAccumulator[acc_index] = 0;
        });
      });
    }


    size_t team_shmem_size (int team_size) const {
      return team_size * sizeof (nnz_lno_t) * 8;
    }
  };


  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
            typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
            typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t,
            typename pool_memory_type>
  struct PortableNumericCHASH{

    nnz_lno_t numrows;

    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;
    a_scalar_view_t valuesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesB;
    b_scalar_view_t valuesB;

    c_row_view_t rowmapC;
    c_nnz_view_t entriesC;
    c_scalar_view_t valuesC;


    nnz_lno_t *pEntriesC;
    scalar_t *pvaluesC;
    size_t shared_memory_size;
    int vector_size;
    pool_memory_type memory_space;

    nnz_lno_t max_nnz;
    nnz_lno_t pow2_hash_size;
    nnz_lno_t pow2_hash_func;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
    nnz_lno_t team_work_size;

    PortableNumericCHASH(
        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,
        a_scalar_view_t valuesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesB_,
        b_scalar_view_t valuesB_,

        c_row_view_t rowmapC_,
        c_nnz_view_t entriesC_,
        c_scalar_view_t valuesC_,

        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):
          numrows(m_),
          row_mapA (row_mapA_),
          entriesA(entriesA_),
          valuesA(valuesA_),

          row_mapB(row_mapB_),
          entriesB(entriesB_),
          valuesB(valuesB_),

          rowmapC(rowmapC_),
          entriesC(entriesC_),
          valuesC(valuesC_),
          pEntriesC(entriesC_.ptr_on_device()), pvaluesC(valuesC_.ptr_on_device()),
          shared_memory_size(0), vector_size (),
          memory_space(),
          max_nnz(),
          pow2_hash_size(),
          pow2_hash_func(),
          my_exec_space(my_exec_space_),
          team_work_size(1)
          {
          }
    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }

    }

    //assumes that the vector lane is 1, as in cpus
    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {

      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm2(pow2_hash_size, pow2_hash_size,NULL, NULL, NULL, NULL);

      volatile nnz_lno_t * tmp = NULL;
      size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
      while (tmp == NULL){
        tmp = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
      }
      nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *) tmp;
      tmp += pow2_hash_size ;

      hm2.hash_begins = (nnz_lno_t *) (tmp);
      tmp += pow2_hash_size;
      hm2.hash_nexts = (nnz_lno_t *) (tmp);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
        nnz_lno_t globally_used_hash_count = 0;
        nnz_lno_t used_hash_sizes = 0;

        const size_type c_row_begin = rowmapC[row_index];
        const size_type c_row_end = rowmapC[row_index + 1];

        const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);
        hm2.max_value_size = global_memory_hash_size;
        hm2.keys = pEntriesC + c_row_begin;
        hm2.values = pvaluesC + c_row_begin;

        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
        for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
          size_type a_col = col_begin + ii;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_workB = row_mapB(rowB + 1) - rowBegin;

          for ( nnz_lno_t i = 0; i < left_workB; ++i){
            const size_type adjind = i + rowBegin;
            nnz_lno_t b_col_ind = entriesB[adjind];
            scalar_t b_val = valuesB[adjind] * valA;
            nnz_lno_t hash = b_col_ind & pow2_hash_func;

            //this has to be a success, we do not need to check for the success.
            int insertion = hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(
                hash, b_col_ind, b_val,
                &used_hash_sizes, hm2.max_value_size
                ,&globally_used_hash_count,
                globally_used_hash_indices
            );
          }
        }
        for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i){
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        }
      });
      memory_space.release_chunk(globally_used_hash_indices);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric2Tag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;

      nnz_lno_t globally_used_hash_count = 0;
      nnz_lno_t used_hash_sizes = 0;


      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];

      const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
        hm2(pow2_hash_size, global_memory_hash_size,
            NULL, NULL, pEntriesC + c_row_begin, pvaluesC + c_row_begin);

      volatile nnz_lno_t * tmp = NULL;
      size_t tid = get_thread_id(row_index);
      while (tmp == NULL){
        tmp = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
      }

      nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *) tmp;
      tmp += pow2_hash_size ;

      hm2.hash_begins = (nnz_lno_t *) (tmp);
      tmp += pow2_hash_size;
      hm2.hash_nexts = (nnz_lno_t *) (tmp);


      {
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
        for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
          size_type a_col = col_begin + ii;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          for ( nnz_lno_t i = 0; i < left_work; ++i){
            const size_type adjind = i + rowBegin;
            nnz_lno_t b_col_ind = entriesB[adjind];
            scalar_t b_val = valuesB[adjind] * valA;
            nnz_lno_t hash = b_col_ind & pow2_hash_func;

            //this has to be a success, we do not need to check for the success.
            int insertion = hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(
                hash,b_col_ind,b_val,
                &used_hash_sizes, hm2.max_value_size
                ,&globally_used_hash_count, globally_used_hash_indices
            );
          }
        }
        for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i){
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        }
      }
      memory_space.release_chunk(globally_used_hash_indices);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const GPUTag&, const team_member_t & teamMember) const {

      nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);


      int thread_memory = (shared_memory_size / 8 / teamMember.team_size()) * 8;
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      nnz_lno_t *globally_used_hash_count = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t) ; //begins, nexts, keys and vals .
      nnz_lno_t shared_memory_hash_size = (thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory;
      if (shared_memory_hash_size & 1) shared_memory_hash_size -= 1;

      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      scalar_t * vals = (scalar_t *) (all_shared_memory);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm2(pow2_hash_size, pow2_hash_size,
          NULL, NULL, NULL, NULL);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
        const size_type c_row_begin = rowmapC[row_index];
        const size_type c_row_end = rowmapC[row_index + 1];

        const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);
        hm2.max_value_size = global_memory_hash_size;
        hm2.keys = pEntriesC + c_row_begin;
        hm2.values = pvaluesC + c_row_begin;

        //initialize begins.
        Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size), [&] (nnz_lno_t i) {
            begins[i] = -1; });

        //initialize hash usage sizes
        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          used_hash_sizes[0] = 0;
          used_hash_sizes[1] = 0;
          globally_used_hash_count[0] = 0;
        });

        bool is_global_alloced = false;
        nnz_lno_t *globally_used_hash_indices = NULL;

        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;

        for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
          size_type a_col = col_begin + ii;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          while (left_work){
            nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);
            nnz_lno_t b_col_ind = -1;
            scalar_t b_val = -1;
            nnz_lno_t hash = -1;
            Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, work_to_handle),
                [&] (nnz_lno_t i) {
              const size_type adjind = i + rowBegin;
              b_col_ind = entriesB[adjind];
              b_val = valuesB[adjind] * valA;
              hash = b_col_ind % shared_memory_hash_size;
            });

            int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeAdd(
                teamMember, vector_size,
                hash, b_col_ind, b_val,
                used_hash_sizes,
                shared_memory_hash_size);

            int overall_num_unsuccess = 0;

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
                [&] (const int threadid, int &overall_num_unsuccess) {
              overall_num_unsuccess += num_unsuccess;
            }, overall_num_unsuccess);

            if (overall_num_unsuccess){
              if (!is_global_alloced){
                volatile nnz_lno_t * tmp = NULL;
                size_t tid = get_thread_id(row_index);
                while (tmp == NULL){
                  Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile nnz_lno_t * &memptr) {
                    memptr = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
                  }, tmp);
                }
                is_global_alloced = true;
                globally_used_hash_indices = (nnz_lno_t *) tmp;
                tmp += pow2_hash_size ;

                hm2.hash_begins = (nnz_lno_t *) (tmp);
                tmp += pow2_hash_size ;
                hm2.hash_nexts = (nnz_lno_t *) (tmp);
              }

              nnz_lno_t hash = -1;
              if (num_unsuccess) {
                hash = b_col_ind & pow2_hash_func;
              }

              //this has to be a success, we do not need to check for the success.
              int insertion = hm2.vector_atomic_insert_into_hash_mergeAdd_TrackHashes(
                  teamMember, vector_size,
                  hash,b_col_ind,b_val,
                  used_hash_sizes + 1, hm2.max_value_size
                  ,globally_used_hash_count, globally_used_hash_indices
              );
            }
            left_work -= work_to_handle;
            rowBegin += work_to_handle;
          }
        }

        if (is_global_alloced){

          nnz_lno_t dirty_hashes = globally_used_hash_count[0];
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, dirty_hashes),
              [&] (nnz_lno_t i) {
            nnz_lno_t dirty_hash = globally_used_hash_indices[i];
            hm2.hash_begins[dirty_hash] = -1;
          });

          Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
            memory_space.release_chunk(globally_used_hash_indices);
          });
        }

        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
        });

        nnz_lno_t num_elements = used_hash_sizes[0];

        nnz_lno_t written_index = used_hash_sizes[1];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, num_elements),
            [&] (nnz_lno_t i) {
          pEntriesC[c_row_begin + written_index + i] = keys[i];
          pvaluesC[c_row_begin + written_index + i] = vals[i];
        });
      });
    }
    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }
  };

  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
            typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
            typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t,
            typename pool_memory_type>
  struct NumericCHASH{
    nnz_lno_t numrows;

    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;
    a_scalar_view_t valuesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesB;
    b_scalar_view_t valuesB;

    c_row_view_t rowmapC;
    c_nnz_view_t entriesC;
    c_scalar_view_t valuesC;


    nnz_lno_t *pEntriesC;
    scalar_t *pvaluesC;

    const size_t shared_memory_size;
    int vector_size;
    pool_memory_type memory_space;

    nnz_lno_t max_nnz;
    nnz_lno_t pow2_hash_size;
    nnz_lno_t pow2_hash_func;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;


    NumericCHASH(
        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,
        a_scalar_view_t valuesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesB_,
        b_scalar_view_t valuesB_,

        c_row_view_t rowmapC_,
        c_nnz_view_t entriesC_,
        c_scalar_view_t valuesC_,

        const size_type sharedMemorySize_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):
          numrows(m_),
          row_mapA (row_mapA_),
          entriesA(entriesA_),
          valuesA(valuesA_),

          row_mapB(row_mapB_),
          entriesB(entriesB_),
          valuesB(valuesB_),

          rowmapC(rowmapC_),
          entriesC(entriesC_),
          valuesC(valuesC_),
          pEntriesC(entriesC_.ptr_on_device()), pvaluesC(valuesC_.ptr_on_device()),
          shared_memory_size(sharedMemorySize_),
          vector_size (), memory_space(),
          max_nnz(),
          pow2_hash_size(),
          pow2_hash_func(),
          my_exec_space(my_exec_space_)
          {
          }
    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }

    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric2Tag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;

      int thread_memory = (shared_memory_size /8 / teamMember.team_size()) * 8;
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));


      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      nnz_lno_t *globally_used_hash_count = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t) ; //begins, nexts, keys and vals .
      nnz_lno_t shared_memory_hash_size = (thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory;
      if (shared_memory_hash_size & 1) shared_memory_hash_size -= 1;

      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      scalar_t * vals = (scalar_t *) (all_shared_memory);

      //printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts, keys, vals);
      //return;
      //first level hashmap


      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];

      const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
        hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
        hm2(pow2_hash_size, global_memory_hash_size,
            NULL, NULL, pEntriesC + c_row_begin, pvaluesC + c_row_begin);

      //initialize begins.
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
          [&] (nnz_lno_t i) {
        begins[i] = -1;
      });
      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
        globally_used_hash_count[0] = 0;
      });
      bool is_global_alloced = false;
      nnz_lno_t *globally_used_hash_indices = NULL;

      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
      for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
        size_type a_col = col_begin + ii;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];

        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

        while (left_work){
          nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);


          nnz_lno_t b_col_ind = -1;
          scalar_t b_val = -1;
          nnz_lno_t hash = -1;
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, work_to_handle),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            b_col_ind = entriesB[adjind];
            b_val = valuesB[adjind] * valA;
            hash = b_col_ind % shared_memory_hash_size;
          });

          int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeAdd(
              teamMember, vector_size,
              hash, b_col_ind, b_val,
              used_hash_sizes,
              shared_memory_hash_size);



          int overall_num_unsuccess = 0;

          Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
              [&] (const int threadid, int &overall_num_unsuccess) {
            overall_num_unsuccess += num_unsuccess;
          }, overall_num_unsuccess);


          if (overall_num_unsuccess){


            if (!is_global_alloced){
              volatile nnz_lno_t * tmp = NULL;
              size_t tid = get_thread_id(row_index);
              while (tmp == NULL){
                Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile nnz_lno_t * &memptr) {
                  memptr = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
                }, tmp);
              }
              is_global_alloced = true;
              globally_used_hash_indices = (nnz_lno_t *) tmp;
              tmp += pow2_hash_size ;

              hm2.hash_begins = (nnz_lno_t *) (tmp);
              tmp += pow2_hash_size ;
              hm2.hash_nexts = (nnz_lno_t *) (tmp);
            }

            nnz_lno_t hash = -1;
            if (num_unsuccess) {
              hash = b_col_ind & pow2_hash_func;
            }

            //this has to be a success, we do not need to check for the success.
            int insertion = hm2.vector_atomic_insert_into_hash_mergeAdd_TrackHashes(
                teamMember, vector_size,
                hash,b_col_ind,b_val,
                used_hash_sizes + 1, hm2.max_value_size
                ,globally_used_hash_count, globally_used_hash_indices
                );
          }
          left_work -= work_to_handle;
          rowBegin += work_to_handle;
        }
      }

      if (is_global_alloced){

        nnz_lno_t dirty_hashes = globally_used_hash_count[0];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, dirty_hashes),
            [&] (nnz_lno_t i) {
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        });

        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          memory_space.release_chunk(globally_used_hash_indices);
        });
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
      });

      nnz_lno_t num_elements = used_hash_sizes[0];

      nnz_lno_t written_index = used_hash_sizes[1];
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, num_elements),
          [&] (nnz_lno_t i) {
        pEntriesC[c_row_begin + written_index + i] = keys[i];
        pvaluesC[c_row_begin + written_index + i] = vals[i];
      });
    }
//#endif
    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }
  };




  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t,
            typename pool_memory_space>
  struct NonzeroesC_CPU{
    nnz_lno_t numrows;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
    c_nnz_view_t entriesSetIndicesC;
    c_nnz_view_t entriesSetsC;


    const nnz_lno_t pow2_hash_size;
    const nnz_lno_t pow2_hash_func;
    const nnz_lno_t MaxRoughNonZero;

    pool_memory_space m_space;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
    NonzeroesC_CPU(

        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,

        c_row_view_t rowmapC_,
        c_nnz_view_t entriesSetIndicesC_,

        const nnz_lno_t hash_size_,
        const nnz_lno_t MaxRoughNonZero_,
        const size_t sharedMemorySize_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):

          numrows(m_),

          row_mapA (row_mapA_),
          entriesA(entriesA_),

          row_mapB(row_mapB_),
          entriesSetIndicesB(entriesSetIndicesB_),
          entriesSetsB(entriesSetsB_),

          rowmapC(rowmapC_),
          entriesSetIndicesC(entriesSetIndicesC_),
          entriesSetsC(),


          pow2_hash_size(hash_size_),
          pow2_hash_func(hash_size_ - 1),
          MaxRoughNonZero(MaxRoughNonZero_),
          m_space(),
          my_exec_space(my_exec_space_)
          {}

    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;


      //printf("row:%d\n", row_index);

      nnz_lno_t *globally_used_hash_indices = NULL;
      nnz_lno_t globally_used_hash_count = 0;
      nnz_lno_t used_hash_size = 0;
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;



      volatile nnz_lno_t * tmp = NULL;
      size_t tid = get_thread_id(row_index);
      while (tmp == NULL){
        tmp = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
      }



      globally_used_hash_indices = (nnz_lno_t *) tmp;
      tmp += pow2_hash_size ;

      hm2.hash_begins = (nnz_lno_t *) (tmp);
      tmp += pow2_hash_size ;

      //poins to the next elements
      hm2.hash_nexts = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;

      //holds the keys
      hm2.keys = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;
      hm2.values = (nnz_lno_t *) (tmp);

      hm2.hash_key_size = pow2_hash_size;
      hm2.max_value_size = MaxRoughNonZero;


      {
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t col_size = row_mapA[row_index + 1] - col_begin;

        for (nnz_lno_t colind = 0; colind < col_size; ++colind){
          size_type a_col = colind + col_begin;

          nnz_lno_t rowB = entriesA[a_col];
          size_type rowBegin = row_mapB(rowB);

          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          for (nnz_lno_t i = 0; i < left_work; ++i){

            const size_type adjind = i + rowBegin;
            nnz_lno_t b_set_ind = entriesSetIndicesB[adjind];
            nnz_lno_t b_set = entriesSetsB[adjind];
            nnz_lno_t hash = b_set_ind & pow2_hash_func;


            int insertion = hm2.sequential_insert_into_hash_mergeOr_TrackHashes(
                hash,b_set_ind,b_set,
                &used_hash_size, hm2.max_value_size
                ,&globally_used_hash_count, globally_used_hash_indices
            );
          }
        }




        int set_size = sizeof(nnz_lno_t) * 8;
        nnz_lno_t num_el = rowmapC(row_index);
        for (nnz_lno_t ii = 0; ii < used_hash_size; ++ii){
          nnz_lno_t c_rows_setind = hm2.keys[ii];
          nnz_lno_t c_rows = hm2.values[ii];

          int current_row = 0;
          nnz_lno_t unit = 1;

          while (c_rows){
            if (c_rows & unit){
              //std::cout << "num_el:" << num_el << " set_size * c_rows_setind + current_row:" << set_size * c_rows_setind + current_row << std::endl;
              entriesSetIndicesC(num_el++) = set_size * c_rows_setind + current_row;
            }
            current_row++;
            c_rows = c_rows & ~unit;
            unit = unit << 1;
          }
        }
        for (int i = 0; i < globally_used_hash_count; ++i){
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        }

      }

      m_space.release_chunk(globally_used_hash_indices);
    }


  };


  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t,
            typename pool_memory_space>
  struct StructureC_CPU{
    nnz_lno_t numrows;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
    c_nnz_view_t entriesSetIndicesC;
    c_nnz_view_t entriesSetsC;


    const nnz_lno_t pow2_hash_size;
    const nnz_lno_t pow2_hash_func;
    const nnz_lno_t MaxRoughNonZero;

    pool_memory_space m_space;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
    StructureC_CPU(

        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,

        c_row_view_t rowmapC_,

        const nnz_lno_t hash_size_,
        const nnz_lno_t MaxRoughNonZero_,
        const size_t sharedMemorySize_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):

          numrows(m_),

          row_mapA (row_mapA_),
          entriesA(entriesA_),

          row_mapB(row_mapB_),
          entriesSetIndicesB(entriesSetIndicesB_),
          entriesSetsB(entriesSetsB_),

          rowmapC(rowmapC_),
          entriesSetIndicesC(),
          entriesSetsC(),


          pow2_hash_size(hash_size_),
          pow2_hash_func(hash_size_ - 1),
          MaxRoughNonZero(MaxRoughNonZero_),
          m_space(),
          my_exec_space(my_exec_space_)
          {}

    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;


      //printf("row:%d\n", row_index);

      nnz_lno_t *globally_used_hash_indices = NULL;
      nnz_lno_t globally_used_hash_count = 0;
      nnz_lno_t used_hash_size = 0;
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;



      volatile nnz_lno_t * tmp = NULL;
      size_t tid = get_thread_id(row_index);
      while (tmp == NULL){
        tmp = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
      }



      globally_used_hash_indices = (nnz_lno_t *) tmp;
      tmp += pow2_hash_size ;

      hm2.hash_begins = (nnz_lno_t *) (tmp);
      tmp += pow2_hash_size ;

      //poins to the next elements
      hm2.hash_nexts = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;

      //holds the keys
      hm2.keys = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;
      hm2.values = (nnz_lno_t *) (tmp);

      hm2.hash_key_size = pow2_hash_size;
      hm2.max_value_size = MaxRoughNonZero;


      {
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t col_size = row_mapA[row_index + 1] - col_begin;

        for (nnz_lno_t colind = 0; colind < col_size; ++colind){
          size_type a_col = colind + col_begin;

          nnz_lno_t rowB = entriesA[a_col];
          size_type rowBegin = row_mapB(rowB);

          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          for (nnz_lno_t i = 0; i < left_work; ++i){

            const size_type adjind = i + rowBegin;
            nnz_lno_t b_set_ind = entriesSetIndicesB[adjind];
            nnz_lno_t b_set = entriesSetsB[adjind];
            nnz_lno_t hash = b_set_ind & pow2_hash_func;


            int insertion = hm2.sequential_insert_into_hash_mergeOr_TrackHashes(
                hash,b_set_ind,b_set,
                &used_hash_size, hm2.max_value_size
                ,&globally_used_hash_count, globally_used_hash_indices
            );
          }
        }




        nnz_lno_t num_el = 0;
        for (nnz_lno_t ii = 0; ii < used_hash_size; ++ii){
          nnz_lno_t c_rows = hm2.values[ii];

          for (; c_rows; num_el++) {
            c_rows &= c_rows - 1; // clear the least significant bit set
          }
        }
        for (int i = 0; i < globally_used_hash_count; ++i){
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        }

        rowmapC(row_index) = num_el;
      }

      m_space.release_chunk(globally_used_hash_indices);
    }


  };



  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t,
            typename pool_memory_space>
  struct StructureC{
    nnz_lno_t numrows;

    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
    c_nnz_view_t entriesSetIndicesC;
    c_nnz_view_t entriesSetsC;


    const nnz_lno_t pow2_hash_size;
    const nnz_lno_t pow2_hash_func;
    const nnz_lno_t MaxRoughNonZero;

    const size_t shared_memory_size;
    int vector_size;
    pool_memory_space m_space;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
    StructureC(

        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,

        c_row_view_t rowmapC_,

        const nnz_lno_t hash_size_,
        const nnz_lno_t MaxRoughNonZero_,
        const size_t sharedMemorySize_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):
          numrows(m_),

          row_mapA (row_mapA_),
          entriesA(entriesA_),

          row_mapB(row_mapB_),
          entriesSetIndicesB(entriesSetIndicesB_),
          entriesSetsB(entriesSetsB_),

          rowmapC(rowmapC_),
          entriesSetIndicesC(),
          entriesSetsC(),


          pow2_hash_size(hash_size_),
          pow2_hash_func(hash_size_ - 1),
          MaxRoughNonZero(MaxRoughNonZero_),

          shared_memory_size(sharedMemorySize_),
          vector_size (), m_space(), my_exec_space(my_exec_space_)
          {}

    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;


      //printf("row:%d\n", row_index);

      int thread_memory = ((shared_memory_size/ 4 / teamMember.team_size())) * 4;
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

      //nnz_lno_t *alloc_global_memory = NULL;
      nnz_lno_t *globally_used_hash_indices = NULL;

      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      nnz_lno_t *globally_used_hash_count = (nnz_lno_t *) (all_shared_memory);

      all_shared_memory += sizeof(nnz_lno_t) ;
      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
      nnz_lno_t shared_memory_hash_size = (thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory;

      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      nnz_lno_t * vals = (nnz_lno_t *) (all_shared_memory);

      //printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts, keys, vals);
      //return;
      //first level hashmap
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;

      //initialize begins.
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
          [&] (int i) {
        begins[i] = -1;
      });

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
        globally_used_hash_count[0] = 0;
      });

      bool is_global_alloced = false;

      const size_type col_end = row_mapA[row_index + 1];
      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t col_size = col_end - col_begin;

      for (nnz_lno_t colind = 0; colind < col_size; ++colind){
        size_type a_col = colind + col_begin;

        nnz_lno_t rowB = entriesA[a_col];
        size_type rowBegin = row_mapB(rowB);

        nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

        while (left_work){
          nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

          nnz_lno_t b_set_ind = -1, b_set = -1;
          nnz_lno_t hash = -1;
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, work_to_handle),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            b_set_ind = entriesSetIndicesB[adjind];
            b_set = entriesSetsB[adjind];
            hash = b_set_ind % shared_memory_hash_size;
          });


          int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(
              teamMember, vector_size,
              hash, b_set_ind, b_set,
              used_hash_sizes,
              shared_memory_hash_size);


          int overall_num_unsuccess = 0;

          Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
              [&] (const int threadid, int &overall_num_unsuccess) {
            overall_num_unsuccess += num_unsuccess;
          }, overall_num_unsuccess);


          if (overall_num_unsuccess){

            //printf("row:%d\n", row_index);
            if (!is_global_alloced){
              volatile nnz_lno_t * tmp = NULL;
              size_t tid = get_thread_id(row_index);
              while (tmp == NULL){
                Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile nnz_lno_t * &memptr) {
                  memptr = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
                }, tmp);
              }
              is_global_alloced = true;

              globally_used_hash_indices = (nnz_lno_t *) tmp;
              tmp += pow2_hash_size ;

              hm2.hash_begins = (nnz_lno_t *) (tmp);
              tmp += pow2_hash_size ;

              //poins to the next elements
              hm2.hash_nexts = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;

              //holds the keys
              hm2.keys = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;
              hm2.values = (nnz_lno_t *) (tmp);

              hm2.hash_key_size = pow2_hash_size;
              hm2.max_value_size = MaxRoughNonZero;
            }

            nnz_lno_t hash = -1;
            if (num_unsuccess) hash = b_set_ind & pow2_hash_func;

            int insertion = hm2.vector_atomic_insert_into_hash_mergeOr_TrackHashes(
                teamMember, vector_size,
                hash,b_set_ind,b_set,
                used_hash_sizes + 1, hm2.max_value_size
                ,globally_used_hash_count, globally_used_hash_indices
                );

          }
          left_work -= work_to_handle;
          rowBegin += work_to_handle;
        }
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
      });

      /*
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[1] > hm2.max_value_size) used_hash_sizes[1] = hm2.max_value_size;
      });
      */

      nnz_lno_t num_elements = 0;

      nnz_lno_t num_compressed_elements = used_hash_sizes[0];

      Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
          [&] (const nnz_lno_t ii, nnz_lno_t &num_nnz_in_row) {
        nnz_lno_t c_rows = hm.values[ii];
        nnz_lno_t num_el = 0;
        for (; c_rows; num_el++) {
          c_rows &= c_rows - 1; // clear the least significant bit set
        }
        num_nnz_in_row += num_el;
      }, num_elements);


      if (is_global_alloced){
        nnz_lno_t num_global_elements = 0;
        nnz_lno_t num_compressed_elements = used_hash_sizes[1];
        Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
            [&] (const nnz_lno_t ii, nnz_lno_t &num_nnz_in_row) {
          nnz_lno_t c_rows = hm2.values[ii];
          nnz_lno_t num_el = 0;
          for (; c_rows; num_el++) {
            c_rows &= c_rows - 1; // clear the least significant bit set
          }
          num_nnz_in_row += num_el;
        }, num_global_elements);


        //now thread leaves the memory as it finds. so there is no need to initialize the hash begins
        nnz_lno_t dirty_hashes = globally_used_hash_count[0];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, dirty_hashes),
            [&] (nnz_lno_t i) {
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        });


        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          m_space.release_chunk(globally_used_hash_indices);
        });
        num_elements += num_global_elements;
      }

      rowmapC(row_index) = num_elements;
    }

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }

  };



  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t,
            typename pool_memory_space>
  struct NonzeroesC{
    nnz_lno_t numrows;

    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
    c_nnz_view_t entriesSetIndicesC;
    c_nnz_view_t entriesSetsC;


    const nnz_lno_t pow2_hash_size;
    const nnz_lno_t pow2_hash_func;
    const nnz_lno_t MaxRoughNonZero;

    const size_t shared_memory_size;
    int vector_size;
    pool_memory_space m_space;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
    NonzeroesC(

        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,

        c_row_view_t rowmapC_,
        c_nnz_view_t entriesSetIndicesC_,

        const nnz_lno_t hash_size_,
        const nnz_lno_t MaxRoughNonZero_,
        const size_t sharedMemorySize_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):
          numrows(m_),

          row_mapA (row_mapA_),
          entriesA(entriesA_),

          row_mapB(row_mapB_),
          entriesSetIndicesB(entriesSetIndicesB_),
          entriesSetsB(entriesSetsB_),

          rowmapC(rowmapC_),
          entriesSetIndicesC(entriesSetIndicesC_),
          entriesSetsC(),


          pow2_hash_size(hash_size_),
          pow2_hash_func(hash_size_ - 1),
          MaxRoughNonZero(MaxRoughNonZero_),

          shared_memory_size(sharedMemorySize_),
          vector_size (), m_space(), my_exec_space(my_exec_space_)
          {}

    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;


      //printf("row:%d\n", row_index);

      int thread_memory = ((shared_memory_size/ 4 / teamMember.team_size())) * 4;
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

      //nnz_lno_t *alloc_global_memory = NULL;
      nnz_lno_t *globally_used_hash_indices = NULL;

      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      nnz_lno_t *globally_used_hash_count = (nnz_lno_t *) (all_shared_memory);

      all_shared_memory += sizeof(nnz_lno_t) ;
      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
      nnz_lno_t shared_memory_hash_size = (thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory;

      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      nnz_lno_t * vals = (nnz_lno_t *) (all_shared_memory);

      //printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts, keys, vals);
      //return;
      //first level hashmap
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;

      //initialize begins.
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
          [&] (int i) {
        begins[i] = -1;
      });

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
        globally_used_hash_count[0] = 0;
      });

      bool is_global_alloced = false;

      const size_type col_end = row_mapA[row_index + 1];
      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t col_size = col_end - col_begin;

      for (nnz_lno_t colind = 0; colind < col_size; ++colind){
        size_type a_col = colind + col_begin;

        nnz_lno_t rowB = entriesA[a_col];
        size_type rowBegin = row_mapB(rowB);

        nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

        while (left_work){
          nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

          nnz_lno_t b_set_ind = -1, b_set = -1;
          nnz_lno_t hash = -1;
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, work_to_handle),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            b_set_ind = entriesSetIndicesB[adjind];
            b_set = entriesSetsB[adjind];
            hash = b_set_ind % shared_memory_hash_size;
          });


          int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(
              teamMember, vector_size,
              hash, b_set_ind, b_set,
              used_hash_sizes,
              shared_memory_hash_size);


          int overall_num_unsuccess = 0;

          Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
              [&] (const int threadid, int &overall_num_unsuccess) {
            overall_num_unsuccess += num_unsuccess;
          }, overall_num_unsuccess);


          if (overall_num_unsuccess){

            //printf("row:%d\n", row_index);
            if (!is_global_alloced){
              volatile nnz_lno_t * tmp = NULL;
              size_t tid = get_thread_id(row_index);
              while (tmp == NULL){
                Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile nnz_lno_t * &memptr) {
                  memptr = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
                }, tmp);
              }
              is_global_alloced = true;

              globally_used_hash_indices = (nnz_lno_t *) tmp;
              tmp += pow2_hash_size ;

              hm2.hash_begins = (nnz_lno_t *) (tmp);
              tmp += pow2_hash_size ;

              //poins to the next elements
              hm2.hash_nexts = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;

              //holds the keys
              hm2.keys = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;
              hm2.values = (nnz_lno_t *) (tmp);

              hm2.hash_key_size = pow2_hash_size;
              hm2.max_value_size = MaxRoughNonZero;
            }

            nnz_lno_t hash = -1;
            if (num_unsuccess) hash = b_set_ind & pow2_hash_func;

            int insertion = hm2.vector_atomic_insert_into_hash_mergeOr_TrackHashes(
                teamMember, vector_size,
                hash,b_set_ind,b_set,
                used_hash_sizes + 1, hm2.max_value_size
                ,globally_used_hash_count, globally_used_hash_indices
                );

          }
          left_work -= work_to_handle;
          rowBegin += work_to_handle;
        }
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
      });

      /*
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[1] > hm2.max_value_size) used_hash_sizes[1] = hm2.max_value_size;
      });
      */

      nnz_lno_t num_compressed_elements = used_hash_sizes[0];
      used_hash_sizes[0] = 0;
      size_type row_begin = rowmapC(row_index);
      int set_size = sizeof(nnz_lno_t) * 8;
      Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
          [&] (const nnz_lno_t ii) {
        nnz_lno_t c_rows_setind = hm.keys[ii];
        nnz_lno_t c_rows = hm.values[ii];

        int current_row = 0;
        nnz_lno_t unit = 1;

        while (c_rows){
          if (c_rows & unit){

            size_type wind = Kokkos::atomic_fetch_add(used_hash_sizes, 1);
            entriesSetIndicesC(wind + row_begin) = set_size * c_rows_setind + current_row;
          }
          current_row++;
          c_rows = c_rows & ~unit;
          unit = unit << 1;
        }

      });


      if (is_global_alloced){

        nnz_lno_t num_compressed_elements = used_hash_sizes[1];
        Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
            [&] (const nnz_lno_t ii) {
          nnz_lno_t c_rows_setind = hm2.keys[ii];
          nnz_lno_t c_rows = hm2.values[ii];

          int current_row = 0;
          nnz_lno_t unit = 1;

          while (c_rows){
            if (c_rows & unit){

              size_type wind = Kokkos::atomic_fetch_add(used_hash_sizes, 1);
              entriesSetIndicesC(wind + row_begin) = set_size * c_rows_setind + current_row;
            }
            current_row++;
            c_rows = c_rows & ~unit;
            unit = unit << 1;
          }
        });

        //now thread leaves the memory as it finds. so there is no need to initialize the hash begins
        nnz_lno_t dirty_hashes = globally_used_hash_count[0];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, dirty_hashes),
            [&] (nnz_lno_t i) {
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        });


        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          m_space.release_chunk(globally_used_hash_indices);
        });
      }

    }

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }

  };

  template <typename a_row_view_t, typename a_nnz_view_t,
              typename b_row_view_t, typename b_nnz_view_t,
              typename c_row_view_t, typename c_nnz_view_t>
  void symbolic_c(
      nnz_lno_t m,
      a_row_view_t row_mapA,
      a_nnz_view_t entriesA,

      b_row_view_t row_mapB,
      b_nnz_view_t entriesSetIndex,
      b_nnz_view_t entriesSets,

      c_row_view_t &rowmapC,
      c_nnz_view_t &entryIndicesC_,
      c_nnz_view_t &entriesSetsC_,
      nnz_lno_t maxNumRoughNonzeros
      ){

    size_t concurrency = MyExecSpace::concurrency();
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size();
    //maxNumRoughNonzeros = maxNumRoughNonzerohans ;
    if (verbose){
      std::cout << "\tSymbolic Concurrency:" << concurrency << std::endl;
    }

    //typedef Kokkos::Experimental::MemoryPool< MyExecSpace> pool_memory_space;
    typedef KokkosKernels::Experimental::Util::UniformMemoryPool< MyExecSpace, nnz_lno_t> pool_memory_space;

    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  entriesSetIndex.dimension_0();
    int teamSizeMax = 0;
    int vector_size = 0;
    size_t shared_memory_size = KOKKOSKERNELS_SPGEMM_SHMEMSIZE;

    nnz_lno_t min_hash_size = 1;
    while (maxNumRoughNonzeros > min_hash_size){
      min_hash_size *= 2;
    }
    //min_hash_size *= 4;


    SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();
    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = KokkosKernels::Experimental::Util::get_exec_space_type<MyExecSpace>();

    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA) {
      StructureC<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t, pool_memory_space>
      sc(
          m,
          row_mapA,
          entriesA,
          row_mapB,
          entriesSetIndex,
          entriesSets,
          rowmapC,
          min_hash_size,
          maxNumRoughNonzeros,
          shared_memory_size,
          my_exec_space);


      int max_allowed_team_size = team_policy_t::team_size_max(sc);
      KokkosKernels::Experimental::Util::get_suggested_vector_team_size<size_type, MyExecSpace>(
          max_allowed_team_size, vector_size, teamSizeMax, brows,bnnz);

      //teamSizeMax = 1;

      size_t chunksize = min_hash_size ; //this is for used hash indices
      chunksize += min_hash_size ; //this is for the hash begins
      chunksize += maxNumRoughNonzeros ; //this is for hash nexts
      chunksize += maxNumRoughNonzeros ; //this is for hash keys
      chunksize += maxNumRoughNonzeros ; //this is for hash values

      //TODO: HERE I would like to check the memory limits, and reduce the number of chunks based on
      //the available memory.
      size_t num_chunks = concurrency / vector_size;

      KokkosKernels::Experimental::Util::PoolType my_pool_type = KokkosKernels::Experimental::Util::ManyThread2OneChunk;


      Kokkos::Impl::Timer timer1;
      MyExecSpace::fence();
      pool_memory_space m_space(num_chunks, chunksize, -1,  my_pool_type);
      MyExecSpace::fence();

      if (verbose){
        std::cout << "\tPool Allocation Time:" << timer1.seconds() << std::endl;
        std::cout << "\tPool Alloc MB:" << (num_chunks * chunksize) / 1024. / 1024.  << std::endl;
      }

      sc.vector_size = vector_size;
      sc.m_space = m_space;

      if (verbose){
        std::cout << "Running symbolic-count with vs:" << vector_size
            << " teamSizeMax:" << teamSizeMax
            << " max_allowed_team_size:" << max_allowed_team_size << std::endl;
      }

      timer1.reset();
      Kokkos::parallel_for( team_count_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), sc);
      MyExecSpace::fence();

      if (verbose){
        std::cout << "\tNNZ Count TIME:" << timer1.seconds() << std::endl;
      }
      //if we need to find the max nnz in a row.
      {
        Kokkos::Impl::Timer timer1;
        size_type c_max_nnz = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<c_row_view_t, MyExecSpace>(m, rowmapC, c_max_nnz);
        MyExecSpace::fence();
        this->handle->get_spgemm_handle()->set_max_result_nnz(c_max_nnz);

        if (verbose){
          std::cout << "\tMax Reduce Count TIME:" << timer1.seconds() << std::endl;
        }
      }


      //perform parallel prefixsum for rowmapC


      KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<row_lno_temp_work_view_t, MyExecSpace>(m+1, rowmapC);
      MyExecSpace::fence();

      //KokkosKernels::Experimental::Util::print_1Dview(rowmapC);

      auto d_c_nnz_size = Kokkos::subview(rowmapC, m);
      auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
      Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
      typename c_row_view_t::non_const_value_type c_nnz_size = h_c_nnz_size();
      if (verbose){
        std::cout << "c_nnz_size:" << c_nnz_size << std::endl;
      }

      if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR ||
          spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR ||
          spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){

        entryIndicesC_ = c_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("entryIndicesC_"), c_nnz_size);

        NonzeroesC<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t, pool_memory_space>
        sc( m,
            row_mapA,
            entriesA,
            row_mapB,
            entriesSetIndex,
            entriesSets,
            rowmapC,
            entryIndicesC_,
            min_hash_size,
            maxNumRoughNonzeros,
            shared_memory_size,
            my_exec_space);

        sc.vector_size = vector_size;
        sc.m_space = m_space;


        timer1.reset();
        Kokkos::parallel_for( team_count_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), sc);
        MyExecSpace::fence();


        if (verbose){
          std::cout << "\tNNZ FILL TIME:" << timer1.seconds() <<  std::endl;
        }

        nnz_lno_t num_colors, num_multi_colors, num_used_colors;
        nnz_lno_persistent_work_host_view_t color_xadj;
        nnz_lno_persistent_work_view_t color_adj, vertex_colors;

        this->d2_color_c_matrix(rowmapC, entryIndicesC_, num_colors, color_xadj, color_adj , vertex_colors, num_multi_colors, num_used_colors, spgemm_algorithm);

        std::cout << "Done with coloring" << std::endl;


        timer1.reset();
        for (nnz_lno_t i = 0; i < num_used_colors; ++i){
          if (color_xadj(i+1) - color_xadj(i) <= 32) continue;
          auto sv = Kokkos::subview(color_adj,Kokkos::pair<nnz_lno_t, nnz_lno_t> (color_xadj(i), color_xadj(i+1)));
          Kokkos::sort(sv);
        }


        if (verbose){
          std::cout << "\tCOLOR SORT TIME:" << timer1.seconds() <<  std::endl;
        }

        //this->write_matrix_to_plot(num_colors, color_xadj, color_adj, rowmapC,entryIndicesC_);

        this->handle->get_spgemm_handle()->set_color_xadj(num_colors, color_xadj, color_adj, vertex_colors, num_multi_colors, num_used_colors);


      }
    }
    else {
      KokkosKernels::Experimental::Util::PoolType my_pool_type = KokkosKernels::Experimental::Util::OneThread2OneChunk;


      StructureC_CPU<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t, pool_memory_space>
      sc(
          m,
          row_mapA,
          entriesA,
          row_mapB,
          entriesSetIndex,
          entriesSets,
          rowmapC,
          min_hash_size,
          maxNumRoughNonzeros,
          shared_memory_size,
          my_exec_space);


      int max_allowed_team_size = team_policy_t::team_size_max(sc);
      KokkosKernels::Experimental::Util::get_suggested_vector_team_size<size_type, MyExecSpace>(
          max_allowed_team_size, vector_size, teamSizeMax, brows,bnnz);

      //teamSizeMax = 1;

      size_t chunksize = min_hash_size ; //this is for used hash indices
      chunksize += min_hash_size ; //this is for the hash begins
      chunksize += maxNumRoughNonzeros ; //this is for hash nexts
      chunksize += maxNumRoughNonzeros ; //this is for hash keys
      chunksize += maxNumRoughNonzeros ; //this is for hash values

      //TODO: HERE I would like to check the memory limits, and reduce the number of chunks based on
      //the available memory.
      size_t num_chunks = concurrency / vector_size;



      Kokkos::Impl::Timer timer1;
      MyExecSpace::fence();
      pool_memory_space m_space(num_chunks, chunksize, -1,  my_pool_type);
      MyExecSpace::fence();

      if (verbose){
        std::cout << "\tPool Allocation Time:" << timer1.seconds() << std::endl;
        std::cout << "\tPool Alloc MB:" << (num_chunks * chunksize) / 1024. / 1024.  << std::endl;
      }

      //sc.vector_size = vector_size;
      sc.m_space = m_space;

      if (verbose){
        std::cout << "Running symbolic-count with vs:" << vector_size
            << " teamSizeMax:" << teamSizeMax
            << " max_allowed_team_size:" << max_allowed_team_size << std::endl;
      }

      timer1.reset();
      Kokkos::parallel_for( team_count_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), sc);
      MyExecSpace::fence();

      if (verbose){
        std::cout << "\tNNZ Count TIME:" << timer1.seconds() << std::endl;
      }
      //if we need to find the max nnz in a row.
      {
        Kokkos::Impl::Timer timer1;
        size_type c_max_nnz = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<c_row_view_t, MyExecSpace>(m, rowmapC, c_max_nnz);
        MyExecSpace::fence();
        this->handle->get_spgemm_handle()->set_max_result_nnz(c_max_nnz);

        if (verbose){
          std::cout << "\tMax Reduce Count TIME:" << timer1.seconds() << std::endl;
        }
      }


      //perform parallel prefixsum for rowmapC


      KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<row_lno_temp_work_view_t, MyExecSpace>(m+1, rowmapC);
      MyExecSpace::fence();

      //KokkosKernels::Experimental::Util::print_1Dview(rowmapC);

      auto d_c_nnz_size = Kokkos::subview(rowmapC, m);
      auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
      Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
      typename c_row_view_t::non_const_value_type c_nnz_size = h_c_nnz_size();
      if (verbose){
        std::cout << "c_nnz_size:" << c_nnz_size << std::endl;
      }
      if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR ||
          spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR ||
          spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){


        KokkosKernels::Experimental::Util::print_1Dview(rowmapC);
        entryIndicesC_ = c_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("entryIndicesC_"), c_nnz_size);

        NonzeroesC_CPU<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t, pool_memory_space>
        sc( m,
            row_mapA,
            entriesA,
            row_mapB,
            entriesSetIndex,
            entriesSets,
            rowmapC,
            entryIndicesC_,
            min_hash_size,
            maxNumRoughNonzeros,
            shared_memory_size,
            my_exec_space);

        sc.m_space = m_space;


        timer1.reset();
        Kokkos::parallel_for( team_count_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), sc);
        MyExecSpace::fence();


        if (verbose){
          std::cout << "\tNNZ FILL TIME:" << timer1.seconds() << std::endl;
        }

        //KokkosKernels::Experimental::Util::print_1Dview(entryIndicesC_);

        nnz_lno_t num_colors, num_multi_colors, num_used_colors;
        nnz_lno_persistent_work_host_view_t color_xadj;
        nnz_lno_persistent_work_view_t color_adj, vertex_colors;
        this->d2_color_c_matrix(rowmapC, entryIndicesC_, num_colors, color_xadj, color_adj, vertex_colors, num_multi_colors, num_used_colors, spgemm_algorithm);

        timer1.reset();
        for (nnz_lno_t i = 0; i < num_used_colors; ++i){
          if (color_xadj(i+1) - color_xadj(i) <= 32) continue;
          auto sv = Kokkos::subview(color_adj,Kokkos::pair<nnz_lno_t, nnz_lno_t> (color_xadj(i), color_xadj(i+1)));
          Kokkos::sort(sv);
        }

        if (verbose){
          std::cout << "\tCOLOR SORT TIME:" << timer1.seconds() <<  std::endl;
        }
        this->handle->get_spgemm_handle()->set_color_xadj(num_colors, color_xadj, color_adj, vertex_colors, num_multi_colors, num_used_colors);

        KokkosKernels::Experimental::Util::print_1Dview(color_xadj);
        KokkosKernels::Experimental::Util::print_1Dview(color_adj);
      }

    }
  }


  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t>
  void Fill_C_Graph(
      size_type m,
      a_row_view_t row_mapA,
      a_nnz_view_t entriesA,

      b_row_view_t row_mapB,
      b_nnz_view_t entriesSetIndex,
      b_nnz_view_t entriesSets,

      c_row_view_t &rowmapC,
      c_nnz_view_t &entriesC_,
      typename c_row_view_t::const_value_type maxNumRoughNonzeros
      ){

    Kokkos::Impl::Timer timer1;
    typedef typename c_row_view_t::non_const_type non_const_c_lno_row_view_t;
    typedef typename const_b_lno_nnz_view_t::non_const_value_type c_nnz_lno_t;

    /*
    int teamSizeMax = 0;
    int vector_size = 0;

    SymbolicC<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t> fnnz(
            m,
            row_mapA,
            entriesA,
            row_mapB,
            entriesSetIndex,
            entriesSets,
            rowmapC);



    int max_allowed_team_size = 1;//team_policy_t::team_size_max(cnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n, entriesB.dimension_0());
    Kokkos::parallel_for( team_count_policy_t((m / (teamSizeMax * KOKKOSKERNELS_SPGEMM_BLOCKSIZE)) + 1 , teamSizeMax, vector_size), fnnz);
    MyExecSpace::fence();
    */
    int num_teams = 1;

#if defined( KOKKOS_HAVE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, MyExecSpace >::value){
#pragma omp parallel
        #pragma omp master
        {
          num_teams = omp_get_num_threads();
        }
    }
#endif
    size_t min_hash_size = 1;
    while (maxNumRoughNonzeros > min_hash_size){
      min_hash_size *= 2;
    }

    if (min_hash_size < 16384){
      min_hash_size = 16384;
    }


    if (verbose){
      std::cout << "\tSYMBOLIC HASH SIZE:" << min_hash_size << std::endl;
    }

    if (this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK4){
      symbolic_sorted<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t> fnnz(
          m,
          row_mapA,
          entriesA,
          row_mapB,
          entriesSetIndex,
          entriesSets,
          rowmapC,
          maxNumRoughNonzeros,
          (m / num_teams + 1));
      int teamSizeMax = 1;
      int vector_size = 1;

      row_lno_temp_work_view_t compressed_rowmap_c (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
      fnnz.compressedrowmapC = compressed_rowmap_c;
      Kokkos::parallel_for(
          team_fillextra_policy_t(num_teams, teamSizeMax, vector_size), fnnz);
      MyExecSpace::fence();

      this->handle->get_spgemm_handle()->set_compressed_c(compressed_rowmap_c);


      size_type c_nnz = 0;
      KokkosKernels::Experimental::Util::view_reduce_max<non_const_c_lno_row_view_t, MyExecSpace>(m, compressed_rowmap_c, c_nnz);

      if (verbose){
        std::cout << "\tMAX NNZ for compressed C:" << c_nnz << std::endl;
      }

      this->handle->get_spgemm_handle()->set_max_compresed_result_nnz(c_nnz);

      KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<row_lno_temp_work_view_t, MyExecSpace>(m+1, compressed_rowmap_c);
      MyExecSpace::fence();
    }
    else {
      SymbolicC<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t> fnnz(
          m,
          row_mapA,
          entriesA,
          row_mapB,
          entriesSetIndex,
          entriesSets,
          rowmapC,
          min_hash_size,
          maxNumRoughNonzeros,
          (m / num_teams + 1));
      int teamSizeMax = 1;
      int vector_size = 1;


      if (this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK3){
        row_lno_temp_work_view_t compressed_rowmap_c (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
        fnnz.compressedrowmapC = compressed_rowmap_c;
        Kokkos::parallel_for(
            team_fillextra_policy_t(num_teams, teamSizeMax, vector_size), fnnz);
        MyExecSpace::fence();

        this->handle->get_spgemm_handle()->set_compressed_c(compressed_rowmap_c);


        size_type c_nnz = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<non_const_c_lno_row_view_t, MyExecSpace>(m, compressed_rowmap_c, c_nnz);

        if (verbose){
          std::cout << "\tMAX NNZ for compressed C:" << c_nnz << std::endl;
        }

        this->handle->get_spgemm_handle()->set_max_compresed_result_nnz(c_nnz);

        KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<row_lno_temp_work_view_t, MyExecSpace>(m+1, compressed_rowmap_c);
        MyExecSpace::fence();


      }
      else {
        Kokkos::parallel_for(
            multicore_team_count_policy_t(num_teams, teamSizeMax, vector_size), fnnz);
        MyExecSpace::fence();
      }
    }


    if (verbose){
      std::cout << "\tNNZ Count TIME:" << timer1.seconds() << std::endl;
    }


    timer1.reset();
    size_type c_nnz = 0;
    KokkosKernels::Experimental::Util::view_reduce_max<non_const_c_lno_row_view_t, MyExecSpace>(m, rowmapC, c_nnz);

    if (verbose){
      std::cout << "\tmax uncomppressed is:" << c_nnz << std::endl;
    }


    this->handle->get_spgemm_handle()->set_max_result_nnz(c_nnz);

    timer1.reset();
    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_row_view_t, MyExecSpace>(m+1, rowmapC);
    MyExecSpace::fence();

    //rowmapC = row_map_c_copy;
    /*
    auto d_c_nnz_size = Kokkos::subview(rowmapC, m);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    typename c_row_view_t::non_const_value_type c_nnz_size = h_c_nnz_size();

    std::cout << "\tCOMPRESSED C NNZ SIZE:" << c_nnz_size << std::endl;



    row_lno_temp_work_view_t c_set_index_entries(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"),c_nnz_size);
    row_lno_temp_work_view_t c_set_entries(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), c_nnz_size);

    fnnz.entriesSetIndicesC = c_set_index_entries;
    fnnz.entriesSetsC = c_set_entries;

    timer1.reset();
    Kokkos::parallel_for( team_fill_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), fnnz);

    MyExecSpace::fence();
    std::cout << "\tActual NNZ FILL TIME:" << timer1.seconds() << std::endl;

    timer1.reset();
    this->uncompressMatrix(rowmapC, entriesC_, row_map_c_copy, c_set_index_entries, c_set_entries);
    std::cout << "\tUncompress TIME:" << timer1.seconds() << std::endl << std::endl << std::endl;
    */
  }


  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_apply_color(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_, SPGEMMAlgorithm spgemm_algorithm){

    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space =
        KokkosKernels::Experimental::Util::get_exec_space_type<MyExecSpace>();
    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  valsB.dimension_0();
    int suggested_vector_size = get_suggested_vector__size(brows, bnnz, my_exec_space);

    nnz_lno_t num_colors, num_multi_colors, num_used_colors;
    nnz_lno_persistent_work_host_view_t color_xadj;
    nnz_lno_persistent_work_view_t color_adj, vertex_colors;
    this->handle->get_spgemm_handle()->get_color_xadj(num_colors, color_xadj, color_adj, vertex_colors, num_multi_colors, num_used_colors);

    const nnz_lno_t block_size = 64;
    const nnz_lno_t shift_divisor = 6;
    scalar_temp_work_view_t denseAccumulator ("Scalar Accumulator", ( (this->k + block_size) * num_multi_colors));


    //bool_temp_view_t denseAccumulatorFlags ("Accumulator flags", ((this->k* 1.5)  * num_multi_colors));
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size();

    NumericCCOLOR<
      const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
      const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
      c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t>
      sc(
        m, k,
        row_mapA,
        entriesA,
        valsA,

        row_mapB,
        entriesB,
        valsB,

        rowmapC_,
        entriesC_,
        valuesC_,

        denseAccumulator,
        //denseAccumulatorFlags,
        team_row_chunk_size);




    if (verbose){
      std::cout << "num_multi_colors:" << num_multi_colors << " num_used_colors:" << num_used_colors << std::endl;
    }
    sc.color_adj = color_adj;
    sc.vertex_colors = vertex_colors;
    sc.chunk_divison = shift_divisor;
    sc.chunk_and = block_size - 1;
    sc.consecutive_chunk_size = block_size;
    sc.consecutive_all_color_chunk_size = sc.consecutive_chunk_size * num_multi_colors;

    Kokkos::Impl::Timer timer1;
    for (nnz_lno_t i = 0; i < num_used_colors; ){
      nnz_lno_t color_begin = color_xadj(i);
      nnz_lno_t lastcolor = i + 1;
      if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){
        lastcolor = KOKKOSKERNELS_MACRO_MIN(i + num_multi_colors, num_used_colors );
        i += num_multi_colors;
      }
      else {
        ++i;
      }

      nnz_lno_t color_end = color_xadj(lastcolor);
      sc.color_begin = color_begin;
      sc.color_end = color_end;

      switch (spgemm_algorithm){
      case KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR:
        Kokkos::parallel_for( team_numeric1_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
            Kokkos::AUTO_t(), suggested_vector_size), sc);
        break;
      case KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2:
        Kokkos::parallel_for( team_numeric2_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
            Kokkos::AUTO_t(), suggested_vector_size), sc);
        break;
      case KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR:
        Kokkos::parallel_for( team_numeric3_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
            Kokkos::AUTO_t(), suggested_vector_size), sc);
        break;
      }
      MyExecSpace::fence();

    }

    if (verbose){
      std::cout << "\tNumeric TIME:" << timer1.seconds() << std::endl;
    }
  }

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_apply_speed(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_,
      KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space){

    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  valsB.dimension_0();
    int suggested_vector_size = get_suggested_vector__size(brows, bnnz, my_exec_space);
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size();

    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      c_lno_nnz_view_t beginsC (Kokkos::ViewAllocateWithoutInitializing("C keys"), valuesC_.dimension_0());
      c_lno_nnz_view_t nextsC (Kokkos::ViewAllocateWithoutInitializing("C nexts"), valuesC_.dimension_0());
      Kokkos::deep_copy(beginsC, -1);
      int shared_memory_size = KOKKOSKERNELS_SPGEMM_SHMEMSIZE;
      NumericCMEM<
      const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
      const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
      c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t>
      sc(
          m,
          row_mapA,
          entriesA,
          valsA,

          row_mapB,
          entriesB,
          valsB,

          rowmapC_,
          entriesC_,
          valuesC_,

          beginsC, nextsC,
          shared_memory_size);

      Kokkos::Impl::Timer timer1;
      MyExecSpace::fence();
      sc.vector_size = suggested_vector_size;
      sc.team_work_size = team_row_chunk_size;
      if (verbose){
        std::cout << "Running GPU NUMERIC MEM numeric with vs:" <<
            suggested_vector_size << " team_row_chunk_size:" << team_row_chunk_size << std::endl;
      }

      timer1.reset();
      Kokkos::parallel_for( gpu_team_policy_t(m / team_row_chunk_size + 1 , 256 / suggested_vector_size, suggested_vector_size), sc);
      MyExecSpace::fence();

      if (verbose){
        std::cout << "\tNumeric TIME:" << timer1.seconds() << std::endl;
      }
    }
    else {
      typedef KokkosKernels::Experimental::Util::UniformMemoryPool< MyExecSpace, scalar_t> pool_memory_space;
      NumericCMEM_CPU<
        const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
        const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
        c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t,
        pool_memory_space>
        sc(
          m,
          k,
          row_mapA,
          entriesA,
          valsA,

          row_mapB,
          entriesB,
          valsB,

          rowmapC_,
          entriesC_,
          valuesC_, my_exec_space);

      Kokkos::Impl::Timer timer1;
      MyExecSpace::fence();

      KokkosKernels::Experimental::Util::PoolType my_pool_type = KokkosKernels::Experimental::Util::OneThread2OneChunk;
      int num_chunks = MyExecSpace::concurrency();
      pool_memory_space m_space(num_chunks, this->k + (this->k) / sizeof(scalar_t) + 1, 0,  my_pool_type);
      MyExecSpace::fence();

      if (verbose){
        std::cout << "\tPool Allocation Time:" << timer1.seconds() << std::endl;
        std::cout << "\tPool Alloc MB:" << (num_chunks * (this->k + (this->k) / sizeof(scalar_t) + 1)) / 1024. / 1024.  << std::endl;
      }
      sc.memory_space = m_space;
      sc.team_work_size = team_row_chunk_size;

      if (verbose){
        std::cout << "Running CPU NUMERIC MEM numeric with vs:" <<
            suggested_vector_size << " team_row_chunk_size:" << team_row_chunk_size << std::endl;
      }
      timer1.reset();
      Kokkos::parallel_for( multicore_team_policy_t(m / team_row_chunk_size + 1 , Kokkos::AUTO_t(), suggested_vector_size), sc);
      MyExecSpace::fence();

      if (verbose){
        std::cout << "\tNumeric TIME:" << timer1.seconds() << std::endl;
      }
    }
  }

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_apply_hash(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_,
      KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space){

    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  valsB.dimension_0();
    int suggested_vector_size = get_suggested_vector__size(brows, bnnz, my_exec_space);

    typedef KokkosKernels::Experimental::Util::UniformMemoryPool< MyExecSpace, nnz_lno_t> pool_memory_space;

    int shared_memory_size = KOKKOSKERNELS_SPGEMM_SHMEMSIZE;
    PortableNumericCHASH<
    const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
    const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
    c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t,
    pool_memory_space>
    sc(
        m,
        row_mapA,
        entriesA,
        valsA,

        row_mapB,
        entriesB,
        valsB,

        rowmapC_,
        entriesC_,
        valuesC_,
        my_exec_space);

    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size();
    Kokkos::Impl::Timer timer1;

    nnz_lno_t max_nnz = this->handle->get_spgemm_handle()->get_max_result_nnz();
    nnz_lno_t min_hash_size = 1;
    while (max_nnz > min_hash_size){
      min_hash_size *= 2;
    }

    size_t chunksize = min_hash_size; //this is for used hash indices
    chunksize += min_hash_size ; //this is for the hash begins
    chunksize += max_nnz; //this is for hash nexts
    int num_chunks = MyExecSpace::concurrency() / suggested_vector_size;

    if (verbose){
      std::cout << "\tNUMERIC:maxNonzeros: " << max_nnz << " chunksize:" << chunksize << " numchunks:" << num_chunks << std::endl;
    }

    KokkosKernels::Experimental::Util::PoolType my_pool_type = KokkosKernels::Experimental::Util::OneThread2OneChunk;
    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      my_pool_type = KokkosKernels::Experimental::Util::ManyThread2OneChunk;
    }

    pool_memory_space m_space(num_chunks, chunksize, -1,  my_pool_type);
    MyExecSpace::fence();

    if (verbose){
      std::cout << "\tPool Allocation Time:" << timer1.seconds() << std::endl;
      std::cout << "\tPool Alloc MB:" << (num_chunks * chunksize) / 1024. / 1024.  << std::endl;
    }


    sc.memory_space = m_space;
    sc.max_nnz = max_nnz;
    sc.pow2_hash_size = min_hash_size;
    sc.pow2_hash_func = min_hash_size - 1;
    sc.vector_size = suggested_vector_size;
    sc.shared_memory_size = shared_memory_size;
    sc.team_work_size = team_row_chunk_size;

    if (verbose){
      std::cout << "\t\tRunning GPU NUMERIC HASH with vs:" << suggested_vector_size  << " team_row_chunk_size:" << team_row_chunk_size << std::endl;
    }
    timer1.reset();

    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      Kokkos::parallel_for( gpu_team_policy_t(m / team_row_chunk_size + 1 , Kokkos::AUTO_t(), suggested_vector_size), sc);
      MyExecSpace::fence();
    }
    else {

      Kokkos::parallel_for( multicore_team_policy_t(m / team_row_chunk_size + 1 , Kokkos::AUTO_t(), suggested_vector_size), sc);
      MyExecSpace::fence();
    }

    if (verbose){
      std::cout << "\tNumeric TIME:" << timer1.seconds() << std::endl;
    }

  }

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_apply(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_){

    Kokkos::Impl::Timer timer1;
    auto d_c_nnz_size = Kokkos::subview(rowmapC_, rowmapC_.dimension_0() - 1);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    size_t c_nnz_size = h_c_nnz_size();

    if (verbose){
      std::cout << "\tC NNZ SIZE:" << c_nnz_size << std::endl;
    }


    valuesC_ = c_scalar_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("C values"), c_nnz_size);
    SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();
    if (!
        (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR ||
        spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR ||
        spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2)){
      entriesC_ = c_lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("C entries"), c_nnz_size);
    }
    else {
      KokkosKernels::Experimental::Util::print_1Dview(entriesC_);
    }
    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = KokkosKernels::Experimental::Util::get_exec_space_type<MyExecSpace>();

    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  valsB.dimension_0();
    int suggested_vector_size = get_suggested_vector__size(brows, bnnz, my_exec_space);

    //SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();
    if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED)
    {
      this->KokkosSPGEMM_apply_speed(rowmapC_, entriesC_, valuesC_, my_exec_space);

    }
    else if ( spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR ||
              spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR ||
              spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){
      this->KokkosSPGEMM_apply_color(rowmapC_, entriesC_, valuesC_, spgemm_algorithm);
    }
    else {
      this->KokkosSPGEMM_apply_hash(rowmapC_, entriesC_, valuesC_, my_exec_space);
    }

    return;
#if 0
    /*
    calculateC <const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
    const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
    c_lno_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t> cc(
        m,
        row_mapA,
        entriesA,
        valsA,
        row_mapB,
        entriesB,
        valsB,
        rowmapC_, entriesC_, valuesC
    );


    int vector_size = 1, teamSizeMax = 1;
    int max_allowed_team_size = 1;//team_policy_t::team_size_max(cnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n, entriesB.dimension_0());
    Kokkos::parallel_for( team_fill_policy_t (m / (teamSizeMax * block_size) + 1 , teamSizeMax, vector_size), cc);
    MyExecSpace::fence();
    */

    int num_teams = 1;
#if defined( KOKKOS_HAVE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, MyExecSpace >::value){
#pragma omp parallel
#pragma omp master
      {
        num_teams = omp_get_num_threads();
      }
    }
#endif
    //bool used_dense_accumulator = true;
    size_t min_pow2_hash_size = 1;

    if (this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK2){
      const size_type max_non_zero = k;
      calculateC <const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
      const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
      c_lno_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t> cc(
          m,
          row_mapA,
          entriesA,
          valsA,
          row_mapB,
          entriesB,
          valsB,
          rowmapC_, entriesC_, valuesC,
          min_pow2_hash_size,
          max_non_zero,
          (m / num_teams + 1)
          //,max_flops
      );


      int vector_size = 1, teamSizeMax = 1;
      Kokkos::parallel_for( multicore_dense_team_count_policy_t(num_teams, teamSizeMax, vector_size), cc);
      MyExecSpace::fence();
    }
    else if (this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK3
        || this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK4){

      row_lno_temp_work_view_t compressed_b_rowmap_, compressed_b_set_indices_, compressed_b_sets_, compressed_b_set_begins_, compressed_b_set_nexts_;

      this->handle->get_spgemm_handle()->get_compressed_b(
          compressed_b_rowmap_, compressed_b_set_indices_,
          compressed_b_sets_, compressed_b_set_begins_,
          compressed_b_set_nexts_);


      const size_type max_non_zero = this->handle->get_spgemm_handle()->get_max_result_nnz();
      const size_type max_compressed_non_zero = this->handle->get_spgemm_handle()->get_max_compresed_result_nnz();

      while (max_non_zero > min_pow2_hash_size){
        min_pow2_hash_size *= 2;
      }

      //min_pow2_hash_size *= 2;
      if (min_pow2_hash_size < 16384){
        min_pow2_hash_size = 16384;
      }
      //min_pow2_hash_size *= 16;
      if (verbose){
        std::cout << "\tNUMERIC HASH SIZE:" << min_pow2_hash_size << std::endl;
      }
      calculateC <const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
      const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
      c_lno_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t> cc(
          m,
          row_mapA,
          entriesA,
          valsA,
          row_mapB,
          entriesB,
          valsB,

          compressed_b_rowmap_,
          compressed_b_set_indices_,
          compressed_b_sets_,
          compressed_b_set_begins_,
          compressed_b_set_nexts_,

          rowmapC_, entriesC_, valuesC,
          min_pow2_hash_size,
          max_non_zero,
          max_compressed_non_zero,
          (m / num_teams + 1)
          //,max_flops
      );


      int vector_size = 1, teamSizeMax = 1;
      Kokkos::parallel_for( team_fillextra_policy_t(num_teams, teamSizeMax, vector_size), cc);
      MyExecSpace::fence();


    }
    else {

      const size_type max_non_zero = this->handle->get_spgemm_handle()->get_max_result_nnz();

      //const row_lno_t max_non_zero == k;
      while (max_non_zero > min_pow2_hash_size){
        min_pow2_hash_size *= 2;
      }

      if (min_pow2_hash_size < 16384){
        min_pow2_hash_size = 16384;
      }
      //min_pow2_hash_size = 2 *k;
      //min_pow2_hash_size *= 16;

      if (verbose){
        std::cout << "NUMERIC HASH SIZE:" << min_pow2_hash_size << std::endl;
      }
      calculateC <const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
      const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
      c_lno_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t> cc(
          m,
          row_mapA,
          entriesA,
          valsA,
          row_mapB,
          entriesB,
          valsB,
          rowmapC_, entriesC_, valuesC,
          min_pow2_hash_size,
          max_non_zero,
          (m / num_teams + 1)
          //,max_flops
      );


      int vector_size = 1, teamSizeMax = 1;
      Kokkos::parallel_for( multicore_team_count_policy_t(num_teams, teamSizeMax, vector_size), cc);
      MyExecSpace::fence();

    }


    if (verbose){
      std::cout << "\tApply TIME:" << timer1.seconds() << std::endl;
    }


#ifdef TRACK_INSERTS
        std::cout << "cc.inserts:" << cc.overal_inserts(0)
            << " cc.hashop:" << cc.overall_hash_op(0)
            << " cc.wasted:" << cc.wasted_ops(0)
            << " hashop/inserts:" << cc.overall_hash_op(0) / (double) cc.overal_inserts(0)
            << " wasted/hashop:" << cc.wasted_ops(0) / (double) cc.overall_hash_op(0)
            << " (overall_hash_op - wasted_ops) / overal_inserts:"

            << (cc.overall_hash_op(0) - cc.wasted_ops(0)) / (double) cc.overal_inserts(0)
            << " used_hashes:" << cc.used_hashes(0) << " TOTALhash:" << 16 * m
            << " hash_use_ratio:" << cc.used_hashes(0) / (16.0 * m)
            << std::endl;

        /*
        typename c_lno_row_view_t::non_const_value_type rough_size = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<c_lno_row_view_t, MyExecSpace>(m, row_mapC, rough_size);

        MyExecSpace::fence();
        */

      //std::cout << "row:" << row_index << " hash:" << hash_ops << " inserts:" << inserts << std::endl;
#endif
#endif
  }

  template <typename c_row_view_t, typename c_nnz_view_t>
  void KokkosSPGEMM_symbolic(c_row_view_t &rowmapC_, c_nnz_view_t &entriesC_){

    row_lno_temp_work_view_t new_row_mapB;
    nnz_lno_temp_work_view_t set_index_entries;
    nnz_lno_temp_work_view_t set_entries;

    //incase the mapping from compressed to uncompressed is needed, we the sets below.
    row_lno_temp_work_view_t set_entries_begins, set_entries_next;

    //First Compress B.
    Kokkos::Impl::Timer timer1;

    this->compressMatrix(this->row_mapB, this->entriesB, new_row_mapB, set_index_entries, set_entries,set_entries_begins,set_entries_next);

    if (verbose){
      std::cout << "\n\n\tCOMPRESSION TIME:" << timer1.seconds() << std::endl;
      std::cout <<  "\tOld NNZ:" << this->entriesB.dimension_0() <<
                    " NEW NNZ:" << set_entries.dimension_0() <<
                    " compression ratio:" << set_entries.dimension_0() / (double) this->entriesB.dimension_0() << std::endl;
    }


    //if KK3 or KK4 save compression info.
    if (this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK3 ||
        this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK4){
      this->handle->get_spgemm_handle()->set_compressed_b(
          new_row_mapB,
          set_index_entries,
          set_entries,
          set_entries_begins,
          set_entries_next);
    }

    //CHECK if the compression is correct.
    /*
    {
      non_const_b_lno_row_view_t row_mapB2;
      non_const_b_lno_nnz_view_t entriesB2;
      this->uncompressMatrix(row_mapB2, entriesB2, new_row_mapB, set_index_entries, set_entries);
      if (!KokkosKernels::Experimental::Util::isSame<const_b_lno_row_view_t, non_const_b_lno_row_view_t, MyExecSpace>(this->n, this->row_mapB, row_mapB2)){
        std::cout << "rowmapB and rowmapB2 differ!!!!" << std::endl;
      }

      if (!KokkosKernels::Experimental::Util::isSame<const_b_lno_nnz_view_t, non_const_b_lno_nnz_view_t, MyExecSpace>
      (entriesB2.dimension_0(), this->entriesB, entriesB2)){
        std::cout << "entriesB and entriesB2 differ!!!!" << std::endl;
      }

      KokkosKernels::Experimental::Util::print_1Dview(row_mapB);
      KokkosKernels::Experimental::Util::print_1Dview(row_mapB2);

      KokkosKernels::Experimental::Util::print_1Dview(entriesB);
      KokkosKernels::Experimental::Util::print_1Dview(entriesB2);
    }
    */

    timer1.reset();
    typedef typename c_row_view_t::non_const_type nonconst_c_row_view_t;
    typedef typename c_nnz_view_t::non_const_type nonconst_c_nnz_view_t;

    nonconst_c_row_view_t rowmapC (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
    //incase structure to be captured as well.
    nonconst_c_nnz_view_t entryIndicesC;
    nonconst_c_nnz_view_t entrySetsC;

    //First get a rough count.
    nnz_lno_t maxNumRoughZeros = this->getMaxRoughRowNNZ(m, row_mapA, entriesA, new_row_mapB, rowmapC);
    if (verbose){
      std::cout << "\tUpper Bound Max NNZ in a Row:" << maxNumRoughZeros  << std::endl;
      std::cout << "\tCalculate Upper Bound Time:" << timer1.seconds()  << std::endl;
    }

    this->symbolic_c(m, row_mapA, entriesA,
                        new_row_mapB, set_index_entries, set_entries,
                        rowmapC, entryIndicesC, entrySetsC,
                        maxNumRoughZeros);

    rowmapC_ = rowmapC;
    entriesC_ = entryIndicesC;

  }


  template <typename in_row_view_t, typename in_nnz_view_t, typename out_rowmap_view_t, typename out_nnz_view_t>
  void compressMatrix(
      in_row_view_t in_row_map,
      in_nnz_view_t in_entries,

      out_rowmap_view_t &out_row_map,
      out_nnz_view_t &out_nnz_indices,
      out_nnz_view_t &out_nnz_sets,

      out_rowmap_view_t &out_column_set_begins,
      out_rowmap_view_t &out_column_set_nexts){

    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = KokkosKernels::Experimental::Util::get_exec_space_type<MyExecSpace>();

    //number of rows
    nnz_lno_t n = in_row_map.dimension_0() - 1;
    size_type nnz = in_entries.dimension_0();

    //size of the lno_t, how many bits it can hold.
    int lnot_size = sizeof(nnz_lno_t) * 8;
    int compression_bit_divide_shift_ = 0;
    int val = lnot_size;
    while (val > 1) {
      ++compression_bit_divide_shift_;
      val = val >> 1;
    }

    nnz_lno_t compression_bit_mask_ = 0;
    for (int i = 0; i < compression_bit_divide_shift_; ++i){
      compression_bit_mask_ = (compression_bit_mask_ << 1) | 1;
    }

    Kokkos::Impl::Timer timer1;

    size_type new_nnz_cnt = 0;
    out_row_map = out_rowmap_view_t (Kokkos::ViewAllocateWithoutInitializing("new row map"), n+1);

    int teamSizeMax = 0;
    int vector_size = 0;


    out_nnz_view_t set_entries_ (Kokkos::ViewAllocateWithoutInitializing("set_entries_"), nnz);
    out_nnz_view_t set_indices_ (Kokkos::ViewAllocateWithoutInitializing("set_indices_"), nnz);

    out_nnz_view_t set_begins_ (Kokkos::ViewAllocateWithoutInitializing("set_begins_"), nnz);
    out_nnz_view_t set_nexts_ (Kokkos::ViewAllocateWithoutInitializing("set_nexts_"), nnz);

    Kokkos::deep_copy (set_begins_, -1);

    MyExecSpace::fence();
    if (this->verbose){
      std::cout << "\tCompression Allocations :" <<  timer1.seconds() << std::endl;
    }

    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      SingleStepZipMatrix <in_row_view_t, in_nnz_view_t, out_rowmap_view_t, out_nnz_view_t>
      sszm_compressMatrix(
          in_row_map, in_entries, //input symbolic matrix
          compression_bit_mask_, compression_bit_divide_shift_, vector_size, //divisor, shifter, vector size.
          out_row_map, //output row map
          set_begins_, set_nexts_, set_indices_, set_entries_ //global hash memory space.
      );

      int max_allowed_team_size = team_policy_t::team_size_max(sszm_compressMatrix);

      KokkosKernels::Experimental::Util::get_suggested_vector_team_size<size_type, MyExecSpace>(
          max_allowed_team_size, vector_size, teamSizeMax, n,nnz);
      sszm_compressMatrix.vector_size = vector_size;


      if (this->verbose){
        std::cout << "\tRunning compression with VECTORSIZE:" << vector_size
            << " teamSizeMax:" << teamSizeMax
            << " max_allowed_team_size:" << max_allowed_team_size << std::endl;
      }

      timer1.reset();
      Kokkos::parallel_for( team_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), sszm_compressMatrix);
      MyExecSpace::fence();
    }
    else {

      SingleStepZipMatrixSingleVector <in_row_view_t, in_nnz_view_t, out_rowmap_view_t, out_nnz_view_t>
      sszm_compressMatrix(
          in_row_map, in_entries, //input symbolic matrix
          compression_bit_mask_, compression_bit_divide_shift_, vector_size, //divisor, shifter, vector size.
          out_row_map, //output row map
          set_begins_, set_nexts_, set_indices_, set_entries_ //global hash memory space.
      );

      int max_allowed_team_size = team_policy_t::team_size_max(sszm_compressMatrix);

      KokkosKernels::Experimental::Util::get_suggested_vector_team_size<size_type, MyExecSpace>(
          max_allowed_team_size, vector_size, teamSizeMax, n,nnz);
      sszm_compressMatrix.vector_size = vector_size;


      if (this->verbose){
        std::cout << "\tRunning compression with VECTORSIZE:" << vector_size
            << " teamSizeMax:" << teamSizeMax
            << " max_allowed_team_size:" << max_allowed_team_size << std::endl;
      }

      timer1.reset();
      Kokkos::parallel_for( team_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), sszm_compressMatrix);
      MyExecSpace::fence();
    }


    if (this->verbose){
      std::cout << "\tFirst Step of Zip Matrix:" <<  timer1.seconds() << std::endl;
    }


    //KokkosKernels::Experimental::Util::print_1Dview(out_row_map, true);
    timer1.reset();
    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum <out_rowmap_view_t, MyExecSpace> (n+1, out_row_map);
    MyExecSpace::fence();

    if (this->verbose){
      std::cout << "\tPPS Zip Matrix:" <<  timer1.seconds() << std::endl;
    }

    auto d_c_nnz_size = Kokkos::subview(out_row_map, n);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    new_nnz_cnt = h_c_nnz_size();

    timer1.reset();
    set_nexts_ = out_nnz_view_t ();
    set_begins_ = out_nnz_view_t ();
    out_nnz_indices = out_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("set index entries"), new_nnz_cnt);
    out_nnz_sets = out_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("set entries"), new_nnz_cnt);
    MyExecSpace::fence();

    if (this->verbose){
      std::cout << "\tCompression Second Allocations:" <<  timer1.seconds() << std::endl;
    }


    copyMatrix <in_row_view_t, out_rowmap_view_t, out_nnz_view_t>
              matrixCopy(
                  in_row_map, set_indices_, set_entries_,
                  out_row_map, out_nnz_indices, out_nnz_sets);

    timer1.reset();
    Kokkos::parallel_for( team_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), matrixCopy/*, new_nnz_cnt*/);
    MyExecSpace::fence();

    if (this->verbose){
      std::cout << "\tCompression Copy Compressed Matrix:" <<  timer1.seconds() << std::endl;
    }

    /*
    cmc.set_index_entries = out_nnz_indices;
    cmc.set_entries = out_nnz_sets;
    if (!(this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK3 ||
        this->handle->get_spgemm_handle()->get_algorithm_type() == KokkosKernels::Experimental::Graph::SPGEMM_KK4)){
      MyExecSpace::fence();
      Kokkos::parallel_for( my_fill_exec_space(0, n), cmc);
      MyExecSpace::fence();
    }
    else {
      out_column_set_begins = out_view_t (Kokkos::ViewAllocateWithoutInitializing("set index entries"), new_nnz_cnt + 1);
      out_column_set_nexts = out_view_t (Kokkos::ViewAllocateWithoutInitializing("set entries"), in_entries.dimension_0());
      cmc.set_entries_begin = out_column_set_begins;
      cmc.set_entries_next = out_column_set_nexts;

      MyExecSpace::fence();
      Kokkos::parallel_for( my_fillextra_exec_space(0, n), cmc);
      MyExecSpace::fence();

    }
    */

  }

  template <typename out_row_view_t, typename out_nnz_view_t, typename in_row_view_t, typename in_nnz_view_t>
  void uncompressMatrix(
      out_row_view_t &out_row_map, out_nnz_view_t &out_entries,
      in_row_view_t in_row_map,
      in_nnz_view_t in_nnz_indices,
      in_nnz_view_t in_nnz_sets){

    size_type n = in_row_map.dimension_0() - 1;
    size_type nnz = in_nnz_indices.dimension_0();
    out_row_map = out_row_view_t(Kokkos::ViewAllocateWithoutInitializing("out_row_view_t out_row_map") , n + 1);
    MyExecSpace::fence();
    int teamSizeMax = 0;
    int vector_size = 0;
    size_t shared_memory_size = KOKKOSKERNELS_SPGEMM_SHMEMSIZE;


    unzipMatrix<out_row_view_t, out_nnz_view_t, in_row_view_t, in_nnz_view_t>
        ucc(n, in_row_map, in_nnz_indices, in_nnz_sets, out_row_map, out_entries, shared_memory_size);

    KokkosKernels::Experimental::Util::print_1Dview(out_row_map);

    int max_allowed_team_size = 256;//team_policy_t::team_size_max(cnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size
        <size_type, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n,nnz);
    //Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), cnnnz);




    Kokkos::parallel_for( team_count_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), ucc);
    MyExecSpace::fence();



    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<out_row_view_t, MyExecSpace>(n+1, out_row_map);
    MyExecSpace::fence();
    auto d_nnz_size = Kokkos::subview(out_row_map, n);
    auto h_nnz_size = Kokkos::create_mirror_view (d_nnz_size);
    Kokkos::deep_copy (h_nnz_size, d_nnz_size);
    MyExecSpace::fence();


    size_type nnz_size = h_nnz_size();
    out_entries = out_row_view_t(Kokkos::ViewAllocateWithoutInitializing("out_entries ") , nnz_size);
    MyExecSpace::fence();
    ucc.out_entries = out_entries;

    Kokkos::parallel_for( team_fill_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), ucc);
    MyExecSpace::fence();

  }

  template <typename c_row_view_t, typename c_nnz_view_t>
  void d2_color_c_matrix(c_row_view_t &rowmapC, c_nnz_view_t &entryIndicesC_,
      nnz_lno_t &num_colors,
      nnz_lno_persistent_work_host_view_t &h_color_xadj,
      nnz_lno_persistent_work_view_t &color_adj,
      nnz_lno_persistent_work_view_t &vertex_colors,
      nnz_lno_t &num_multi_colors,
      nnz_lno_t &num_used_colors,
      SPGEMMAlgorithm spgemm_algorithm){

    nnz_lno_persistent_work_view_t color_xadj;

    size_type c_nnz_size = entryIndicesC_.dimension_0();
    row_lno_temp_work_view_t transpose_col_xadj ("transpose_col_xadj", k + 1);
    nnz_lno_temp_work_view_t transpose_col_adj (Kokkos::ViewAllocateWithoutInitializing("tmp_row_view"), c_nnz_size);
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size();


    Kokkos::Impl::Timer timer1;
    KokkosKernels::Experimental::Util::transpose_graph<
        c_row_view_t, c_nnz_view_t,
        row_lno_temp_work_view_t, nnz_lno_temp_work_view_t, row_lno_temp_work_view_t, MyExecSpace>
            (m, k, rowmapC, entryIndicesC_, transpose_col_xadj,transpose_col_adj, team_row_chunk_size);
    MyExecSpace::fence();
    if (verbose){
      std::cout << "\tTranspose Time:" << timer1.seconds() << std::endl;
    }

    /*
    {

      timer1.reset();
      this->handle->create_graph_coloring_handle();
      //handle->get_graph_coloring_handle()->set_coloring_type(KokkosKernels::Experimental::Graph::Distance2);
      KokkosKernels::Experimental::Graph::d2_graph_color
        <HandleType,c_row_view_t,c_nnz_view_t, row_lno_temp_work_view_t,nnz_lno_temp_work_view_t>
        (this->handle, m, k, rowmapC, entryIndicesC_, transpose_col_xadj, transpose_col_adj);

      num_colors = handle->get_graph_coloring_handle()->get_num_colors();
      std::cout << "Num colors:" << num_colors <<  " coloring time:" << timer1.seconds()  << std::endl;
      typename HandleType::GraphColoringHandleType::color_view_t color_view = handle->get_graph_coloring_handle()->get_vertex_colors();


      nnz_lno_temp_work_view_t histogram ("histogram", handle->get_graph_coloring_handle()->get_num_colors() + 1);
      MyExecSpace::fence();
      timer1.reset();
      KokkosKernels::Experimental::Util::get_histogram
        <typename HandleType::GraphColoringHandleType::color_view_t, nnz_lno_temp_work_view_t, MyExecSpace>(m, color_view, histogram);
      std::cout << "Histogram" << " time:" << timer1.seconds()  << std::endl;
      KokkosKernels::Experimental::Util::print_1Dview(histogram);


      this->handle->destroy_graph_coloring_handle();
    }
    */


    {

      timer1.reset();
      this->handle->create_graph_coloring_handle();
      //handle->get_graph_coloring_handle()->set_coloring_type(KokkosKernels::Experimental::Graph::Distance2);
      handle->get_graph_coloring_handle()->set_algorithm(KokkosKernels::Experimental::Graph::COLORING_SERIAL2);
      KokkosKernels::Experimental::Graph::d2_graph_color
        <HandleType,c_row_view_t,c_nnz_view_t, row_lno_temp_work_view_t,nnz_lno_temp_work_view_t>
        (this->handle, m, k, rowmapC, entryIndicesC_, transpose_col_xadj, transpose_col_adj);

      num_colors = handle->get_graph_coloring_handle()->get_num_colors();
      num_used_colors =  num_colors;
      num_multi_colors = 1;

      std::cout << "Num colors:" << handle->get_graph_coloring_handle()->get_num_colors() <<  " coloring time:" << timer1.seconds()  << std::endl;
      typename HandleType::GraphColoringHandleType::color_view_t color_view = handle->get_graph_coloring_handle()->get_vertex_colors();

      vertex_colors = nnz_lno_persistent_work_view_t(
          Kokkos::ViewAllocateWithoutInitializing("persistent_color_view"), m);



      nnz_lno_temp_work_view_t histogram ("histogram", handle->get_graph_coloring_handle()->get_num_colors() + 1);

      MyExecSpace::fence();
      timer1.reset();
      KokkosKernels::Experimental::Util::get_histogram
        <typename HandleType::GraphColoringHandleType::color_view_t, nnz_lno_temp_work_view_t, MyExecSpace>(m, color_view, histogram);



      std::cout << "Histogram" << " time:" << timer1.seconds()  << std::endl;
      KokkosKernels::Experimental::Util::print_1Dview(histogram);

      {
        nnz_lno_temp_work_view_t tmp_color_view = color_view;
        if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR){
          tmp_color_view = nnz_lno_temp_work_view_t( Kokkos::ViewAllocateWithoutInitializing("tmp_color_view"), m);
          num_multi_colors = c_nnz_size / this->k;

          double scale_ = this->handle->get_spgemm_handle()->get_multi_color_scale();
          num_multi_colors = scale_ * num_multi_colors;

          if (num_multi_colors > 1){
            float scale_factor = 1.0 / num_multi_colors;
            KokkosKernels::Experimental::Util::a_times_x_plus_b
                <nnz_lno_temp_work_view_t, typename HandleType::GraphColoringHandleType::color_view_t, float, float, MyExecSpace>(
                    m, tmp_color_view, color_view, scale_factor, 1);
            num_used_colors = num_colors / num_multi_colors;
            if (num_colors % num_multi_colors) ++num_used_colors;
          }
          else {
            num_multi_colors = 1;
          }
        }

        /*
        std::cout << "tmp_color_view" << std::endl;
        KokkosKernels::Experimental::Util::print_1Dview(tmp_color_view,true);
        std::cout << "tmp_color_view" << std::endl;
        */

        KokkosKernels::Experimental::Util::modular_view
          <nnz_lno_persistent_work_view_t, typename HandleType::GraphColoringHandleType::color_view_t, MyExecSpace>
          (m, vertex_colors, color_view, num_multi_colors);

        timer1.reset();
        KokkosKernels::Experimental::Util::create_reverse_map
          <typename HandleType::GraphColoringHandleType::color_view_t,
          nnz_lno_persistent_work_view_t, MyExecSpace>
              (m, num_used_colors, tmp_color_view, color_xadj, color_adj);
        MyExecSpace::fence();

        if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){
          num_used_colors = num_colors;
          num_multi_colors = c_nnz_size / this->k;
          double scale_ = this->handle->get_spgemm_handle()->get_multi_color_scale();
          num_multi_colors = scale_ * num_multi_colors;
        }

        if (verbose){
          std::cout << "\tColor Reverse Map Create Time:" << timer1.seconds() << std::endl;
        }

        h_color_xadj = Kokkos::create_mirror_view (color_xadj);
        Kokkos::deep_copy (h_color_xadj, color_xadj);
        MyExecSpace::fence();

      }

      this->handle->destroy_graph_coloring_handle();
    }

  }

  template <typename c_row_view_t, typename c_nnz_view_t>
  void write_matrix_to_plot(
      nnz_lno_t &num_colors,
      nnz_lno_persistent_work_host_view_t &h_color_xadj,
      nnz_lno_persistent_work_view_t &color_adj,
      c_row_view_t &rowmapC, c_nnz_view_t &entryIndicesC_){
    std::cout << "writing to plot" << std::endl;

    nnz_lno_persistent_work_host_view_t h_color_adj = Kokkos::create_mirror_view (color_adj);
    Kokkos::deep_copy (h_color_adj, color_adj);
    auto h_rowmapC = Kokkos::create_mirror_view (rowmapC);
    Kokkos::deep_copy (h_rowmapC, rowmapC);
    auto h_entryIndicesC = Kokkos::create_mirror_view (entryIndicesC_);
    Kokkos::deep_copy (h_entryIndicesC, entryIndicesC_);

    for (nnz_lno_t i = 0; i < num_colors; ++i){
      nnz_lno_t color_begin = h_color_xadj(i);
      nnz_lno_t color_end = h_color_xadj(i + 1);

      std::string colorind = "";
      std::stringstream ss;
      ss << i;


      ss >> colorind;
      colorind += ".coords";
      std::fstream fs;
      fs.open(colorind.c_str(), std::fstream::out);

      std::cout << "COLOR:" << i << " colorbegin:" << color_begin << " colorend:" << color_end << " size:" << color_end - color_begin << std::endl;
      for (nnz_lno_t j = color_begin; j < color_end; ++j){
        nnz_lno_t row = h_color_adj(j);
        for (size_type k = h_rowmapC(row); k < h_rowmapC(row + 1); ++k){
          nnz_lno_t column = h_entryIndicesC(k);
          //std::cout << row << " " << column << std::endl;
          fs << row << " " << column << std::endl;
        }
      }
      fs.close();
    }



    std::fstream fs;
    fs.open("plot1.gnuplot", std::fstream::out);
    for (nnz_lno_t i = 0; i < num_colors; ++i){
      std::string colorind = "\"";
      std::stringstream ss;
      ss << i;

      ss >> colorind;
      colorind += ".coords\"";
      if (i > 0) fs << "re";
      fs << "plot " << colorind << std::endl;

    }
    fs << "pause -1" << std::endl;
    fs.close();
  }

};
}
}
}
}

#endif
