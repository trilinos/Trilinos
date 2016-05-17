
#ifndef _KOKKOSSPGEMMIMPL_HPP
#define _KOKKOSSPGEMMIMPL_HPP
#define HASHTRACK
//#define TRACK_INSERTS

#include <KokkosKernels_Utils.hpp>

#define KOKKOSKERNELS_SPGEMM_BLOCKSIZE 256
#define KOKKOSKERNELS_SPGEMM_HASHSIZE 16
#define KOKKOSKERNELS_SPGEMM_ANDHASHFUNCTION 15

namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{

template <typename HandleType,
  typename a_lno_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
  typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_
  >
class KokkosSPGEMM{
public:

  typedef a_lno_row_view_t_ a_in_lno_row_view_t;
  typedef a_lno_nnz_view_t_ a_in_lno_nnz_view_t;
  typedef a_scalar_nnz_view_t_ a_in_scalar_nnz_view_t;

  typedef b_lno_row_view_t_ b_in_lno_row_view_t;
  typedef b_lno_nnz_view_t_ b_in_lno_nnz_view_t;
  typedef b_scalar_nnz_view_t_ b_in_scalar_nnz_view_t;



  typedef typename a_in_lno_row_view_t::non_const_value_type row_lno_t;
  typedef typename a_in_lno_row_view_t::const_value_type const_row_lno_t;


  typedef typename a_in_lno_nnz_view_t::non_const_value_type nnz_lno_t;
  typedef typename a_in_lno_nnz_view_t::const_value_type const_nnz_lno_t;

  typedef typename a_in_scalar_nnz_view_t::non_const_value_type nnz_scalar_t;
  typedef typename a_in_scalar_nnz_view_t::const_value_type const_nnz_scalar_t;


  typedef typename a_in_lno_row_view_t::const_type const_a_lno_row_view_t;
  typedef typename a_in_lno_row_view_t::non_const_type non_const_a_lno_row_view_t;

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


  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t ;
  typedef typename team_policy_t::member_type team_member_t ;

  struct CountTag{};
  struct FillTag{};
  struct MultiCoreTagCount{};


  //typedef Kokkos::RangePolicy<MultiCoreTag, MyExecSpace> my_count_multicore_space;
  typedef Kokkos::TeamPolicy<MultiCoreTagCount, MyExecSpace> multicore_team_count_policy_t ;

  typedef Kokkos::RangePolicy<CountTag, MyExecSpace> my_count_exec_space;
  typedef Kokkos::RangePolicy<FillTag, MyExecSpace> my_fill_exec_space;
  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace> team_count_policy_t ;
  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace> team_fill_policy_t ;

private:
  HandleType *handle;
  row_lno_t m;
  row_lno_t n;
  row_lno_t k;


  const_a_lno_row_view_t row_mapA;
  const_a_lno_nnz_view_t entriesA;
  const_a_scalar_nnz_view_t valsA;
  bool transposeA;

  const_b_lno_row_view_t row_mapB;
  const_b_lno_nnz_view_t entriesB;
  const_b_scalar_nnz_view_t valsB;
  bool transposeB;

public:
  KokkosSPGEMM(
      HandleType *handle_,
      row_lno_t m_,
      row_lno_t n_,
      row_lno_t k_,
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
        row_lno_t m_,
        row_lno_t n_,
        row_lno_t k_,
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
            typename b_row_view_t, typename c_lno_row_view_t>
  struct Predict_C_NUM_NNZ{
    row_lno_t m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;
    b_row_view_t row_mapB;
    c_lno_row_view_t rough_row_mapC;
    Predict_C_NUM_NNZ(
        row_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_row_view_t row_mapB_,
        c_lno_row_view_t rough_row_mapC_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_mapB(row_mapB_), rough_row_mapC(rough_row_mapC_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {
      row_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //check ii is out of range. if it is, just return.
      if (row_index >= m) return;
      const row_lno_t col_begin = row_mapA[row_index];
      const row_lno_t col_end = row_mapA[row_index + 1];

      row_lno_t max_num_results_in_row = 0;

      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
          [&] (row_lno_t i, typename c_lno_row_view_t::non_const_value_type & valueToUpdate) {
        const row_lno_t adjind = i + col_begin;
        const row_lno_t colIndex = entriesA[adjind];
        valueToUpdate += row_mapB [colIndex + 1] - row_mapB[colIndex];
      },
      max_num_results_in_row);
      Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
        rough_row_mapC[row_index] = max_num_results_in_row;
      });
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

    zipMatrix(
        row_view_t row_map_,nnz_view_t entries_,
        lno_t compression_bit_mask_, int compression_bit_divide_shift_,
        new_row_view_t new_row_map_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_ ):
      row_map(row_map_), entries(entries_),
      compression_bit_mask(compression_bit_mask_), compression_bit_divide_shift(compression_bit_divide_shift_),
      new_row_map(new_row_map_), set_index_entries(set_index_entries_), set_entries(set_entries_){}

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

  };

  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t,typename c_lno_row_view_t>
  size_t getMaxRoughRowLen(
      row_lno_t m,
      a_row_view_t row_mapA,
      a_nnz_view_t entriesA,
      b_row_view_t row_mapB,
      c_lno_row_view_t row_mapC){
    int teamSizeMax = 0;
    int vector_size = 0;

    Predict_C_NUM_NNZ<a_row_view_t, a_nnz_view_t, b_row_view_t, c_lno_row_view_t> pcnnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        row_mapC);

    int max_allowed_team_size = team_policy_t::team_size_max(pcnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(
        max_allowed_team_size, vector_size, teamSizeMax, m, entriesA.dimension_0());

    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), pcnnnz);
    MyExecSpace::fence();

    typename c_lno_row_view_t::non_const_value_type rough_size = 0;
    KokkosKernels::Experimental::Util::view_reduce_max<c_lno_row_view_t, MyExecSpace>(m, row_mapC, rough_size);

    MyExecSpace::fence();

    //KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_lno_row_view_t, MyExecSpace>(m+1, row_mapC);
    //MyExecSpace::fence();
    return rough_size;
  }

  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
            typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
            typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t>
  struct calculateC{
    row_lno_t m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;
    a_scalar_view_t valuesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesB;
    b_scalar_view_t valuesB;

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
    const c_nnz_lno_t MaxRoughNonZero;

    const size_t requiredMem;
    const size_t multicore_chunk_size;

#ifdef TRACK_INSERTS
    //c_row_view_t overal_inserts, overall_hash_op, wasted_ops;
    Kokkos::View<long, Kokkos::OpenMP> overal_inserts, overall_hash_op, wasted_ops, used_hashes;
    //c_scalar_view_t hash_op_per_insert;
#endif

    calculateC(
        row_lno_t m_,
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
          MaxRoughNonZero(16384), requiredMem(16384),
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
        row_lno_t m_,
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
        size_t multicore_chunk_size_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_), valuesA(valuesA_),
          row_mapB(row_mapB_), entriesB(entriesB_), valuesB(valuesB_),
          row_mapC(row_mapC_), entriesC(entriesC_), valuesC(valuesC_),
          hash_size(min_pow2_hash_size), pow2_hash_func(min_pow2_hash_size -1),
          MaxRoughNonZero(MaxRoughNonZero_),
          requiredMem(MaxRoughNonZero_ * 2 * sizeof(c_nnz_lno_t) +
              MaxRoughNonZero_ * sizeof(c_scalar_t) +
              40 + hash_size * 2 * sizeof(c_nnz_lno_t) ),
          multicore_chunk_size(multicore_chunk_size_)

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
    void operator()(const MultiCoreTagCount&, const team_member_t & teamMember) const {
      row_lno_t row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);


      const size_t hash_index_array_size = MaxRoughNonZero * sizeof(c_nnz_lno_t);
      const size_t hash_begin_array_size = hash_size * sizeof(c_nnz_lno_t);
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);
      c_nnz_lno_t *hash_nexts = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_scalar_t *hash_values = (c_scalar_t *) teamMember.team_shmem().get_shmem(MaxRoughNonZero * sizeof(c_scalar_t));
      c_nnz_lno_t *used_hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);

      for (size_t i = 0; i < hash_size; ++i){
        hash_nexts[i] = hash_begins[i] = -1;
      }
#ifdef HASHTRACK
      size_t numhashes = 0;
      size_t numpoints = 0;
#endif
      for (; row_index < row_max; ++row_index){
        c_nnz_lno_t used_hash_counts = 0;
        c_nnz_lno_t next_available_entry = 0;

        //check ii is out of range. if it is, just return.
        row_lno_t col_begin = row_mapA[row_index];
        const row_lno_t col_end = row_mapA[row_index + 1];

        /*
        if (col_begin < col_end){
          const nnz_lno_t rowb = entriesA[col_begin];
          const c_scalar_t aval = valuesA[col_begin];

          //std::cout << "  for column:" << rowb << std::endl;
          row_lno_t colb_ind = row_mapB[rowb];
          const row_lno_t col_b_end = row_mapB[rowb + 1];
          for (; colb_ind < col_b_end; ++colb_ind){
            //printf("colb_ind:%d\n", colb_ind);

            const c_nnz_lno_t b_col_ind = entriesB[colb_ind];
            //const c_scalar_t bmul = valuesB[colb_ind] * aval;
            {
              //printf("b_col_set_index:%d", b_col_set_index)
              const c_nnz_lno_t hash = b_col_ind & pow2_hash_func;
              //std::cout << "\tinserting:" << b_col_set_index << " hash:" << hash << std::endl;
              c_nnz_lno_t hash_ind = hash_begins[hash];
              //std::cout << "\t\thash begin:" << hash_ind << std::endl;
              hash_begins[hash] = next_available_entry;
              hash_indices[next_available_entry] = b_col_ind;
              hash_values[next_available_entry] = valuesB[colb_ind] * aval;
              if (hash_ind == -1){
                used_hash_indices[used_hash_counts++] = hash;
              }
              else {
                hash_nexts[next_available_entry] = hash_ind;
              }
              ++next_available_entry;
            }


          }

        }
        ++col_begin;
        */
        for (; col_begin < col_end; ++col_begin){
          //printf("col_begin:%d\n", col_begin);
          const nnz_lno_t rowb = entriesA[col_begin];
          const c_scalar_t aval = valuesA[col_begin];

          //std::cout << "  for column:" << rowb << std::endl;
          row_lno_t colb_ind = row_mapB[rowb];
          const row_lno_t col_b_end = row_mapB[rowb + 1];
          for (; colb_ind < col_b_end; ){
            //printf("colb_ind:%d\n", colb_ind);

            const c_nnz_lno_t b_col_ind = entriesB[colb_ind];
            //const c_scalar_t bmul = valuesB[colb_ind] * aval;
            {
              //printf("b_col_set_index:%d", b_col_set_index)
              const c_nnz_lno_t hash = b_col_ind & pow2_hash_func;
              //std::cout << "\tinserting:" << b_col_set_index << " hash:" << hash << std::endl;
              c_nnz_lno_t hash_ind = hash_begins[hash];
              //std::cout << "\t\thash begin:" << hash_ind << std::endl;
              if (hash_ind == -1){
                //std::cout << "\t\tinserting for the first time" << std::endl;
                hash_begins[hash] = next_available_entry;
                hash_indices[next_available_entry] = b_col_ind;
                hash_values[next_available_entry++] = valuesB[colb_ind] * aval;
                used_hash_indices[used_hash_counts++] = hash;
              }
              else {
                while (hash_nexts[hash_ind] ^ -1){
                  if(b_col_ind ^ hash_indices[hash_ind]){ //if they are not equal
                    hash_ind = hash_nexts[hash_ind];
                  }
                  else {
                    hash_values[hash_ind] = hash_values[hash_ind] +  valuesB[colb_ind] * aval;
                    goto endloop; //break;
                  }
                }
                if(b_col_ind ^ hash_indices[hash_ind]){ //if they are not equal
                  hash_nexts[hash_ind] = next_available_entry;
                  hash_indices[next_available_entry] = b_col_ind;
                  hash_values[next_available_entry++] = valuesB[colb_ind] * aval;
                }
                else {
                  hash_values[hash_ind] = hash_values[hash_ind] + valuesB[colb_ind] * aval;
                }
              }
            }
            endloop:
            ++colb_ind;
          }
        }
        c_nnz_lno_t  my_nonzeros = row_mapC(row_index);
        for (c_nnz_lno_t i = 0; i < used_hash_counts; ++i){
#ifdef HASHTRACK
          ++numhashes;
#endif
          const int hash = used_hash_indices[i];
          c_nnz_lno_t hash_ind = hash_begins[hash];

          hash_begins[hash] = -1;
          do {
#ifdef HASHTRACK
            ++numpoints;
#endif
            entriesC(my_nonzeros) = hash_indices[hash_ind];
            valuesC(my_nonzeros++) = hash_values[hash_ind];
            c_nnz_lno_t newhash_ind = hash_nexts[hash_ind];
            hash_nexts[hash_ind] = -1;
            hash_ind = newhash_ind;
          }
          while (hash_ind != -1);
        }

      }
#ifdef HASHTRACK
      printf("numhash:%ld np:%ld\n", numhashes, numpoints);
#endif
    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const FillTag&, const team_member_t & teamMember) const {

      const int team_size = teamMember.team_size();
      const int team_rank = teamMember.team_rank();

      row_lno_t row_index = (teamMember.league_rank()  * team_size + team_rank) * KOKKOSKERNELS_SPGEMM_BLOCKSIZE;
      //if (row_index >= m) return;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + KOKKOSKERNELS_SPGEMM_BLOCKSIZE, m);
      //int row_max = row_index + KOKKOSKERNELS_SPGEMM_BLOCKSIZE;
      //if(row_max > m) row_max = m;



      const int hash_array_bytes = KOKKOSKERNELS_SPGEMM_HASHSIZE * sizeof(c_nnz_lno_t) * team_size;
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_bytes);
      c_nnz_lno_t *hash_ends = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_bytes);

      const int overall_set_size = ((16384 - hash_array_bytes * 2) - 16) / (sizeof(c_nnz_lno_t) + sizeof(c_scalar_t));


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
        row_lno_t col_a_index = row_mapA[row_index];
        const row_lno_t col_end = row_mapA[row_index + 1];

#ifdef TRACK_INSERTS
        long hash_ops = 0, inserts = 0;
        long wasted_op_count = 0;
#endif
        //c_nnz_lno_t  current_heap_size = 0;

        if (col_a_index < col_end){
          const nnz_lno_t rowb = entriesA[col_a_index];
          const c_scalar_t aval = valuesA[col_a_index];
          row_lno_t col_b_begin = row_mapB[rowb];
          const row_lno_t col_b_end = row_mapB[rowb + 1];

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

          row_lno_t col_b_begin = row_mapB[rowb];
          const row_lno_t col_b_end = row_mapB[rowb + 1];

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
    row_lno_t m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
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
        row_lno_t m_,
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
          MaxRoughNonZero(16384), requiredMem(16384),
          multicore_chunk_size(0)
          //,KOKKOSKERNELS_SPGEMM_HASHSIZE(16),BLOCKSIZE (block_size)
          {}

    SymbolicC(
        row_lno_t m_,
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
          multicore_chunk_size(multicore_chunk_size_)
          {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTagCount&, const team_member_t & teamMember) const {




      row_lno_t row_index = teamMember.league_rank() * multicore_chunk_size;
      const int row_max = KOKKOSKERNELS_MACRO_MIN(row_index + multicore_chunk_size, m);
      //printf("myRowMin:%d rowMax:%d\n", row_index, row_max);


      const size_t hash_index_array_size = MaxRoughNonZero * sizeof(c_nnz_lno_t);
      const size_t hash_begin_array_size = hash_size * sizeof(c_nnz_lno_t);
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);
      c_nnz_lno_t *hash_nexts = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *hash_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_index_array_size);
      c_nnz_lno_t *used_hash_indices = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_begin_array_size);

      for (size_t i = 0; i < hash_size; ++i){
        hash_nexts[i] = hash_begins[i] = -1;
      }
#ifdef HASHTRACK
      size_t numhashes = 0;
      size_t numpoints = 0;
#endif
      for (; row_index < row_max; ++row_index){
        c_nnz_lno_t used_hash_counts = 0;
        c_nnz_lno_t next_available_entry = 0;

        //check ii is out of range. if it is, just return.
        row_lno_t col_begin = row_mapA[row_index];
        const row_lno_t col_end = row_mapA[row_index + 1];
        for (; col_begin < col_end; ++col_begin){
          //printf("col_begin:%d\n", col_begin);
          const nnz_lno_t rowb = entriesA[col_begin];

          //std::cout << "  for column:" << rowb << std::endl;
          row_lno_t colb_ind = row_mapB[rowb];
          const row_lno_t col_b_end = row_mapB[rowb + 1];
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
              if (hash_ind == -1){
                //std::cout << "\t\tinserting for the first time" << std::endl;
                hash_begins[hash] = next_available_entry;
                hash_indices[next_available_entry] = b_col_set_index;
                hash_values[next_available_entry++] = b_col_set;
                used_hash_indices[used_hash_counts++] = hash;
              }
              else {
                while (hash_nexts[hash_ind] != -1){
                  if(b_col_set_index ^ hash_indices[hash_ind]){ //if they are not equal
                    hash_ind = hash_nexts[hash_ind];
                  }
                  else {
                    hash_values[hash_ind] = hash_values[hash_ind] | b_col_set;
                    goto endloop; //break;
                  }
                }
                if(b_col_set_index ^ hash_indices[hash_ind]){ //if they are not equal
                  hash_nexts[hash_ind] = next_available_entry;
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
#ifdef HASHTRACK
          ++numhashes;
#endif
          const int hash = used_hash_indices[i];
          c_nnz_lno_t hash_ind = hash_begins[hash];

          hash_begins[hash] = -1;
          do {
#ifdef HASHTRACK
            ++numpoints;
#endif
            c_nnz_lno_t c_rows = hash_values[hash_ind];
            c_nnz_lno_t newhash_ind = hash_nexts[hash_ind];
            hash_nexts[hash_ind] = -1;
            hash_ind = newhash_ind;
            for (; c_rows; my_nonzeros++)
            {
              c_rows &= c_rows - 1; // clear the least significant bit set
            }
          }
          while (hash_ind != -1);
        }
        rowmapC(row_index) = my_nonzeros;
      }
#ifdef HASHTRACK
      printf("numhash:%ld np:%ld\n", numhashes, numpoints);
#endif

    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {
      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      row_lno_t row_index = (teamMember.league_rank()  * team_size + team_rank) * KOKKOSKERNELS_SPGEMM_BLOCKSIZE;

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
        const row_lno_t col_begin = row_mapA[row_index];
        const row_lno_t col_end = row_mapA[row_index + 1];

        for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
          const nnz_lno_t rowb = entriesA[col_a_index];
          row_lno_t colb_ind = row_mapB[rowb];
          const row_lno_t col_b_end = row_mapB[rowb + 1];
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
      row_lno_t row_index = teamMember.league_rank()  * team_size + team_rank;
      if (row_index >= m) return;


      int hash_array_size = KOKKOSKERNELS_SPGEMM_HASHSIZE * sizeof(c_nnz_lno_t) * team_size;
      c_nnz_lno_t *hash_begins = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_size);
      c_nnz_lno_t *hash_ends = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(hash_array_size);

      int set_size = (16384 - hash_array_size * 2) / 2 - 8;
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
      const row_lno_t col_begin = row_mapA[row_index];
      const row_lno_t col_end = row_mapA[row_index + 1];

      c_nnz_lno_t  current_heap_size = 0;
      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t rowb = entriesA[col_a_index];
        const row_lno_t col_b_begin = row_mapB[rowb];
        const row_lno_t col_b_end = row_mapB[rowb + 1];
        for (row_lno_t colb_ind = col_b_begin; colb_ind < col_b_end; ++colb_ind){
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



  template <
      typename out_row_view_t,
      typename out_nnz_view_t,
      typename in_row_map_t,
      typename in_nnz_t>
  struct unzipMatrix{
    typedef typename out_row_view_t::non_const_type non_const_c_lno_row_view_t;
    row_lno_t m;
    in_row_map_t &in_row_map;
    in_nnz_t &in_set_index_entries;
    in_nnz_t &in_set_entries;
    out_row_view_t &out_rowmap;
    out_nnz_view_t &out_entries;
    int set_size;

    unzipMatrix(row_lno_t m_,
                in_row_map_t row_map_c_copy_,
                in_nnz_t c_set_index_entries_,
                in_nnz_t c_set_entries_,
                out_row_view_t rowmapC_,
                out_row_view_t entriesC_): m(m_), in_row_map(row_map_c_copy_),
                    in_set_index_entries(c_set_index_entries_), in_set_entries(c_set_entries_),
                    out_rowmap(rowmapC_), out_entries(entriesC_), set_size(sizeof(typename in_nnz_t::value_type) * 8) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {
      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      row_lno_t row_index = teamMember.league_rank()  * team_size + team_rank;

      if (row_index >= m) return;
      //check ii is out of range. if it is, just return.
      const row_lno_t col_begin = in_row_map[row_index];
      const row_lno_t col_end = in_row_map[row_index + 1];
      row_lno_t c = 0;
      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t c_rows = in_set_entries[col_a_index];
        while (c_rows){
          c_rows &= (c_rows-1) ;
          c++;
        }
      }
      out_rowmap(row_index) = c;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const FillTag&, const team_member_t & teamMember) const {

      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      row_lno_t row_index = teamMember.league_rank()  * team_size + team_rank;

      if (row_index >= m) return;
      //check ii is out of range. if it is, just return.
      const row_lno_t col_begin = in_row_map[row_index];
      const row_lno_t col_end = in_row_map[row_index + 1];
      row_lno_t c = 0;
      row_lno_t out_row_map_index = out_rowmap(row_index);

      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t c_rows = in_set_entries[col_a_index];
        nnz_lno_t c_rows_set_index = in_set_index_entries[col_a_index];
        int current_row = 0;
        nnz_lno_t unit = 1;
        /*
        if (row_index == 0)
        std::cout << "col_begin:" << col_begin <<" col_end:" << col_end
                  <<" c_rows_set_index:" << c_rows_set_index << " c_rows:" << c_rows << std::endl;
        */
        while (c_rows){
          if (c_rows & unit){
            out_entries(out_row_map_index + c++) = set_size * c_rows_set_index + current_row;
            /*
            if (row_index == 0){
                std::cout << "\tcol:" << set_size * c_rows_set_index + current_row << std::endl;
            }
            */
          }
          current_row++;
          c_rows = c_rows & ~unit;
          unit = unit << 1;
        }
      }
    }
  };


  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename c_nnz_view_t>
  void Fill_C_Graph(
      row_lno_t m,
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

    std::cout << "chunksize:" << m / num_teams + 1 << std::endl;

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

    Kokkos::parallel_for(
        multicore_team_count_policy_t(num_teams, teamSizeMax, vector_size), fnnz);
    MyExecSpace::fence();

    std::cout << "\tActual NNZ Count TIME:" << timer1.seconds() << std::endl;


    timer1.reset();
    row_lno_t c_nnz = 0;
    KokkosKernels::Experimental::Util::view_reduce_max<non_const_c_lno_row_view_t, MyExecSpace>(m, rowmapC, c_nnz);
    std::cout << "\t\tmax is:" << c_nnz << std::endl;
    std::cout << "\tmax reduction:" << timer1.seconds() << std::endl;

    this->handle->get_spgemm_handle()->set_num_result_nnz(c_nnz);

    timer1.reset();
    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_row_view_t, MyExecSpace>(m+1, rowmapC);
    MyExecSpace::fence();
    std::cout << "\tpps reduction:" << timer1.seconds() << std::endl;
    //rowmapC = row_map_c_copy;
    /*
    auto d_c_nnz_size = Kokkos::subview(row_map_c_copy, m);
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

  template <typename c_lno_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_apply(c_lno_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC){

    Kokkos::Impl::Timer timer1;
    auto d_c_nnz_size = Kokkos::subview(rowmapC_, rowmapC_.dimension_0() - 1);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    size_t c_nnz_size = h_c_nnz_size();

    std::cout << "\tBefore Apply:C NNZ SIZE:" << c_nnz_size << std::endl;

    valuesC = c_scalar_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("C values"), c_nnz_size);
    entriesC_ = c_lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("C entries"), c_nnz_size);

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
    size_t min_pow2_hash_size = 1;
    const row_lno_t max_non_zero = this->handle->get_spgemm_handle()->get_num_result_nnz();
    while (max_non_zero > min_pow2_hash_size){
      min_pow2_hash_size *= 2;
    }

    min_pow2_hash_size *= 2;
    std::cout << "max_non_zero:" << max_non_zero << std::endl;

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
    );


    int vector_size = 1, teamSizeMax = 1;
    Kokkos::parallel_for( multicore_team_count_policy_t(num_teams, teamSizeMax, vector_size), cc);
    MyExecSpace::fence();

    std::cout << "\tApply TIME:" << timer1.seconds() << std::endl;


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

  }

  template <typename c_lno_row_view_t, typename c_lno_nnz_view_t>
  void KokkosSPGEMM_symbolic(c_lno_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_){

    row_lno_temp_work_view_t new_row_mapB;
    row_lno_temp_work_view_t set_index_entries;
    row_lno_temp_work_view_t set_entries;
    //First Compress B.
    Kokkos::Impl::Timer timer1;

    this->compressMatrix(this->row_mapB, this->entriesB, new_row_mapB, set_index_entries, set_entries);
    std::cout << "\n\n\tCOMPRESSION TIME:" << timer1.seconds() << std::endl;
    std::cout << "Old NNZ:" << this->entriesB.dimension_0() << " NEW NNZ:" << set_entries.dimension_0() << std::endl;

    /*
    //CHECK if the compression is correct.
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
    typedef typename c_lno_row_view_t::non_const_type non_const_c_lno_row_view_t;
    typedef typename c_lno_nnz_view_t::non_const_type non_const_c_lno_nnz_view_t;

    non_const_c_lno_row_view_t rowmapC (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
    non_const_c_lno_nnz_view_t entriesC;

    //First get a rough count.
    typedef typename c_lno_row_view_t::const_value_type const_c_row_lno_t;
    const_c_row_lno_t maxNumRoughZeros = this->getMaxRoughRowLen(m, row_mapA, entriesA, new_row_mapB, rowmapC);
    std::cout << "\tmaxNumRoughZeros:" << maxNumRoughZeros  << std::endl;
    std::cout << "\tmaxTime:" << timer1.seconds()  << std::endl;

    this->Fill_C_Graph(m, row_mapA, entriesA,
                        new_row_mapB, set_index_entries, set_entries,
                        rowmapC, entriesC, maxNumRoughZeros);

    rowmapC_ = rowmapC;
    entriesC_ = entriesC;
  }


  template <typename in_row_view_t, typename in_nnz_view_t, typename out_view_t>
  void compressMatrix(
      in_row_view_t in_row_map,
      in_nnz_view_t in_entries,

      out_view_t &out_row_map,
      out_view_t &out_nnz_indices,
      out_view_t &out_nnz_sets){

    //number of rows
    row_lno_t n = in_row_map.dimension_0() - 1;
    typedef typename out_view_t::non_const_value_type lno_t;

    //size of the lno_t, how many bits it can hold.
    int lnot_size = sizeof(lno_t) * 8;
    int compression_bit_divide_shift_ = 0;
    int val = lnot_size;
    while (val > 1) {
      ++compression_bit_divide_shift_;
      val = val >> 1;
    }

    lno_t compression_bit_mask_ = 0;
    for (int i = 0; i < compression_bit_divide_shift_; ++i){
      compression_bit_mask_ = (compression_bit_mask_ << 1) | 1;
    }
    /*
    std::cout << "lnot_size:" << lnot_size << " compression_bit_divide_shift_:" << compression_bit_divide_shift_
                  << " compression_bit_mask_:" << compression_bit_mask_ << std::endl;
                  */

    lno_t new_nnz_cnt = 0;
    out_row_map = out_view_t (Kokkos::ViewAllocateWithoutInitializing("new row map"), n+1);
    MyExecSpace::fence();

    zipMatrix <const_b_lno_row_view_t, const_b_lno_nnz_view_t, out_view_t, out_view_t>
        cmc(in_row_map, in_entries, compression_bit_mask_, compression_bit_divide_shift_,out_row_map, out_nnz_indices, out_nnz_sets);


    Kokkos::parallel_reduce( my_count_exec_space(0, n), cmc, new_nnz_cnt);

    MyExecSpace::fence();
    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum <out_view_t, MyExecSpace> (n+1, out_row_map);
    MyExecSpace::fence();


    out_nnz_indices = out_view_t (Kokkos::ViewAllocateWithoutInitializing("set index entries"), new_nnz_cnt);
    out_nnz_sets = out_view_t (Kokkos::ViewAllocateWithoutInitializing("set entries"), new_nnz_cnt);

    cmc.set_index_entries = out_nnz_indices;
    cmc.set_entries = out_nnz_sets;
    MyExecSpace::fence();
    Kokkos::parallel_for( my_fill_exec_space(0, n), cmc);
    MyExecSpace::fence();

  }

  template <typename out_row_view_t, typename out_nnz_view_t, typename in_row_view_t, typename in_nnz_view_t>
  void uncompressMatrix(
      out_row_view_t out_row_map, out_nnz_view_t &out_entries,
      in_row_view_t in_row_map,
      in_nnz_view_t in_nnz_indices,
      in_nnz_view_t in_nnz_sets){

    row_lno_t n = in_row_map.dimension_0() - 1;
    row_lno_t nnz = in_nnz_indices.dimension_0() - 1;
    //out_row_map = out_row_view_t(Kokkos::ViewAllocateWithoutInitializing("out_row_view_t out_row_map") , n + 1);
    MyExecSpace::fence();
    int teamSizeMax = 0;
    int vector_size = 0;

    unzipMatrix<out_row_view_t, out_nnz_view_t, in_row_view_t, in_nnz_view_t>
        ucc(n, in_row_map, in_nnz_indices, in_nnz_sets, out_row_map, out_entries);

    int max_allowed_team_size = 1;//team_policy_t::team_size_max(cnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size
        <row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n,nnz);
    //Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), cnnnz);


    Kokkos::parallel_for( team_count_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), ucc);
    MyExecSpace::fence();
    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<out_row_view_t, MyExecSpace>(n+1, out_row_map);
    MyExecSpace::fence();

    auto d_nnz_size = Kokkos::subview(out_row_map, n);
    auto h_nnz_size = Kokkos::create_mirror_view (d_nnz_size);
    Kokkos::deep_copy (h_nnz_size, d_nnz_size);
    MyExecSpace::fence();
    row_lno_t nnz_size = h_nnz_size();
    out_entries = out_row_view_t(Kokkos::ViewAllocateWithoutInitializing("out_entries ") , nnz_size);
    MyExecSpace::fence();
    Kokkos::parallel_for( team_fill_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), ucc);
    MyExecSpace::fence();
  }





};
}
}
}
}

#endif
