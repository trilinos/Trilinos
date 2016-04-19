
#ifndef _KOKKOSSPGEMMIMPL_HPP
#define _KOKKOSSPGEMMIMPL_HPP

#include <KokkosKernels_Utils.hpp>
#include <Kokkos_UnorderedMap.hpp>

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

  template <
            typename a_row_view_t, typename a_nnz_view_t,
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
      row_lno_t col_begin = row_mapA[row_index];
      row_lno_t col_end = row_mapA[row_index + 1];

      typename c_lno_row_view_t::value_type max_num_results_in_row = 0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
          [&] (row_lno_t i, typename c_lno_row_view_t::non_const_value_type & valueToUpdate) {
        row_lno_t adjind = i + col_begin;
        row_lno_t colIndex = entriesA[adjind];
        valueToUpdate += row_mapB [colIndex + 1] - row_mapB[colIndex];
      },
      max_num_results_in_row);

      Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
        rough_row_mapC[row_index] = max_num_results_in_row;
      });
    }
  };

  template <typename c_lno_row_view_t>
  size_t getNumRoughNonZeros(c_lno_row_view_t row_mapC){

    int teamSizeMax = 0;
    int vector_size = 0;



    Predict_C_NUM_NNZ<c_lno_row_view_t> pcnnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        row_mapC);

    int max_allowed_team_size = team_policy_t::team_size_max(pcnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, m, entriesA.dimension_0());

    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), pcnnnz);
    MyExecSpace::fence();

    typename c_lno_row_view_t::non_const_value_type rough_size = 0;
    KokkosKernels::Experimental::Util::view_reduce_sum<c_lno_row_view_t, MyExecSpace>(m+1, row_mapC, rough_size);
    return rough_size;
    /*
    {
    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_lno_row_view_t, MyExecSpace>(m+1, row_mapC);
    MyExecSpace::fence();
    auto d_rough_size = Kokkos::subview(row_mapC, m);
    auto h_rough_size = Kokkos::create_mirror_view (d_rough_size);
    Kokkos::deep_copy (h_rough_size, d_rough_size);
    typename c_lno_row_view_t::non_const_value_type rough_size = h_rough_size();
    return rough_size;
    }
    */
  }


  template <typename c_lno_row_view_t, typename hashmap_t>
  struct CorrectNumNNZ{

    row_lno_t m;
    const_a_lno_row_view_t row_mapA;
    const_a_lno_nnz_view_t entriesA;

    const_b_lno_row_view_t row_mapB;
    const_b_lno_nnz_view_t entriesB;

    c_lno_row_view_t rowmapC;
    hashmap_t hashmap;

    typedef typename const_a_lno_row_view_t::non_const_value_type c_row_lno_t;
    typedef typename const_b_lno_nnz_view_t::non_const_value_type c_nnz_lno_t;

    CorrectNumNNZ(
        row_lno_t m_,
        const_a_lno_row_view_t row_mapA_,
        const_a_lno_nnz_view_t entriesA_,

        const_b_lno_row_view_t row_mapB_,
        const_b_lno_nnz_view_t entriesB_,
        c_lno_row_view_t rowmapC_,
        hashmap_t hashmap_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_mapB(row_mapB_), entriesB(entriesB_),
          rowmapC(rowmapC_),
          hashmap(hashmap_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {
      row_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //check ii is out of range. if it is, just return.
      if (row_index >= m) return;
      row_lno_t col_begin = row_mapA[row_index];
      row_lno_t col_end = row_mapA[row_index + 1];

      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t rowb = entriesA[col_a_index];
        row_lno_t col_b_begin = row_mapB[rowb];
        row_lno_t col_b_end = row_mapB[rowb + 1];

        //std::cout << "col_b_begin:" << col_b_begin<< " col_b_end:" << col_b_end << std::endl;
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, col_b_end - col_b_begin),
            [&] (row_lno_t col_b_index) {
          row_lno_t b_adjind = col_b_index + col_b_begin;

          //std::cout << "b_adjind:" << b_adjind<< " col_b_end:" << col_b_end << " entriesB.dimension_0():" << entriesB.dimension_0() << std::endl;
          row_lno_t b_col = entriesB[b_adjind];

          Kokkos::UnorderedMapInsertResult r = hashmap.insert (Kokkos::pair<row_lno_t, c_nnz_lno_t>(row_index, b_col));
          if (r.success()){
            Kokkos::atomic_fetch_add(&(rowmapC(row_index)),1);
          }
        });
      }
    }
  };

  template <typename c_lno_row_view_t, typename c_lno_nnz_view_t, typename hashmap_t>
  struct FillCNNZ{
    row_lno_t m;
    const_a_lno_row_view_t row_mapA;
    const_a_lno_nnz_view_t entriesA;

    const_b_lno_row_view_t row_mapB;
    const_b_lno_nnz_view_t entriesB;

    c_lno_row_view_t rowmapC;
    c_lno_nnz_view_t entriesC;
    hashmap_t hashmap;

    typedef typename c_lno_row_view_t::non_const_value_type c_row_lno_t;
    typedef typename c_lno_nnz_view_t::non_const_value_type c_nnz_lno_t;

    FillCNNZ(
        row_lno_t m_,
        const_a_lno_row_view_t row_mapA_,
        const_a_lno_nnz_view_t entriesA_,

        const_b_lno_row_view_t row_mapB_,
        const_b_lno_nnz_view_t entriesB_,
        c_lno_row_view_t rowmapC_,
        c_lno_nnz_view_t entriesC_,
        hashmap_t hashmap_): m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_mapB(row_mapB_), entriesB(entriesB_),
          rowmapC(rowmapC_), entriesC(entriesC_),
          hashmap(hashmap_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {
      row_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //check ii is out of range. if it is, just return.
      if (row_index >= m) return;
      row_lno_t col_begin = row_mapA[row_index];
      row_lno_t col_end = row_mapA[row_index + 1];

      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t rowb = entriesA[col_a_index];
        row_lno_t col_b_begin = row_mapB[rowb];
        row_lno_t col_b_end = row_mapB[rowb + 1];


        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, col_b_end - col_b_begin),
            [&] (row_lno_t col_b_index) {
          row_lno_t b_adjind = col_b_index + col_b_begin;
          row_lno_t b_col = entriesB[b_adjind];
          Kokkos::UnorderedMapInsertResult r = hashmap.insert (Kokkos::pair<row_lno_t, c_nnz_lno_t>(row_index, b_col));
          if (r.success()){
            c_row_lno_t index_on_nnz = Kokkos::atomic_fetch_add(&(rowmapC(row_index)),1);
            entriesC[index_on_nnz] = b_col;
          }
        });
      }
    }
  };

  template <typename c_lno_row_view_t, typename c_lno_nnz_view_t>
  void Fill_C_Graph(
      typename c_lno_row_view_t::const_value_type numRoughNonzeros,
      c_lno_row_view_t &rowmapC,
      c_lno_nnz_view_t &entriesC_){

    typedef typename c_lno_row_view_t::non_const_type non_const_c_lno_row_view_t;
    non_const_c_lno_row_view_t row_map_c_copy("non_const_lnow_row", m + 1);
    //Kokkos::deep_copy (rough_row_map_c, rowmapC);

    typedef typename const_b_lno_nnz_view_t::non_const_value_type c_nnz_lno_t;

    typedef Kokkos::UnorderedMap< Kokkos::pair<row_lno_t, c_nnz_lno_t> , void , MyExecSpace> hashmap_t;
    hashmap_t umap(numRoughNonzeros);

    int teamSizeMax = 0;
    int vector_size = 0;

    CorrectNumNNZ<c_lno_row_view_t, hashmap_t> cnnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        entriesB,
        row_map_c_copy,
        umap);

    int max_allowed_team_size = team_policy_t::team_size_max(cnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n, entriesB.dimension_0());
    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), cnnnz);
    MyExecSpace::fence();

    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_lno_row_view_t, MyExecSpace>(m+1, row_map_c_copy);
    MyExecSpace::fence();
    auto d_c_nnz_size = Kokkos::subview(row_map_c_copy, m);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    typename c_lno_row_view_t::non_const_value_type c_nnz_size = h_c_nnz_size();


    typedef typename c_lno_nnz_view_t::non_const_type non_const_c_lno_nnz_view_t;

    umap = hashmap_t (c_nnz_size);
    non_const_c_lno_nnz_view_t entriesC(Kokkos::ViewAllocateWithoutInitializing("EntriesC"), c_nnz_size);
    Kokkos::deep_copy(rowmapC, row_map_c_copy);
    MyExecSpace::fence();



    FillCNNZ<c_lno_row_view_t, c_lno_nnz_view_t, hashmap_t> fcnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        entriesB,
        row_map_c_copy,
        entriesC,
        umap);
    max_allowed_team_size = team_policy_t::team_size_max(fcnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n, entriesB.dimension_0());
    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), fcnnz);
    MyExecSpace::fence();
    entriesC_ = entriesC;
  }

  /*
  template <typename c_lno_row_view_t, typename c_lno_nnz_view_t>
  void KokkosSPGEMM_symbolic(c_lno_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_){
    typedef typename c_lno_row_view_t::non_const_type non_const_c_lno_row_view_t;
    //typedef typename c_lno_nnz_view_t::non_const_type non_const_c_lno_nnz_view_t;
    non_const_c_lno_row_view_t rowmapC (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
    //First get a rough count.
    //typedef typename c_lno_row_view_t::non_const_value_type c_row_lno_t;
    typedef typename c_lno_row_view_t::const_value_type const_c_row_lno_t;
    //typedef typename c_lno_nnz_view_t::non_const_value_type c_nnz_lno_t;
    //typedef typename c_lno_nnz_view_t::const_value_type const_c_nnz_lno_t;
    const_c_row_lno_t numRoughNonzeros = this->getNumRoughNonZeros(rowmapC);

    this->Fill_C_Graph(numRoughNonzeros,rowmapC,entriesC_);
    rowmapC_ = rowmapC;
  }
  */

  template <typename row_view_t, typename nnz_view_t, typename new_row_view_t>
  struct CompressMatrixCount{
    const row_view_t row_map;
    const nnz_view_t entries;
    new_row_view_t new_row_map;
    const int compression_bit_divide_shift;
    typedef typename row_view_t::non_const_value_type index_t;
    typedef typename nnz_view_t::non_const_value_type lno_t;

    CompressMatrixCount(row_view_t row_map_,nnz_view_t entries_,
        int compression_bit_divide_shift_,
        new_row_view_t new_row_map_):
      row_map(row_map_), entries(entries_), compression_bit_divide_shift(compression_bit_divide_shift_),
      new_row_map(new_row_map_){
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const lno_t & i, lno_t & overall_size) const {
      const index_t rowBegin = row_map(i);
      const index_t rowEnd = row_map(i + 1);
      lno_t prev_n_set = -1;
      lno_t neighbor_set_count = 0;
      //std::cout << "row:" << i << std::endl;
      for (index_t ii = rowBegin; ii < rowEnd; ++ii){

        lno_t n = entries(ii);
        lno_t n_set = n >> compression_bit_divide_shift;
        neighbor_set_count += (prev_n_set != n_set);
        //std::cout << "\t" << " col:" << n << " n_set:" << n_set << std::endl;
        prev_n_set = n_set;
      }

      //std::cout << "\t neighbor_set_count:" << neighbor_set_count << std::endl << std::endl;
      new_row_map(i) = neighbor_set_count;
      overall_size += neighbor_set_count;
    }
  };

  template <typename row_view_t, typename nnz_view_t, typename new_row_view_t, typename new_nnz_view_t>
  struct CompressMatrixFill{
    typedef typename row_view_t::non_const_value_type index_t;
    typedef typename nnz_view_t::non_const_value_type lno_t;
    row_view_t row_map;
    nnz_view_t entries;
    const lno_t compression_bit_mask;
    const int compression_bit_divide_shift;
    new_row_view_t new_row_map;
    new_nnz_view_t set_index_entries;
    new_nnz_view_t set_entries;



    CompressMatrixFill(
        row_view_t row_map_,nnz_view_t entries_,
        lno_t compression_bit_mask_, int compression_bit_divide_shift_,
        new_row_view_t new_row_map_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_ ):
      row_map(row_map_), entries(entries_),
      compression_bit_mask(compression_bit_mask_), compression_bit_divide_shift(compression_bit_divide_shift_),
      new_row_map(new_row_map_), set_index_entries(set_index_entries_), set_entries(set_entries_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t & i) const {

      index_t newrowBegin = new_row_map(i);

      const index_t rowBegin = row_map(i);
      const index_t rowEnd = row_map(i + 1);

      index_t ii = rowBegin;

      lno_t col_set_index = 0;

      lno_t n = entries[ii++];
      lno_t prev_n_set_index = n >> compression_bit_divide_shift;

      lno_t prev_n_set = 1; lno_t n_mask = n & compression_bit_mask;  prev_n_set = prev_n_set << n_mask;

      //std::cout << "--- brow:" << i << " col:" << n << " prev_n_set:" << prev_n_set << std::endl;


      for (; ii < rowEnd; ++ii){
        lno_t n = entries[ii];




        lno_t n_set_index = n >> compression_bit_divide_shift;

        /*
        if (i == 2){
          lno_t set_mark = 1; lno_t n_mask = n & compression_bit_mask;  set_mark = set_mark << n_mask;
          std::cout
          << "- brow:" << i
          << "  col:" << n
          << "  prev_n_set:" << prev_n_set
          << "  prev_n_set_index:" << prev_n_set_index
          << " n_set_index:" << n_set_index
          << " set_mark:" << set_mark
          << std::endl;

        }
        */



        if (n_set_index != prev_n_set_index){
          set_entries(newrowBegin) = prev_n_set;
          set_index_entries(newrowBegin++) = prev_n_set_index;
          /*
          if (i == 2){
            std::cout << "--- brow:" << i
                      << " prev_n_set:" << prev_n_set
                      << " prev_n_set_index:" << prev_n_set_index
                      << std::endl;
          }
          */
          prev_n_set_index = n_set_index;
          prev_n_set = 0;
        }
        lno_t set_mark = 1; lno_t n_mask = n & compression_bit_mask;  set_mark = set_mark << n_mask;
        prev_n_set = prev_n_set | set_mark;

        /*
        std::cout
            << "-- brow:" << i
            << " col:" << n
            << " prev_n_set:" << prev_n_set << " set_mark:" << set_mark
            << " prev_n_set_index:" << prev_n_set_index << std::endl;
            */


      }
      set_entries(newrowBegin) = prev_n_set;
      set_index_entries(newrowBegin) = prev_n_set_index;
      /*
      if (i == 2){
        std::cout << "--- brow:" << i
            << " prev_n_set:" << prev_n_set
            << " prev_n_set_index:" << prev_n_set_index
            << std::endl;
      }
      */
    }
  };


  void compressB(
      row_lno_temp_work_view_t &new_row_mapB,
      row_lno_temp_work_view_t &set_index_entries,
      row_lno_temp_work_view_t &set_entries){

    typedef typename const_b_lno_nnz_view_t::non_const_value_type lno_t;

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

    //compression_bit_mask_ = ((size_t) compression_bit_mask_) >> (lnot_size - compression_bit_divide_shift_);

    std::cout << "lnot_size:" << lnot_size << " compression_bit_divide_shift_:" << compression_bit_divide_shift_
              << " compression_bit_mask_:" << compression_bit_mask_ << std::endl;

    lno_t new_nnz_cnt = 0;

    new_row_mapB = row_lno_temp_work_view_t (Kokkos::ViewAllocateWithoutInitializing("new row map"), n+1);
    CompressMatrixCount <const_b_lno_row_view_t, const_b_lno_nnz_view_t, row_lno_temp_work_view_t>
          cmc(row_mapB, entriesB, compression_bit_divide_shift_,new_row_mapB);
    Kokkos::parallel_reduce( my_exec_space(0, n), cmc, new_nnz_cnt);
    MyExecSpace::fence();

    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum <row_lno_temp_work_view_t, MyExecSpace> (n+1, new_row_mapB);
    MyExecSpace::fence();
    KokkosKernels::Experimental::Util::print_1Dview(new_row_mapB);

    set_index_entries = row_lno_temp_work_view_t (Kokkos::ViewAllocateWithoutInitializing("set index entries"), new_nnz_cnt);
    set_entries = row_lno_temp_work_view_t (Kokkos::ViewAllocateWithoutInitializing("set entries"), new_nnz_cnt);

    CompressMatrixFill <const_b_lno_row_view_t, const_b_lno_nnz_view_t, row_lno_temp_work_view_t, row_lno_temp_work_view_t>
          cmf(row_mapB, entriesB, compression_bit_mask_, compression_bit_divide_shift_,new_row_mapB, set_index_entries, set_entries);
    Kokkos::parallel_for( my_exec_space(0, n), cmf);
    std::cout << "nnz:" << entriesB.dimension_0() << " new:" << new_nnz_cnt << std::endl;
  }


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
    KokkosKernels::Experimental::Util::view_reduce_max<c_lno_row_view_t, MyExecSpace>(m+1, row_mapC, rough_size);

    MyExecSpace::fence();
    //KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_lno_row_view_t, MyExecSpace>(m+1, row_mapC);
    //MyExecSpace::fence();
    return rough_size;
  }




  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t>
  struct DetermineCSize{
    row_lno_t m;
    a_row_view_t row_mapA;
    a_nnz_view_t entriesA;

    b_row_view_t row_mapB;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;


    typedef typename a_row_view_t::non_const_value_type c_row_lno_t;
    typedef typename b_nnz_view_t::non_const_value_type c_nnz_lno_t;


    DetermineCSize(
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
          rowmapC(rowmapC_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {
      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      row_lno_t row_index = teamMember.league_rank()  * team_size + team_rank;

      if (row_index >= m) return;
      int my_local_work_array_size = 16384 / team_size / 2;
      size_t double_size = team_size * sizeof(c_nnz_lno_t);

      c_nnz_lno_t * team_shared_keys = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(16384 / 2);
      c_nnz_lno_t * team_shared_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(16384 / 2);
      team_shared_keys += team_rank;
      team_shared_values += team_rank;

      //check ii is out of range. if it is, just return.
      const row_lno_t col_begin = row_mapA[row_index];
      const row_lno_t col_end = row_mapA[row_index + 1];

      c_nnz_lno_t  current_heap_size = 0;
      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t rowb = entriesA[col_a_index];
        const row_lno_t col_b_begin = row_mapB[rowb];
        const row_lno_t col_b_end = row_mapB[rowb + 1];

        //std::cout << "row:" << row_index << " rowb:" << rowb << " col_b_begin:" << col_b_begin << " col_b_end:" << col_b_end  << std::endl;
        for (row_lno_t colb_ind = col_b_begin; colb_ind < col_b_end; ++colb_ind){
          c_nnz_lno_t b_col_set_index = entriesSetIndicesB[colb_ind];
          c_nnz_lno_t b_col_set = entriesSetsB[colb_ind];
          //std::cout << "row:" << row_index << " col:" << b_col_set_index << " current_heap_size:" << current_heap_size << std::endl;
          insert_to_heap(b_col_set_index, b_col_set, team_shared_keys, team_shared_values, current_heap_size);
        }
      }
      rowmapC(row_index) = current_heap_size;
      //if (current_heap_size < 8) std::cout << "row_index:" << row_index << " current_heap_size:" << current_heap_size<< std::endl;
    }

    KOKKOS_INLINE_FUNCTION
    void insert_to_heap(
        c_nnz_lno_t &key, c_nnz_lno_t &val,
        c_nnz_lno_t *keys, c_nnz_lno_t *vals,
        c_nnz_lno_t &current_size) const {
      for (int i = 0; i < current_size; ++i){
        if (key == keys[i]){
          vals[i] = vals[i] | val;
          return;
        }
      }
      keys[current_size] = key;
      vals[current_size++] = val;
    }

    // Provide the shared memory capacity.
    // This function takes the team_size as an argument ,
    // which allows team_size dependent allocations.
    size_t team_shmem_size (int team_size) const {
      return 16384;
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

      /*
      if (row_index == 2){
        std::cout << "row:" << row_index << " colsize:" << col_end - col_begin << std::endl;
      }
      */
      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t c_rows = in_set_entries[col_a_index];
        nnz_lno_t c_rows_set_index = in_set_index_entries[col_a_index];
        int current_row = 0;
        nnz_lno_t unit = 1;
        while (c_rows){
          if (c_rows & unit){

            /*
            if (row_index == 2){
              std::cout << "row:" << row_index << " colsize:" << col_end - col_begin
                        << " col_a_index:" << col_a_index
                        << " c_rows:" << c_rows
                        << " out_row_map_index:" << out_row_map_index
                        << " c:" << c
                        << " c_rows_set_index:" << c_rows_set_index
                        << " current_row:" << current_row
                        << " val:" << set_size * c_rows_set_index + current_row
                        << std::endl;
            }
            */

            out_entries(out_row_map_index + c++) = set_size * c_rows_set_index + current_row;
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
  struct FillC{
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


    FillC(
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
          rowmapC(rowmapC_){}



    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {

      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      row_lno_t row_index = teamMember.league_rank()  * team_size + team_rank;
      if (row_index >= m) return;

      int my_local_work_array_size = 16384 / team_size / 2;
      size_t double_size = team_size * sizeof(c_nnz_lno_t);

      c_nnz_lno_t * team_shared_keys = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(16384 / 2);
      c_nnz_lno_t * team_shared_values = (c_nnz_lno_t *) teamMember.team_shmem().get_shmem(16384 / 2);
      team_shared_keys += team_rank;
      team_shared_values += team_rank;

      //check ii is out of range. if it is, just return.
      const row_lno_t col_begin = row_mapA[row_index];
      const row_lno_t col_end = row_mapA[row_index + 1];

      c_nnz_lno_t  current_heap_size = 0;
      for (row_lno_t col_a_index = col_begin; col_a_index < col_end; ++col_a_index){
        nnz_lno_t rowb = entriesA[col_a_index];
        const row_lno_t col_b_begin = row_mapB[rowb];
        const row_lno_t col_b_end = row_mapB[rowb + 1];

        //std::cout << "row:" << row_index << " rowb:" << rowb << " col_b_begin:" << col_b_begin << " col_b_end:" << col_b_end  << std::endl;
        for (row_lno_t colb_ind = col_b_begin; colb_ind < col_b_end; ++colb_ind){
          c_nnz_lno_t b_col_set_index = entriesSetIndicesB[colb_ind];
          c_nnz_lno_t b_col_set = entriesSetsB[colb_ind];
          /*
          std::cout << "row:" << row_index
                    << " col:" << b_col_set_index << " b_col_set:" << b_col_set
                    << " current_heap_size:" << current_heap_size << std::endl;
                    */
          insert_to_heap(b_col_set_index, b_col_set, team_shared_keys, team_shared_values, current_heap_size);
        }
      }



      c_nnz_lno_t  current_heap_ind = 0;
      c_nnz_lno_t adjInd = rowmapC(row_index);

      c_nnz_lno_t my_nonzeros = 0;

      for (c_nnz_lno_t current_heap_ind = 0; current_heap_ind < current_heap_size; ++current_heap_ind){
        entriesSetIndicesC(current_heap_ind + adjInd) = team_shared_keys[current_heap_ind];
        entriesSetsC(current_heap_ind + adjInd) = team_shared_values[current_heap_ind];

        /*
        if (current_heap_ind + adjInd < 18 && current_heap_ind + adjInd >= 12){
          std::cout << "current_heap_ind + adjInd:" << current_heap_ind + adjInd
                    << "row_index:" << row_index
                    << std::endl;
        }
        */
        /*
        if (row_index == 2){
              std::cout << "row_index:" << row_index
                        << " adjInd:" << adjInd
                        << " setindex:" << team_shared_keys[current_heap_ind]
                        << " set:" << team_shared_values[current_heap_ind]
                        << std::endl;
        }
        */

        c_nnz_lno_t c_rows = team_shared_values[current_heap_ind];
        nnz_lno_t unit = 1;
        while (c_rows){
          if (c_rows & unit){
            my_nonzeros++;
          }
          c_rows = c_rows & ~unit;
          unit = unit << 1;
        }
      }

      //if (current_heap_size < 8) std::cout << "row_index:" << row_index << " current_heap_size:" << current_heap_size<< std::endl;

    }

    KOKKOS_INLINE_FUNCTION
    void insert_to_heap(c_nnz_lno_t &key, c_nnz_lno_t &val, c_nnz_lno_t *keys, c_nnz_lno_t *vals, c_nnz_lno_t &current_size) const {
      for (int i = 0; i < current_size; ++i){
        if (key == keys[i]){
          vals[i] = vals[i] | val;
          return;
        }
      }
      keys[current_size] = key;
      vals[current_size++] = val;
    }

    // Provide the shared memory capacity.
    // This function takes the team_size as an argument ,
    // which allows team_size dependent allocations.
    size_t team_shmem_size (int team_size) const {
      return 16384;
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
      typename c_row_view_t::const_value_type numRoughNonzeros
      ){

    typedef typename c_row_view_t::non_const_type non_const_c_lno_row_view_t;
    non_const_c_lno_row_view_t row_map_c_copy(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
    //Kokkos::deep_copy (rough_row_map_c, rowmapC);

    typedef typename const_b_lno_nnz_view_t::non_const_value_type c_nnz_lno_t;

    int teamSizeMax = 0;
    int vector_size = 0;

    DetermineCSize<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t> cnnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        entriesSetIndex,
        entriesSets,
        row_map_c_copy);

    int max_allowed_team_size = 1;//team_policy_t::team_size_max(cnnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n, entriesB.dimension_0());
    //Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), cnnnz);

    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), cnnnz);
    MyExecSpace::fence();


    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_row_view_t, MyExecSpace>(m+1, row_map_c_copy);
    MyExecSpace::fence();



    //KokkosKernels::Experimental::Util::print_1Dview(row_map_c_copy);

    auto d_c_nnz_size = Kokkos::subview(row_map_c_copy, m);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    typename c_row_view_t::non_const_value_type c_nnz_size = h_c_nnz_size();


    std::cout << "c_nnz_compressed_size:" << c_nnz_size << std::endl;
    row_lno_temp_work_view_t c_set_index_entries(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"),c_nnz_size);
    row_lno_temp_work_view_t c_set_entries(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), c_nnz_size);
    FillC<a_row_view_t, a_nnz_view_t, b_row_view_t, b_nnz_view_t, c_row_view_t, c_nnz_view_t> fnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        entriesSetIndex,
        entriesSets,
        row_map_c_copy);

    fnnz.entriesSetIndicesC = c_set_index_entries;
    fnnz.entriesSetsC = c_set_entries;
    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), fnnz);
    //Kokkos::parallel_for( team_policy_t(3 , teamSizeMax, vector_size), fnnz);

    MyExecSpace::fence();

    this->uncompressMatrix(rowmapC, entriesC_, row_map_c_copy, c_set_index_entries, c_set_entries);

    //KokkosKernels::Experimental::Util::print_1Dview(c_set_index_entries);

    /*
    UnCompressC<c_row_view_t, c_nnz_view_t>
        ucc(m, row_map_c_copy, c_set_index_entries, c_set_entries, rowmapC, entriesC_);
    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), ucc);
    MyExecSpace::fence();

    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum<c_row_view_t, MyExecSpace>(m+1, rowmapC);
    MyExecSpace::fence();

    KokkosKernels::Experimental::Util::print_1Dview(rowmapC);

    auto d_actual_c_nnz_size = Kokkos::subview(rowmapC, m);
    auto h_actual_c_nnz_size = Kokkos::create_mirror_view (d_actual_c_nnz_size);
    Kokkos::deep_copy (h_actual_c_nnz_size, d_actual_c_nnz_size);
    typename c_row_view_t::non_const_value_type c_actual_nnz_size = h_actual_c_nnz_size();
    std::cout << "c_actual_nnz_size:" << c_actual_nnz_size << std::endl;
    */

    /*
    typedef typename c_nnz_view_t::non_const_type non_const_c_lno_nnz_view_t;

    umap = hashmap_t (c_nnz_size);
    non_const_c_lno_nnz_view_t entriesC(Kokkos::ViewAllocateWithoutInitializing("EntriesC"), c_nnz_size);
    Kokkos::deep_copy(rowmapC, row_map_c_copy);
    MyExecSpace::fence();



    FillCNNZ<c_row_view_t, c_nnz_view_t, hashmap_t> fcnnz(
        m,
        row_mapA,
        entriesA,
        row_mapB,
        entriesB,
        row_map_c_copy,
        entriesC,
        umap);
    max_allowed_team_size = team_policy_t::team_size_max(fcnnz);
    KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, MyExecSpace>(max_allowed_team_size, vector_size, teamSizeMax, n, entriesB.dimension_0());
    Kokkos::parallel_for( team_policy_t(m / teamSizeMax + 1 , teamSizeMax, vector_size), fcnnz);
    MyExecSpace::fence();
    entriesC_ = entriesC;
    */
  }


  template <typename c_lno_row_view_t, typename c_lno_nnz_view_t>
  void KokkosSPGEMM_symbolic(c_lno_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_){


    row_lno_temp_work_view_t new_row_mapB;
    row_lno_temp_work_view_t set_index_entries;
    row_lno_temp_work_view_t set_entries;
    //First Compress B.
    this->compressMatrix(this->row_mapB, this->entriesB, new_row_mapB, set_index_entries, set_entries);



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



    typedef typename c_lno_row_view_t::non_const_type non_const_c_lno_row_view_t;
    //typedef typename c_lno_nnz_view_t::non_const_type non_const_c_lno_nnz_view_t;
    non_const_c_lno_row_view_t rowmapC (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);

    //First get a rough count.
    typedef typename c_lno_row_view_t::const_value_type const_c_row_lno_t;
    const_c_row_lno_t maxNumRoughZeros = this->getMaxRoughRowLen(m, row_mapA, entriesA, new_row_mapB, rowmapC);
    std::cout << "maxNumRoughZeros:" << maxNumRoughZeros  << std::endl;

    this->Fill_C_Graph(m, row_mapA, entriesA, new_row_mapB, set_index_entries, set_entries, rowmapC, entriesC_, maxNumRoughZeros);



  }


  template <typename in_row_view_t, typename in_nnz_view_t, typename out_view_t>
  void compressMatrix(
      in_row_view_t in_row_map, in_nnz_view_t in_entries,
      out_view_t &out_row_map,
      out_view_t &out_nnz_indices,
      out_view_t &out_nnz_sets){

    row_lno_t n = in_row_map.dimension_0() - 1;
    typedef typename in_nnz_view_t::non_const_value_type lno_t;
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
    std::cout << "lnot_size:" << lnot_size << " compression_bit_divide_shift_:" << compression_bit_divide_shift_
                  << " compression_bit_mask_:" << compression_bit_mask_ << std::endl;

    lno_t new_nnz_cnt = 0;
    out_row_map = out_view_t (Kokkos::ViewAllocateWithoutInitializing("new row map"), n+1);
    MyExecSpace::fence();

    CompressMatrixCount <in_row_view_t, in_nnz_view_t, out_view_t>
        cmc(in_row_map, in_entries, compression_bit_divide_shift_,out_row_map);

    Kokkos::parallel_reduce( my_exec_space(0, n), cmc, new_nnz_cnt);
    MyExecSpace::fence();

    KokkosKernels::Experimental::Util::exclusive_parallel_prefix_sum <out_view_t, MyExecSpace> (n+1, out_row_map);
    MyExecSpace::fence();
    KokkosKernels::Experimental::Util::print_1Dview(out_row_map);

    out_nnz_indices = out_view_t (Kokkos::ViewAllocateWithoutInitializing("set index entries"), new_nnz_cnt);
    out_nnz_sets = out_view_t (Kokkos::ViewAllocateWithoutInitializing("set entries"), new_nnz_cnt);
    MyExecSpace::fence();
    CompressMatrixFill <const_b_lno_row_view_t, const_b_lno_nnz_view_t, out_view_t, out_view_t>
          cmf(in_row_map, in_entries, compression_bit_mask_, compression_bit_divide_shift_,out_row_map, out_nnz_indices, out_nnz_sets);
    Kokkos::parallel_for( my_exec_space(0, n), cmf);
    MyExecSpace::fence();

    std::cout << "nnz:" << in_entries.dimension_0() << " new:" << new_nnz_cnt << std::endl;
  }

  template <typename out_row_view_t, typename out_nnz_view_t, typename in_row_view_t, typename in_nnz_view_t>
  void uncompressMatrix(
      out_row_view_t &out_row_map, out_nnz_view_t &out_entries,
      in_row_view_t in_row_map,
      in_nnz_view_t in_nnz_indices,
      in_nnz_view_t in_nnz_sets){

    row_lno_t n = in_row_map.dimension_0() - 1;
    row_lno_t nnz = in_nnz_indices.dimension_0() - 1;
    out_row_map = out_row_view_t(Kokkos::ViewAllocateWithoutInitializing("out_row_view_t out_row_map") , n + 1);
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

    /*
    {
      auto nnz_sets = Kokkos::subview(in_nnz_sets, Kokkos::make_pair(12,18));
      auto nnz_indices = Kokkos::subview(in_nnz_indices, Kokkos::make_pair(12,18));

      KokkosKernels::Experimental::Util::print_1Dview(nnz_indices);
      KokkosKernels::Experimental::Util::print_1Dview(nnz_sets);

    }
    */

    auto d_nnz_size = Kokkos::subview(out_row_map, n);
    auto h_nnz_size = Kokkos::create_mirror_view (d_nnz_size);
    Kokkos::deep_copy (h_nnz_size, d_nnz_size);
    MyExecSpace::fence();
    row_lno_t nnz_size = h_nnz_size();
    std::cout << "nnz_size:" << nnz_size << std::endl;
    out_entries = out_row_view_t(Kokkos::ViewAllocateWithoutInitializing("out_entries ") , nnz_size);
    MyExecSpace::fence();
    Kokkos::parallel_for( team_fill_policy_t(n / teamSizeMax + 1 , teamSizeMax, vector_size), ucc);
    MyExecSpace::fence();
    //KokkosKernels::Experimental::Util::print_1Dview(out_row_map);
  }





};
}
}
}
}

#endif
