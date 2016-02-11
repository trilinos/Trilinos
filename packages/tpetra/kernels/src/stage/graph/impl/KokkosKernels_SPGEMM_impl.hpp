
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

  template <typename c_lno_row_view_t>
  struct Predict_C_NUM_NNZ{

    row_lno_t m;
    const_a_lno_row_view_t row_mapA;
    const_a_lno_nnz_view_t entriesA;

    const_b_lno_row_view_t row_mapB;
    c_lno_row_view_t rough_row_mapC;


    Predict_C_NUM_NNZ(
        row_lno_t m_,
        const_a_lno_row_view_t row_mapA_,
        const_a_lno_nnz_view_t entriesA_,

        const_b_lno_row_view_t row_mapB_,
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
      rough_row_mapC[row_index] = max_num_results_in_row;
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


};
}
}
}
}

#endif
