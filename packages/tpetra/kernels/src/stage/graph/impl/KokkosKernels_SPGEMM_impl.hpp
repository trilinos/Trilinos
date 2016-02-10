
#ifndef _KOKKOSSPGEMMIMPL_HPP
#define _KOKKOSSPGEMMIMPL_HPP


namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{

template <typename HandleType,
  typename a_lno_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
  typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_
  //,typename c_lno_row_view_t_, typename c_lno_nnz_view_t_, typename c_scalar_nnz_view_t_
  >
class KokkosSPGEMM{
public:

  typedef a_lno_row_view_t_ a_in_lno_row_view_t;
  typedef a_lno_nnz_view_t_ a_in_lno_nnz_view_t;
  typedef a_scalar_nnz_view_t_ a_in_scalar_nnz_view_t;

  typedef b_lno_row_view_t_ b_in_lno_row_view_t;
  typedef b_lno_nnz_view_t_ b_in_lno_nnz_view_t;
  typedef b_scalar_nnz_view_t_ b_in_scalar_nnz_view_t;

  /*
  typedef c_lno_row_view_t_ c_in_lno_row_view_t;
  typedef c_lno_nnz_view_t_ c_in_lno_nnz_view_t;
  typedef c_scalar_nnz_view_t_ c_in_scalar_nnz_view_t;
  */


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





  /*
  typedef typename c_in_lno_row_view_t::const_type const_c_lno_row_view_t;
  typedef typename c_in_lno_row_view_t::non_const_type non_const_c_lno_row_view_t;

  typedef typename c_in_lno_nnz_view_t::const_type const_c_lno_nnz_view_t;
  typedef typename c_in_lno_nnz_view_t::non_const_type non_const_c_lno_nnz_view_t;

  typedef typename c_in_scalar_nnz_view_t::const_type const_c_scalar_nnz_view_t;
  typedef typename c_in_scalar_nnz_view_t::non_const_type non_const_c_scalar_nnz_view_t;
  */




  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t; //Host view type


  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef row_lno_t color_t;
  typedef Kokkos::View<row_lno_t *, MyTempMemorySpace> color_view_t;

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
  /*
  non_const_c_lno_row_view_t row_mapC;
  non_const_c_lno_nnz_view_t entriesC;
  non_const_c_scalar_nnz_view_t valsC;
  */

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

  template <typename c_lno_row_view_t_, typename c_lno_nnz_view_t_>
  void KokkosSPGEMM_symbolic(c_lno_row_view_t_ &rowmapC_, c_lno_nnz_view_t_ &entriesC_){
    typedef typename c_lno_row_view_t_::non_const_type non_const_c_lno_row_view_t;
    typedef typename c_lno_nnz_view_t_::non_const_type non_const_c_lno_nnz_view_t;

    non_const_c_lno_row_view_t rowmapC(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);

  }


  template <typename c_lno_row_view_t_, typename c_lno_nnz_view_t_>
  struct SymbolicSPGEMM{
    row_lno_t m,n,k;
    const_a_lno_row_view_t axadj;
    const_b_lno_row_view_t bxadj;
    c_lno_row_view_t_ cxadj;

    const_a_lno_nnz_view_t aadj;
    const_b_lno_nnz_view_t badj;
    c_lno_nnz_view_t_ cadj;

    SymbolicSPGEMM(
        row_lno_t m_,row_lno_t n_,row_lno_t k_,
        const_a_lno_row_view_t axadj_,
        const_b_lno_row_view_t bxadj_,
        c_lno_row_view_t_ cxadj_,
        const_a_lno_nnz_view_t aadj_,
        const_b_lno_nnz_view_t badj_):
          m(m_),
          n(n_),
          k(k_),
          axadj(axadj_),
          bxadj(bxadj_),
          cxadj(cxadj_),
          aadj(aadj_),
          badj(badj_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {
      row_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //check ii is out of range. if it is, just return.
      if (row_index >= m) return;



      row_lno_t col_begin = axadj[row_index];
      row_lno_t col_end = axadj[row_index + 1];

      row_lno_t max_num_results_in_row = 0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
          [&] (row_lno_t i, row_lno_t & valueToUpdate) {
        row_lno_t adjind = i + col_begin;
        row_lno_t colIndex = aadj[adjind];
        valueToUpdate += bxadj [colIndex + 1] - bxadj[colIndex];
      },
      max_num_results_in_row);
     }
  };

};
}
}
}
}

#endif
