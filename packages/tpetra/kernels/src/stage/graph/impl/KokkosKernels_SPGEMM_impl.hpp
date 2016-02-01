
#ifndef _KOKKOSSPGEMMIMPL_HPP
#define _KOKKOSSPGEMMIMPL_HPP


namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{

template <typename KernelHandle>
class KokkosSPGEMM{
public:

  typedef typename HandleType::idx_array_type idx_array_type;
  typedef typename HandleType::idx_edge_array_type idx_edge_array_type;
  typedef typename HandleType::value_array_type value_array_type;

  typedef typename idx_array_type::value_type idx;
  typedef typename idx_array_type::array_layout idx_array_layout;
  typedef typename idx_array_type::device_type idx_device_type;
  typedef typename idx_array_type::memory_traits idx_memory_traits;
  typedef typename idx_array_type::HostMirror host_view_type; //Host view type



  typedef typename idx_edge_array_type::value_type idx_edge;
  typedef typename idx_edge_array_type::array_layout idx_edge_array_layout;
  typedef typename idx_edge_array_type::device_type idx_edge_device_type;
  typedef typename idx_edge_array_type::memory_traits idx_edge_memory_traits;
  typedef typename idx_edge_array_type::HostMirror host_edge_view_type; //Host view type


  typedef typename value_array_type::value_type value_type;
  typedef typename value_array_type::array_layout value_type_array_layout;
  typedef typename value_array_type::device_type value_type_device_type;
  typedef typename value_array_type::memory_traits value_type_memory_traits;
  typedef typename value_array_type::HostMirror host_value_view_type; //Host view type

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


  typedef Kokkos::View<idx *, MyTempMemorySpace> idx_temp_work_array_type;
  typedef Kokkos::View<idx *, MyPersistentMemorySpace> idx_persistent_work_array_type;
  typedef typename idx_persistent_work_array_type::HostMirror host_idx_persistent_view_type; //Host view type


  typedef Kokkos::View<value_type *, MyTempMemorySpace> value_temp_work_array_type;
  typedef Kokkos::View<value_type *, MyPersistentMemorySpace> value_persistent_work_array_type;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy ;
  typedef typename team_policy::member_type team_member ;
private:
  KernelHandle *handle;
  idx m;
  idx n;
  idx k;
  idx_array_type row_mapA;
  idx_edge_array_type entriesA;
  value_array_type valsA;
  bool transposeA;
  idx_array_type row_mapB;
  idx_edge_array_type entriesB;
  value_array_type valsB;
  bool transposeB;
  idx_array_type row_mapC;
  idx_array_type entriesC;
  value_array_type valsC;

  KokkosSPGEMM(
      KernelHandle *handle_,
      idx m_,
      idx n_,
      idx k_,
      idx_array_type row_mapA_,
      idx_edge_array_type entriesA_,
      bool transposeA_,
      idx_array_type row_mapB_,
      idx_edge_array_type entriesB_,
      bool transposeB_):handle (handle_), m(m_), n(n_), k(k_),
          row_mapA(row_mapA_), entriesA(entriesA_), valsA(), transposeA(transposeA_),
          row_mapB(row_mapB_), entriesB(entriesB_), valsB(), transposeB(transposeB_),
          row_mapC(), entriesC(), valsC(){}

  KokkosSPGEMM(
        KernelHandle *handle_,
        idx m_,
        idx n_,
        idx k_,
        idx_array_type row_mapA_,
        idx_edge_array_type entriesA_,
        value_array_type valsA_,
        bool transposeA_,
        idx_array_type row_mapB_,
        idx_edge_array_type entriesB_,
        value_array_type valsB_,
        bool transposeB_):handle (handle_), m(m_), n(n_), k(k_),
            row_mapA(row_mapA_), entriesA(entriesA_), valsA(valsA_), transposeA(transposeA_),
            row_mapB(row_mapB_), entriesB(entriesB_), valsB(valsB_), transposeB(transposeB_),
            row_mapB(), entriesC(), valsC(){}


  void KokkosSPGEMM_symbolic(idx_array_type row_mapC){
    //assign rows of A to warps
  }

public:

  template <typename idxview, typename idxedgeview>
  struct SymbolicSPGEMM{
    idx m,n,k;
    idxview axadj;
    idxview bxadj;
    idxview cxadj;

    idxedgeview aadj;
    idxedgeview badj;

    SymbolicSPGEMM(
        idx m_,idx n_,idx k_,
        idxview axadj_,
        idxview bxadj_,
        idxview cxadj_,
        idxedgeview aadj_,
        idxedgeview badj_):
          m(m_),
          n(n_),
          k(k_),
          axadj(axadj_),
          bxadj(bxadj_),
          cxadj(cxadj_),
          aadj(aadj_),
          badj(badj_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member & teamMember) const {
      idx row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //check ii is out of range. if it is, just return.
      if (row_index >= m) return;



      idx col_begin = axadj[row_index];
      idx col_end = axadj[row_index + 1];

      idx max_num_results_in_row = 0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
          [&] (idx i, value_type & valueToUpdate) {
        idx adjind = i + col_begin;
        idx colIndex = aadj[adjind];
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
