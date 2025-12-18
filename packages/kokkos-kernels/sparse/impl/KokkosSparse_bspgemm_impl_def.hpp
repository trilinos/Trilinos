// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosBSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                   b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::KokkosBSPGEMM_numeric(c_row_view_t &rowmapC_,
                                                                                   c_lno_nnz_view_t &entriesC_,
                                                                                   c_scalar_nnz_view_t &valuesC_) {
  // get the algorithm and execution space.
  // SPGEMMAlgorithm spgemm_algorithm =
  // this->handle->get_spgemm_handle()->get_algorithm_type();
  KokkosKernels::Impl::ExecSpaceType my_exec_space_ = KokkosKernels::Impl::get_exec_space_type<MyExecSpace>();

  if (Base::KOKKOSKERNELS_VERBOSE) {
    std::cout << "Numeric PHASE" << std::endl;
  }

  if (Base::spgemm_algorithm == SPGEMM_KK_SPEED || Base::spgemm_algorithm == SPGEMM_KK_DENSE) {
    this->KokkosBSPGEMM_numeric_speed(rowmapC_, entriesC_, valuesC_, my_exec_space_);
  } else {
    this->KokkosBSPGEMM_numeric_hash(rowmapC_, entriesC_, valuesC_, my_exec_space_);
  }
}

}  // namespace Impl
}  // namespace KokkosSparse
