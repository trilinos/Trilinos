/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_,
          typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_,
          typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t,
          typename c_scalar_nnz_view_t>
void KokkosBSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_,
                   a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
                   b_scalar_nnz_view_t_>::
    KokkosBSPGEMM_numeric(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_,
                          c_scalar_nnz_view_t &valuesC_) {
  // get the algorithm and execution space.
  // SPGEMMAlgorithm spgemm_algorithm =
  // this->handle->get_spgemm_handle()->get_algorithm_type();
  KokkosKernels::Impl::ExecSpaceType my_exec_space_ =
      KokkosKernels::Impl::get_exec_space_type<MyExecSpace>();

  if (Base::KOKKOSKERNELS_VERBOSE) {
    std::cout << "Numeric PHASE" << std::endl;
  }

  if (Base::spgemm_algorithm == SPGEMM_KK_SPEED ||
      Base::spgemm_algorithm == SPGEMM_KK_DENSE) {
    this->KokkosBSPGEMM_numeric_speed(rowmapC_, entriesC_, valuesC_,
                                      my_exec_space_);
  } else {
    this->KokkosBSPGEMM_numeric_hash(rowmapC_, entriesC_, valuesC_,
                                     my_exec_space_);
  }
}

}  // namespace Impl
}  // namespace KokkosSparse
