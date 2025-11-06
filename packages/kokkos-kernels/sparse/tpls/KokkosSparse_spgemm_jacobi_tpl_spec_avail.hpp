// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_SPGEMM_JACOBI_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPGEMM_JACOBI_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class b_size_view_t_,
          class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t, class c_scalar_view_t,
          class dinv_scalar_view_t>
struct spgemm_jacobi_tpl_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif
