// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_GAUSS_SEIDEL_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_GAUSS_SEIDEL_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t>
struct gauss_seidel_symbolic_tpl_spec_avail {
  enum : bool { value = false };
};

// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t>
struct gauss_seidel_numeric_tpl_spec_avail {
  enum : bool { value = false };
};

// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class x_scalar_view_t,
          class y_scalar_view_t>
struct gauss_seidel_apply_tpl_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif
