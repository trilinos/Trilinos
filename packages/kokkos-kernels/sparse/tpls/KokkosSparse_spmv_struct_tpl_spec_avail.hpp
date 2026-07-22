// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_SPMV_STRUCT_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPMV_STRUCT_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class XVector, class YVector>
struct spmv_struct_tpl_spec_avail {
  enum : bool { value = false };
};

// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class XVector, class YVector>
struct spmv_mv_struct_tpl_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSPARSE_SPMV_STRUCT_TPL_SPEC_AVAIL_HPP_
