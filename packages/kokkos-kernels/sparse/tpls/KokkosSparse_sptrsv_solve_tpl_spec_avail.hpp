// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_SPTRSV_SOLVE_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPTRSV_SOLVE_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class KernelHandle, class RowMapType, class EntriesType, class ValuesType, class BType,
          class XType>
struct sptrsv_solve_tpl_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif
