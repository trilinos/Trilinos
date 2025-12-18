// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_TRSV_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_TRSV_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class CrsMatrixType, class DomainMultiVectorType, class RangeMultiVectorType>
struct trsv_tpl_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif
