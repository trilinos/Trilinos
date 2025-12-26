// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_SPILUK_SYMBOLIC_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPILUK_SYMBOLIC_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class ARowMapType, class AEntriesType, class LRowMapType, class LEntriesType,
          class URowMapType, class UEntriesType>
struct spiluk_symbolic_tpl_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif
