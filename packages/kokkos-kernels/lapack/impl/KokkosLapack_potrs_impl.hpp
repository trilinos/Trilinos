// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_POTRS_IMPL_HPP_
#define KOKKOSLAPACK_POTRS_IMPL_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {
namespace Impl {

// Implementation struct for potrs
template <class AViewType, class BViewType>
struct PotrsImpl {
  static void potrs([[maybe_unused]] const char uplo[], [[maybe_unused]] const AViewType& A,
                    [[maybe_unused]] BViewType& B) {
    // TODO: Implement in the future?
    KokkosKernels::Impl::throw_runtime_exception(
        "KokkosLapack::potrs: Not yet implemented. Enable LAPACK or CUSOLVER|ROCSOLVER");
  }
};

}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_POTRS_IMPL_HPP_
