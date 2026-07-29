// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_POTRF_IMPL_HPP_
#define KOKKOSLAPACK_POTRF_IMPL_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {
namespace Impl {

// Implementation struct for potrf
template <class AViewType>
struct PotrfImpl {
  static void potrf([[maybe_unused]] const char uplo[], [[maybe_unused]] AViewType& A) {
    // TODO: Implement?
    KokkosKernels::Impl::throw_runtime_exception(
        "KokkosLapack::potrf: Not yet implemented. Enable LAPACK or CUSOLVER|ROCSOLVER");
  }
};

}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_POTRF_IMPL_HPP_
