// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_LACGV_HPP_
#define KOKKOSBATCHED_LACGV_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Lacgv:
/// Conjugates the elements of a complex vector x.
/// No operation is performed if x is real.
///
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \param x [in]: x is a length n vector, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///
struct SerialLacgv {
  template <typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Lacgv_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_LACGV_HPP_
