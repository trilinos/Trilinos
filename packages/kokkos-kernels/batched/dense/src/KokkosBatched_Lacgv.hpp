//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
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
