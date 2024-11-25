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
#ifndef KOKKOSBATCHED_IAMAX_HPP_
#define KOKKOSBATCHED_IAMAX_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Iamax:
/// Iamax finds the index of the first element having maximum absolute value.
///
/// \tparam XViewType: Input view type, needs to be a 1D view
///
/// \param X [in]: Input view type
///
/// \return The index of the first element having maximum absolute value
/// As well as Blas, this returns 0 (0 in Fortran) for an empty vector
/// No nested parallel_for is used inside of the function.
///

struct SerialIamax {
  template <typename XViewType>
  KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const XViewType &x);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Iamax_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_IAMAX_HPP_
