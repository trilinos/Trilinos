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
#ifndef __KOKKOSBATCHED_EIGENDECOMPOSITION_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_EIGENDECOMPOSITION_SERIAL_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Eigendecomposition_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========
template <typename AViewType, typename EViewType, typename UViewType, typename WViewType>
KOKKOS_INLINE_FUNCTION int SerialEigendecomposition::invoke(const AViewType &A, const EViewType &er,
                                                            const EViewType &ei, const UViewType &UL,
                                                            const UViewType &UR, const WViewType &W) {
  /// view checking
  const int m = A.extent(0);
  assert(m == int(A.extent(1)) && "Eigendecomposition: A is not square");
  assert(m == int(er.extent(0)) && "Eigendecomposition: Length of er does not match to A's dimension");
  assert(m == int(ei.extent(0)) && "Eigendecomposition: Length of ei does not match to A's dimension");
  assert(m == int(UL.extent(0)) && "Eigendecomposition: Length of UL does not match to A's dimension");
  assert(m == int(UL.extent(1)) && "Eigendecomposition: Width of UL does not match to A's dimension");
  assert(m == int(UR.extent(0)) && "Eigendecomposition: Length of UR does not match to A's dimension");
  assert(m == int(UR.extent(1)) && "Eigendecomposition: Width of UR does not match to A's dimension");
  // assert(int(W.extent(0)) >= int(2*m*m+5*m) && "Eigendecomposition: workspace
  // size is too small");
  assert(int(W.stride(0)) == int(1) && "Eigendecomposition: Provided workspace is not contiguous");

  /// static assert A,er,ei,UL,UR,W has the same value_type
  /// static assert all views have the same memory space
  return m ? SerialEigendecompositionInternal ::invoke(
                 A.extent(0), A.data(), A.stride(0), A.stride(1), er.data(), er.stride(0), ei.data(), ei.stride(0),
                 UL.data(), UL.stride(0), UL.stride(1), UR.data(), UR.stride(0), UR.stride(1), W.data(), W.extent(0))
           : 0;
}

}  // namespace KokkosBatched

#endif
