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
#ifndef KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_
#define KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Functor that implements the single-vector, two-argument
///   version of KokkosBlas::dot (dot product of two vectors).
///
/// \tparam XVector Type of the first vector x; 1-D View
/// \tparam YVector Type of the second vector y; 1-D View
/// \tparam SizeType Type of the row index used in the dot product.
///   For best performance, use int instead of size_t here.
template <class execution_space, class AV, class XVector, class YVector, typename SizeType>
struct DotFunctor {
  typedef SizeType size_type;
  typedef typename AV::non_const_value_type avalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<avalue_type> IPT;
  typedef typename IPT::dot_type value_type;

  XVector m_x;
  YVector m_y;

  DotFunctor(const XVector& x, const YVector& y) : m_x(x), m_y(y) {}

  void run(const char* label, const execution_space& space, AV result) {
    Kokkos::RangePolicy<execution_space, size_type> policy(space, 0, m_x.extent(0));
    Kokkos::parallel_reduce(label, policy, *this, result);
  }

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void operator()(const size_type& i, value_type& sum) const {
    Kokkos::Details::updateDot(sum, m_x(i), m_y(i));  // sum += m_x(i) * m_y(i)
  }

  KOKKOS_INLINE_FUNCTION void init(value_type& update) const { update = Kokkos::ArithTraits<value_type>::zero(); }

  KOKKOS_INLINE_FUNCTION void join(value_type& update, const value_type& source) const { update += source; }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_
