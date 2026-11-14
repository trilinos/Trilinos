// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_
#define KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Functor that implements the single-vector, two-argument
///   version of KokkosBlas::dot (dot product of two vectors).
///
/// \tparam AV 0-D output View
/// \tparam XVector Type of the first vector x; 1-D View
/// \tparam YVector Type of the second vector y; 1-D View
/// \tparam SizeType Type of the row index used in the dot product.
///   For best performance, use int instead of size_t here.
template <class Exec, class AV, class XVector, class YVector, typename SizeType>
struct DotFunctor {
  typedef SizeType size_type;
  typedef typename AV::non_const_value_type avalue_type;
  typedef KokkosKernels::Details::InnerProductSpaceTraits<avalue_type> IPT;
  typedef typename IPT::dot_type dot_type;
  // Use a higher precision accumulator type when the element
  // type is float or complex<float>.
  typedef typename KokkosKernels::Details::DotAccumulatingScalar<dot_type>::type accum_type;
  // Define value_type for the functor so that we can run it without a Reducer argument.
  typedef accum_type value_type;

  AV m_result;
  XVector m_x;
  YVector m_y;

  void run(const char* label, const Exec& exec) const {
    Kokkos::RangePolicy<Exec, size_type> policy(exec, 0, m_x.extent(0));
    // Perform the parallel_reduce in higher-precision accumulator type,
    // using a temporary scalar.
    if constexpr (Kokkos::SpaceAccessibility<Exec, typename AV::memory_space>::accessible) {
      // If Exec can access m_result, then have the functor write its result directly into m_result in final()
      Kokkos::parallel_reduce(label, policy, *this);
    } else {
      // Otherwise, compute acc first and then write or copy it to m_result
      accum_type acc;
      Kokkos::parallel_reduce(label, policy, *this, acc);
      avalue_type result = KokkosKernels::Details::CastPossiblyComplex<avalue_type, accum_type>::cast(acc);
      if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace,
                                               typename AV::memory_space>::accessible) {
        m_result() = result;
      } else {
        Kokkos::deep_copy(exec, m_result, result);
      }
    }
  }

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void operator()(const size_type& i, accum_type& sum) const {
    KokkosKernels::Details::updateDot(sum, m_x(i), m_y(i));  // sum += m_x(i) * m_y(i)
  }

  KOKKOS_INLINE_FUNCTION void init(accum_type& update) const {
    update = KokkosKernels::ArithTraits<accum_type>::zero();
  }

  KOKKOS_INLINE_FUNCTION void join(accum_type& update, const accum_type& source) const { update += source; }

  KOKKOS_INLINE_FUNCTION void final(accum_type& result) const {
    if constexpr (Kokkos::SpaceAccessibility<Exec, typename AV::memory_space>::accessible) {
      // Can write the final result directly into m_result on device
      m_result() = KokkosKernels::Details::CastPossiblyComplex<avalue_type, accum_type>::cast(result);
    }
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_
