// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_ROTM_IMPL_HPP_
#define KOKKOSBATCHED_ROTM_IMPL_HPP_

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Rotm_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename XViewType, typename YViewType, typename ParamViewType>
KOKKOS_INLINE_FUNCTION static int checkRotmInput([[maybe_unused]] const XViewType &x,
                                                 [[maybe_unused]] const YViewType &y,
                                                 [[maybe_unused]] const ParamViewType &param) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::rotm: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::rotm: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<ParamViewType>, "KokkosBatched::rotm: ParamViewType is not a Kokkos::View.");
  static_assert(XViewType::rank() == 1, "KokkosBatched::rotm: XViewType must have rank 1.");
  static_assert(YViewType::rank() == 1, "KokkosBatched::rotm: YViewType must have rank 1.");
  static_assert(ParamViewType::rank() == 1, "KokkosBatched::rotm: ParamViewType must have rank 1.");
  static_assert(std::is_same_v<typename XViewType::value_type, typename XViewType::non_const_value_type>,
                "KokkosBatched::rotm: XViewType must have non-const value type.");
  static_assert(std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>,
                "KokkosBatched::rotm: YViewType must have non-const value type.");
  using x_value_type     = typename XViewType::non_const_value_type;
  using y_value_type     = typename YViewType::non_const_value_type;
  using param_value_type = typename ParamViewType::non_const_value_type;

  static_assert(!KokkosKernels::ArithTraits<x_value_type>::is_complex &&
                    !KokkosKernels::ArithTraits<y_value_type>::is_complex &&
                    !KokkosKernels::ArithTraits<param_value_type>::is_complex,
                "KokkosBatched::rotm: Complex types are not supported for Rotm.");

#ifndef NDEBUG
  const int n = x.extent_int(0);

  if (y.extent_int(0) != n) {
    Kokkos::printf(
        "KokkosBatched::rotm: x and y must have the same length: x length "
        "= "
        "%d, y length = %d\n",
        n, y.extent_int(0));
    return 1;
  }

  // param must have length 5: param(0) = flag, param(1) = h11, param(2) = h21, param(3) = h12, param(4) = h22
  if (param.extent_int(0) != 5) {
    Kokkos::printf("KokkosBatched::rotm: param must have length 5: param length = %d\n", param.extent_int(0));
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

///
/// Serial Impl
/// ===========

template <typename XViewType, typename YViewType, typename ParamViewType>
KOKKOS_INLINE_FUNCTION int SerialRotm::invoke(const XViewType &x, const YViewType &y, const ParamViewType &param) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0) return 0;

  auto info = Impl::checkRotmInput(x, y, param);
  if (info) return info;

  using ScalarType      = typename XViewType::non_const_value_type;
  const ScalarType flag = param(0);
  const ScalarType zero = KokkosKernels::ArithTraits<ScalarType>::zero();
  const ScalarType one  = KokkosKernels::ArithTraits<ScalarType>::one();
  const ScalarType two  = one + one;

  if (flag == -two) {
    // flag == -2.0: identity, no need to do anything
  } else if (flag == -one) {
    Impl::SerialRotmInternal<-1>::invoke(n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(),
                                         param.stride(0));
  } else if (flag == zero) {
    Impl::SerialRotmInternal<0>::invoke(n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(), param.stride(0));
  } else if (flag == one) {
    Impl::SerialRotmInternal<1>::invoke(n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(), param.stride(0));
  } else {
#ifndef NDEBUG
    // Invalid flag value
    Kokkos::printf(
        "KokkosBatched::SerialRotm: Invalid flag value. Flag must be one of -1, 0, 1, or -2. Flag value = %f\n",
        static_cast<double>(flag));
#endif
    return 1;
  }
  return 0;
}

///
/// Team Impl
/// ===========

template <typename MemberType>
template <typename XViewType, typename YViewType, typename ParamViewType>
KOKKOS_INLINE_FUNCTION int TeamRotm<MemberType>::invoke(const MemberType &member, const XViewType &x,
                                                        const YViewType &y, const ParamViewType &param) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0) return 0;

  auto info = Impl::checkRotmInput(x, y, param);
  if (info) return info;

  using ScalarType      = typename XViewType::non_const_value_type;
  const ScalarType flag = param(0);
  const ScalarType zero = KokkosKernels::ArithTraits<ScalarType>::zero();
  const ScalarType one  = KokkosKernels::ArithTraits<ScalarType>::one();
  const ScalarType two  = one + one;

  if (flag == -two) {
    // flag == -2.0: identity, no need to do anything
  } else if (flag == -one) {
    Impl::TeamRotmInternal<-1>::invoke(member, n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(),
                                       param.stride(0));
  } else if (flag == zero) {
    Impl::TeamRotmInternal<0>::invoke(member, n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(),
                                      param.stride(0));
  } else if (flag == one) {
    Impl::TeamRotmInternal<1>::invoke(member, n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(),
                                      param.stride(0));
  } else {
#ifndef NDEBUG
    // Invalid flag value
    Kokkos::printf(
        "KokkosBatched::TeamRotm: Invalid flag value. Flag must be one of -1, 0, 1, or -2. Flag value = %f\n",
        static_cast<double>(flag));
#endif
    return 1;
  }
  return 0;
}

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
template <typename XViewType, typename YViewType, typename ParamViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorRotm<MemberType>::invoke(const MemberType &member, const XViewType &x,
                                                              const YViewType &y, const ParamViewType &param) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0) return 0;

  auto info = Impl::checkRotmInput(x, y, param);
  if (info) return info;

  using ScalarType      = typename XViewType::non_const_value_type;
  const ScalarType flag = param(0);
  const ScalarType zero = KokkosKernels::ArithTraits<ScalarType>::zero();
  const ScalarType one  = KokkosKernels::ArithTraits<ScalarType>::one();
  const ScalarType two  = one + one;

  if (flag == -two) {
    // flag == -2.0: identity, no need to do anything
  } else if (flag == -one) {
    Impl::TeamVectorRotmInternal<-1>::invoke(member, n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(),
                                             param.stride(0));
  } else if (flag == zero) {
    Impl::TeamVectorRotmInternal<0>::invoke(member, n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(),
                                            param.stride(0));
  } else if (flag == one) {
    Impl::TeamVectorRotmInternal<1>::invoke(member, n, x.data(), x.stride(0), y.data(), y.stride(0), param.data(),
                                            param.stride(0));
  } else {
#ifndef NDEBUG
    // Invalid flag value
    Kokkos::printf(
        "KokkosBatched::TeamVectorRotm: Invalid flag value. Flag must be one of -1, 0, 1, or -2. Flag value = %f\n",
        static_cast<double>(flag));
#endif
    return 1;
  }
  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_ROTM_IMPL_HPP_
