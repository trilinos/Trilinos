// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_LIGHTWEIGHTVIEW_HPP
#define IFPACK2_DETAILS_LIGHTWEIGHTVIEW_HPP

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <type_traits>

/// \file Ifpack2_Details_LightweightView.hpp
/// \brief A barebones stand-in for Kokkos::View.
///
/// LightweightView represents a minimal nonowning multidimensional array,
/// like an unmanaged Kokkos::View or an mdspan.
/// It uses 32-bit extents and strides, unlike the mdspan version of Kokkos::View which
/// always uses 64 bit. This does not depend on mdspan, so it works without C++23 and
/// with or without Kokkos_ENABLE_IMPL_MDSPAN.
/// TODO: when Kokkos requires C++23 and removes the legacy View impl, replace this
/// with std::mdspan (or a device-friendly version of it, if Kokkos provides that).

namespace Ifpack2 {
namespace Details {

/// \class LightweightView
/// \tparam ValueType Type of each entry.
/// \tparam Rank Number of (dynamic) dimensions.
/// \tparam Layout Memory layout: Kokkos::LayoutLeft, Kokkos::LayoutRight, or
///   Kokkos::LayoutStride.
///
/// The primary template is intentionally left undefined; only the layout
/// specializations below are usable.
template <typename ValueType, int Rank, typename Layout>
class LightweightView;

namespace Impl {

// Number of explicitly stored strides for each layout. LayoutStride stores all
// R strides. LayoutLeft fixes stride(0) == 1 and LayoutRight fixes
// stride(R-1) == 1, so each stores R-1 explicit strides.
template <typename Layout, int Rank>
struct NumStoredStrides;

template <int Rank>
struct NumStoredStrides<Kokkos::LayoutStride, Rank> {
  static constexpr int value = Rank;
};
template <int Rank>
struct NumStoredStrides<Kokkos::LayoutLeft, Rank> {
  static constexpr int value = (Rank > 0) ? (Rank - 1) : 0;
};
template <int Rank>
struct NumStoredStrides<Kokkos::LayoutRight, Rank> {
  static constexpr int value = (Rank > 0) ? (Rank - 1) : 0;
};

}  // namespace Impl

/// LayoutStride specialization: stores all R strides explicitly.
template <typename ValueType, int Rank>
class LightweightView<ValueType, Rank, Kokkos::LayoutStride> {
 public:
  using value_type   = ValueType;
  using array_layout = Kokkos::LayoutStride;
  using size_type    = uint32_t;

  static constexpr int rank_value = Rank;
  KOKKOS_INLINE_FUNCTION static constexpr int rank() { return Rank; }

  LightweightView() = default;

  /// \brief Construct from a Kokkos::View.
  ///
  /// Copies the data pointer, extents and strides out of \c v, narrowing the
  /// extents and strides to 32 bits.
  template <typename View,
            typename = std::enable_if_t<Kokkos::is_view<View>::value>>
  KOKKOS_INLINE_FUNCTION explicit LightweightView(const View& v)
    : data_(v.data()) {
    static_assert(static_cast<int>(View::rank()) == Rank,
                  "LightweightView: rank mismatch with source Kokkos::View");
    for (int i = 0; i < Rank; ++i) {
      extents_[i] = static_cast<size_type>(v.extent(i));
      strides_[i] = static_cast<size_type>(v.stride(i));
    }
  }

  KOKKOS_INLINE_FUNCTION
  value_type* data() const { return data_; }

  KOKKOS_INLINE_FUNCTION
  size_type extent(int i) const { return extents_[i]; }
  KOKKOS_INLINE_FUNCTION
  size_type stride(int i) const { return strides_[i]; }

  template <typename... Indices>
  KOKKOS_INLINE_FUNCTION value_type& operator()(Indices... idxs) const {
    static_assert(sizeof...(Indices) == Rank,
                  "LightweightView: wrong number of indices");
    return data_[offset<0>(idxs...)];
  }

 private:
  template <int Dim, typename Index, typename... Rest>
  KOKKOS_INLINE_FUNCTION size_type offset(Index i, Rest... rest) const {
    return static_cast<size_type>(i) * strides_[Dim] +
           offset<Dim + 1>(rest...);
  }
  template <int Dim>
  KOKKOS_INLINE_FUNCTION size_type offset() const {
    return 0u;
  }

  value_type* data_ = nullptr;
  Kokkos::Array<size_type, Rank> extents_{};
  Kokkos::Array<size_type, Rank> strides_{};
};

/// LayoutLeft specialization: leftmost stride is fixed to 1, so only strides
/// 1..R-1 are stored.
template <typename ValueType, int Rank>
class LightweightView<ValueType, Rank, Kokkos::LayoutLeft> {
 public:
  using value_type   = ValueType;
  using array_layout = Kokkos::LayoutLeft;
  using size_type    = uint32_t;

  static constexpr int num_stored_strides =
      Impl::NumStoredStrides<Kokkos::LayoutLeft, Rank>::value;

  static constexpr int rank_value = Rank;
  KOKKOS_INLINE_FUNCTION static constexpr int rank() { return Rank; }

  LightweightView() = default;

  template <typename View,
            typename = std::enable_if_t<Kokkos::is_view<View>::value>>
  KOKKOS_INLINE_FUNCTION explicit LightweightView(const View& v)
    : data_(v.data()) {
    static_assert(static_cast<int>(View::rank()) == Rank,
                  "LightweightView: rank mismatch with source Kokkos::View");
    for (int i = 0; i < Rank; ++i) {
      extents_[i] = static_cast<size_type>(v.extent(i));
    }
    // stride(0) is implicitly 1 for LayoutLeft; store the remaining strides.
    for (int i = 1; i < Rank; ++i) {
      strides_[i - 1] = static_cast<size_type>(v.stride(i));
    }
  }

  KOKKOS_INLINE_FUNCTION
  value_type* data() const { return data_; }

  KOKKOS_INLINE_FUNCTION
  size_type extent(int i) const { return extents_[i]; }
  KOKKOS_INLINE_FUNCTION
  size_type stride(int i) const {
    return (i == 0) ? size_type(1) : strides_[i - 1];
  }

  template <typename... Indices>
  KOKKOS_INLINE_FUNCTION value_type& operator()(Indices... idxs) const {
    static_assert(sizeof...(Indices) == Rank,
                  "LightweightView: wrong number of indices");
    return data_[offset<0>(idxs...)];
  }

 private:
  template <int Dim, typename Index, typename... Rest>
  KOKKOS_INLINE_FUNCTION size_type offset(Index i, Rest... rest) const {
    const size_type s = (Dim == 0) ? size_type(1) : strides_[Dim - 1];
    return static_cast<size_type>(i) * s + offset<Dim + 1>(rest...);
  }
  template <int Dim>
  KOKKOS_INLINE_FUNCTION size_type offset() const {
    return 0u;
  }

  value_type* data_ = nullptr;
  Kokkos::Array<size_type, Rank> extents_{};
  Kokkos::Array<size_type, (num_stored_strides > 0 ? num_stored_strides : 1)>
      strides_{};
};

/// LayoutRight specialization: rightmost stride is fixed to 1, so only strides
/// 0..R-2 are stored.
template <typename ValueType, int Rank>
class LightweightView<ValueType, Rank, Kokkos::LayoutRight> {
 public:
  using value_type   = ValueType;
  using array_layout = Kokkos::LayoutRight;
  using size_type    = uint32_t;

  static constexpr int num_stored_strides =
      Impl::NumStoredStrides<Kokkos::LayoutRight, Rank>::value;

  static constexpr int rank_value = Rank;
  KOKKOS_INLINE_FUNCTION static constexpr int rank() { return Rank; }

  LightweightView() = default;

  template <typename View,
            typename = std::enable_if_t<Kokkos::is_view<View>::value>>
  KOKKOS_INLINE_FUNCTION explicit LightweightView(const View& v)
    : data_(v.data()) {
    static_assert(static_cast<int>(View::rank()) == Rank,
                  "LightweightView: rank mismatch with source Kokkos::View");
    for (int i = 0; i < Rank; ++i) {
      extents_[i] = static_cast<size_type>(v.extent(i));
    }
    // stride(R-1) is implicitly 1 for LayoutRight; store the leading strides.
    for (int i = 0; i < Rank - 1; ++i) {
      strides_[i] = static_cast<size_type>(v.stride(i));
    }
  }

  /// \brief Construct from a raw pointer and explicit extents (row-major).
  ///
  /// Convenience constructor for wrapping a freshly allocated, contiguous
  /// scratch buffer. The strides are computed assuming a fully packed
  /// (unpadded) row-major LayoutRight layout; do not use this constructor to
  /// wrap a padded allocation.
  template <typename... Extents>
  KOKKOS_INLINE_FUNCTION LightweightView(value_type* ptr, Extents... exts)
    : data_(ptr) {
    static_assert(sizeof...(Extents) == Rank,
                  "LightweightView: wrong number of extents");
    const size_type e[Rank] = {static_cast<size_type>(exts)...};
    for (int i = 0; i < Rank; ++i) extents_[i] = e[i];
    size_type run = 1;
    for (int i = Rank - 1; i >= 1; --i) {
      run *= extents_[i];
      strides_[i - 1] = run;
    }
  }

  KOKKOS_INLINE_FUNCTION
  value_type* data() const { return data_; }

  KOKKOS_INLINE_FUNCTION
  size_type extent(int i) const { return extents_[i]; }
  KOKKOS_INLINE_FUNCTION
  size_type stride(int i) const {
    return (i == Rank - 1) ? size_type(1) : strides_[i];
  }

  template <typename... Indices>
  KOKKOS_INLINE_FUNCTION value_type& operator()(Indices... idxs) const {
    static_assert(sizeof...(Indices) == Rank,
                  "LightweightView: wrong number of indices");
    return data_[offset<0>(idxs...)];
  }

 private:
  template <int Dim, typename Index, typename... Rest>
  KOKKOS_INLINE_FUNCTION size_type offset(Index i, Rest... rest) const {
    const size_type s = (Dim == Rank - 1) ? size_type(1) : strides_[Dim];
    return static_cast<size_type>(i) * s + offset<Dim + 1>(rest...);
  }
  template <int Dim>
  KOKKOS_INLINE_FUNCTION size_type offset() const {
    return 0u;
  }

  value_type* data_ = nullptr;
  Kokkos::Array<size_type, Rank> extents_{};
  Kokkos::Array<size_type, (num_stored_strides > 0 ? num_stored_strides : 1)>
      strides_{};
};

/// \brief Helper that converts a Kokkos::View to a matching LightweightView.
///
/// The resulting LightweightView uses the same layout as the source View, with
/// 32-bit extents and strides.
template <typename View>
KOKKOS_INLINE_FUNCTION auto view_to_lightweight(const View& v) {
  using value_type = typename View::value_type;
  using src_layout = typename View::array_layout;
  constexpr int R  = static_cast<int>(View::rank());
  using layout =
      std::conditional_t<std::is_same<src_layout, Kokkos::LayoutLeft>::value,
                         Kokkos::LayoutLeft,
                         std::conditional_t<std::is_same<src_layout, Kokkos::LayoutRight>::value,
                                            Kokkos::LayoutRight, Kokkos::LayoutStride>>;
  return LightweightView<value_type, R, layout>(v);
}

}  // namespace Details
}  // namespace Ifpack2

#endif  // IFPACK2_DETAILS_LIGHTWEIGHTVIEW_HPP
