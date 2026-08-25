// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_LAYOUT_CONTIGUOUS_HPP
#define SACADO_FAD_KOKKOS_LAYOUT_CONTIGUOUS_HPP

// We are hooking into Kokkos Core internals here
// Need to define this macro since we include non-public headers
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Layout.hpp"
#include "Kokkos_DynRankView.hpp"
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

namespace Sacado {

// Contiguous layout for scalar types -- equivalent to the wrapped
// layout type
template <typename Layout, unsigned Stride = 1>
struct LayoutContiguous : public Layout {

  enum { scalar_stride = Stride };

  //! Tag this class as a kokkos array layout
  typedef LayoutContiguous array_layout ;

  // Pull in Layout's constructors
  using Layout::Layout;

  KOKKOS_INLINE_FUNCTION
  constexpr LayoutContiguous( Layout const & layout ) : Layout(layout) {}
  KOKKOS_INLINE_FUNCTION
  constexpr LayoutContiguous( Layout && layout ) : Layout(layout) {}

  KOKKOS_INLINE_FUNCTION
  Layout base_layout() const { return *this; }
};

// Is Layout == LayoutContiguous<L> for some L
template <class Layout>
struct is_layout_contiguous {
  static const bool value = false;
};

template <class Layout, unsigned Stride>
struct is_layout_contiguous< LayoutContiguous<Layout, Stride> > {
  static const bool value = true;
};

// Extract inner layout from LayoutContiguous
template <class Layout>
struct inner_layout {
  typedef Layout type;
};

template <class Layout>
struct inner_layout< LayoutContiguous<Layout> > {
  typedef Layout type;
};

template <class Layout, unsigned Stride>
struct inner_layout< LayoutContiguous<Layout, Stride> > {
  typedef Layout type;
};

namespace Impl {

template <typename Layout>
struct LayoutScalarStride {
  static const unsigned stride = 1;
  static const bool is_unit_stride = true;
};

template <typename Layout, unsigned Stride>
struct LayoutScalarStride< LayoutContiguous<Layout,Stride> > {
  static const unsigned stride = Stride;
  static const bool is_unit_stride = (Stride == 1);
};

} // namespace Impl

} // namespace Sacado

namespace Kokkos {

namespace Impl {

template <class Layout, unsigned Stride>
struct LayoutFromArrayLayout<Sacado::LayoutContiguous<Layout, Stride>> {
  using type = typename LayoutFromArrayLayout<Layout>::type;
};

// Specialize DynRankDimTraits for LayoutContiguous.
// Weirdly, we need a full specialization on bool for the non-legacy View impl case
// because of this code in Kokkos_DynRankView.hpp (around line 429):
// #ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
//   using drdtraits = Impl::DynRankDimTraits<typename view_type::specialize>;
// #else
//   using drdtraits = Impl::DynRankDimTraits<
//       std::conditional_t<view_type::traits::impl_is_customized, bool, void>>;
// #endif
template <>
struct DynRankDimTraits<bool> {
  using drdtraits = DynRankDimTraits<void>;
  enum : size_t { unspecified = drdtraits::unspecified };

  // Compute the rank of the view from the nonzero dimension arguments.
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank(const size_t N0, const size_t N1, const size_t N2,
                            const size_t N3, const size_t N4, const size_t N5,
                            const size_t N6, const size_t N7) {
    return drdtraits::computeRank(N0,N1,N2,N3,N4,N5,N6,N7);
  }

  // Compute the rank of the view from the nonzero layout arguments.
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION static size_t computeRank(const Layout& layout) {
    return drdtraits::computeRank(layout);
  }

  // Extra overload to match that for specialize types v2
  template <typename Layout, typename... P>
  KOKKOS_INLINE_FUNCTION static size_t computeRank(
      const Kokkos::Impl::ViewCtorProp<P...>& prop ,
      const Layout& layout) {
    return drdtraits::computeRank(prop, layout);
  }

  // Create the layout for the rank-7 view.
  // Because the underlying View is rank-7, preserve "unspecified" for
  // dimension 8.

  // Non-contiguous Layout
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<
      !(Sacado::is_layout_contiguous<Layout>::value),
      Layout>
  createLayout(const Layout& layout,
              size_t new_rank = unspecified) {
    return drdtraits::createLayout(layout, new_rank);
  }

  // Contiguous Layout
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<
      (Sacado::is_layout_contiguous<Layout>::value),
      Layout >
  createLayout(const Layout& layout,
              size_t new_rank = unspecified) {
    return Layout(drdtraits::createLayout(layout.base_layout(), new_rank));
  }

  // Extra overload to match that for specialize types
  template <typename Traits, typename... P>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<
      !(Sacado::is_layout_contiguous<typename Traits::array_layout>::value),
      typename Traits::array_layout>
  createLayout(const Kokkos::Impl::ViewCtorProp<P...>& prop,
              typename Traits::array_layout layout) {
    if constexpr (Traits::impl_is_customized &&
                  !Kokkos::Impl::ViewCtorProp<P...>::has_accessor_arg) {
      auto rank              = computeRank(prop, layout) - 1;
      layout.dimension[rank] = unspecified;
    }
    return createLayout(layout);
  }

  // Extra overload to match that for specialize types
  template <typename Traits, typename... P>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<
      (Sacado::is_layout_contiguous<typename Traits::array_layout>::value),
      typename Traits::array_layout>
  createLayout(const Kokkos::Impl::ViewCtorProp<P...>& prop,
              typename Traits::array_layout layout) {
    if constexpr (Traits::impl_is_customized &&
                  !Kokkos::Impl::ViewCtorProp<P...>::has_accessor_arg) {
      auto rank              = computeRank(prop, layout) - 1;
      layout.dimension[rank] = unspecified;
    }
    return createLayout(layout);
  }

  // Create a view from the given dimension arguments.
  // This is only necessary because the shmem constructor doesn't take a layout.
  //   NDE shmem View's are not compatible with the added view_alloc value_type
  //   / fad_dim deduction functionality
  template <typename ViewType, typename ViewArg>
  static ViewType createView(const ViewArg& arg, const size_t N0,
                            const size_t N1, const size_t N2, const size_t N3,
                            const size_t N4, const size_t N5, const size_t N6,
                            const size_t N7) {
    return drdtraits::createView(arg, N0, N1, N2, N3, N4, N5, N6, N7);
  }
};

// Overload reconstructLayout for LayoutContiguous
template <typename Layout, unsigned Stride, typename iType>
KOKKOS_INLINE_FUNCTION Sacado::LayoutContiguous<Layout,Stride>
reconstructLayout(const Sacado::LayoutContiguous<Layout,Stride>& layout, iType dynrank) {
  return Sacado::LayoutContiguous<Layout,Stride>(reconstructLayout(layout.base_layout(), dynrank));
}

} // namespace Impl

} // namespace Kokkos

// Inject some things into Kokkos for backwards compatibility
#ifdef SACADO_ENABLE_DEPRECATED_CODE
namespace Kokkos {

using Sacado::LayoutContiguous;
using Sacado::is_layout_contiguous;
using Sacado::inner_layout;
namespace Impl {
using Sacado::Impl::LayoutScalarStride;
}

}
#endif

#endif // #ifndef SACADO_FAD_KOKKOS_LAYOUT_CONTIGUOUS_HPP
