// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_LAYOUT_CONTIGUOUS_HPP
#define KOKKOS_LAYOUT_CONTIGUOUS_HPP

// We are hooking into Kokkos Core internals here
// Need to define this macro since we include non-public headers
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Layout.hpp"
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

namespace Kokkos {

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
};

// Is Layout == LayoutContiguous<L> for some L
template <class Layout>
struct is_layout_contiguous {
  static const bool value = false;
};

template <class Layout>
struct is_layout_contiguous< LayoutContiguous<Layout> > {
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

} // namespace Kokkos

// FIXME This is evil and needs refactoring urgently.
// Make LayoutContiguous<Layout> equivalent to Layout
namespace std {

  template <class Layout, unsigned Stride>
  struct is_same< Kokkos::LayoutContiguous<Layout,Stride>, Layout> {
    static const bool value = true;
  };

  template <class Layout, unsigned Stride>
#if defined(KOKKOS_COMPILER_INTEL)
  inline constexpr bool is_same_v< Kokkos::LayoutContiguous<Layout,Stride>, Layout> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
#else
  static constexpr bool is_same_v< Kokkos::LayoutContiguous<Layout,Stride>, Layout> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
#endif

  template <class Layout, unsigned Stride>
  struct is_same< Layout, Kokkos::LayoutContiguous<Layout,Stride> > {
    static const bool value = true;
  };

  template <class Layout, unsigned Stride>
#if defined(KOKKOS_COMPILER_INTEL)
  inline constexpr bool is_same_v< Layout, Kokkos::LayoutContiguous<Layout,Stride>> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
#else
  static constexpr bool is_same_v< Layout, Kokkos::LayoutContiguous<Layout,Stride>> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
#endif
}

#if KOKKOS_VERSION >= 40499
#include "View/Kokkos_ViewMapping.hpp"
#else
#include "impl/Kokkos_ViewMapping.hpp"
#endif

namespace Kokkos {
namespace Impl {

// Implement ViewOffset for LayoutContiguous
template < class Dimension , class Layout , unsigned Stride >
struct ViewOffset<Dimension, LayoutContiguous<Layout,Stride>, void>
  : public ViewOffset<Dimension,Layout> {
public:

  // Would like to use inherited constructors, but gcc 4.7 doesn't support it
  //using ViewOffset<Dimension,Layout>::ViewOffset;

  typedef ViewOffset<Dimension,Layout> Base;

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  // All constructors take one or two arguments

  template <typename Arg1>
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset(const Arg1& arg1) : Base(arg1) {}

  template <typename Arg1, typename Arg2>
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset(const Arg1& arg1, const Arg2& arg2) : Base(arg1,arg2) {}
};

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
} // namespace Kokkos

#endif // #ifndef KOKKOS_LAYOUT_CONTIGUOUS_HPP
