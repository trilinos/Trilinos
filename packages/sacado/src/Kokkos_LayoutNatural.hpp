// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_LAYOUT_NATURAL_HPP
#define KOKKOS_LAYOUT_NATURAL_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_Macros.hpp"
#include "Kokkos_Layout.hpp"
#include "Kokkos_LayoutContiguous.hpp" // for inner_layout<>

namespace Kokkos {

// Natural layout for scalar types -- equivalent to the wrapped
// layout type
template <typename Layout>
struct LayoutNatural : public Layout {

  //! Tag this class as a kokkos array layout
  typedef LayoutNatural array_layout ;

  LayoutNatural( LayoutNatural const & ) = default ;
  LayoutNatural( LayoutNatural && ) = default ;
  LayoutNatural & operator = ( LayoutNatural const & ) = default ;
  LayoutNatural & operator = ( LayoutNatural && ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr LayoutNatural(
    size_t N0 = 0 , size_t N1 = 0 , size_t N2 = 0 , size_t N3 = 0
  , size_t N4 = 0 , size_t N5 = 0 , size_t N6 = 0 , size_t N7 = 0 )
    : Layout( N0 , N1 , N2 , N3 , N4 , N5 , N6 , N7 ) {}
};

// Is Layout == LayoutNatural<L> for some L
template <class Layout>
struct is_layout_natural {
  static const bool value = false;
};

template <class Layout>
struct is_layout_natural< LayoutNatural<Layout> > {
  static const bool value = true;
};

template <class Layout>
struct inner_layout< LayoutNatural<Layout> > {
  typedef Layout type;
};

} // namespace Kokkos

// Make LayoutNatural<Layout> equivalent to Layout
namespace std {

  template <class Layout>
  struct is_same< Kokkos::LayoutNatural<Layout>, Layout> {
    static const bool value = true;
  };

  template <class Layout>
  struct is_same< Layout, Kokkos::LayoutNatural<Layout> > {
    static const bool value = true;
  };

}

#include "impl/Kokkos_ViewMapping.hpp"

namespace Kokkos {
namespace Impl {

// Implement ViewOffset for LayoutNatural
template < class Dimension , class Layout >
struct ViewOffset<Dimension, LayoutNatural<Layout>, void>
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

} // namespace Impl
} // namespace Kokkos

#endif // #ifndef KOKKOS_LAYOUT_NATURAL_HPP
