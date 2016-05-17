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

#ifndef KOKKOS_EXPERIMENTAL_LAYOUT_CONTIGUOUS_HPP
#define KOKKOS_EXPERIMENTAL_LAYOUT_CONTIGUOUS_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_Macros.hpp"
#include "Kokkos_Layout.hpp"

namespace Kokkos {

// Contiguous layout for scalar types -- equivalent to the wrapped
// layout type
template <typename Layout>
struct LayoutContiguous : public Layout {

  //! Tag this class as a kokkos array layout
  typedef LayoutContiguous array_layout ;

  LayoutContiguous( LayoutContiguous const & ) = default ;
  LayoutContiguous( LayoutContiguous && ) = default ;
  LayoutContiguous & operator = ( LayoutContiguous const & ) = default ;
  LayoutContiguous & operator = ( LayoutContiguous && ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr LayoutContiguous(
    size_t N0 = 0 , size_t N1 = 0 , size_t N2 = 0 , size_t N3 = 0
  , size_t N4 = 0 , size_t N5 = 0 , size_t N6 = 0 , size_t N7 = 0 )
    : Layout( N0 , N1 , N2 , N3 , N4 , N5 , N6 , N7 ) {}
};

} // namespace Kokkos

// Make LayoutContiguous<Layout> equivalent to Layout
namespace std {

  template <class Layout>
  struct is_same< Kokkos::LayoutContiguous<Layout>, Layout> {
    static const bool value = true;
  };

  template <class Layout>
  struct is_same< Layout, Kokkos::LayoutContiguous<Layout> > {
    static const bool value = true;
  };

}

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#include "impl/KokkosExp_ViewMapping.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

// Implement ViewOffset for LayoutContiguous
template < class Dimension , class Layout >
struct ViewOffset<Dimension, LayoutContiguous<Layout>, void>
  : public ViewOffset<Dimension,Layout> {
public:
  using ViewOffset<Dimension,Layout>::ViewOffset;
};

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

#endif // defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#endif // #ifndef KOKKOS_EXPERIMENTAL_LAYOUT_CONTIGUOUS_HPP
