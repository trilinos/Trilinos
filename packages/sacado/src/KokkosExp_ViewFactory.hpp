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

#ifndef KOKKOS_EXPERIMENTAL_VIEW_FACTORY_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_FACTORY_HPP

// This only works with the experimental view enabled
#include "Sacado_ConfigDefs.h"
#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#include "Sacado_Traits.hpp"
#include "KokkosExp_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"

namespace Kokkos {

namespace Impl {

// Class to determine the value_type for a view as a function of one or more
// input views
template <class ... ViewPack>
struct ViewFactoryType {};

template <class View>
struct ViewFactoryType<View> {
  typedef typename View::value_type type;
};

template <class View, class ... ViewPack>
struct ViewFactoryType<View,ViewPack...> {
  typedef typename Sacado::Promote<
    typename View::value_type,
    typename ViewFactoryType<ViewPack...>::type
    >::type type;
};

}

// Function to compute the scalar dimension (e.g., Fad dimesion) from one or
// more views.  It relies on the overload for a single view provided by Sacado
template <class View, class ... ViewPack>
unsigned dimension_scalar(const View& v, const ViewPack&... views) {
  const unsigned dim0 = dimension_scalar(v);
  const unsigned dim1 = dimension_scalar(views...);
  return dim0 >= dim1 ? dim0 : dim1 ;
}

// Traits class used to create a view for a given rank and dimension as a
// function of one or more views.  The value_type for the view is determined
// by value_type, and the view is created through the create_view() function.
// The calling code must determine the rank and dimensions of the view to create
// however internal Sacado dimension will be determined automatically.
template <class ... ViewPack>
struct ViewFactory {

  typedef typename Impl::ViewFactoryType<ViewPack...>::type value_type;

  template <class ResultView, class CtorProp, class ... Dims>
  static ResultView
  create_view(const ViewPack& ... views,
              const CtorProp& prop,
              const Dims ... dims) {
    typename ResultView::array_layout layout(dims...);
    const unsigned rank = computeRank(layout);
    layout.dimension[rank] = dimension_scalar(views...);
    return ResultView(prop, layout);
  }

  // Compute the view rank from the dimension arguments
  // This allows the code to work for both static and dynamic-rank views
  template <typename Layout>
  static size_t
  computeRank(const Layout& layout) {
    return computeRank( layout.dimension[0], layout.dimension[1],
                        layout.dimension[2], layout.dimension[3],
                        layout.dimension[4], layout.dimension[5],
                        layout.dimension[6], layout.dimension[7] );
  }

  // Compute the view rank from the dimension arguments
  // This allows the code to work for both static and dynamic-rank views
  static size_t
  computeRank(
    const size_t N0, const size_t N1, const size_t N2, const size_t N3,
    const size_t N4, const size_t N5, const size_t N6, const size_t N7 ) {
    return  ( (N7 == 0) ?
            ( (N6 == 0) ?
            ( (N5 == 0) ?
            ( (N4 == 0) ?
            ( (N3 == 0) ?
            ( (N2 == 0) ?
            ( (N1 == 0) ?
            ( (N0 == 0) ? 0 : 1 ) : 2 ) : 3 ) : 4 ) : 5 ) : 6 ) : 7 ) : 8 );
  }

};

}

#endif

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_FACTORY_HPP */
