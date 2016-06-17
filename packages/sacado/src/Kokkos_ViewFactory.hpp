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

#ifndef KOKKOS_VIEW_FACTORY_HPP
#define KOKKOS_VIEW_FACTORY_HPP

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

    using value_type = typename ResultView::non_const_value_type;
    constexpr bool is_scalar = Sacado::IsScalarType<value_type>::value;
    constexpr bool is_dyn_rank = is_dyn_rank_view<ResultView>::value;

    // rank == number of arguments
    constexpr unsigned rank = sizeof...(Dims);

    // Check rank is valid
    static_assert( rank <= 7, "Invalid rank...too many dimension arguments" );

    // Create layout from our dimension arguments
    typename ResultView::array_layout layout(dims...);

    // Set scalar dimension
    layout.dimension[rank] = dimension_scalar(views...);

    // Handle the case where all of the input view's are scalar's, but the
    // result isn't (e.g., a Fad), in which case we have to specify a valid
    // scalar dimension
    if (!is_scalar && layout.dimension[rank] == 0)
      layout.dimension[rank] = 1;

    // Reconstruct layout for dynamic rank
    if (is_dyn_rank) {
      constexpr unsigned r = is_scalar ? rank : rank + 1;
      layout = Experimental::Impl::reconstructLayout(layout, r);
    }

    return ResultView(prop, layout);
  }

};

//! Wrapper to simplify use of Sacado ViewFactory
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
typename std::enable_if<
  is_view<InputViewType>::value || is_dyn_rank_view<InputViewType>::value,
  ResultViewType>::type
createDynRankViewWithType(const InputViewType& a,
                          const CtorProp& prop,
                          const Dims... dims)
{
  using view_factory = Kokkos::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory
template <typename InputViewType, typename CtorProp, typename ... Dims >
typename std::enable_if<
  is_view<InputViewType>::value || is_dyn_rank_view<InputViewType>::value,
  Kokkos::DynRankView<typename InputViewType::non_const_value_type,
                      typename InputViewType::array_layout,
                      typename InputViewType::device_type> >::type
createDynRankView(const InputViewType& a,
                  const CtorProp& prop,
                  const Dims... dims)
{
  using ResultViewType =
    Kokkos::DynRankView<typename InputViewType::non_const_value_type,
                        typename InputViewType::array_layout,
                        typename InputViewType::device_type>;
  return createDynRankViewWithType<ResultViewType>(a, prop, dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
typename std::enable_if<
  is_view<InputViewType>::value || is_dyn_rank_view<InputViewType>::value,
  ResultViewType>::type
createViewWithType(const InputViewType& a,
                   const CtorProp& prop,
                   const Dims... dims)
{
  using view_factory = Kokkos::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

}

#endif

#endif /* #ifndef KOKKOS_VIEW_FACTORY_HPP */
