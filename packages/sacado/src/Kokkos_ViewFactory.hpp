// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_FACTORY_HPP
#define KOKKOS_VIEW_FACTORY_HPP

#include <type_traits>

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

    using nc_value_type = typename ResultView::non_const_value_type;
    constexpr bool is_scalar = Sacado::IsScalarType<nc_value_type>::value;
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
      layout = Impl::reconstructLayout(layout, r);
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

namespace Impl {
  // Helper type trait to determine type of resulting DynRankView from
  // createDynRankView below
  template <typename InputView>
  struct ResultDynRankView {
    // Allow for use of LayoutStride in InputViewType.  We don't want to create
    // a new view with LayoutStride, so replace it with the default layout
    // instead.
    using input_value    = typename InputView::non_const_value_type;
    using input_layout   = typename InputView::array_layout;
    using input_device   = typename InputView::device_type;
    using default_layout = typename input_device::execution_space::array_layout;
    using result_layout  =
      typename std::conditional<
        std::is_same< input_layout, Kokkos::LayoutStride >::value,
        default_layout,
        input_layout >::type;
    using type =
      Kokkos::DynRankView<input_value, result_layout, input_device>;
  };

}

//! Wrapper to simplify use of Sacado ViewFactory
template <typename InputViewType, typename CtorProp, typename ... Dims >
typename std::enable_if<
  is_view<InputViewType>::value || is_dyn_rank_view<InputViewType>::value,
  typename Impl::ResultDynRankView<InputViewType>::type
  >::type
createDynRankView(const InputViewType& a,
                  const CtorProp& prop,
                  const Dims... dims)
{
  using ResultViewType = typename Impl::ResultDynRankView<InputViewType>::type;
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

#endif /* #ifndef KOKKOS_VIEW_FACTORY_HPP */
